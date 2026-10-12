"""0.2.0 commit 2b: a claim is about one case, not about the whole upload.

A cohort of unrelated people used to be intersected as if it were one family, so 12
samples that were each 93% callable produced a BLOCKING "exploratory" at 38%. That
number answers "the same variant in all 12", which is check 5's pooled sweep, not
"was each sample searched". The claim is now made per analysis unit (a PED family,
or one sample), the headline is the worst unit, and the pooled models stay for check 5.
"""

from __future__ import annotations

import json
import random

import pytest

from admissible.checks.callability import CallabilityConfig, check_callability, load_coverage
from admissible.cli import main
from admissible.model import Status
from admissible.ped import read_ped

T = 100_000


def _cov(tmp_path, label, ivs):
    p = tmp_path / f"{label}.bed"
    p.write_text("".join(f"chr1\t{s}\t{e}\n" for s, e in ivs))
    return p


def _cohort_cov(tmp_path, labels, seed=1):
    """Each sample ~92% callable, with gaps at different places."""
    rng = random.Random(seed)
    spec = {}
    for lab in labels:
        gaps = sorted({rng.randrange(0, T - 3000) for _ in range(3)})
        ivs, pos = [], 0
        for g in gaps:
            ivs.append((pos, g))
            pos = g + 2500
        ivs.append((pos, T))
        spec[lab] = _cov(tmp_path, lab, [(s, e) for s, e in ivs if e > s])
    return spec


def _run(tmp_path, spec, ped_text=None, **cfg):
    (tmp_path / "t.bed").write_text(f"chr1\t0\t{T}\n")
    ped = None
    if ped_text is not None:
        (tmp_path / "x.ped").write_text(ped_text)
        ped = read_ped(tmp_path / "x.ped")
    return check_callability(
        load_coverage(spec), ped, target_bed=tmp_path / "t.bed",
        cfg=CallabilityConfig(**cfg) if cfg else None,
    )


def _one_family_each(labels):
    return "#FID\tIID\tPAT\tMAT\tSEX\tPHENO\n" + "".join(
        f"F{i}\t{lab}\t0\t0\t0\t2\n" for i, lab in enumerate(labels)
    )


LABELS = [f"C{i}" for i in range(12)]


# --- a cohort -----------------------------------------------------------------


def test_a_cohort_of_unrelated_cases_is_judged_per_case_not_as_one_intersection(tmp_path):
    res = _run(tmp_path, _cohort_cov(tmp_path, LABELS), _one_family_each(LABELS))
    m = res.metrics
    assert m["analysis"]["mode"] == "family" and m["analysis"]["n_units"] == 12
    # every case is searched adequately, so the run passes ...
    assert res.status is Status.PASS
    assert m["gating_fraction"] > 0.9
    assert m["cohort"]["units_by_rung"] == {"SEARCHED_MOST": 12}
    # ... although what is callable in all 12 at once is under half (the old verdict)
    assert m["models"]["recessive"]["fraction_of_target"] < 0.5
    assert m["joint_callable_fraction"] < 0.5


def test_the_pooled_models_check_5_reads_are_unchanged_by_units(tmp_path):
    res = _run(tmp_path, _cohort_cov(tmp_path, LABELS), _one_family_each(LABELS))
    rec = res.metrics["models"]["recessive"]
    assert sorted(rec["required_samples"]) == sorted(LABELS)
    assert rec["status"] == "computed" and rec["fraction_of_target"] is not None


def test_independent_makes_each_sample_its_own_case_without_a_pedigree(tmp_path):
    res = _run(tmp_path, _cohort_cov(tmp_path, LABELS), None, independent=True)
    assert res.metrics["analysis"]["mode"] == "independent"
    assert [u["unit"] for u in res.metrics["units"]] == LABELS
    assert all(u["affected"] == [u["unit"]] for u in res.metrics["units"])
    assert res.status is Status.PASS
    assert any("--independent" in n for n in res.notes)


def test_without_a_pedigree_the_default_stays_one_group_and_says_how_to_change_it(tmp_path):
    """The safe reading for a real family whose PED was not supplied."""
    res = _run(tmp_path, _cohort_cov(tmp_path, LABELS), None)
    assert res.metrics["analysis"]["n_units"] == 1
    assert res.status is Status.FAIL
    assert any("--independent" in n for n in res.notes)


def test_the_worst_case_drives_the_headline_and_is_named(tmp_path):
    spec = _cohort_cov(tmp_path, LABELS)
    spec["C7"] = _cov(tmp_path, "C7", [(0, 30_000)])  # one sample is 30% callable
    res = _run(tmp_path, spec, None, independent=True)
    m = res.metrics
    assert res.status is Status.FAIL
    assert m["gating_unit"] == "C7" and m["gating_fraction"] == pytest.approx(0.30)
    assert m["cohort"]["units_by_rung"] == {"SEARCHED_MOST": 11, "EXPLORATORY": 1}
    assert m["cohort"]["worst_unit"] == "C7"
    f = next(x for x in res.findings if x.code == "CALLABILITY_EXPLORATORY")
    assert "C7" in f.message and "11 SEARCHED_MOST" in f.message
    assert res.summary.startswith("0.30 worst of 12 units (C7)")


# --- families -----------------------------------------------------------------

TWO_FAMILIES = (
    "#FID\tIID\tPAT\tMAT\tSEX\tPHENO\n"
    "A\ta_dad\t0\t0\t1\t1\nA\ta_mum\t0\t0\t2\t1\nA\ta_kid\ta_dad\ta_mum\t1\t2\n"
    "B\tb_dad\t0\t0\t1\t1\nB\tb_mum\t0\t0\t2\t1\nB\tb_kid\tb_dad\tb_mum\t2\t2\n"
)


def test_each_family_is_held_to_its_own_parents(tmp_path):
    spec = {
        "a_kid": _cov(tmp_path, "a_kid", [(0, 95_000)]),
        "a_dad": _cov(tmp_path, "a_dad", [(0, 90_000)]),
        "a_mum": _cov(tmp_path, "a_mum", [(0, 90_000)]),
        "b_kid": _cov(tmp_path, "b_kid", [(0, 95_000)]),
        "b_dad": _cov(tmp_path, "b_dad", [(0, 40_000)]),   # family B's father is poor
        "b_mum": _cov(tmp_path, "b_mum", [(0, 90_000)]),
    }
    res = _run(tmp_path, spec, TWO_FAMILIES)
    units = {u["unit"]: u for u in res.metrics["units"]}
    assert units["A"]["models"]["recessive"]["fraction_of_target"] == pytest.approx(0.90)
    assert units["B"]["models"]["recessive"]["fraction_of_target"] == pytest.approx(0.40)
    assert units["A"]["models"]["recessive"]["required_samples"] == ["a_kid", "a_dad", "a_mum"]
    assert res.metrics["gating_unit"] == "B"
    assert res.metrics["gating_fraction"] == pytest.approx(0.40)


def test_one_family_missing_a_parent_is_unknown_without_hiding_the_other(tmp_path):
    spec = {
        "a_kid": _cov(tmp_path, "a_kid", [(0, 95_000)]),
        "a_dad": _cov(tmp_path, "a_dad", [(0, 95_000)]),
        "a_mum": _cov(tmp_path, "a_mum", [(0, 95_000)]),
        "b_kid": _cov(tmp_path, "b_kid", [(0, 95_000)]),   # B's parents have no files
    }
    res = _run(tmp_path, spec, TWO_FAMILIES)
    units = {u["unit"]: u for u in res.metrics["units"]}
    assert units["A"]["models"]["recessive"]["status"] == "computed"
    assert units["B"]["models"]["recessive"]["status"] == "unknown"
    assert res.status is Status.WARN  # 0.95 would be a PASS were B not unassessed
    f = next(x for x in res.findings if x.code == "MODEL_REQUIRED_SAMPLE_WITHOUT_COVERAGE")
    assert set(f.evidence["units"]) == {"B"}
    assert f.evidence["units"]["B"]["recessive"] == ["b_dad", "b_mum"]
    assert "--affected-only" in f.evidence["next_steps"][0]


def test_a_single_family_is_identical_to_the_pooled_models(tmp_path):
    spec = {
        "KID": _cov(tmp_path, "KID", [(0, 95_000)]),
        "DAD": _cov(tmp_path, "DAD", [(0, 80_000)]),
        "MUM": _cov(tmp_path, "MUM", [(0, 70_000)]),
    }
    ped = ("#FID\tIID\tPAT\tMAT\tSEX\tPHENO\nF\tDAD\t0\t0\t1\t1\nF\tMUM\t0\t0\t2\t1\n"
           "F\tKID\tDAD\tMUM\t1\t2\n")
    res = _run(tmp_path, spec, ped)
    assert len(res.metrics["units"]) == 1
    assert res.metrics["units"][0]["models"] == res.metrics["models"]


# --- proband only ---------------------------------------------------------------

TRIO = (
    "#FID\tIID\tPAT\tMAT\tSEX\tPHENO\nF\tDAD\t0\t0\t1\t1\nF\tMUM\t0\t0\t2\t1\n"
    "F\tKID\tDAD\tMUM\t1\t2\n"
)


def test_a_proband_whose_parents_are_not_sequenced_can_still_get_a_result(tmp_path):
    spec = {"KID": _cov(tmp_path, "KID", [(0, 95_000)])}
    strict = _run(tmp_path, spec, TRIO)
    assert strict.metrics["models"]["recessive"]["status"] == "unknown"
    assert strict.status is Status.WARN

    res = _run(tmp_path, spec, TRIO, require_relatives=False)
    rec = res.metrics["models"]["recessive"]
    assert rec["status"] == "computed" and rec["required_samples"] == ["KID"]
    assert rec["relatives_required"] is False
    assert any("affected-only" in n for n in rec["notes"])
    assert res.metrics["analysis"]["relatives_required"] is False
    assert any("affected-only mode" in n for n in res.notes)
    # de novo needs parents, so it is not applicable here, and not "missing a file"
    dn = res.metrics["models"]["de-novo"]
    assert dn["status"] == "unknown" and dn["missing_coverage"] == []
    assert not any(f.code == "MODEL_REQUIRED_SAMPLE_WITHOUT_COVERAGE" for f in res.findings)
    assert res.status is Status.PASS


def test_affected_only_ignores_parents_even_when_their_coverage_is_supplied(tmp_path):
    spec = {
        "KID": _cov(tmp_path, "KID", [(0, 95_000)]),
        "DAD": _cov(tmp_path, "DAD", [(0, 10_000)]),
        "MUM": _cov(tmp_path, "MUM", [(0, 10_000)]),
    }
    res = _run(tmp_path, spec, TRIO, require_relatives=False)
    assert res.metrics["models"]["recessive"]["fraction_of_target"] == pytest.approx(0.95)


def test_independent_ignores_the_pedigree_structure(tmp_path):
    spec = {"KID": _cov(tmp_path, "KID", [(0, 95_000)])}
    res = _run(tmp_path, spec, TRIO, independent=True)
    assert [u["unit"] for u in res.metrics["units"]] == ["KID"]
    assert res.metrics["models"]["recessive"]["required_samples"] == ["KID"]
    assert res.status is Status.PASS


# --- a unit nothing can be said about --------------------------------------------


def test_a_case_with_no_coverage_blocks_a_pass_and_is_counted(tmp_path):
    ped = _one_family_each(["C0", "C1", "C2"])
    spec = {"C0": _cov(tmp_path, "C0", [(0, 95_000)]), "C1": _cov(tmp_path, "C1", [(0, 95_000)])}
    res = _run(tmp_path, spec, ped)  # C2 is affected and has no file
    units = {u["unit"]: u for u in res.metrics["units"]}
    assert units["F2"]["status"] == "unknown" and units["F2"]["claim_rung"] == "UNKNOWN"
    assert res.metrics["cohort"]["n_units_unknown"] == 1
    assert res.status is Status.WARN
    assert any(f.code == "AFFECTED_WITHOUT_COVERAGE" for f in res.findings)


# --- the CLI ------------------------------------------------------------------------


def test_cli_flags_reach_the_check_and_the_json(tmp_path, capsys):
    spec = _cohort_cov(tmp_path, LABELS[:4])
    (tmp_path / "t.bed").write_text(f"chr1\t0\t{T}\n")
    out = tmp_path / "r.json"
    flat = tmp_path / "flat.ped"
    flat.write_text(
        "".join(f"{lab}\t{lab}\t0\t0\t0\t2\n" for lab in spec)
    )
    args = [
        "callability", "--target", str(tmp_path / "t.bed"), "--ped", str(flat),
        "--json", str(out),
    ]
    for lab, p in spec.items():
        args += ["--coverage", f"{lab}={p}"]
    assert main(args + ["--independent"]) == 0
    c = next(x for x in json.loads(out.read_text())["checks"] if x["check"] == "callability")
    assert c["metrics"]["analysis"]["mode"] == "independent"
    assert c["metrics"]["analysis"]["n_units"] == 4
    assert c["status"] == "PASS"

    (tmp_path / "x.ped").write_text(TRIO)
    out2 = tmp_path / "r2.json"
    main([
        "callability", "--target", str(tmp_path / "t.bed"), "--ped", str(tmp_path / "x.ped"),
        "--coverage", f"KID={spec['C0']}", "--affected-only", "--json", str(out2),
    ])
    c2 = next(x for x in json.loads(out2.read_text())["checks"] if x["check"] == "callability")
    assert c2["metrics"]["analysis"]["relatives_required"] is False
    assert c2["metrics"]["models"]["recessive"]["status"] == "computed"


# --- fuzz ------------------------------------------------------------------------------


def _random_ped(rng):
    ids = [f"S{i}" for i in range(rng.randrange(1, 9))]
    rows = ["#FID\tIID\tPAT\tMAT\tSEX\tPHENO"]
    for i in ids:
        def parent(i=i):
            r = rng.random()
            return "0" if r < 0.5 else i if r < 0.6 else "GHOST" if r < 0.7 else rng.choice(ids)

        rows.append(
            f"{rng.choice(['F', 'F', 'G', 'H'])}\t{i}\t{parent()}\t{parent()}\t"
            f"{rng.choice(['0', '1', '2'])}\t{rng.choice(['1', '2', '2', '-9'])}"
        )
    return "\n".join(rows) + "\n", ids


def test_fuzz_modes_never_raise_and_the_headline_is_the_worst_unit(tmp_path):
    rng = random.Random(20261013)
    reached = set()
    for trial in range(150):
        d = tmp_path / f"t{trial}"
        d.mkdir()
        text, ids = _random_ped(rng)
        labels = [i for i in ids if rng.random() < 0.75] or ids[:1]
        spec = {}
        good = rng.random() < 0.5
        for lab in labels:
            end = rng.randrange(int(T * 0.85), T + 1) if good else rng.randrange(1, T)
            spec[lab] = _cov(d, lab, [(0, end)])
        res = _run(
            d, spec, text if rng.random() < 0.8 else None,
            independent=rng.random() < 0.4, require_relatives=rng.random() < 0.6,
        )
        reached.add(res.status.value)
        units = res.metrics.get("units", [])
        if not units:
            continue
        assert res.metrics["analysis"]["n_units"] == len(units)
        c = res.metrics.get("cohort")
        if c:
            assert sum(c["units_by_rung"].values()) == len(units)
        assessable = [u for u in units if u["gating_fraction"] is not None]
        if assessable and res.metrics.get("gating_fraction") is not None:
            assert res.metrics["gating_fraction"] == min(u["gating_fraction"] for u in assessable)
        if res.status is Status.PASS:
            assert res.metrics["gating_fraction"] >= CallabilityConfig().strong
            assert all(u["status"] == "computed" for u in units)
    assert {"PASS", "FAIL", "WARN"} <= reached, reached
