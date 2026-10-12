"""0.2.0 commit 2: a model's callable fraction is over every sample it needs.

Before this, only affected samples were intersected, so a recessive or de novo
fraction could read 95% while a parent was uncallable at most of those positions.
Check 5 drops any site at which a relative it relies on was not observed, so the
fraction for such a model has to be held to the same samples.
"""

from __future__ import annotations

import random

import pytest

from admissible.checks.callability import (
    CallabilityConfig,
    check_callability,
    claim_for,
    claim_rung_for,
    load_coverage,
)
from admissible.model import Status
from admissible.ped import read_ped

TRIO = (
    "#FID\tIID\tPAT\tMAT\tSEX\tPHENO\n"
    "F\tDAD\t0\t0\t1\t1\nF\tMUM\t0\t0\t2\t1\nF\tKID\tDAD\tMUM\t1\t2\n"
)
QUARTET = TRIO + "F\tSIB\tDAD\tMUM\t2\t2\n"  # two affected children, same parents
T = 1000  # target length


def run(tmp_path, ped_text, cov, cfg=None, family=None):
    """``cov`` maps label -> list of (start, end) callable intervals on chr1."""
    ped = None
    if ped_text is not None:
        pf = tmp_path / "x.ped"
        pf.write_text(ped_text)
        ped = read_ped(pf)
    (tmp_path / "target.bed").write_text(f"chr1\t0\t{T}\n")
    spec = {}
    for label, ivs in cov.items():
        p = tmp_path / f"{label}.bed"
        p.write_text("".join(f"chr1\t{s}\t{e}\n" for s, e in ivs))
        spec[label] = p
    return check_callability(
        load_coverage(spec), ped, target_bed=tmp_path / "target.bed", cfg=cfg, family_id=family
    )


def frac(res, model):
    return res.metrics["models"][model]["fraction_of_target"]


# --- the claim never says "excluded" and follows the configured thresholds ----


@pytest.mark.parametrize("f", [0.0, 0.2, 0.49, 0.5, 0.79, 0.8, 0.95, 1.0])
def test_no_rung_says_excluded_and_each_states_what_was_not_searched(f):
    cfg = CallabilityConfig()
    claim = claim_for(f, cfg)
    assert "exclud" not in claim.lower()
    assert f"{1 - f:.1%} of the target was not searched" in claim


def test_claim_follows_the_configured_thresholds_not_a_fixed_ladder():
    """0.1.1 printed "monogenic cause excluded" under a BLOCKING finding when the
    thresholds were configured above the fraction: the claim ignored the config."""
    strict = CallabilityConfig(strong=0.95, weak=0.90)
    assert claim_rung_for(0.85, strict) == "EXPLORATORY"
    assert claim_rung_for(0.85, CallabilityConfig()) == "SEARCHED_MOST"


def test_the_rung_always_agrees_with_the_finding_it_sits_in(tmp_path):
    rng = random.Random(3)
    for _ in range(60):
        weak = rng.uniform(0.1, 0.9)
        strong = rng.uniform(weak, 1.0)
        cfg = CallabilityConfig(strong=strong, weak=weak)
        end = rng.randrange(1, T + 1)
        res = run(tmp_path, None, {"A": [(0, end)], "B": [(0, end)]}, cfg=cfg)
        f = next(x for x in res.findings if x.code.startswith("CALLABILITY_"))
        want = {
            "CALLABILITY_ADEQUATE": "SEARCHED_MOST",
            "CALLABILITY_LIMITED": "SEARCHED_PARTLY",
            "CALLABILITY_EXPLORATORY": "EXPLORATORY",
        }[f.code]
        assert f.evidence["claim_rung"] == want
        assert f.evidence["permitted_claim"] == claim_for(f.evidence["gating_fraction"], cfg)


# --- parents are required where the model relies on them ----------------------


def test_a_model_that_needs_an_uncovered_parent_is_unknown_not_dropped(tmp_path):
    res = run(tmp_path, TRIO, {"KID": [(0, 950)]})
    m = res.metrics["models"]
    for name in ("recessive", "compound-het", "dominant", "x-linked-recessive",
                 "x-linked-dominant"):
        assert m[name]["status"] == "unknown", name
        assert m[name]["fraction_of_target"] is None
        assert m[name]["permitted_claim"] is None
        assert m[name]["claim_rung"] == "UNKNOWN"
        assert m[name]["missing_coverage"] == ["DAD", "MUM"]
        assert "DAD" in m[name]["reason"]
    # models that do not lean on unaffected relatives are unaffected
    for name in ("dominant-reduced-penetrance", "mitochondrial"):
        assert m[name]["status"] == "computed"
        assert m[name]["fraction_of_target"] == pytest.approx(0.95)
        assert m[name]["required_samples"] == ["KID"]
    assert m["de-novo"]["status"] == "unknown"
    assert m["de-novo"]["missing_coverage"] == ["DAD", "MUM"]


def test_unknown_models_are_reported_and_stop_a_pass(tmp_path):
    res = run(tmp_path, TRIO, {"KID": [(0, 950)]})
    assert res.status is Status.WARN  # 0.95 would otherwise be a PASS
    f = next(x for x in res.findings if x.code == "MODEL_REQUIRED_SAMPLE_WITHOUT_COVERAGE")
    assert f.evidence["models"]["recessive"] == ["DAD", "MUM"]
    assert "overstate" in f.message


def test_with_the_parents_covered_the_fraction_is_the_intersection_of_all_three(tmp_path):
    res = run(tmp_path, TRIO, {"KID": [(0, 950)], "DAD": [(0, 900)], "MUM": [(0, 700)]})
    rec = res.metrics["models"]["recessive"]
    assert rec["fraction_of_target"] == pytest.approx(0.70)
    assert rec["required_samples"] == ["KID", "DAD", "MUM"]
    assert rec["samples_required"] == 3
    # the affected-only number is kept, and is the larger, misleading one
    assert res.metrics["joint_callable_fraction"] == pytest.approx(0.95)
    # the affected-only models do not feel the parents
    assert frac(res, "mitochondrial") == pytest.approx(0.95)
    # the headline is the strictest assessable model
    assert res.metrics["gating_model"] in {"recessive", "compound-het", "dominant",
                                           "x-linked-recessive", "x-linked-dominant",
                                           "de-novo"}
    assert res.metrics["gating_fraction"] == pytest.approx(0.70)
    assert res.status is Status.WARN  # 0.70 is between the thresholds


def test_a_parent_who_is_also_affected_is_counted_once(tmp_path):
    ped = (
        "#FID\tIID\tPAT\tMAT\tSEX\tPHENO\n"
        "F\tMUM\t0\t0\t2\t2\nF\tDAD\t0\t0\t1\t1\nF\tKID\tDAD\tMUM\t1\t2\n"
    )
    res = run(tmp_path, ped, {s: [(0, 900)] for s in ("KID", "DAD", "MUM")})
    assert res.metrics["models"]["recessive"]["required_samples"] == ["MUM", "KID", "DAD"]


# --- de novo is a union over trios ------------------------------------------


def test_de_novo_needs_only_one_complete_trio_to_be_callable_at_a_position(tmp_path):
    cov = {"KID": [(0, 500)], "SIB": [(500, 1000)], "DAD": [(0, 1000)], "MUM": [(0, 1000)]}
    res = run(tmp_path, QUARTET, cov)
    # one trio covers each half: together the whole target
    assert frac(res, "de-novo") == pytest.approx(1.0)
    # while recessive, which needs both affected children callable together, gets none
    assert frac(res, "recessive") == pytest.approx(0.0)


def test_a_skipped_trio_makes_de_novo_a_floor_and_says_so(tmp_path):
    cov = {"KID": [(0, 600)], "DAD": [(0, 1000)], "MUM": [(0, 1000)]}  # SIB has no file
    res = run(tmp_path, QUARTET, cov)
    dn = res.metrics["models"]["de-novo"]
    assert dn["status"] == "computed" and dn["lower_bound"] is True
    assert dn["fraction_of_target"] == pytest.approx(0.60)
    assert dn["missing_coverage"] == ["SIB"]
    assert any("de-novo" in n and "floor" in n for n in res.notes)


def test_de_novo_is_unknown_when_no_affected_sample_has_both_parents(tmp_path):
    one_parent = "#FID\tIID\tPAT\tMAT\tSEX\tPHENO\nF\tDAD\t0\t0\t1\t1\nF\tKID\tDAD\t0\t1\t2\n"
    res = run(tmp_path, one_parent, {"KID": [(0, 900)], "DAD": [(0, 900)]})
    dn = res.metrics["models"]["de-novo"]
    assert dn["status"] == "unknown" and dn["missing_coverage"] == []
    assert "both parents" in dn["reason"]


# --- a pedigree gap is not a missing file -----------------------------------


def test_a_parent_the_pedigree_does_not_list_is_noted_and_does_not_block_a_pass(tmp_path):
    one_parent = "#FID\tIID\tPAT\tMAT\tSEX\tPHENO\nF\tDAD\t0\t0\t1\t1\nF\tKID\tDAD\t0\t1\t2\n"
    res = run(tmp_path, one_parent, {"KID": [(0, 950)], "DAD": [(0, 960)]})
    rec = res.metrics["models"]["recessive"]
    assert rec["required_samples"] == ["KID", "DAD"]
    assert any("no mother in the pedigree" in n for n in rec["notes"])
    assert frac(res, "recessive") == pytest.approx(0.95)
    assert any(f.code == "MODEL_REQUIREMENTS_UNDEFINED" for f in res.findings)
    assert not any(f.code == "MODEL_REQUIRED_SAMPLE_WITHOUT_COVERAGE" for f in res.findings)
    assert res.status is Status.PASS


def test_without_a_pedigree_everyone_is_affected_and_de_novo_is_undefined(tmp_path):
    res = run(tmp_path, None, {"A": [(0, 900)], "B": [(0, 800)]})
    assert frac(res, "recessive") == pytest.approx(0.80)
    dn = res.metrics["models"]["de-novo"]
    assert dn["status"] == "unknown" and "no pedigree" in dn["reason"]


# --- phenocopy(n-1): the exception budget ------------------------------------

THREE = (
    "#FID\tIID\tPAT\tMAT\tSEX\tPHENO\n"
    + "".join(f"F\tK{i}\t0\t0\t1\t2\n" for i in (1, 2, 3))
)


def test_phenocopy_with_all_covered_is_two_of_three(tmp_path):
    res = run(tmp_path, THREE, {"K1": [(0, 800)], "K2": [(0, 800)], "K3": [(0, 200)]})
    assert frac(res, "phenocopy(n-1)") == pytest.approx(0.80)
    assert res.metrics["models"]["phenocopy(n-1)"]["samples_required"] == 2


def test_one_affected_sample_without_coverage_spends_the_phenocopy_exception(tmp_path):
    res = run(tmp_path, THREE, {"K1": [(0, 800)], "K2": [(0, 700)]})
    ph = res.metrics["models"]["phenocopy(n-1)"]
    assert ph["status"] == "computed" and ph["lower_bound"] is True
    assert ph["fraction_of_target"] == pytest.approx(0.70)
    assert ph["missing_coverage"] == ["K3"]
    # everything that needs all three is unknown, not computed over two
    assert res.metrics["models"]["recessive"]["status"] == "unknown"


def test_two_affected_samples_without_coverage_make_phenocopy_unknown(tmp_path):
    res = run(tmp_path, THREE, {"K1": [(0, 800)]})
    assert res.metrics["models"]["phenocopy(n-1)"]["status"] == "unknown"
    # with two of three affected samples uncovered, every model is unassessable
    assert all(e["status"] == "unknown" for e in res.metrics["models"].values())
    assert res.status is Status.UNKNOWN
    assert res.metrics["gating_model"] is None and res.metrics["gating_fraction"] is None


def test_every_model_unknown_is_status_unknown_not_a_fraction(tmp_path):
    res = run(tmp_path, TRIO, {"DAD": [(0, 900)], "MUM": [(0, 900)], "KID": []})
    # KID has an empty file: it is covered nowhere, which IS an answer (0%)
    assert res.metrics["models"]["mitochondrial"]["fraction_of_target"] == 0.0
    res2 = run(tmp_path, QUARTET, {"DAD": [(0, 900)], "MUM": [(0, 900)]})
    assert res2.status is Status.UNKNOWN  # no affected sample has coverage at all


# --- consumers: check 5 reads these keys --------------------------------------


def test_the_keys_check_5_reads_survive_for_every_model(tmp_path):
    res = run(tmp_path, TRIO, {"KID": [(0, 950)], "DAD": [(0, 900)], "MUM": [(0, 900)]})
    for name, entry in res.metrics["models"].items():
        for key in ("fraction_of_target", "permitted_claim", "samples_required", "status",
                    "required_samples", "claim_rung", "missing_coverage", "lower_bound"):
            assert key in entry, (name, key)
    assert "dominant-reduced-penetrance" in res.metrics["models"]  # check 5 reports it too


# --- fuzz: nothing raises and the invariants hold on any pedigree -------------


def _random_ped(rng):
    ids = [f"S{i}" for i in range(rng.randrange(1, 8))]
    rows = ["#FID\tIID\tPAT\tMAT\tSEX\tPHENO"]
    for i in ids:
        def parent(i=i):
            r = rng.random()
            if r < 0.5:
                return "0"
            if r < 0.6:
                return i  # its own parent
            if r < 0.7:
                return "GHOST"  # not in the pedigree
            return rng.choice(ids)
        rows.append(
            f"{rng.choice(['F', 'F', 'G'])}\t{i}\t{parent()}\t{parent()}\t"
            f"{rng.choice(['0', '1', '2', 'x'])}\t{rng.choice(['0', '1', '2', '2', '-9', 'q'])}"
        )
    if rng.random() < 0.15:
        rows.append("garbage line")
    return "\n".join(rows) + "\n", ids


def test_fuzz_random_pedigrees_and_coverage_never_raise_and_stay_consistent(tmp_path):
    rng = random.Random(20261012)
    reached: dict[str, int] = {}
    for trial in range(120):
        text, ids = _random_ped(rng)
        labels = [i for i in ids if rng.random() < 0.7]
        if rng.random() < 0.2:
            labels.append("NOT_IN_PED")
        cov = {}
        good_trial = rng.random() < 0.4  # most samples well covered: reaches PASS/WARN
        for lab in labels:
            ivs = []
            if good_trial:
                ivs.append((0, rng.randrange(int(T * 0.85), T + 1)))
            else:
                for _ in range(rng.randrange(0, 4)):
                    s = rng.randrange(0, T)
                    ivs.append((s, min(T, s + rng.randrange(1, 400))))
            cov[lab] = ivs
        d = tmp_path / f"t{trial}"
        d.mkdir()
        res = run(d, text if rng.random() < 0.9 else None, cov)
        assert res.status in set(Status)
        reached[res.status.value] = reached.get(res.status.value, 0) + 1
        models = res.metrics.get("models", {})
        affected_only = res.metrics.get("joint_callable_fraction")
        for name, e in models.items():
            assert e["status"] in ("computed", "unknown")
            if e["status"] == "unknown":
                assert e["fraction_of_target"] is None and e["permitted_claim"] is None
                assert e["claim_rung"] == "UNKNOWN" and e["reason"]
            else:
                f = e["fraction_of_target"]
                assert f is not None and 0.0 <= f <= 1.0
                assert e["claim_rung"] in ("SEARCHED_MOST", "SEARCHED_PARTLY", "EXPLORATORY")
                assert "exclud" not in e["permitted_claim"].lower()
                # a model never searches more than the affected-only intersection,
                # except those that deliberately need fewer affected samples
                if (
                    affected_only is not None
                    and not e["lower_bound"]
                    and name not in ("phenocopy(n-1)", "de-novo")
                ):
                    assert f <= affected_only + 1e-9
        if res.status is Status.PASS:
            assert not any(e["missing_coverage"] and e["status"] == "unknown"
                           for e in models.values())
    # the fuzz must actually reach every outcome, or it proves little
    assert {"UNKNOWN", "PASS", "WARN", "FAIL"} <= set(reached), reached
