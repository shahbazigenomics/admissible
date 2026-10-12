"""0.2.0 commit 3: named regions, the exclusion BED, scope and kit sanity.

One number per run averages a dead exon away.  These tests pin the behaviour that
replaces it: one row per named region, the worst region per gene, a gene claim
that a wholly unsearched exon caps, NOT_REPORTED kept apart from "uncallable", an
optional exclusion BED reported with and without, and two plausibility warnings.
Everything is also fuzzed with malformed input, because none of it may raise.
"""

from __future__ import annotations

import gzip
import json
import random

import pytest

from admissible.checks.callability import (
    CallabilityConfig,
    check_callability,
    load_coverage,
)
from admissible.cli import main
from admissible.intervals import IntervalIndex, intersect, merge, total_bp
from admissible.model import Status

# Four exons of GENE_A (the second has no callable base) and two of GENE_B.
TARGET = (
    "chr1\t0\t1000\tGENE_A_exon1\n"
    "chr1\t2000\t3000\tGENE_A_exon2\n"
    "chr1\t4000\t5000\tGENE_A_exon3\n"
    "chr1\t6000\t7000\tGENE_A_exon4\n"
    "chr1\t8000\t9000\tGENE_B_1\n"
    "chr1\t10000\t11000\tGENE_B_2\n"
)


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


def _run(tmp_path, cov_lines, target=TARGET, exclude=None, **cfg):
    cov = _write(tmp_path, "S1.bed", cov_lines)
    tgt = _write(tmp_path, "t.bed", target)
    return check_callability(
        load_coverage({"S1": cov}),
        None,
        target_bed=tgt,
        cfg=CallabilityConfig(**cfg),
        exclude_bed=exclude,
    )


# Everything callable except GENE_A_exon2 (a gap), but with a quantize label so the
# file states the gap rather than just omitting it.
FULL_BUT_EXON2 = (
    "chr1\t0\t2000\t10:inf\n"
    "chr1\t2000\t3000\t0:10\n"
    "chr1\t3000\t12000\t10:inf\n"
)


def _block(res):
    return res.metrics["by_region"]["units"]["all"]


def _row(res, name):
    return next(r for r in _block(res)["regions"] if r["region"] == name)


def _gene(res, name):
    return next(g for g in _block(res)["genes"] if g["gene"] == name)


# --- problem 1: named-region output ---------------------------------------------------


def test_default_output_has_no_region_block_and_is_otherwise_unchanged(tmp_path):
    plain = _run(tmp_path, FULL_BUT_EXON2)
    with_regions = _run(tmp_path, FULL_BUT_EXON2, by_region=True)
    assert "by_region" not in plain.metrics and "exclusion" not in plain.metrics
    # asking for regions adds a key and changes no existing number
    a, b = dict(plain.metrics), dict(with_regions.metrics)
    b.pop("by_region")
    assert a == b
    assert plain.status == with_regions.status and plain.summary == with_regions.summary


def test_one_row_per_region_and_the_dead_exon_is_visible(tmp_path):
    res = _run(tmp_path, FULL_BUT_EXON2, by_region=True)
    blk = _block(res)
    assert blk["n_regions"] == 6 and blk["n_genes"] == 2
    dead = _row(res, "GENE_A_exon2")
    assert dead["fraction"] == 0.0 and dead["fully_uncallable"] is True
    assert dead["status"] == "ok"
    assert dead["uncallable_intervals"] == ["1:2000-3000"]
    # worst first
    assert blk["regions"][0]["region"] == "GENE_A_exon2"
    assert _row(res, "GENE_B_1")["fraction"] == 1.0


def test_worst_region_per_gene_and_gene_share_hide_nothing(tmp_path):
    res = _run(tmp_path, FULL_BUT_EXON2, by_region=True)
    a = _gene(res, "GENE_A")
    assert a["n_regions"] == 4
    assert a["worst_region"] == "GENE_A_exon2" and a["worst_fraction"] == 0.0
    assert a["fraction"] == pytest.approx(0.75)  # 3 of 4 equal-sized exons
    assert _gene(res, "GENE_B")["worst_fraction"] == 1.0


# --- problem 2: a wholly unsearched exon caps the claim --------------------------------


def test_an_uncovered_exon_caps_a_gene_that_the_fraction_alone_would_rate_higher(tmp_path):
    # 0.75 is SEARCHED_PARTLY at the default 0.80, so tighten the threshold to make
    # the fraction alone reach SEARCHED_MOST and show the cap is what holds it back.
    res = _run(tmp_path, FULL_BUT_EXON2, by_region=True, strong=0.70, weak=0.40)
    a = _gene(res, "GENE_A")
    assert a["fraction"] == pytest.approx(0.75)
    assert a["claim_rung"] == "SEARCHED_PARTLY" and a["capped"] is True
    assert "GENE_A_exon2" in a["permitted_claim"] and "capped" in a["permitted_claim"]
    assert "excluded" not in a["permitted_claim"].lower()
    # an uncapped gene is left alone
    assert _gene(res, "GENE_B")["claim_rung"] == "SEARCHED_MOST"
    assert _gene(res, "GENE_B")["capped"] is False


def test_the_cap_level_is_configurable(tmp_path):
    res = _run(
        tmp_path, FULL_BUT_EXON2, by_region=True, strong=0.70, weak=0.40,
        blind_region_cap="EXPLORATORY",
    )
    assert _gene(res, "GENE_A")["claim_rung"] == "EXPLORATORY"


def test_a_flagged_region_is_a_finding(tmp_path):
    res = _run(tmp_path, FULL_BUT_EXON2, by_region=True)
    f = next(f for f in res.findings if f.code == "TARGET_REGIONS_NOT_SEARCHED")
    assert f.evidence["units"]["all"]["n_fully_uncallable_regions"] == 1
    assert "GENE_A" in f.evidence["units"]["all"]["genes"]


# --- NOT_REPORTED is not uncallable ----------------------------------------------------


def test_a_contig_missing_from_the_coverage_file_is_not_reported_not_zero(tmp_path):
    target = TARGET + "chr2\t0\t1000\tGENE_C_exon1\n"
    res = _run(tmp_path, FULL_BUT_EXON2, target=target, by_region=True)
    r = _row(res, "GENE_C_exon1")
    assert r["status"] == "NOT_REPORTED"
    assert r["fraction"] is None and r["callable_bp"] is None
    assert r["fully_uncallable"] is False
    assert r["claim_rung"] == "UNKNOWN"
    g = _gene(res, "GENE_C")
    assert g["fraction"] is None and g["claim_rung"] == "UNKNOWN"
    assert any(f.code == "TARGET_CONTIG_NOT_IN_COVERAGE" for f in res.findings)


def test_a_contig_reported_with_nothing_callable_is_uncallable_not_unreported(tmp_path):
    target = TARGET + "chr2\t0\t1000\tGENE_C_exon1\n"
    cov = FULL_BUT_EXON2 + "chr2\t0\t5000\t0:10\n"
    res = _run(tmp_path, cov, target=target, by_region=True)
    r = _row(res, "GENE_C_exon1")
    assert r["status"] == "ok" and r["fraction"] == 0.0 and r["fully_uncallable"] is True
    assert not any(f.code == "TARGET_CONTIG_NOT_IN_COVERAGE" for f in res.findings)


def test_not_reported_regions_do_not_inflate_the_gene_share(tmp_path):
    target = "chr1\t0\t1000\tG_exon1\nchr2\t0\t1000\tG_exon2\n"
    res = _run(tmp_path, "chr1\t0\t1000\t10:inf\n", target=target, by_region=True)
    g = _gene(res, "G")
    assert g["fraction"] == 1.0 and g["n_not_reported"] == 1
    assert g["claim_rung"] == "SEARCHED_PARTLY" and g["capped"] is True


# --- region mechanics -------------------------------------------------------------------


def test_one_name_on_several_lines_is_one_region(tmp_path):
    target = "chr1\t0\t500\tX_exon1\nchr1\t500\t1000\tX_exon1\nchr1\t2000\t2500\tX_exon2\n"
    res = _run(tmp_path, "chr1\t0\t1000\t10:inf\n", target=target, by_region=True)
    assert _block(res)["n_regions"] == 2
    assert _row(res, "X_exon1")["bp"] == 1000


def test_an_unnamed_target_gets_one_row_per_line_and_a_note(tmp_path):
    target = "chr1\t0\t1000\nchr1\t2000\t3000\n"
    res = _run(tmp_path, "chr1\t0\t1000\t10:inf\n", target=target, by_region=True)
    blk = res.metrics["by_region"]
    assert blk["has_names"] is False and blk["n_unnamed"] == 2
    assert _block(res)["n_regions"] == 2
    assert any("no column-4 names" in n for n in res.notes)


def test_gene_pattern_is_configurable_and_a_bad_one_cannot_raise(tmp_path):
    target = "chr1\t0\t1000\tFOO_part1\nchr1\t2000\t3000\tFOO_part2\n"
    cov = "chr1\t0\t3000\t10:inf\n"
    default = _run(tmp_path, cov, target=target, by_region=True)
    assert _block(default)["n_genes"] == 2  # "_part1" is not a marker the default knows
    custom = _run(tmp_path, cov, target=target, by_region=True, gene_regex=r"_part\d+$")
    assert _block(custom)["n_genes"] == 1
    bad = _run(tmp_path, cov, target=target, by_region=True, gene_regex="(")
    assert bad.status is not Status.UNKNOWN or bad.metrics.get("by_region")
    assert any("not a valid regular expression" in n for n in bad.notes)
    assert _block(bad)["n_genes"] == 2


def test_long_tables_are_cut_worst_first_and_say_so(tmp_path):
    lines = "".join(f"chr1\t{i * 100}\t{i * 100 + 50}\tG{i}_exon1\n" for i in range(20))
    cov = "chr1\t0\t100\t10:inf\n"  # only the first region is covered
    res = _run(tmp_path, cov, target=lines, by_region=True, max_region_rows=3)
    blk = _block(res)
    assert blk["n_regions"] == 20 and len(blk["regions"]) == 3 and blk["truncated"] is True
    assert blk["regions"][0]["fraction"] == 0.0
    assert blk["n_fully_uncallable_regions"] == 19  # counted, not just those listed


def test_region_fractions_match_a_brute_force_bitmap(tmp_path):
    rng = random.Random(7)
    for trial in range(40):
        cov_iv = [(s := rng.randrange(0, 900), s + rng.randrange(1, 120)) for _ in range(6)]
        regions = []
        for i in range(8):
            s = rng.randrange(0, 950)
            regions.append((f"R{i}_exon1", s, s + rng.randrange(1, 100)))
        cov = "".join(f"chr1\t{s}\t{e}\t10:inf\n" for s, e in cov_iv)
        cov += "chr1\t2000\t2001\t0:10\n"
        target = "".join(f"chr1\t{s}\t{e}\t{n}\n" for n, s, e in regions)
        res = _run(tmp_path, cov, target=target, by_region=True)
        callable_pos = {p for s, e in cov_iv for p in range(s, e)}
        for n, s, e in regions:
            expect = len(callable_pos & set(range(s, e))) / (e - s)
            assert _row(res, n)["fraction"] == pytest.approx(expect), (trial, n)


# --- problem 4: exclusion BED ------------------------------------------------------------


def test_exclusion_lowers_the_fraction_and_the_unfiltered_one_is_kept(tmp_path):
    cov = "chr1\t0\t12000\t10:inf\n"
    base = _run(tmp_path, cov)
    ex = _write(tmp_path, "ex.bed", "chr1\t0\t3000\n")
    res = _run(tmp_path, cov, exclude=ex)
    m = res.metrics["models"]["recessive"]
    assert m["without_exclusion"]["fraction_of_target"] == pytest.approx(
        base.metrics["models"]["recessive"]["fraction_of_target"]
    )
    # exons 1 and 2 (2,000 of 6,000 target bp) sit inside the exclusion; exon 3.. do not
    assert m["fraction_of_target"] == pytest.approx(4 / 6)
    assert res.metrics["exclusion"]["bp_in_target"] == 2000
    assert res.metrics["gating_fraction"] == pytest.approx(4 / 6)
    # the claim follows the stricter figure
    assert res.metrics["units"][0]["claim_rung"] == "SEARCHED_PARTLY"


def test_exclusion_applies_to_region_rows_too(tmp_path):
    cov = "chr1\t0\t12000\t10:inf\n"
    ex = _write(tmp_path, "ex.bed", "chr1\t0\t500\n")
    res = _run(tmp_path, cov, exclude=ex, by_region=True)
    assert _row(res, "GENE_A_exon1")["fraction"] == pytest.approx(0.5)


def test_an_unusable_exclusion_is_a_finding_and_changes_nothing(tmp_path):
    cov = "chr1\t0\t12000\t10:inf\n"
    ex = _write(tmp_path, "ex.bed", "not\ta\tbed\n")
    res = _run(tmp_path, cov, exclude=ex)
    assert any(f.code == "EXCLUDE_BED_UNUSABLE" for f in res.findings)
    assert "exclusion" not in res.metrics
    assert "without_exclusion" not in res.metrics["models"]["recessive"]
    missing = _run(tmp_path, cov, exclude=tmp_path / "nope.bed")
    assert any(f.code == "EXCLUDE_BED_UNUSABLE" for f in missing.findings)


def test_an_exclusion_off_the_target_is_noted(tmp_path):
    ex = _write(tmp_path, "ex.bed", "chr9\t0\t1000\n")
    res = _run(tmp_path, "chr1\t0\t12000\t10:inf\n", exclude=ex)
    assert any("does not overlap the target" in n for n in res.notes)


def test_subtracting_once_from_the_joint_set_equals_subtracting_from_each_sample(tmp_path):
    from admissible.intervals import coverage_at_least, subtract

    rng = random.Random(11)

    def rand_set():
        return {"1": merge([(s := rng.randrange(0, 300), s + rng.randrange(1, 60))
                            for _ in range(rng.randrange(0, 7))])}

    for _ in range(300):
        sets = [rand_set() for _ in range(rng.randrange(1, 5))]
        ex = rand_set()
        k = rng.randrange(1, len(sets) + 1)
        once = subtract(coverage_at_least(sets, k), ex)
        each = coverage_at_least([subtract(s, ex) for s in sets], k)
        assert once == each


# --- problem 5: scope and kit sanity ------------------------------------------------------


def test_kit_mismatch_warns_for_a_kit_sized_target_and_never_for_a_panel(tmp_path):
    # Coverage lives elsewhere on the genome; the "kit" target holds none of it.
    cov = "chr1\t100000\t900000\t10:inf\nchr1\t0\t1000\t10:inf\n"
    target = "chr1\t0\t2000\tK_exon1\nchr1\t2000\t4000\tK_exon2\n"
    gated = _run(tmp_path, cov, target=target, kit_min_bp=1000)
    f = next(f for f in gated.findings if f.code == "TARGET_KIT_MISMATCH_SUSPECTED")
    assert f.evidence["calibrated"] is False and f.evidence["share_callable_in_target"]["S1"] < 0.01
    # the same data against a panel-sized target (below the gate) is not flagged
    panel = _run(tmp_path, cov, target=target, kit_min_bp=5_000_000)
    assert not any(f.code == "TARGET_KIT_MISMATCH_SUSPECTED" for f in panel.findings)
    # coverage that is the target itself is not flagged
    same = _run(tmp_path, "chr1\t0\t4000\t10:inf\n", target=target, kit_min_bp=1000)
    assert not any(f.code == "TARGET_KIT_MISMATCH_SUSPECTED" for f in same.findings)


def test_the_sanity_warning_does_not_change_status(tmp_path):
    cov = "chr1\t100000\t900000\t10:inf\nchr1\t0\t1000\t10:inf\n"
    target = "chr1\t0\t2000\tK_exon1\nchr1\t2000\t4000\tK_exon2\n"
    a = _run(tmp_path, cov, target=target, kit_min_bp=1000)
    b = _run(tmp_path, cov, target=target, kit_min_bp=10**9)
    assert a.status == b.status and a.summary == b.summary


# --- interval index ------------------------------------------------------------------------


def test_interval_index_matches_intersect():
    rng = random.Random(3)
    for _ in range(2000):
        a = {"1": [(s := rng.randrange(0, 200), s + rng.randrange(1, 30))
                   for _ in range(rng.randrange(0, 8))]}
        b = {"1": merge([(s := rng.randrange(0, 200), s + rng.randrange(1, 40))
                         for _ in range(rng.randrange(0, 5))])}
        assert IntervalIndex(a).overlap_set_bp(b) == total_bp(
            intersect({"1": merge(list(a["1"]))}, b)
        )


# --- fuzz: nothing here may raise ------------------------------------------------------------


def _junk_line(rng):
    """One BED line, usually malformed in some way a real file could be."""
    lo, hi = rng.randrange(0, 40), rng.randrange(40, 90)
    kind = rng.randrange(10)
    if kind == 0:  # wrong number of columns, odd contigs
        cols = ["chr1", "1", "chrX", "", "chrUn", "MT"]
        return "\t".join(rng.choice(cols) for _ in range(rng.randrange(1, 5)))
    if kind == 1:  # negative or inverted coordinates, blank-ish names
        a, b = rng.randrange(-5, 50), rng.randrange(-5, 50)
        return f"chr1\t{a}\t{b}\t" + rng.choice(["", " ", "G_exon1", "G_1"])
    if kind == 2:
        return "# comment"
    if kind == 3:
        return "chr1\tabc\t10\tG_exon1"
    if kind == 4:  # a very long name
        return f"chr1\t{lo}\t{hi}\t" + "x" * rng.randrange(0, 300)
    if kind == 5:  # a contig the coverage never mentions
        return f"chr2\t{lo}\t{hi}\tH_{rng.randrange(5)}"
    if kind == 6:  # a name that is only a marker
        return f"chr1\t{lo}\t{hi}\t_exon1"
    if kind == 7:  # non-ASCII name
        return f"chr1\t{lo}\t{hi}\t\u0416_{rng.randrange(9)}"
    g = f"G{rng.randrange(4)}_exon{rng.randrange(3)}"
    return f"chr1\t{rng.randrange(0, 80)}\t{rng.randrange(0, 90)}\t{g}"


def test_fuzz_region_and_exclusion_inputs_never_raise(tmp_path):
    rng = random.Random(20261014)
    statuses: set = set()
    reached: set[str] = set()
    for trial in range(250):
        tgt = "\n".join(_junk_line(rng) for _ in range(rng.randrange(0, 12)))
        ex = "\n".join(_junk_line(rng) for _ in range(rng.randrange(0, 6)))
        cov = "\n".join(
            rng.choice(["chr1\t%d\t%d\t10:inf", "chr1\t%d\t%d\t0:10", "chr1\t%d\t%d"])
            % (s := rng.randrange(0, 70), s + rng.randrange(1, 40))
            for _ in range(rng.randrange(0, 8))
        )
        tpath = _write(tmp_path, "t.bed", tgt + "\n")
        cpath = _write(tmp_path, "S1.bed", cov + "\n")
        epath = _write(tmp_path, "ex.bed", ex + "\n")
        if trial % 7 == 0:  # a corrupt gzip as the exclusion
            good = gzip.compress(b"chr1\t0\t50\n" * 30)
            epath = tmp_path / "ex.bed.gz"
            epath.write_bytes(good[:-12] + bytes(rng.randrange(256) for _ in range(12)))
        res = check_callability(
            load_coverage({"S1": cpath}), None, target_bed=tpath,
            cfg=CallabilityConfig(
                by_region=rng.random() < 0.8,
                gene_regex=rng.choice([r"[_.\-:](?:exon|ex|e)?\d+$", "(", "", ".*", "[", r"\d+"]),
                max_region_rows=rng.choice([0, 1, 5, 5000]),
                kit_min_bp=rng.choice([1, 10**9]),
            ),
            exclude_bed=rng.choice([epath, None, tmp_path / "missing.bed"]),
        )
        statuses.add(res.status)
        json.dumps(res.to_dict())  # the output must always be serialisable
        reached.update(f.code for f in res.findings)
        if "exclusion" in res.metrics:
            reached.add("exclusion-applied")
        blk = res.metrics.get("by_region")
        if blk:
            for u in blk["units"].values():
                if u.get("truncated"):
                    reached.add("truncated")
                if u.get("n_not_reported_regions"):
                    reached.add("not-reported")
                if u.get("n_fully_uncallable_regions"):
                    reached.add("fully-uncallable")
                for r in u.get("regions", []):
                    if r["fraction"] is not None:
                        assert 0.0 <= r["fraction"] <= 1.0
                        assert r["callable_bp"] <= r["bp"]
                    else:
                        assert r["status"] != "ok" or r["bp"] <= 0
    assert len(statuses) >= 2
    # a fuzz test that never reaches the interesting paths proves nothing
    assert {
        "exclusion-applied", "EXCLUDE_BED_UNUSABLE", "truncated", "not-reported",
        "fully-uncallable", "TARGET_KIT_MISMATCH_SUSPECTED",
        "TARGET_CONTIG_NOT_IN_COVERAGE",
    } <= reached, reached


# --- CLI ----------------------------------------------------------------------------------------


def test_cli_flags_reach_the_check_and_the_json(tmp_path):
    cov = _write(tmp_path, "S1.bed", FULL_BUT_EXON2)
    tgt = _write(tmp_path, "t.bed", TARGET)
    ex = _write(tmp_path, "ex.bed", "chr1\t0\t500\n")
    ped = _write(tmp_path, "x.ped", "F\tS1\t0\t0\t0\t2\n")
    out = tmp_path / "r.json"
    main([
        "callability", "--ped", str(ped), "--coverage", f"S1={cov}", "--target", str(tgt),
        "--by-region", "--exclude-bed", str(ex), "--gene-regex", r"_exon\d+$|_\d+$",
        "--json", str(out),
    ])
    c = next(x for x in json.loads(out.read_text())["checks"] if x["check"] == "callability")
    assert "by_region" in c["metrics"] and "exclusion" in c["metrics"]
    unit = next(iter(c["metrics"]["by_region"]["units"].values()))
    assert unit["n_genes"] == 2
