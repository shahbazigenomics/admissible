"""Check 4 - callability.

The central test is ``test_product_understates_the_intersection``: it constructs
coverage with the correlation structure real exomes have and shows the product of
per-sample fractions landing far below the true joint fraction. That is the error
this check exists to prevent.
"""

from __future__ import annotations

import pytest

from admissible.checks.callability import (
    CallabilityConfig,
    check_callability,
    load_coverage,
)
from admissible.intervals import coverage_at_least, intersect, merge, read_bed, total_bp
from admissible.model import Status
from admissible.ped import read_ped


def write_bed(path, rows, label=None):
    with open(path, "w") as fh:
        for chrom, start, end in rows:
            fh.write(f"{chrom}\t{start}\t{end}" + (f"\t{label}" if label else "") + "\n")
    return path


# --- interval primitives --------------------------------------------------


def test_merge_coalesces_touching_intervals():
    assert merge([(0, 10), (10, 20), (30, 40), (5, 8)]) == [(0, 20), (30, 40)]


def test_intersection_is_not_union():
    a = {"1": [(0, 100)]}
    b = {"1": [(50, 200)]}
    assert intersect(a, b) == {"1": [(50, 100)]}
    assert total_bp(intersect(a, b)) == 50


def test_coverage_at_least_k_of_n():
    s1 = {"1": [(0, 100)]}
    s2 = {"1": [(0, 60)]}
    s3 = {"1": [(0, 30)]}
    assert total_bp(coverage_at_least([s1, s2, s3], 3)) == 30
    assert total_bp(coverage_at_least([s1, s2, s3], 2)) == 60
    assert total_bp(coverage_at_least([s1, s2, s3], 1)) == 100


def test_quantize_label_filtering(tmp_path):
    p = tmp_path / "q.bed"
    with open(p, "w") as fh:
        fh.write("chr1\t0\t100\t0:10\nchr1\t100\t300\t10:inf\n")
    assert read_bed(p, keep_label="10:inf").total_bp == 200
    assert read_bed(p).total_bp == 300


def test_bed_reader_never_raises(tmp_path):
    p = tmp_path / "bad.bed"
    p.write_text("chr1\tx\ty\nchr1\t100\n\x00\x01\nchr1\t50\t10\nchr1\t0\t10\n")
    bed = read_bed(p)
    assert bed.total_bp == 10
    # non-numeric coords, a two-column row, a binary row, and end <= start
    assert bed.n_malformed == 4


# --- the check ------------------------------------------------------------


@pytest.fixture
def ped(tmp_path):
    p = tmp_path / "f.ped"
    p.write_text(
        "#FID\tIID\tPAT\tMAT\tSEX\tPHENO\n"
        "F\tDAD\t0\t0\t1\t1\nF\tMUM\t0\t0\t2\t1\n"
        "F\tKID1\tDAD\tMUM\t1\t2\nF\tKID2\tDAD\tMUM\t2\t2\nF\tKID3\tDAD\tMUM\t1\t2\n"
    )
    return read_ped(p)


def test_no_coverage_is_unknown_not_an_estimate():
    res = check_callability()
    assert res.status is Status.UNKNOWN
    assert "VCF cannot answer" in res.summary


def test_missing_target_refuses_to_invent_a_denominator(tmp_path, ped):
    spec = {}
    for s in ("KID1", "KID2", "KID3"):
        spec[s] = write_bed(tmp_path / f"{s}.bed", [("chr1", 0, 1000)])
    res = check_callability(load_coverage(spec), ped)
    assert res.status is Status.UNKNOWN
    assert any(f.code == "NO_TARGET_INTERVALS" for f in res.findings)
    assert res.metrics["joint_callable_bp"] == 1000


def test_product_understates_the_intersection(tmp_path, ped):
    """Correlated dropout: every sample fails on the same hard region.

    Target is 10,000 bp. Each affected sample covers 7,000 bp, and 6,000 bp of
    that is the *same* 6,000 bp in all three - which is how exome coverage
    actually behaves. True joint fraction 0.60; the product says 0.343.
    """
    target = write_bed(tmp_path / "target.bed", [("chr1", 0, 10_000)])
    shared = (0, 6_000)
    private = {"KID1": (6_000, 7_000), "KID2": (7_000, 8_000), "KID3": (8_000, 9_000)}
    spec = {}
    for s, (ps, pe) in private.items():
        spec[s] = write_bed(
            tmp_path / f"{s}.bed", [("chr1", *shared), ("chr1", ps, pe)]
        )

    res = check_callability(load_coverage(spec), ped, target_bed=target)
    m = res.metrics
    assert m["per_sample"]["KID1"]["fraction_of_target"] == pytest.approx(0.70)
    assert m["joint_callable_fraction"] == pytest.approx(0.60)
    assert m["product_of_per_sample_fractions"] == pytest.approx(0.343, abs=1e-3)
    # The whole point: the product is far below the truth.
    assert m["product_of_per_sample_fractions"] < m["joint_callable_fraction"] * 0.6
    assert any("must never" in n for n in res.notes)


def test_callability_is_computed_per_model(tmp_path, ped):
    target = write_bed(tmp_path / "target.bed", [("chr1", 0, 10_000)])
    spec = {
        "KID1": write_bed(tmp_path / "a.bed", [("chr1", 0, 8_000)]),
        "KID2": write_bed(tmp_path / "b.bed", [("chr1", 0, 8_000)]),
        "KID3": write_bed(tmp_path / "c.bed", [("chr1", 0, 2_000)]),
    }
    res = check_callability(load_coverage(spec), ped, target_bed=target)
    models = res.metrics["models"]
    # All three callable: only the 2,000 bp KID3 has.
    assert models["recessive"]["fraction_of_target"] == pytest.approx(0.20)
    # Two of three: the 8,000 bp KID1 and KID2 share.
    assert models["phenocopy(n-1)"]["fraction_of_target"] == pytest.approx(0.80)
    assert models["recessive"]["samples_required"] == 3
    assert models["phenocopy(n-1)"]["samples_required"] == 2


def test_low_callability_blocks_and_forbids_the_negative_claim(tmp_path, ped):
    target = write_bed(tmp_path / "target.bed", [("chr1", 0, 10_000)])
    spec = {
        s: write_bed(tmp_path / f"{s}.bed", [("chr1", 0, 2_000)])
        for s in ("KID1", "KID2", "KID3")
    }
    res = check_callability(load_coverage(spec), ped, target_bed=target)
    assert res.status is Status.FAIL
    f = next(x for x in res.findings if x.code == "CALLABILITY_EXPLORATORY")
    assert f.evidence["do_not_conclude"] == "no monogenic cause"
    assert "cannot be excluded" in f.evidence["permitted_claim"]


def test_high_callability_permits_the_strong_claim(tmp_path, ped):
    target = write_bed(tmp_path / "target.bed", [("chr1", 0, 10_000)])
    spec = {
        s: write_bed(tmp_path / f"{s}.bed", [("chr1", 0, 9_500)])
        for s in ("KID1", "KID2", "KID3")
    }
    res = check_callability(load_coverage(spec), ped, target_bed=target)
    assert res.status is Status.PASS
    f = next(x for x in res.findings if x.code == "CALLABILITY_ADEQUATE")
    assert "excluded for the models tested" in f.evidence["permitted_claim"]
    assert "do_not_conclude" not in f.evidence


def test_only_affected_samples_count(tmp_path, ped):
    target = write_bed(tmp_path / "target.bed", [("chr1", 0, 10_000)])
    spec = {
        "KID1": write_bed(tmp_path / "k1.bed", [("chr1", 0, 9_000)]),
        "KID2": write_bed(tmp_path / "k2.bed", [("chr1", 0, 9_000)]),
        "KID3": write_bed(tmp_path / "k3.bed", [("chr1", 0, 9_000)]),
        "DAD": write_bed(tmp_path / "dad.bed", [("chr1", 0, 100)]),
    }
    res = check_callability(load_coverage(spec), ped, target_bed=target)
    assert res.metrics["n_affected_with_coverage"] == 3
    assert "DAD" not in res.metrics["per_sample"]
    assert res.metrics["joint_callable_fraction"] == pytest.approx(0.90)


def test_thresholds_are_labelled_as_conventions(tmp_path, ped):
    target = write_bed(tmp_path / "target.bed", [("chr1", 0, 1_000)])
    spec = {
        s: write_bed(tmp_path / f"{s}.bed", [("chr1", 0, 900)])
        for s in ("KID1", "KID2", "KID3")
    }
    res = check_callability(load_coverage(spec), ped, target_bed=target)
    assert any("conventions" in n for n in res.notes)


def test_config_thresholds_are_honoured(tmp_path, ped):
    target = write_bed(tmp_path / "target.bed", [("chr1", 0, 1_000)])
    spec = {
        s: write_bed(tmp_path / f"{s}.bed", [("chr1", 0, 700)])
        for s in ("KID1", "KID2", "KID3")
    }
    strict = CallabilityConfig(strong=0.95, weak=0.90)
    res = check_callability(load_coverage(spec), ped, target_bed=target, cfg=strict)
    assert res.status is Status.FAIL  # 0.70 is below the stricter weak threshold


# --- check 3 off-target detection (needs the target BED) -------------------


def test_unrestricted_calling_is_detected(tmp_path):
    """80% off-target means the caller was never restricted to the kit."""
    from admissible.checks.provenance import check_provenance
    from admissible.vcfio import SiteIndex, scan_vcf

    target = write_bed(tmp_path / "t.bed", [("chr1", 0, 1000)])
    vcf = tmp_path / "wide.vcf"
    rows = "".join(
        f"chr1\t{p}\t.\tA\tG\t50\tPASS\tAC=1\tGT:DP\t0/1:30\n"
        for p in list(range(100, 300, 100)) + list(range(50_000, 50_800, 100))
    )
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=249250621>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" + rows
    )
    res = check_provenance([scan_vcf(vcf, SiteIndex())], target_bed=str(target))
    f = next(x for x in res.findings if x.code == "NOT_INTERVAL_RESTRICTED")
    assert f.evidence["on_target_fraction"] < 0.3
    assert "not restricted to the kit" in f.message


def test_interval_restricted_calling_is_not_flagged(tmp_path):
    from admissible.checks.provenance import check_provenance
    from admissible.vcfio import SiteIndex, scan_vcf

    target = write_bed(tmp_path / "t.bed", [("chr1", 0, 1000)])
    vcf = tmp_path / "tight.vcf"
    rows = "".join(
        f"chr1\t{p}\t.\tA\tG\t50\tPASS\tAC=1\tGT:DP\t0/1:30\n" for p in range(100, 900, 100)
    )
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=249250621>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" + rows
    )
    res = check_provenance([scan_vcf(vcf, SiteIndex())], target_bed=str(target))
    assert not [x for x in res.findings if x.code == "NOT_INTERVAL_RESTRICTED"]
