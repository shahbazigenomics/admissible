"""Regression tests for three gaps found auditing a real multi-family
WES cohort mixing native VCFs and ANNOVAR multianno tables (see
``admissible_improvement_prompt.md``):

1. The duplicate-pair AND-gate (``agreement >= dup_agreement AND jaccard >=
   dup_jaccard``) silently dropped a real cross-format duplicate pair
   (agreement=1.000, jaccard=0.181, >5x difference in non-ref count) because
   jaccard collapses whenever one sample's call set is far smaller than the
   other's - exactly what a native VCF vs. an AC/AN-reconstructed multianno
   table looks like.
2. The full pairwise evidence table was written to the JSON report only; the
   text report (what most users actually read) never showed the numbers that
   would have let a human catch (1) even before the gate itself was fixed.
3. A cohort mixing native-VCF and multianno-reconstructed genotypes was
   calibrated as one pool with no indication that some samples' genotypes
   carry no per-site depth/GQ/allele-balance evidence at all.
"""

from __future__ import annotations

from admissible.checks.identity import IdentityConfig, check_identity
from admissible.model import Report, Severity
from admissible.report.text import render_text
from admissible.vcfio import load_cohort

HEADER_TWO = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=1,length=249250621>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
)


def _build_high_agreement_low_jaccard_vcf(tmp_path):
    """A is het at 100 sites; B is het at the first 20 of those (a subset)
    and hom-ref at the rest. Every shared site matches (agreement=1.0), the
    site sets only overlap 20/100 (jaccard=0.20), and non-ref counts differ
    5x (100 vs 20) - the same shape as the real IPK-5/IPM-BJ case (agreement
    1.000, jaccard 0.181, ~5.5x non-ref-count difference) that produced zero
    duplicate-related findings under the old AND-gate.
    """
    lines = [HEADER_TWO]
    for i in range(100):
        pos = 1000 + i
        b_gt = "0/1" if i < 20 else "0/0"
        lines.append(f"1\t{pos}\t.\tA\tG\t99\tPASS\t.\tGT\t0/1\t{b_gt}\n")
    vcf = tmp_path / "pair.vcf"
    vcf.write_text("".join(lines))
    return vcf


def test_high_agreement_low_jaccard_pair_is_flagged_not_silently_dropped(tmp_path):
    vcf = _build_high_agreement_low_jaccard_vcf(tmp_path)
    matrix = load_cohort([vcf])
    cfg = IdentityConfig(min_shared_sites=15)
    res = check_identity(matrix, ped=None, cfg=cfg)

    hits = [f for f in res.findings if f.code == "HIGH_AGREEMENT_LOW_JACCARD"]
    assert len(hits) == 1
    f = hits[0]
    assert set(f.subjects) == {"A", "B"}
    assert f.severity is Severity.WARN
    assert f.evidence["agreement"] == 1.0
    assert f.evidence["jaccard"] == 0.20
    assert f.evidence["n_nonref_a"] == 100
    assert f.evidence["n_nonref_b"] == 20

    # Not auto-classified as a duplicate - a human decides - but this must
    # never again be the ONLY thing the report says about this pair.
    assert res.metrics["n_duplicate_pairs"] == 0


def test_similar_sized_call_sets_with_low_jaccard_are_flagged_too(tmp_path):
    """The same person called by two pipelines: each file has about the same
    number of non-ref sites but they only partly overlap. Agreement on shared
    sites is perfect, Jaccard is low, and the non-ref counts are within the
    2x completeness ratio. This used to produce no finding at all.
    """
    lines = [HEADER_TWO]
    # 60 sites both call, 20 only A calls, 20 only B calls -> jaccard 60/100.
    # Make it clearly low: 50 shared, 40 A-only, 40 B-only -> 50/130 = 0.38.
    for i in range(130):
        pos = 1000 + i
        if i < 50:
            a_gt = b_gt = "0/1"
        elif i < 90:
            a_gt, b_gt = "0/1", "0/0"
        else:
            a_gt, b_gt = "0/0", "0/1"
        lines.append(f"1\t{pos}\t.\tA\tG\t99\tPASS\t.\tGT\t{a_gt}\t{b_gt}\n")
    vcf = tmp_path / "pair.vcf"
    vcf.write_text("".join(lines))
    matrix = load_cohort([vcf])
    cfg = IdentityConfig(min_shared_sites=15)
    res = check_identity(matrix, ped=None, cfg=cfg)

    hits = [f for f in res.findings if f.code == "HIGH_AGREEMENT_LOW_JACCARD"]
    assert len(hits) == 1
    assert hits[0].severity is Severity.WARN
    assert hits[0].evidence["jaccard"] < 0.6
    assert hits[0].evidence["n_nonref_a"] == hits[0].evidence["n_nonref_b"] == 90
    assert "called separately" in hits[0].message
    assert res.metrics["n_duplicate_pairs"] == 0


def test_low_agreement_pair_is_still_not_flagged(tmp_path):
    """Low Jaccard alone (agreement below the bar) must stay silent."""
    lines = [HEADER_TWO]
    for i in range(60):
        pos = 1000 + i
        # shared on all 60 sites but genotypes differ on half -> agreement 0.5
        a_gt, b_gt = ("0/1", "0/1") if i % 2 else ("0/1", "1/1")
        lines.append(f"1\t{pos}\t.\tA\tG\t99\tPASS\t.\tGT\t{a_gt}\t{b_gt}\n")
    vcf = tmp_path / "pair.vcf"
    vcf.write_text("".join(lines))
    matrix = load_cohort([vcf])
    res = check_identity(matrix, ped=None, cfg=IdentityConfig(min_shared_sites=15))
    assert not [f for f in res.findings if f.code == "HIGH_AGREEMENT_LOW_JACCARD"]


def test_pairwise_evidence_is_printed_in_verbose_text_report(tmp_path):
    vcf = _build_high_agreement_low_jaccard_vcf(tmp_path)
    matrix = load_cohort([vcf])
    cfg = IdentityConfig(min_shared_sites=15)
    res = check_identity(matrix, ped=None, cfg=cfg)
    report = Report(family_id="TEST", checks=[res])

    quiet = render_text(report, verbose=False)
    assert "PAIRWISE EVIDENCE" not in quiet

    text = render_text(report, verbose=True)
    assert "PAIRWISE EVIDENCE" in text
    assert "A / B" in text
    assert "agreement=1.0000" in text
    assert "jaccard=0.2000" in text


def _build_mixed_source_cohort(tmp_path):
    native = tmp_path / "native.vcf"
    native.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1,length=249250621>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNATIVE1\n"
        "1\t1000\t.\tA\tG\t99\tPASS\t.\tGT\t0/1\n"
        "1\t2000\t.\tA\tG\t99\tPASS\t.\tGT\t0/0\n"
    )
    multianno_dir = tmp_path / "FAM1" / "S3"
    multianno_dir.mkdir(parents=True)
    multianno = multianno_dir / "sample.hg19_multianno.txt"
    multianno.write_text(
        "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\t"
        "Otherinfo1\tOtherinfo2\tOtherinfo3\tOtherinfo4\tOtherinfo5\t"
        "Otherinfo6\tOtherinfo7\tOtherinfo8\n"
        "chr1\t1000\t1000\tA\tG\texonic\tGENEX\t"
        "chr1\t1000\t.\tA\tG\t50\tPASS\tAC=1;AN=2\n"
    )
    return native, multianno


def test_mixed_genotype_source_cohort_is_flagged(tmp_path):
    native, multianno = _build_mixed_source_cohort(tmp_path)
    matrix = load_cohort([native, multianno])
    assert set(matrix.samples) == {"NATIVE1", "FAM1-S3"}

    res = check_identity(matrix, ped=None)
    hits = [f for f in res.findings if f.code == "MIXED_GENOTYPE_SOURCE"]
    assert len(hits) == 1
    f = hits[0]
    assert f.severity is Severity.WARN
    assert f.subjects == ["FAM1-S3"]
    assert f.evidence["native_format_samples"] == ["NATIVE1"]
    assert f.evidence["reconstructed_ac_an_samples"] == ["FAM1-S3"]
    assert res.metrics["genotype_source"] == {
        "NATIVE1": "native/format",
        "FAM1-S3": "reconstructed-AC/AN",
    }


def test_single_source_cohort_is_not_flagged(tmp_path):
    """Two native VCFs, no multianno involved - the finding must not fire
    just because there are two samples; it fires only on an actual mix.
    """
    a = tmp_path / "a.vcf"
    a.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n"
        "1\t1000\t.\tA\tG\t99\tPASS\t.\tGT\t0/1\n"
    )
    b = tmp_path / "b.vcf"
    b.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tB\n"
        "1\t1000\t.\tA\tG\t99\tPASS\t.\tGT\t0/1\n"
    )
    matrix = load_cohort([a, b])
    res = check_identity(matrix, ped=None)
    assert not [f for f in res.findings if f.code == "MIXED_GENOTYPE_SOURCE"]
