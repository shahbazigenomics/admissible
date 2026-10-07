"""0.1.1: skipped-pair note and honest wording when no pedigree is given."""

from __future__ import annotations

from admissible.checks.identity import IdentityConfig, check_identity
from admissible.model import Severity
from admissible.vcfio import load_cohort

HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=1,length=249250621>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
)


def _vcf(tmp_path, n, a_gt="0/1", b_gt="0/1"):
    rows = [HEADER]
    for i in range(n):
        rows.append(f"1\t{1000 + i}\t.\tA\tG\t99\tPASS\t.\tGT\t{a_gt}\t{b_gt}\n")
    p = tmp_path / "c.vcf"
    p.write_text("".join(rows))
    return p


def test_sparse_pairs_are_reported_not_silently_skipped(tmp_path):
    matrix = load_cohort([_vcf(tmp_path, 40)])
    res = check_identity(matrix, ped=None, cfg=IdentityConfig(min_shared_sites=500))
    hits = [f for f in res.findings if f.code == "DUPLICATE_CHECK_SKIPPED"]
    assert len(hits) == 1
    assert hits[0].severity is Severity.INFO
    assert hits[0].evidence["n_pairs_skipped"] == 1
    assert hits[0].evidence["n_pairs_total"] == 1


def test_no_skipped_note_when_all_pairs_are_checked(tmp_path):
    matrix = load_cohort([_vcf(tmp_path, 40)])
    res = check_identity(matrix, ped=None, cfg=IdentityConfig(min_shared_sites=15))
    assert not [f for f in res.findings if f.code == "DUPLICATE_CHECK_SKIPPED"]


def test_duplicate_without_pedigree_does_not_claim_across_families(tmp_path):
    matrix = load_cohort([_vcf(tmp_path, 40)])
    res = check_identity(matrix, ped=None, cfg=IdentityConfig(min_shared_sites=15))
    dup = [f for f in res.findings if f.severity is Severity.BLOCKING]
    assert len(dup) == 1
    assert "ACROSS FAMILIES" not in dup[0].message
    assert "declares them" not in dup[0].message
    assert "does not list both" in dup[0].message
    assert dup[0].evidence["cross_family"] is False
    assert dup[0].evidence["declared_relationship"] is None
