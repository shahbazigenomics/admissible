"""The contract: no input makes this tool raise.

A tool whose purpose is to describe broken files must not die on one.  These
cases are all real things that arrive in a shared drive: truncated downloads,
files gzipped twice, sites-only VCFs, VCFs with no header, headers with no
records, and files that are not VCFs at all.
"""

from __future__ import annotations

import gzip
import random

import pytest

from admissible.checks.identity import check_identity
from admissible.checks.provenance import check_provenance
from admissible.model import Status
from admissible.ped import read_ped
from admissible.vcfio import SiteIndex, encode_gt, load_cohort, scan_vcf

GOOD = """##fileformat=VCFv4.2
##contig=<ID=1,length=249250621>
##contig=<ID=X,length=155270560>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t1000\t.\tA\tG\t50\tPASS\tAC=1\tGT:DP\t0/1:30
1\t2000\t.\tC\tT\t50\tPASS\tAC=2\tGT:DP\t1/1:30
"""

BROKEN = {
    "empty.vcf": "",
    "header_only.vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    "no_header.vcf": "1\t1000\t.\tA\tG\t50\tPASS\tAC=1\tGT\t0/1\n",
    "truncated.vcf": GOOD[: len(GOOD) // 2],
    "short_columns.vcf": GOOD + "1\t3000\t.\tA\n",
    "bad_pos.vcf": GOOD + "1\tNOT_A_NUMBER\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n",
    "no_format.vcf": GOOD.replace("\tGT:DP\t0/1:30", "").replace("\tGT:DP\t1/1:30", ""),
    "weird_gt.vcf": GOOD + "1\t4000\t.\tA\tG\t50\tPASS\t.\tGT\t2|3|4\n",
    "nul_bytes.vcf": GOOD + "\x00\x01\x02\x03\n",
    "not_a_vcf.vcf": "this is a PDF, honestly\n" * 50,
    "sites_only.vcf": "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t1\t.\tA\tG\t1\tPASS\t.\n",
}


@pytest.fixture
def broken_dir(tmp_path):
    for name, body in BROKEN.items():
        (tmp_path / name).write_text(body)
    (tmp_path / "truncated.vcf.gz").write_bytes(gzip.compress(GOOD.encode())[:40])
    (tmp_path / "binary.vcf.gz").write_bytes(bytes(range(256)) * 20)
    (tmp_path / "sample.ped").write_text("FAM\tS1\t0\t0\t1\t2\n")
    return tmp_path


@pytest.mark.parametrize("name", sorted(BROKEN) + ["truncated.vcf.gz", "binary.vcf.gz"])
def test_scan_never_raises(broken_dir, name):
    scan = scan_vcf(broken_dir / name, SiteIndex())
    assert scan.path.endswith(name)
    assert isinstance(scan.problems, list)


def test_missing_file_is_a_problem_not_an_exception(tmp_path):
    scan = scan_vcf(tmp_path / "nope.vcf", SiteIndex())
    assert scan.problems and "cannot open" in scan.problems[0]


def test_checks_never_raise_on_broken_inputs(broken_dir):
    matrix = load_cohort(sorted(broken_dir.glob("*.vcf")) + sorted(broken_dir.glob("*.gz")))
    ped = read_ped(broken_dir / "sample.ped")
    for res in (check_identity(matrix, ped), check_provenance(list(matrix.scans.values()))):
        assert res.status in tuple(Status)
        assert not any(f.code == "CHECK_CRASHED" for f in res.findings), res.summary


def test_ped_parser_never_raises(tmp_path):
    p = tmp_path / "bad.ped"
    p.write_text("FAM\tA\n\nFAM\tB\t0\t0\t9\t9\nFAM\tB\t0\t0\t1\t1\nFAM\tC\tGHOST\t0\t1\t2\n")
    ped = read_ped(p)
    assert ped.warnings
    assert "GHOST" in " ".join(ped.warnings)


def test_missing_ped_file_is_a_warning(tmp_path):
    ped = read_ped(tmp_path / "absent.ped")
    assert ped.warnings and not ped.individuals


@pytest.mark.parametrize(
    "gt,expected",
    [
        ("0/0", 0), ("0|0", 0), ("0/1", 1), ("1|0", 1), ("1/1", 2), ("2/2", 2),
        ("1/2", 1), ("./.", 3), (".", 3), ("", 3), ("0", 0), ("1", 2),
        ("0/.", 3), ("A/B", 3), ("1/1/1", 2), ("0/1/2", 1),
    ],
)
def test_gt_encoding(gt, expected):
    assert encode_gt(gt) == expected


def test_random_byte_soup_never_raises(tmp_path):
    rng = random.Random(1)
    for i in range(25):
        blob = bytes(rng.randrange(256) for _ in range(400))
        path = tmp_path / f"soup{i}.vcf"
        path.write_bytes(blob)
        scan_vcf(path, SiteIndex())


def test_a_gvcf_and_a_vcf_of_the_same_variant_are_the_same_site(tmp_path):
    """``<NON_REF>`` must not become part of a variant's identity.

    GATK writes every gVCF ALT as ``G,<NON_REF>``; a plain VCF writes ``G``.
    Carrying the symbolic allele into the site key gave one variant two
    identities, so a gVCF and a VCF *of the same person* shared no sites at all
    and the pair came out unrelated. Measured on a real GATK gVCF and its own
    call set before the fix: 0 shared sites out of 65 that should have matched.

    Failing towards "unrelated" is the one direction an identity check must
    never fail in - it turns a duplicate into two strangers.
    """
    head = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    )
    gvcf = tmp_path / "a.g.vcf"
    gvcf.write_text(
        head
        + "chr1\t100\t.\tG\t<NON_REF>\t.\t.\tEND=150\tGT:DP\t0/0:30\n"
        + "chr1\t200\t.\tC\tT,<NON_REF>\t50\t.\tDP=30\tGT:DP\t0/1:30\n"
    )
    vcf = tmp_path / "b.vcf"
    vcf.write_text(head + "chr1\t200\t.\tC\tT\t50\tPASS\tDP=30\tGT:DP\t0/1:30\n")

    index = SiteIndex()
    g = scan_vcf(gvcf, index)
    v = scan_vcf(vcf, index)
    # The reference block is counted, never indexed as a variant site.
    assert g.stats.n_reference_blocks == 1
    assert g.stats.n_records == 1
    assert len(index) == 1, index.keys          # one variant, one identity
    assert index.keys[0] == ("1", 200, "C", "T")
    assert set(g.genotypes["S1"]) == set(v.genotypes["S1"])
