"""The duplicate finding offers 'monozygotic twins' as a compatible cause only
when the pedigree makes that plausible.

Identical twins reach the same agreement (~0.999) and site-set overlap as a
re-run of one sample, so ``DUPLICATE_SAME_INDIVIDUAL`` fires and blocks the
audit.  A message that never mentions twins sends the user hunting for a swap
that does not exist.  But a hint shown on every duplicate would be read as
"expected" and wave a real swap through, so it is shown only when:

* the pedigree puts the pair in one family (full sibs, or same family with no
  stated link), and
* the pedigree records the same, known sex for both.

It is never shown for a cross-family pair - the failure this check exists for.
The finding stays BLOCKING in every case; this changes the explanation only.
"""

from __future__ import annotations

from admissible.checks.identity import IdentityConfig, check_identity
from admissible.model import Severity
from admissible.ped import read_ped
from admissible.vcfio import load_cohort

N_SITES = 600


def _identical_pair_vcf(tmp_path):
    """Two samples with identical genotypes at 600 sites: agreement 1.0, jaccard 1.0."""
    lines = [
        "##fileformat=VCFv4.2\n",
        "##contig=<ID=1,length=249250621>\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT1\tT2\n",
    ]
    for i in range(N_SITES):
        gt = "0/1" if i % 3 else "1/1"
        lines.append(f"1\t{1000 + i}\t.\tA\tG\t99\tPASS\t.\tGT\t{gt}\t{gt}\n")
    vcf = tmp_path / "twins.vcf"
    vcf.write_text("".join(lines))
    return vcf


def _duplicate_finding(tmp_path, ped_text):
    vcf = _identical_pair_vcf(tmp_path)
    ped_path = tmp_path / "t.ped"
    ped_path.write_text(ped_text)
    res = check_identity(load_cohort([vcf]), ped=read_ped(ped_path), cfg=IdentityConfig())
    hits = [f for f in res.findings if f.code.startswith("DUPLICATE_")]
    assert len(hits) == 1
    return hits[0]


SAME_SEX_SIBS = (
    "F1 DAD 0 0 1 1\nF1 MOM 0 0 2 1\nF1 T1 DAD MOM 1 2\nF1 T2 DAD MOM 1 2\n"
)


def test_same_sex_full_sibs_get_the_twin_hint_and_stay_blocking(tmp_path):
    f = _duplicate_finding(tmp_path, SAME_SEX_SIBS)
    assert f.severity is Severity.BLOCKING
    assert f.evidence["possible_mz_twins"] is True
    assert "monozygotic twins" in f.message
    # the hint must not read as "this is fine"
    assert "only if" in f.message and "swap" in f.message


def test_opposite_sex_sibs_never_get_the_twin_hint(tmp_path):
    ped = "F1 DAD 0 0 1 1\nF1 MOM 0 0 2 1\nF1 T1 DAD MOM 1 2\nF1 T2 DAD MOM 2 2\n"
    f = _duplicate_finding(tmp_path, ped)
    assert f.severity is Severity.BLOCKING
    assert f.evidence["possible_mz_twins"] is False
    assert "monozygotic" not in f.message


def test_unknown_sex_never_gets_the_twin_hint(tmp_path):
    ped = "F1 DAD 0 0 1 1\nF1 MOM 0 0 2 1\nF1 T1 DAD MOM 0 2\nF1 T2 DAD MOM 0 2\n"
    f = _duplicate_finding(tmp_path, ped)
    assert f.evidence["possible_mz_twins"] is False
    assert "monozygotic" not in f.message


def test_cross_family_duplicate_never_gets_the_twin_hint(tmp_path):
    """The original failure: two files filed under different families are one person."""
    ped = "F1 T1 0 0 1 2\nF2 T2 0 0 1 2\n"
    f = _duplicate_finding(tmp_path, ped)
    assert f.severity is Severity.BLOCKING
    assert f.evidence["cross_family"] is True
    assert f.evidence["possible_mz_twins"] is False
    assert "monozygotic" not in f.message


def test_same_family_no_stated_link_same_sex_gets_the_hint(tmp_path):
    ped = "F1 T1 0 0 2 2\nF1 T2 0 0 2 2\n"
    f = _duplicate_finding(tmp_path, ped)
    assert f.evidence["possible_mz_twins"] is True
    assert "monozygotic twins" in f.message


def test_half_sibs_and_parent_offspring_never_get_the_hint(tmp_path):
    half = "F1 DAD 0 0 1 1\nF1 MOM 0 0 2 1\nF1 T1 DAD MOM 1 2\nF1 T2 DAD 0 1 2\n"
    f = _duplicate_finding(tmp_path, half)
    assert f.evidence["possible_mz_twins"] is False
    po = "F1 T1 0 0 1 1\nF1 T2 T1 0 1 2\n"
    f = _duplicate_finding(tmp_path, po)
    assert f.evidence["possible_mz_twins"] is False
