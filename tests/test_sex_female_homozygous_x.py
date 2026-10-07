"""A female with a homozygous X looks male by chrX het; the report must say what
separates her from a sample swap, without changing the call.

The cohort is simulated, so this checks the plumbing and the direction of the
signal (runs of homozygosity lower autosomal het, a swap does not). It is not a
calibration: no consanguineous data was available.
"""

from __future__ import annotations

import random

from admissible.checks.identity import IdentityConfig, check_identity
from admissible.model import Severity
from admissible.ped import read_ped
from admissible.vcfio import load_cohort

SAMPLES = ["F1", "F2", "F3", "M1", "M2", "M3", "ROHF", "SWAP"]
HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=1,length=249250621>\n##contig=<ID=X,length=155270560>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(SAMPLES) + "\n"
)


def _vcf(tmp_path, seed=3):
    rng = random.Random(seed)
    rows = [HEADER]
    # Autosome: het/(het+homalt) about 0.8 normally; the ROH female has the first
    # 45% of the chromosome in a run of homozygosity (no het calls there).
    n_auto = 3000
    for i in range(n_auto):
        gts = []
        for s in SAMPLES:
            roh = s == "ROHF" and i < 0.45 * n_auto
            if rng.random() < 0.7:  # site is non-ref in this sample
                het = (not roh) and rng.random() < 0.8
                gts.append("0/1" if het else "1/1")
            else:
                gts.append("0/0")
        if all(g == "0/0" for g in gts):
            gts[0] = "0/1"
        rows.append(f"1\t{10_000 + i * 50}\t.\tA\tG\t99\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")
    # chrX (outside PAR): females het ~0.6, males and the homozygous-X female ~0.1.
    for i in range(700):
        gts = []
        for s in SAMPLES:
            p_het = 0.6 if s.startswith("F") else 0.1
            gts.append("0/1" if rng.random() < p_het else "1/1")
        pos = 60_000_000 + i * 1000
        rows.append(f"X\t{pos}\t.\tA\tG\t99\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")
    p = tmp_path / "c.vcf"
    p.write_text("".join(rows))
    return p


def _ped(tmp_path):
    p = tmp_path / "c.ped"
    # ROHF and SWAP are both declared female; SWAP's genotypes are male-like on X
    # and normal on the autosomes.
    lines = []
    for s in SAMPLES:
        sex = "2" if s.startswith("F") or s in ("ROHF", "SWAP") else "1"
        lines.append(f"FAM\t{s}\t0\t0\t{sex}\t1")
    p.write_text("\n".join(lines) + "\n")
    return read_ped(p)


def _result(tmp_path):
    # SWAP needs male-like X: give it the male het rate by rewriting its column.
    vcf = _vcf(tmp_path)
    lines = vcf.read_text().split("\n")
    out, rng = [], random.Random(11)
    swap_col = 9 + SAMPLES.index("SWAP")
    for line in lines:
        if line.startswith("X\t"):
            f = line.split("\t")
            f[swap_col] = "0/1" if rng.random() < 0.1 else "1/1"
            line = "\t".join(f)
        out.append(line)
    vcf.write_text("\n".join(out))
    return check_identity(load_cohort([vcf]), ped=_ped(tmp_path), cfg=IdentityConfig())


def test_both_are_still_called_mismatches_and_still_blocking(tmp_path):
    res = _result(tmp_path)
    hits = {f.subjects[0]: f for f in res.findings if f.code == "SEX_MISMATCH"}
    assert set(hits) == {"ROHF", "SWAP"}
    assert all(f.severity is Severity.BLOCKING for f in hits.values())


def test_the_report_gives_the_evidence_that_separates_them(tmp_path):
    res = _result(tmp_path)
    hits = {f.subjects[0]: f for f in res.findings if f.code == "SEX_MISMATCH"}
    roh, swap = hits["ROHF"].evidence, hits["SWAP"].evidence
    med = roh["cohort_median_autosomal_het_frac"]
    assert roh["autosomal_het_frac"] < med - 0.15  # homozygosity lowers it
    assert abs(swap["autosomal_het_frac"] - med) < 0.05  # a swap does not
    assert "autosomal het fraction" in hits["ROHF"].message
    assert "homozygous X" in hits["ROHF"].message
    assert roh["female_homozygous_x_possible"] is True


def test_context_is_only_added_for_a_declared_female_called_male(tmp_path):
    res = _result(tmp_path)
    for f in res.findings:
        if f.code != "SEX_MISMATCH":
            continue
        assert f.evidence["declared_sex"] == "female"
        assert f.evidence["inferred_sex"] == "male"
