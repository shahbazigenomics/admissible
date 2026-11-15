"""Sex is called female above the cohort boundary, male only when clearly male-like,
and left undetermined (a WARN, never blocking) in between."""

from __future__ import annotations

import random

from admissible.checks.identity import IdentityConfig, check_identity
from admissible.model import Severity
from admissible.ped import read_ped
from admissible.vcfio import load_cohort

SAMPLES = ["F1", "F2", "F3", "F4", "M1", "M2", "M3", "M4", "MIDF", "MIDM", "SWAP"]
# declared sex of each sample; MIDF/MIDM have chrX het between ceiling and boundary,
# SWAP is declared female with clearly male-like genotypes.
FEMALE = {"F1", "F2", "F3", "F4", "MIDF", "SWAP"}
RATE = {"MIDF": 0.38, "MIDM": 0.38, "SWAP": 0.03}


def _build(tmp_path, seed=5):
    rng = random.Random(seed)
    head = (
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n##contig=<ID=X,length=155270560>\n"
        "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t".replace("\\t", "\t")
        + "\t".join(SAMPLES) + "\n"
    )
    rows = [head]
    for i in range(1500):
        g = ["0/1" if rng.random() < 0.7 else "1/1" for _ in SAMPLES]
        rows.append(f"1\t{10_000 + i * 50}\t.\tA\tG\t99\tPASS\t.\tGT\t" + "\t".join(g) + "\n")
    for i in range(800):
        g = []
        for s in SAMPLES:
            p = RATE.get(s, 0.6 if s in FEMALE else 0.03)
            g.append("0/1" if rng.random() < p else "1/1")
        rows.append(f"X\t{60_000_000 + i * 1000}\t.\tA\tG\t99\tPASS\t.\tGT\t" + "\t".join(g) + "\n")
    vcf = tmp_path / "c.vcf"
    vcf.write_text("".join(rows))
    ped = tmp_path / "c.ped"
    ped.write_text("".join(f"FAM\t{s}\t0\t0\t{'2' if s in FEMALE else '1'}\t1\n" for s in SAMPLES))
    return check_identity(load_cohort([vcf]), ped=read_ped(ped), cfg=IdentityConfig())


def _by_code(res, code):
    return {f.subjects[0]: f for f in res.findings if f.code == code}


def test_undetermined_zone_is_a_warning_not_a_mismatch(tmp_path):
    res = _build(tmp_path)
    mism = _by_code(res, "SEX_MISMATCH")
    unc = _by_code(res, "SEX_UNCERTAIN")
    assert "MIDF" not in mism and "MIDM" not in mism
    assert set(unc) == {"MIDF", "MIDM"}
    assert all(f.severity is Severity.WARN for f in unc.values())


def test_clearly_male_genotypes_declared_female_still_block(tmp_path):
    res = _build(tmp_path)
    mism = _by_code(res, "SEX_MISMATCH")
    assert set(mism) == {"SWAP"}
    assert mism["SWAP"].severity is Severity.BLOCKING


def test_clear_calls_give_no_sex_finding(tmp_path):
    res = _build(tmp_path)
    flagged = set(_by_code(res, "SEX_MISMATCH")) | set(_by_code(res, "SEX_UNCERTAIN"))
    assert flagged.isdisjoint({"F1", "F2", "F3", "F4", "M1", "M2", "M3", "M4"})


def test_ceiling_is_configurable(tmp_path):
    res = _build(tmp_path)  # default 0.30
    assert _by_code(res, "SEX_UNCERTAIN")
    # a ceiling above the boundary reduces to the old two-way rule
    cfg = IdentityConfig(male_het_ceiling=1.0)
    from admissible.checks import identity as ident
    ev, boundary, _ = ident.sex_evidence(
        load_cohort([tmp_path / "c.vcf"]),
        ident.build_sets(load_cohort([tmp_path / "c.vcf"])),
        cfg,
        None,
    )
    assert ev["MIDF"].inferred.value == "male"
