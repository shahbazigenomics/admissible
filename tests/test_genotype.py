"""Check 2 - genotype confidence.

The fixture is written inline because each record exists to exercise exactly one
rule, and a record whose purpose has to be inferred from a generator is a worse
test than one you can read.
"""

from __future__ import annotations

import pytest

from admissible.checks.genotype import GenotypeConfig, check_genotype, evaluate_genotype
from admissible.model import Status

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr15,length=102531392>
##contig=<ID=chr1,length=249250621>
##reference=GRCh37
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tpatient
"""

# (id, pos, INFO, sample cell) - one per rule under test.
ROWS = [
    # The real-world shape: hom-alt from two reads, GATK's own MLE says het.
    ("false_hom_low_depth", 1000000, "AC=2;MLEAC=1;MLEAF=0.5;DP=2", "1/1:0,2:2:6:90,6,0"),
    # Adequate depth, reference reads still present: a different failure entirely.
    ("imbalanced_hom", 1100000, "AC=2;MLEAC=2;DP=50", "1/1:8,42:50:99:900,80,0"),
    # A clean homozygote: nothing should fire.
    ("clean_hom", 1200000, "AC=2;MLEAC=2;DP=40", "1/1:0,40:40:99:1200,120,0"),
    # MLEAC disagrees but every independent rule is satisfied.
    ("mleac_only", 1300000, "AC=2;MLEAC=1;DP=40", "1/1:0,40:40:99:1200,120,0"),
    # Clean heterozygote.
    ("clean_het", 1400000, "AC=1;MLEAC=1;DP=40", "0/1:20,20:40:99:600,0,600"),
    # Allele balance 0.20 at DP=10: a fixed 0.30-0.70 window flags this; the
    # binomial does not, because 2 of 10 is unremarkable sampling.
    ("skew_shallow", 1500000, "AC=1;MLEAC=1;DP=10", "0/1:8,2:10:45:120,0,300"),
    # Allele balance 0.20 at DP=200: same ratio, real signal.
    ("skew_deep", 1600000, "AC=1;MLEAC=1;DP=200", "0/1:160,40:200:99:900,0,4000"),
]


@pytest.fixture
def vcf(tmp_path):
    body = HEADER
    for _name, pos, info, cell in ROWS:
        body += f"chr15\t{pos}\t.\tA\tG\t500\tPASS\t{info}\tGT:AD:DP:GQ:PL\t{cell}\n"
    p = tmp_path / "patient.vcf"
    p.write_text(body)
    return p


def flags_for(name: str) -> set[str]:
    _n, _pos, info, cell = next(r for r in ROWS if r[0] == name)
    keys = "GT:AD:DP:GQ:PL".split(":")
    cell_d = dict(zip(keys, cell.split(":"), strict=False))
    info_d = dict(kv.split("=", 1) for kv in info.split(";"))
    got, _ev = evaluate_genotype(cell_d["GT"], cell_d, info_d, GenotypeConfig(), True)
    return set(got)


def test_the_real_shape_is_caught_by_depth_and_likelihood_together():
    f = flags_for("false_hom_low_depth")
    assert "FALSE_HOM_SUSPECT" in f
    assert "LOW_DEPTH_HOM" in f
    assert "HOM_CONTRADICTED_BY_LIKELIHOOD" in f


def test_imbalanced_hom_is_not_labelled_low_depth():
    f = flags_for("imbalanced_hom")
    assert "ALLELE_IMBALANCE_HOM" in f
    assert "LOW_DEPTH_HOM" not in f, "adequate depth must not be reported as shallow"
    assert "FALSE_HOM_SUSPECT" in f


def test_clean_homozygote_raises_nothing():
    assert flags_for("clean_hom") == set()


def test_mleac_does_not_condemn_a_genotype_that_has_gq_and_pl():
    """Where GQ/PL exist, MLEAC is redundant with them and must not act.

    For a biallelic hom-alt call PL[het] is GQ by construction, so acting on
    MLEAC as well would double-count the same evidence and bury the real signal
    under a flag that fires on a large fraction of all homozygous calls.
    """
    f = flags_for("mleac_only")
    assert "MLEAC_DISCORDANT" in f
    assert "MLEAC_CONTRADICTS_GT" not in f
    assert "FALSE_HOM_SUSPECT" not in f


def test_mleac_does_condemn_when_gq_and_pl_are_gone():
    """With the FORMAT block stripped, MLEAC is the only witness left."""
    cell = {"GT": "1/1"}
    info = {"AC": "2", "MLEAC": "1", "AN": "2", "DP": "40"}
    got, _ = evaluate_genotype("1/1", cell, info, GenotypeConfig(), True)
    assert "MLEAC_CONTRADICTS_GT" in got
    assert "FALSE_HOM_SUSPECT" in got


def test_allele_balance_is_depth_scaled_not_a_fixed_window():
    shallow, deep = flags_for("skew_shallow"), flags_for("skew_deep")
    # Identical allele balance (0.20), opposite verdicts.
    assert "AB_SKEW" not in shallow
    assert "AB_SKEW" in deep


def test_clean_het_raises_nothing():
    assert flags_for("clean_het") == set()


def test_end_to_end_counts_and_overlap(vcf):
    res = check_genotype([vcf])
    assert res.status is Status.WARN
    ev = next(f for f in res.findings if f.code == "FALSE_HOM_SUSPECT").evidence
    assert ev["n_false_hom"] == 2
    assert ev["n_hom_alt"] == 4
    assert ev["breakdown"]["LOW_DEPTH_HOM"] == 1
    assert ev["breakdown"]["ALLELE_IMBALANCE_HOM"] == 1
    loci = {e["locus"] for e in ev["examples"]}
    assert "15:1000000" in loci

    ov = res.metrics["rule_overlap"]
    assert ov["mleac_flagged"] == 2
    assert ov["mleac_also_caught_by_rules"] == 1  # the other is MLEAC-only


def test_files_without_format_report_unknown(tmp_path):
    p = tmp_path / "sites_only.vcf"
    p.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tA\tG\t50\tPASS\tAC=1\n"
    )
    res = check_genotype([p])
    assert res.status is Status.UNKNOWN
    assert "cannot be assessed" in res.summary


def test_never_raises_on_garbage(tmp_path):
    p = tmp_path / "junk.vcf"
    p.write_text("not a vcf\n\x00\x01\n#CHROM\tPOS\n1\t\n")
    assert check_genotype([p]).status in tuple(Status)
