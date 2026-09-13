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


# --- allele depth is not always spelled AD ---------------------------------


def _one_sample_vcf(tmp_path, fmt, cell, info="AC=2"):
    p = tmp_path / "one.vcf"
    p.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        f"chr1\t1000\t.\tA\tG\t900\tPASS\t{info}\t{fmt}\t{cell}\n"
    )
    return p


def test_freebayes_ro_ao_is_read_as_allele_depth(tmp_path):
    """freebayes writes RO/AO, never AD.

    Reading only AD left the whole allele-balance arm of this check inert on
    freebayes output while the report still looked complete - no AB_SKEW on
    heterozygotes, no ALLELE_IMBALANCE_HOM on homozygotes. Found by running the
    check on the public CEPH 1463 call set.
    """
    # 40 reads, 4 of them alt, called hom-alt: contradicted by its own evidence.
    vcf = _one_sample_vcf(tmp_path, "GT:DP:RO:AO", "1/1:40:36:4")
    res = check_genotype([vcf])
    flags = res.metrics["per_sample"]["S1"]["flags"]
    assert flags.get("ALLELE_IMBALANCE_HOM") == 1
    assert flags.get("FALSE_HOM_SUSPECT") == 1


def test_dp4_is_read_as_allele_depth(tmp_path):
    """Older samtools/bcftools pipelines write DP4 and nothing else."""
    vcf = _one_sample_vcf(tmp_path, "GT:DP:DP4", "0/1:40:18,18,2,2")
    res = check_genotype([vcf])
    assert res.metrics["per_sample"]["S1"]["flags"].get("AB_SKEW") == 1


def test_varscan_single_valued_ad_is_not_mistaken_for_gatk_ad(tmp_path):
    """VarScan's AD is the ALT count alone; reading it as GATK's would invert it."""
    vcf = _one_sample_vcf(tmp_path, "GT:DP:RD:AD", "1/1:40:2:38")
    res = check_genotype([vcf])
    # No allele-balance flag either way: the field was correctly not trusted.
    assert "ALLELE_IMBALANCE_HOM" not in res.metrics["per_sample"]["S1"]["flags"]


# --- three bugs found by asking "what does this do when the field is absent?" ---


def test_info_dp_is_not_used_as_a_per_sample_depth_in_a_multisample_vcf(tmp_path):
    """INFO/DP is the cohort total; using it per sample multiplies depth by N.

    A ten-sample VCF with INFO DP=60 and no FORMAT/DP gave every sample an
    apparent depth of 60 - six times the truth - so ten 6x homozygotes were
    tallied as "at adequate depth" and the check returned PASS with no findings.
    That is the exact failure this check exists to catch, committed by the check.
    """
    samples = [f"S{i}" for i in range(1, 11)]
    p = tmp_path / "multi.vcf"
    p.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n"
        "chr1\t1000\t.\tA\tG\t900\tPASS\tDP=60\tGT:GQ\t"
        + "\t".join(["1/1:99"] * len(samples)) + "\n"
    )
    res = check_genotype([p])
    assert res.metrics["hom_at_adequate_depth"] == {
        "n": 0, "flagged": 0, "depth_unknown": 10
    }
    assert res.metrics["per_sample"]["S1"]["flags"].get("HOM_NOT_ASSESSABLE") == 1
    assert any(f.code == "HOM_NOT_ASSESSABLE" for f in res.findings)


def test_info_dp_is_still_used_when_there_is_only_one_sample(tmp_path):
    """Single-sample: INFO/DP and FORMAT/DP are the same number, so use it."""
    p = _one_sample_vcf(tmp_path, "GT", "1/1", info="DP=4")
    res = check_genotype([p])
    assert res.metrics["per_sample"]["S1"]["flags"].get("LOW_DEPTH_HOM") == 1


def test_a_multiallelic_hom_alt_is_judged_on_the_allele_it_called(tmp_path):
    """Summing every ALT calls a 5-of-30 homozygote perfectly balanced."""
    from admissible.checks.genotype import GenotypeConfig, evaluate_genotype

    cfg = GenotypeConfig()
    flags, ev = evaluate_genotype("1/1", {"DP": "30", "AD": "0,5,25"}, {}, cfg, True)
    assert ev["allele_balance"] == pytest.approx(5 / 30)
    assert "ALLELE_IMBALANCE_HOM" in flags and "FALSE_HOM_SUSPECT" in flags

    # The same site with the reads actually behind the called allele is clean.
    flags, ev = evaluate_genotype("1/1", {"DP": "30", "AD": "0,30,0"}, {}, cfg, True)
    assert ev["allele_balance"] == pytest.approx(1.0)
    assert flags == []


def test_a_multiallelic_het_is_tested_on_its_two_called_alleles(tmp_path):
    """Reads supporting a third allele say nothing about a 1/2 call's balance."""
    from admissible.checks.genotype import GenotypeConfig, evaluate_genotype

    cfg = GenotypeConfig()
    flags, _ = evaluate_genotype("1/2", {"DP": "30", "AD": "0,15,15"}, {}, cfg, True)
    assert "AB_SKEW" not in flags
    flags, _ = evaluate_genotype("1/2", {"DP": "32", "AD": "2,28,2"}, {}, cfg, True)
    assert "AB_SKEW" in flags
