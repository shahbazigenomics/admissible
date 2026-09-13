"""Check 5 - the inheritance model sweep.

Each fixture is a tiny joint VCF carrying one planted variant that satisfies
exactly one model, so a test failure names the model that broke.
"""

from __future__ import annotations

import pytest

from admissible.checks.models import ModelConfig, check_models
from admissible.model import Status
from admissible.ped import read_ped
from admissible.vcfio import load_cohort

SAMPLES = ["DAD", "MUM", "AFF1", "AFF2", "WELL"]
HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=1,length=249250621>\n##contig=<ID=X,length=155270560>\n"
    "##contig=<ID=MT,length=16569>\n"
    "##INFO=<ID=GENEAF,Number=1,Type=Float,Description=\"af\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(SAMPLES) + "\n"
)

PED = (
    "#FID\tIID\tPAT\tMAT\tSEX\tPHENO\n"
    "F\tDAD\t0\t0\t1\t1\n"
    "F\tMUM\t0\t0\t2\t1\n"
    "F\tAFF1\tDAD\tMUM\t1\t2\n"
    "F\tAFF2\tDAD\tMUM\t2\t2\n"
    "F\tWELL\tDAD\tMUM\t2\t1\n"
)

# name, chrom, pos, gene, af, genotypes for DAD MUM AFF1 AFF2 WELL
ROWS = [
    ("recessive",     "1",  1000, "GENEA", 0.001, ["0/1", "0/1", "1/1", "1/1", "0/1"]),
    ("dominant",      "1",  2000, "GENEB", 0.001, ["0/1", "0/0", "0/1", "0/1", "0/0"]),
    ("reduced_pen",   "1",  3000, "GENEC", 0.001, ["0/1", "0/0", "0/1", "0/1", "0/1"]),
    ("de_novo",       "1",  4000, "GENED", 0.001, ["0/0", "0/0", "0/1", "0/0", "0/0"]),
    ("phenocopy",     "1",  5000, "GENEE", 0.001, ["0/1", "0/1", "1/1", "0/1", "0/1"]),
    ("comphet_a",     "1",  6000, "GENEF", 0.001, ["0/1", "0/0", "0/1", "0/1", "0/0"]),
    ("comphet_b",     "1",  6500, "GENEF", 0.001, ["0/0", "0/1", "0/1", "0/1", "0/0"]),
    ("x_recessive",   "X", 10000, "GENEX", 0.001, ["0/0", "0/1", "1/1", "1/1", "0/1"]),
    ("mito",          "MT",  300, "MTGENE", 0.001, ["0/0", "0/1", "0/1", "0/1", "0/0"]),
    # common, segregates recessively: only a rarity filter removes it
    ("common_rec",    "1",  9000, "GENEG", 0.400, ["0/1", "0/1", "1/1", "1/1", "0/1"]),
]


@pytest.fixture
def cohort(tmp_path):
    body = HEADER
    for _n, chrom, pos, gene, af, gts in ROWS:
        cells = "\t".join(f"{g}:30" for g in gts)
        body += (
            f"{chrom}\t{pos}\t.\tA\tG\t500\tPASS\t"
            f"Gene.refGene={gene};GENEAF={af}\tGT:DP\t{cells}\n"
        )
    vcf = tmp_path / "family.vcf"
    vcf.write_text(body)
    ped = tmp_path / "family.ped"
    ped.write_text(PED)
    return vcf, ped


def sweep(cohort, af_key=None, max_af=None, cfg=None):
    vcf, ped = cohort
    matrix = load_cohort([vcf], af_key=af_key)
    return check_models(matrix, read_ped(ped), cfg=cfg or ModelConfig(max_af=max_af))


def counts(res) -> dict[str, int | None]:
    return {m["model"]: m["n_candidates"] for m in res.metrics["profile"]}


def status_of(res, model) -> str:
    return next(m["status"] for m in res.metrics["profile"] if m["model"] == model)


def test_each_model_finds_its_planted_variant(cohort):
    c = counts(sweep(cohort))
    assert c["recessive"] == 2  # the rare one and the common one both segregate
    assert c["de-novo"] == 1
    assert c["phenocopy(n-1)"] == 1
    assert c["mitochondrial"] == 1
    assert c["x-linked-recessive"] == 1
    assert c["compound-het"] == 1  # one gene with two qualifying hets


def test_dominant_is_impossible_with_two_unaffected_genotyped_parents(cohort):
    """Not a limitation - the right answer.

    A dominant allele carried by an unaffected, genotyped parent is not a
    dominant model; it is a reduced-penetrance model. The strict filter must
    return nothing here, and the reduced-penetrance filter must not.
    """
    c = counts(sweep(cohort))
    assert c["dominant"] == 0
    assert c["dominant-reduced-penetrance"] > 0


def test_dominant_fires_when_the_transmitting_parent_is_affected(tmp_path):
    vcf = tmp_path / "dom.vcf"
    vcf.write_text(
        HEADER
        + "1\t2000\t.\tA\tG\t500\tPASS\tGene.refGene=B\tGT:DP\t"
        + "\t".join(f"{g}:30" for g in ["0/1", "0/0", "0/1", "0/1", "0/0"])
        + "\n"
    )
    ped = tmp_path / "dom.ped"
    ped.write_text(PED.replace("F\tDAD\t0\t0\t1\t1", "F\tDAD\t0\t0\t1\t2"))
    res = check_models(load_cohort([vcf]), read_ped(ped))
    assert counts(res)["dominant"] == 1


def test_rarity_filter_removes_the_common_segregating_variant(cohort):
    without = counts(sweep(cohort))["recessive"]
    with_af = counts(sweep(cohort, af_key="GENEAF", max_af=0.01))["recessive"]
    assert without == 2
    assert with_af == 1, "a common variant that happens to segregate must be dropped"


def test_unfiltered_counts_say_so(cohort):
    res = sweep(cohort)
    assert res.metrics["rarity_filtered"] is False
    assert any("NOT rarity-filtered" in n for n in res.notes)
    assert "not rarity-filtered" in res.summary


def test_inapplicable_models_are_not_reported_as_zero(cohort):
    """The dangerous failure: 0 candidates looks like evidence of absence."""
    res = sweep(cohort)
    assert status_of(res, "two-locus") == "not-applicable"
    assert counts(res)["two-locus"] is None


def test_de_novo_needs_both_parents(tmp_path):
    vcf = tmp_path / "solo.vcf"
    hdr = HEADER.replace("\t".join(SAMPLES), "AFF1")
    vcf.write_text(hdr + "1\t4000\t.\tA\tG\t500\tPASS\tGene.refGene=X\tGT:DP\t0/1:30\n")
    ped = tmp_path / "solo.ped"
    ped.write_text("#FID\tIID\tPAT\tMAT\tSEX\tPHENO\nF\tAFF1\t0\t0\t1\t2\n")
    res = check_models(load_cohort([vcf]), read_ped(ped))
    assert status_of(res, "de-novo") == "not-applicable"
    assert "parents" in next(
        m["reason"] for m in res.metrics["profile"] if m["model"] == "de-novo"
    )


def test_compound_het_without_gene_annotation_is_unknown(tmp_path):
    vcf = tmp_path / "nogene.vcf"
    body = HEADER
    for _n, chrom, pos, _g, _af, gts in ROWS[:3]:
        body += f"{chrom}\t{pos}\t.\tA\tG\t500\tPASS\t.\tGT:DP\t" + "\t".join(
            f"{g}:30" for g in gts
        ) + "\n"
    vcf.write_text(body)
    ped = tmp_path / "p.ped"
    ped.write_text(PED)
    res = check_models(load_cohort([vcf]), read_ped(ped))
    assert status_of(res, "compound-het") == "not-applicable"


def test_zero_candidate_models_forbid_the_negative_claim(tmp_path):
    vcf = tmp_path / "empty.vcf"
    vcf.write_text(
        HEADER + "1\t1000\t.\tA\tG\t500\tPASS\tGene.refGene=G\tGT:DP\t"
        + "\t".join("0/0:30" for _ in SAMPLES) + "\n"
    )
    ped = tmp_path / "p.ped"
    ped.write_text(PED)
    res = check_models(load_cohort([vcf]), read_ped(ped))
    f = next(x for x in res.findings if x.code == "NO_CANDIDATES")
    assert f.evidence["do_not_conclude"] == "no monogenic cause"
    assert "not evidence of absence" in f.message


def test_each_model_carries_its_own_callable_fraction(cohort):
    vcf, ped = cohort
    coverage = {
        "recessive": {"fraction_of_target": 0.20, "permitted_claim": "exploratory"},
        "phenocopy(n-1)": {"fraction_of_target": 0.80, "permitted_claim": "excluded"},
    }
    res = check_models(load_cohort([vcf]), read_ped(ped), callability=coverage)
    prof = {m["model"]: m for m in res.metrics["profile"]}
    assert prof["recessive"]["callable_fraction"] == 0.20
    assert prof["phenocopy(n-1)"]["callable_fraction"] == 0.80
    assert prof["dominant"]["callable_fraction"] is None


def test_no_pedigree_is_unknown_not_empty():
    from admissible.ped import Pedigree
    from admissible.vcfio import GenotypeMatrix, SiteIndex

    res = check_models(GenotypeMatrix(index=SiteIndex()), Pedigree())
    assert res.status is Status.UNKNOWN


def test_never_raises_on_garbage(tmp_path):
    vcf = tmp_path / "junk.vcf"
    vcf.write_text("not a vcf\n\x00\x01\n")
    ped = tmp_path / "p.ped"
    ped.write_text(PED)
    assert check_models(load_cohort([vcf]), read_ped(ped)).status in tuple(Status)


def test_vacuous_exclusion_is_called_out(tmp_path):
    """With no genotyped unaffected sample, the dominant filter excludes nothing."""
    vcf = tmp_path / "noctrl.vcf"
    hdr = HEADER.replace("\t".join(SAMPLES), "AFF1\tAFF2")
    vcf.write_text(
        hdr + "1\t2000\t.\tA\tG\t500\tPASS\tGene.refGene=B\tGT:DP\t0/1:30\t0/1:30\n"
    )
    ped = tmp_path / "noctrl.ped"
    ped.write_text(
        "#FID\tIID\tPAT\tMAT\tSEX\tPHENO\nF\tAFF1\t0\t0\t1\t2\nF\tAFF2\t0\t0\t2\t2\n"
    )
    res = check_models(load_cohort([vcf]), read_ped(ped))
    assert any("vacuous" in n for n in res.notes)
    assert any("no genotyped unaffected carrier" in n for n in res.notes)


def test_multisample_file_without_homref_is_flagged(tmp_path):
    """A file that looks jointly called but carries no hom-reference calls.

    Real merges do this: every non-carrier becomes './.', so a single string now
    means both "hom-reference" and "never covered". IBS0 depends on telling those
    apart, so the tool must refuse it AND say why.
    """
    from admissible.checks.identity import check_identity
    from admissible.checks.provenance import check_provenance
    from admissible.vcfio import SiteIndex, load_cohort, scan_vcf

    vcf = tmp_path / "merged.vcf"
    rows = "".join(
        f"1\t{p}\t.\tA\tG\t500\tPASS\tAC=2;AN=4\tGT:DP\t"
        + "\t".join(["0/1:30", "./.:.", "0/1:30", "1/1:30"]) + "\n"
        for p in range(1000, 6000, 1000)
    )
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD\n" + rows
    )
    res = check_provenance([scan_vcf(vcf, SiteIndex())])
    f = next(x for x in res.findings if x.code == "HOMREF_STRIPPED")
    assert "looks jointly called and is not" in f.message

    ped = tmp_path / "p.ped"
    ped.write_text(
        "#FID\tIID\tPAT\tMAT\tSEX\tPHENO\nF\tA\t0\t0\t1\t1\nF\tB\t0\t0\t2\t1\n"
        "F\tC\tA\tB\t1\t2\nF\tD\tA\tB\t2\t2\n"
    )
    ident = check_identity(load_cohort([vcf]), read_ped(ped))
    assert any("looks jointly called and is not" in n for n in ident.notes)
    assert all(p["mode"] == "sparse-nonref" for p in ident.metrics["pairs"])


def test_a_missing_genotype_does_not_satisfy_an_exclusion_clause(tmp_path):
    """The bug this test exists for, found on a real merged family VCF.

    Every model says some version of "and the unaffected do not carry it". A
    sample with no genotype does not satisfy that - it says nothing. Counting it
    as a non-carrier silently inflates every candidate list on any call set with
    a meaningful missing rate, and on a merge that writes non-carriers as './.'
    it inflates them enormously.
    """
    vcf = tmp_path / "gappy.vcf"
    # Both sites look recessive for the affected pair. At the first, the healthy
    # sibling is observed hom-ref; at the second, it has no genotype at all.
    body = HEADER
    body += (
        "1\t1000\t.\tA\tG\t500\tPASS\tGene.refGene=A\tGT:DP\t"
        + "\t".join(f"{g}:30" for g in ["0/1", "0/1", "1/1", "1/1", "0/0"]) + "\n"
    )
    body += (
        "1\t2000\t.\tA\tG\t500\tPASS\tGene.refGene=B\tGT:DP\t"
        + "\t".join(["0/1:30", "0/1:30", "1/1:30", "1/1:30", "./.:."]) + "\n"
    )
    vcf.write_text(body)
    ped = tmp_path / "p.ped"
    ped.write_text(PED)
    res = check_models(load_cohort([vcf]), read_ped(ped))
    prof = {m["model"]: m for m in res.metrics["profile"]}

    assert prof["recessive"]["n_candidates"] == 1, "only the verified site counts"
    assert prof["recessive"]["n_exclusion_unverified"] == 1
    assert res.metrics["n_exclusion_unverified_total"] >= 1
    assert any("not a negative result" in n for n in res.notes)
    assert "unverifiable" in res.summary


def test_models_without_an_exclusion_clause_are_unaffected(tmp_path):
    """Reduced penetrance tolerates unaffected carriers, so nothing to verify."""
    vcf = tmp_path / "rp.vcf"
    vcf.write_text(
        HEADER + "1\t3000\t.\tA\tG\t500\tPASS\tGene.refGene=C\tGT:DP\t"
        + "\t".join(["0/1:30", "0/1:30", "0/1:30", "0/1:30", "./.:."]) + "\n"
    )
    ped = tmp_path / "p.ped"
    ped.write_text(PED)
    res = check_models(load_cohort([vcf]), read_ped(ped))
    prof = {m["model"]: m for m in res.metrics["profile"]}
    assert prof["dominant-reduced-penetrance"]["n_candidates"] == 1
    assert prof["dominant-reduced-penetrance"]["n_exclusion_unverified"] == 0


# --- gene assignment from VEP / SnpEff --------------------------------------


def test_vep_csq_supplies_the_gene_for_compound_het(tmp_path):
    """Compound-het was unevaluable on VEP-annotated VCFs, i.e. on most of them.

    VEP declares its own pipe-delimited layout in the CSQ header line, so the
    gene column has to be located rather than assumed. Found by running check 5
    on the public CEPH 1463 call set, where compound-het reported
    not-applicable despite the file carrying full VEP annotation.
    """
    from admissible.vcfio import SiteIndex, scan_vcf

    p = tmp_path / "vep.vcf"
    p.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type as '
        'predicted by VEP. Format: Consequence|Codons|Gene|SYMBOL|Feature">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t100\t.\tA\tG\t50\tPASS\tCSQ=missense|gTt/gCt|ENSG001|NOD2|ENST1"
        "\tGT\t0/1\n"
        "chr1\t200\t.\tC\tT\t50\tPASS\tCSQ=missense|cGg/cAg|ENSG001|NOD2|ENST1"
        "\tGT\t0/1\n"
    )
    scan = scan_vcf(p, SiteIndex())
    assert scan.header.csq_gene == ("CSQ", 3)      # the SYMBOL column
    assert sorted(scan.site_gene.values()) == ["NOD2", "NOD2"]


def test_snpeff_ann_supplies_the_gene(tmp_path):
    from admissible.vcfio import SiteIndex, scan_vcf

    p = tmp_path / "snpeff.vcf"
    p.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: '
        "'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID'\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t100\t.\tA\tG\t50\tPASS\tANN=G|missense_variant|MODERATE|IL10RA|ENSG9"
        "\tGT\t0/1\n"
    )
    scan = scan_vcf(p, SiteIndex())
    assert scan.header.csq_gene == ("ANN", 3)      # Gene_Name
    assert list(scan.site_gene.values()) == ["IL10RA"]


def test_a_de_novo_count_carries_its_error_caveat(joint):
    """A raw de-novo count at exome scale is error, not mutation."""
    from admissible.checks.models import check_models
    from admissible.ped import read_ped
    from admissible.vcfio import load_cohort

    matrix = load_cohort([joint["vcf"]])
    res = check_models(matrix, read_ped(joint["ped"]))
    entry = next(e for e in res.metrics["profile"] if e["model"] == "de-novo")
    if entry["status"] == "computed" and entry["n_candidates"]:
        assert "dominated by genotyping error" in entry["caveat"]
        assert any("de-novo" in n and "worklist" in n for n in res.notes)
