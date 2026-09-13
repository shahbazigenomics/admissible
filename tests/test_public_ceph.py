"""Validation against public data with published truth: CEPH pedigree 1463.

Every other test in this suite runs on fixtures written by the same person who
wrote the code, which is a weak form of evidence.  This one runs on the Utah
three-generation CEPH family as distributed with `peddy
<https://github.com/brentp/peddy>`_ - a real joint-called VCF, a real published
pedigree, and, usefully, a deliberately corrupted copy of that pedigree that
peddy ships as its own canonical failure case.

The data is not vendored (6.7 MB), so these tests skip unless it is present::

    git clone --depth 1 https://github.com/brentp/peddy /tmp/peddy

What this cohort tests that the synthetic fixtures cannot:

* a genuine multi-generation pedigree containing first-degree, **second-degree**
  (grandparent-grandchild) and unrelated pairs in one file;
* a different caller (freebayes, ``GL`` rather than ``PL``);
* GRCh37 with no ``##contig`` header lines, so the build cannot be inferred;
* realised IBD variance between true full sibs, which is what makes band-label
  comparison brittle and motivated the conservative reconciliation rule.
"""

from __future__ import annotations

import itertools
from pathlib import Path

import pytest

from admissible.checks.identity import check_identity
from admissible.model import Status
from admissible.ped import read_ped
from admissible.vcfio import load_cohort

DATA = Path("/tmp/peddy/data")
VCF = DATA / "ceph1463.peddy.vcf.gz"
FULL_PED = DATA / "ceph1463.ped"
GOOD_PED = DATA / "ceph1463.good.ped"
BAD_PED = DATA / "ceph1463.bad.ped"

pytestmark = pytest.mark.skipif(
    not VCF.exists(), reason="CEPH 1463 data absent; see this module's docstring"
)


@pytest.fixture(scope="module")
def result():
    return check_identity(load_cohort([VCF]), read_ped(FULL_PED))


def test_the_correct_pedigree_is_accepted(result):
    """No false positives on a real, correct, three-generation pedigree."""
    assert result.status is Status.PASS, result.summary
    assert not [f for f in result.findings if f.code == "PEDIGREE_MISMATCH"]
    assert not [f for f in result.findings if f.code == "SEX_MISMATCH"]


def test_every_sex_call_matches_the_published_pedigree(result):
    ped = read_ped(FULL_PED)
    called = wrong = 0
    for s in result.metrics["sex"]:
        inferred = s["inferred_sex"]
        if inferred == "unknown":
            continue
        called += 1
        if inferred != ped.individuals[s["sample"]].sex.value:
            wrong += 1
    assert called >= 15
    assert wrong == 0


def _by_expected_kinship(result):
    ped = read_ped(FULL_PED)
    out: dict[float, list[float]] = {}
    for p in result.metrics["pairs"]:
        k = ped.kinship(p["a"], p["b"])
        out.setdefault(round(k, 4), []).append(p["king_robust_phi"])
    return out


def test_unrelated_pairs_never_look_related(result):
    """The safety property: no false relatedness. This must hold absolutely."""
    unrelated = _by_expected_kinship(result)[0.0]
    assert len(unrelated) >= 10
    assert max(unrelated) < 0.0442


def test_first_degree_pairs_are_always_recognised(result):
    """The power property, at the degree that matters for a family study."""
    first = _by_expected_kinship(result)[0.25]
    assert len(first) >= 70
    assert min(first) > 0.0884  # comfortably above the unrelated band


def test_second_degree_is_correct_on_average_but_not_per_pair(result):
    """A real limitation, asserted rather than hidden.

    Grandparent-grandchild pairs average the textbook 0.125, but at ~20k sites
    the weakest of them falls below the unrelated band floor. Second-degree
    relatedness cannot be called reliably for an individual pair at this site
    count - which is precisely why the reconciliation rule refuses to report
    first-vs-second-degree discrepancies as pedigree errors.
    """
    second = _by_expected_kinship(result)[0.125]
    assert len(second) >= 20
    assert sum(second) / len(second) == pytest.approx(0.125, abs=0.03)
    assert min(second) < 0.0442  # documents the failure mode
    assert sum(p > 0.0442 for p in second) / len(second) > 0.85


def test_the_corrupted_pedigree_is_caught():
    """peddy's own failure case: a father and daughter swapped.

    Worth noting *why* both engines are needed. Swapping a parent with their
    own child leaves every kinship coefficient unchanged - 0.25 either way - so
    kinship alone cannot see it. The sex check is what catches the swap itself;
    kinship catches the collateral damage to the spouse relationships.
    """
    res = check_identity(load_cohort([VCF]), read_ped(BAD_PED))
    assert res.status is Status.FAIL
    swapped = {s for f in res.findings if f.code == "SEX_MISMATCH" for s in f.subjects}
    assert swapped == {"NA12877", "NA12880"}
    assert [f for f in res.findings if f.code == "PEDIGREE_MISMATCH"]


def test_no_false_duplicates_anywhere_in_a_17_sample_cohort(result):
    assert result.metrics["n_duplicate_pairs"] == 0


def test_kinship_matches_textbook_values():
    ped = read_ped(FULL_PED)
    expected = {
        ("NA12877", "NA12878"): 0.0,  # spouses
        ("NA12877", "NA12879"): 0.25,  # parent-offspring
        ("NA12879", "NA12880"): 0.25,  # full sibs
        ("NA12889", "NA12879"): 0.125,  # grandparent-grandchild
        ("NA12889", "NA12891"): 0.0,  # unrelated founders
        ("NA12877", "NA12891"): 0.0,  # in-laws
    }
    for (a, b), want in expected.items():
        assert ped.kinship(a, b) == pytest.approx(want), f"{a}/{b}"
    assert ped.kinship("NA12879", "NA12879") == pytest.approx(0.5)


def test_all_pairs_are_dense_in_a_joint_vcf(result):
    assert all(p["mode"] == "dense" for p in result.metrics["pairs"])
    assert len(result.metrics["pairs"]) == len(
        list(itertools.combinations(range(17), 2))
    )


def _ped_with_three_affected_sibs(tmp_path):
    """Real genotypes and a real pedigree, with an affection pattern declared.

    CEPH 1463 is a healthy reference family - every phenotype is -9 - so a
    segregation sweep has nothing to segregate. Declaring three of the eleven
    NA12877 x NA12878 children affected keeps every genotype and every
    relationship real while giving the models something to run on.
    """
    with open(FULL_PED) as fh:
        rows = [ln.split() for ln in fh if ln.strip() and not ln.startswith("#")]
    sibs = {r[1] for r in rows if r[2] == "NA12877" and r[3] == "NA12878"}
    affected = set(sorted(sibs)[:3])
    ped_path = tmp_path / "affected.ped"
    with open(ped_path, "w") as fh:
        for r in rows:
            fh.write("\t".join(r[:5] + ["2" if r[1] in affected else "1"]) + "\n")
    return ped_path


# --- checks 2 and 5 on the same real data ----------------------------------
#
# Check 1 was validated here first; checks 2 and 5 rested on synthetic fixtures
# alone, which is the weaker evidence this file exists to avoid.  Running them
# against CEPH found two real defects, and these tests hold the repairs in
# place against real freebayes + VEP output rather than against a fixture
# written to match the code.


def test_allele_balance_is_actually_computed_on_freebayes_output():
    """freebayes writes RO/AO; reading only AD made this arm silently inert."""
    from admissible.checks.genotype import check_genotype

    res = check_genotype([str(VCF)])
    flags: dict[str, int] = {}
    for v in res.metrics["per_sample"].values():
        for k, n in v["flags"].items():
            flags[k] = flags.get(k, 0) + n
    # Before the fix both of these were exactly zero across all 17 samples.
    assert flags.get("AB_SKEW", 0) > 100
    assert flags.get("ALLELE_IMBALANCE_HOM", 0) > 0
    # Real data carries FORMAT/DP, so nothing here is depth-unassessable; if this
    # ever trips, per-sample depth has stopped being read.
    assert res.metrics["hom_at_adequate_depth"]["depth_unknown"] == 0
    # And allele balance finds false homozygotes the likelihood arm does not.
    assert flags["FALSE_HOM_SUSPECT"] > flags["HOM_CONTRADICTED_BY_LIKELIHOOD"]


def test_compound_het_is_evaluable_on_vep_annotation(tmp_path):
    """The CEPH VCF carries full VEP CSQ; compound-het must not be UNKNOWN."""
    from admissible.checks.models import check_models

    ped_path = _ped_with_three_affected_sibs(tmp_path)

    matrix = load_cohort([str(VCF)])
    res = check_models(matrix, read_ped(ped_path))
    entry = next(e for e in res.metrics["profile"] if e["model"] == "compound-het")
    assert entry["status"] == "computed"     # was "not-applicable" before the fix
    assert entry["n_candidates"] > 0


def test_the_de_novo_count_on_real_data_is_qualified(tmp_path):
    """29 apparent de novos over ~20k sites in three sibs is all error.

    The number itself is not wrong - it is what the segregation filter finds -
    but reported bare it reads as a mutation count, which is how it ends up in
    a manuscript. The caveat has to travel with it.
    """
    from admissible.checks.models import check_models

    ped_path = _ped_with_three_affected_sibs(tmp_path)

    res = check_models(load_cohort([str(VCF)]), read_ped(ped_path))
    entry = next(e for e in res.metrics["profile"] if e["model"] == "de-novo")
    assert entry["n_candidates"] > 0
    assert "dominated by genotyping error" in entry["caveat"]
