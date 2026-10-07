"""The likelihood rule must compare a homozygous-alt call with ITS OWN hets.

PL is ordered 0/0, 0/1, 1/1, 0/2, 1/2, 2/2 ... so PL[1] is the 0/1 likelihood. For a
2/2 call that is the likelihood of a genotype without allele 2 in it. Reading it as
"the het margin" let a contradicted 2/2 call pass.
"""

from __future__ import annotations

import pytest

from admissible.checks.genotype import (
    GenotypeConfig,
    evaluate_genotype,
    hom_alt_het_margin,
)

CFG = GenotypeConfig()


def _flags(gt, pl, ad):
    cell = {"GT": gt, "DP": "30", "GQ": "99", "PL": pl, "AD": ad}
    flags, ev = evaluate_genotype(gt, cell, {}, CFG, True)
    return set(flags), ev


def test_a_contradicted_2_2_call_is_flagged():
    # ALT=G,T. 0/2 is only 5 behind 2/2; 0/1 is far away (the old code looked here).
    flags, ev = _flags("2/2", "255,200,255,5,60,0", "0,0,30")
    assert "HOM_CONTRADICTED_BY_LIKELIHOOD" in flags
    assert "FALSE_HOM_SUSPECT" in flags
    assert ev["PL_het"] == 5


def test_a_clean_2_2_call_is_not_flagged():
    flags, ev = _flags("2/2", "255,200,255,150,90,0", "0,0,30")
    assert "HOM_CONTRADICTED_BY_LIKELIHOOD" not in flags
    assert ev["PL_het"] == 90


def test_an_alt_alt_het_counts_as_a_competitor():
    # 1/1 call where 1/2 (index 4) is 6 behind: still a het, still contradicts.
    flags, ev = _flags("1/1", "255,200,0,255,6,255", "0,30,0")
    assert "HOM_CONTRADICTED_BY_LIKELIHOOD" in flags
    assert ev["PL_het"] == 6


def test_biallelic_behaviour_is_unchanged():
    flags, ev = _flags("1/1", "255,12,0", "0,30")
    assert "HOM_CONTRADICTED_BY_LIKELIHOOD" in flags and ev["PL_het"] == 12
    flags, ev = _flags("1/1", "255,60,0", "0,30")
    assert "HOM_CONTRADICTED_BY_LIKELIHOOD" not in flags and ev["PL_het"] == 60


@pytest.mark.parametrize(
    ("pl", "allele", "want"),
    [
        ([90, 40, 0], 1, 40),                      # biallelic
        ([90, 40, 0, 7, 55, 30], 1, 40),           # triallelic, 1/1: 0/1=40, 1/2=55
        ([90, 40, 70, 7, 55, 0], 2, 7),            # triallelic, 2/2: 0/2=7, 1/2=55
        ([90, 40, 0, 7, 55, 30], 2, 7 - 30),       # called 2/2 is not the best: negative
        ([0, 40], 1, None),                        # list too short to hold 1/1
    ],
)
def test_margin_helper(pl, allele, want):
    assert hom_alt_het_margin(pl, allele) == want


# --- haploid calls (male chrX outside PAR, chrY, mitochondrion) -------------


def test_a_clean_haploid_alt_call_is_not_flagged_by_the_likelihood_rule():
    # Haploid PL has one entry per allele: ref 255, alt 0. The diploid indexing read
    # PL[1]=0 as "het is no worse than the call" and flagged every haploid alt call.
    flags, ev = _flags("1", "255,0", "0,30")
    assert "HOM_CONTRADICTED_BY_LIKELIHOOD" not in flags
    assert "FALSE_HOM_SUSPECT" not in flags
    assert ev["PL_het"] == 255


def test_a_haploid_call_whose_reference_is_nearly_as_likely_is_flagged():
    flags, ev = _flags("1", "6,0", "10,20")
    assert "HOM_CONTRADICTED_BY_LIKELIHOOD" in flags
    assert ev["PL_het"] == 6


def test_a_half_call_has_no_likelihood_margin():
    flags, ev = _flags("./1", "255,5,0", "0,30")
    assert "HOM_CONTRADICTED_BY_LIKELIHOOD" not in flags
    assert "PL_het" not in ev


def test_haploid_margin_helper():
    assert hom_alt_het_margin([255, 0], 1, ploidy=1) == 255
    assert hom_alt_het_margin([255, 0, 40], 1, ploidy=1) == 40
    assert hom_alt_het_margin([0, 255], 1, ploidy=1) == -255  # called allele is not best
    assert hom_alt_het_margin([0], 1, ploidy=1) is None
    assert hom_alt_het_margin([0, 5, 9], 1, ploidy=3) is None
