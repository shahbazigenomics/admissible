"""The fixtures must keep reproducing the numbers the checks are tuned against."""

from __future__ import annotations

import pytest


def test_duplicate_pair_is_exact(cohort):
    j, agree, shared = cohort["metrics"]["FAMB_04_vs_FAMC_04"]
    assert j == 1.0
    assert agree == 1.0
    assert shared == 18_205


def test_duplicate_files_differ_at_one_byte(cohort):
    m = cohort["metrics"]
    assert m["file_bytes_FAMB_04"] == m["file_bytes_FAMC_04"]
    assert m["byte_differences"] == 1


def test_resequenced_pair_hits_targets(cohort):
    j, agree = cohort["metrics"]["FAMB_05_vs_FAMC_03"]
    assert j == 0.7508
    assert agree == 0.9726


def test_relatedness_bands_are_narrowly_separated(cohort):
    lo_u, hi_u = cohort["metrics"]["unrelated_jaccard_range"]
    lo_s, hi_s = cohort["metrics"]["sib_jaccard_range"]
    assert 0.33 <= lo_u and hi_u <= 0.39
    assert 0.47 <= lo_s and hi_s <= 0.53
    # The bands must not overlap, but the gap is deliberately thin.
    assert hi_u < lo_s
    assert lo_s - hi_u < 0.15


def test_planted_sex_profiles(cohort):
    assert cohort["metrics"]["chrX_het_FAMB_04"] == 0.297
    assert cohort["metrics"]["chrX_het_FAMB_05"] == 0.376


def test_all_files_in_wes_size_range(cohort):
    counts = cohort["metrics"]["variant_counts"].values()
    assert all(12_000 <= n <= 19_000 for n in counts)


# --- the exact binomial has to survive real depth -------------------------


def test_allele_balance_survives_deep_coverage():
    """math.comb(n, k) overflows a float long before real data runs out of depth.

    A real public VCF (nf-core genmod_compound) carries AD values up to 7,166 -
    an ordinary pileup depth in panel or amplicon data - and the direct
    ``comb(n, k) * p**k`` form raised OverflowError there, taking the whole of
    check 2 down with it for that file. Working in log space fixes it, and the
    normal approximation above 20,000 trials keeps the cost from growing with
    depth.
    """
    from admissible.stats import binom_two_sided_p

    for k, n in [(500, 7166), (3583, 7166), (1, 7166), (10_000, 60_000)]:
        p = binom_two_sided_p(k, n)
        assert 0.0 <= p <= 1.0
    # A balanced deep heterozygote is unremarkable; a skewed one is not.
    assert binom_two_sided_p(3583, 7166) == pytest.approx(1.0, abs=1e-9)
    assert binom_two_sided_p(500, 7166) < 1e-100


def test_binomial_matches_known_exact_values():
    """Log-space arithmetic must not change any answer it could already give."""
    from admissible.stats import binom_two_sided_p

    assert binom_two_sided_p(1, 10) == pytest.approx(0.021484375)
    assert binom_two_sided_p(3, 10) == pytest.approx(0.34375)
    assert binom_two_sided_p(5, 10) == pytest.approx(1.0)
    assert binom_two_sided_p(0, 20) == pytest.approx(2 / 2**20)


def test_the_approximation_agrees_with_the_exact_sum_where_it_matters():
    """Near any usable alpha the two forms must be interchangeable."""
    import math

    from admissible.stats import binom_two_sided_p

    n = 20_000                      # the exact/approximate switch point
    for k in (9800, 9750):
        exact = binom_two_sided_p(k, n)
        sd = math.sqrt(n * 0.25)
        z = (abs(k - n / 2) - 0.5) / sd
        approx = math.erfc(z / math.sqrt(2))
        assert abs(exact - approx) / exact < 0.01
