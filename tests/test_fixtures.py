"""The fixtures must keep reproducing the numbers the checks are tuned against."""

from __future__ import annotations


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
