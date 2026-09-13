"""Check 1 must recover every error planted in the fixture cohort."""

from __future__ import annotations

import pytest

from admissible.checks.identity import check_identity
from admissible.model import Severity, Status
from admissible.ped import read_ped
from admissible.vcfio import load_cohort


@pytest.fixture(scope="module")
def result(cohort):
    matrix = load_cohort(list(cohort["vcfs"]))
    ped = read_ped(cohort["ped"])
    return check_identity(matrix, ped)


def codes(result, code):
    return [f for f in result.findings if f.code == code]


def test_overall_status_is_fail(result):
    assert result.status is Status.FAIL
    assert result.blocking


def test_both_sex_mismatches_found(result):
    found = {f.subjects[0] for f in codes(result, "SEX_MISMATCH")}
    assert found == {"FAMB_04", "FAMB_05"}


def test_sex_call_reports_an_interval_not_a_bare_number(result):
    f = codes(result, "SEX_MISMATCH")[0]
    ev = f.evidence
    assert ev["chrX_sites"] > 300
    lo, hi = ev["chrX_het_ci95"]
    assert lo < ev["chrX_het_frac"] < hi


def test_chry_is_not_used_for_the_sex_call(result):
    # Every male in the fixture has a single-digit chrY count, as a coding-only
    # VCF gives.  The sex call must not depend on it.
    for f in codes(result, "SEX_MISMATCH"):
        assert f.evidence["chrY_sites"] < 10
        assert "chrX" in f.evidence["basis"]


def test_identical_duplicate_is_found_and_typed(result):
    hits = codes(result, "DUPLICATE_CONTENT_IDENTICAL")
    assert len(hits) == 1
    f = hits[0]
    assert set(f.subjects) == {"FAMB_04", "FAMC_04"}
    assert f.evidence["cross_family"] is True
    assert f.evidence["jaccard"] == 1.0
    assert f.evidence["file_byte_size_delta"] == 0
    assert f.evidence["file_sha256_equal"] is False
    assert f.severity is Severity.BLOCKING


def test_resequenced_duplicate_is_found_and_typed_differently(result):
    hits = codes(result, "DUPLICATE_SAME_INDIVIDUAL")
    assert len(hits) == 1
    f = hits[0]
    assert set(f.subjects) == {"FAMB_05", "FAMC_03"}
    assert f.evidence["jaccard"] == pytest.approx(0.7508, abs=5e-4)
    assert f.evidence["agreement"] == pytest.approx(0.9726, abs=5e-4)


def test_duplicate_message_states_hypotheses_not_a_cause(result):
    for code in ("DUPLICATE_CONTENT_IDENTICAL", "DUPLICATE_SAME_INDIVIDUAL"):
        assert "compatible with" in codes(result, code)[0].message


def test_duplicates_are_detected_across_families_not_within(result):
    # Both planted swaps live between two families.  A within-family scan
    # would miss both, which is the point of scanning cohort-wide.
    for code in ("DUPLICATE_CONTENT_IDENTICAL", "DUPLICATE_SAME_INDIVIDUAL"):
        assert codes(result, code)[0].evidence["cross_family"] is True


def test_pedigree_mismatches_implicate_the_swapped_samples(result):
    hits = codes(result, "PEDIGREE_MISMATCH")
    assert hits
    for f in hits:
        assert {"FAMB_04", "FAMB_05"} & set(f.subjects)


def test_clean_families_raise_nothing(result):
    implicated = {s for f in result.findings for s in f.subjects}
    assert not any(s.startswith("FAMA_") for s in implicated)


def test_boundary_is_calibrated_on_the_cohort(result):
    calib = result.metrics["relatedness_calibration"]
    assert calib["n_declared_unrelated_pairs"] >= 3
    assert calib["n_declared_first_degree_pairs"] >= 3
    assert 0.38 < calib["boundary"] < 0.46
    assert calib["separation"] > 2.0


def test_sparse_mode_declines_to_estimate_degree(result):
    joined = " ".join(result.notes)
    assert "KING-robust" in joined and "cannot be computed" in joined
    assert all(p["mode"] == "sparse-nonref" for p in result.metrics["pairs"])


def test_build_was_inferred(result):
    assert result.metrics["genome_build"] == "GRCh37"
