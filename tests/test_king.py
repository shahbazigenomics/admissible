"""KING-robust and IBS0, on a jointly called family with real transmission.

This is the test that earns the claim that KING + IBS0 beats a site-set overlap
statistic: within one family, Jaccard puts parent-offspring and full sibs within
a couple of percent of each other, while IBS0 separates them absolutely.
"""

from __future__ import annotations

import pytest

from admissible.checks.identity import check_identity
from admissible.model import Status
from admissible.ped import read_ped
from admissible.vcfio import load_cohort

PARENTS = {"JF_FATHER", "JF_MOTHER"}


@pytest.fixture(scope="module")
def pairs(joint):
    matrix = load_cohort([joint["vcf"]])
    ped = read_ped(joint["ped"])
    res = check_identity(matrix, ped)
    assert res.status is Status.PASS, res.summary
    return {frozenset((p["a"], p["b"])): p for p in res.metrics["pairs"]}


def _class(pair):
    a, b = tuple(pair)
    if {a, b} == PARENTS:
        return "spouses"
    if a in PARENTS or b in PARENTS:
        return "parent-offspring"
    return "sibs"


def test_joint_vcf_is_recognised_as_dense(pairs):
    assert all(p["mode"] == "dense" for p in pairs.values())


def test_unrelated_founders_land_near_zero(pairs):
    for key, p in pairs.items():
        if _class(key) == "spouses":
            assert abs(p["king_robust_phi"]) < 0.0442
            assert p["king_band"] == "unrelated"


def test_first_degree_relatives_land_near_one_quarter(pairs):
    for key, p in pairs.items():
        if _class(key) in ("parent-offspring", "sibs"):
            assert p["king_robust_phi"] == pytest.approx(0.25, abs=0.03)
            assert p["king_band"] == "first-degree"


def test_ibs0_separates_parent_offspring_from_sibs(pairs):
    po = [p["ibs0_rate"] for k, p in pairs.items() if _class(k) == "parent-offspring"]
    sib = [p["ibs0_rate"] for k, p in pairs.items() if _class(k) == "sibs"]
    assert all(x == 0.0 for x in po), "a parent and child cannot be IBS0 anywhere"
    assert all(x > 0.005 for x in sib)
    assert min(sib) > max(po)


def test_jaccard_alone_cannot_make_that_distinction(pairs):
    """The negative control for the positive result above."""
    po = [p["jaccard"] for k, p in pairs.items() if _class(k) == "parent-offspring"]
    sib = [p["jaccard"] for k, p in pairs.items() if _class(k) == "sibs"]
    assert abs(sum(sib) / len(sib) - sum(po) / len(po)) < 0.05
