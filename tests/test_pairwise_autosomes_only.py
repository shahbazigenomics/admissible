"""0.1.1: pairwise kinship and IBS0 use autosomes only.

Found by Ali's validation, then confirmed on CEPH 1463: father-son pairs were
reported as "sibling-like" (PARENT_OFFSPRING_NOT_SUPPORTED, BLOCKING) because
a son's X comes from his mother, so a father and his son share no X. Written as
hom-ref / hom-alt (hemizygous), they disagree at roughly a third of their X sites
and every disagreement of that kind is an IBS0 site.
"""

from __future__ import annotations

import random

from admissible.checks.identity import IdentityConfig, check_identity
from admissible.model import Severity
from admissible.ped import read_ped
from admissible.vcfio import load_cohort

HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=1,length=249250621>\n"
    "##contig=<ID=X,length=155270560>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tF\tM\tSON\tDAU\n"
)


def _gt(a, b):
    return "/".join(sorted((str(a), str(b))))


def _vcf(tmp_path, n_auto=6000, n_x=3000, seed=11):
    rng = random.Random(seed)
    rows = [HEADER]
    for i in range(n_auto):
        p = rng.uniform(0.15, 0.5)
        f = [int(rng.random() < p) for _ in range(2)]
        m = [int(rng.random() < p) for _ in range(2)]
        gts = [
            _gt(*f),
            _gt(*m),
            _gt(rng.choice(f), rng.choice(m)),
            _gt(rng.choice(f), rng.choice(m)),
        ]
        if all(g == "0/0" for g in gts):
            continue  # not a variant site; forcing a het would bias one sample
        rows.append(f"1\t{1000 + i}\t.\tA\tG\t99\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")
    # non-PAR chrX: the father has one X, written hom. A son takes his mother's X,
    # a daughter takes her father's X and one of her mother's.
    for i in range(n_x):
        p = rng.uniform(0.15, 0.5)
        fx = int(rng.random() < p)
        m = [int(rng.random() < p) for _ in range(2)]
        son = rng.choice(m)
        dau = [fx, rng.choice(m)]
        gts = [_gt(fx, fx), _gt(*m), _gt(son, son), _gt(*dau)]
        if all(g == "0/0" for g in gts):
            continue
        rows.append(f"X\t{80_000_000 + i * 10}\t.\tA\tG\t99\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")
    path = tmp_path / "x.vcf"
    path.write_text("".join(rows))
    return path


def _ped(tmp_path):
    path = tmp_path / "x.ped"
    path.write_text(
        "FAM\tF\t0\t0\t1\t1\n"
        "FAM\tM\t0\t0\t2\t1\n"
        "FAM\tSON\tF\tM\t1\t1\n"
        "FAM\tDAU\tF\tM\t2\t1\n"
    )
    return read_ped(path)


def _run(tmp_path):
    return check_identity(load_cohort([_vcf(tmp_path)]), ped=_ped(tmp_path), cfg=IdentityConfig())


def test_father_and_son_are_not_flagged_for_not_sharing_an_x(tmp_path):
    res = _run(tmp_path)
    assert not [f for f in res.findings if f.code == "PARENT_OFFSPRING_NOT_SUPPORTED"]
    assert not [f for f in res.findings if f.severity is Severity.BLOCKING]


def test_father_son_kinship_and_ibs0_are_those_of_a_parent_and_child(tmp_path):
    res = _run(tmp_path)
    pair = next(p for p in res.metrics["pairs"] if {p["a"], p["b"]} == {"F", "SON"})
    assert abs(pair["king_robust_phi"] - 0.25) < 0.03
    assert pair["n_ibs0"] <= 5


def test_every_parent_child_pair_is_equally_clean(tmp_path):
    """Father-daughter, mother-son and mother-daughter are not special either."""
    res = _run(tmp_path)
    pairs = {frozenset((p["a"], p["b"])): p for p in res.metrics["pairs"]}
    for a, b in (("F", "DAU"), ("M", "SON"), ("M", "DAU"), ("F", "SON")):
        assert abs(pairs[frozenset((a, b))]["king_robust_phi"] - 0.25) < 0.03, (a, b)
