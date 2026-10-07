"""0.1.1: a sibling declared as a parent is caught through IBS0."""

from __future__ import annotations

import random

from admissible.checks.identity import IdentityConfig, check_identity
from admissible.model import Severity
from admissible.ped import Relationship, read_ped
from admissible.vcfio import load_cohort

# --- bug 2: a sibling declared as a parent ---------------------------------

HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=1,length=249250621>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tF\tM\tC1\tC2\n"
)


def _gt(a, b):
    return "/".join(sorted((str(a), str(b))))


def _family_vcf(tmp_path, n_sites=6000, seed=7):
    rng = random.Random(seed)
    rows = [HEADER]
    for i in range(n_sites):
        p = rng.uniform(0.15, 0.5)
        f = [int(rng.random() < p) for _ in range(2)]
        m = [int(rng.random() < p) for _ in range(2)]
        c1 = [rng.choice(f), rng.choice(m)]
        c2 = [rng.choice(f), rng.choice(m)]
        gts = [_gt(*f), _gt(*m), _gt(*c1), _gt(*c2)]
        if all(g == "0/0" for g in gts):
            gts[0] = "0/1"
        rows.append(f"1\t{1000 + i}\t.\tA\tG\t99\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")
    path = tmp_path / "fam.vcf"
    path.write_text("".join(rows))
    return path


def _ped(tmp_path, name, c2_father):
    path = tmp_path / name
    path.write_text(
        "FAM\tF\t0\t0\t1\t1\n"
        "FAM\tM\t0\t0\t2\t1\n"
        "FAM\tC1\tF\tM\t1\t1\n"
        f"FAM\tC2\t{c2_father}\tM\t1\t2\n"
    )
    return read_ped(path)


def _audit(tmp_path, ped):
    matrix = load_cohort([_family_vcf(tmp_path)])
    return check_identity(matrix, ped=ped, cfg=IdentityConfig())


def test_correct_pedigree_has_no_parent_offspring_finding(tmp_path):
    res = _audit(tmp_path, _ped(tmp_path, "ok.ped", "F"))
    assert not [f for f in res.findings if f.code == "PARENT_OFFSPRING_NOT_SUPPORTED"]
    assert not [f for f in res.findings if f.severity is Severity.BLOCKING]


def test_a_sibling_declared_as_a_parent_is_caught(tmp_path):
    ped = _ped(tmp_path, "swap.ped", "C1")  # C2's "father" is really C2's brother
    assert ped.relationship("C1", "C2") is Relationship.PARENT_OFFSPRING
    res = _audit(tmp_path, ped)
    hits = [f for f in res.findings if f.code == "PARENT_OFFSPRING_NOT_SUPPORTED"]
    assert [set(f.subjects) for f in hits] == [{"C1", "C2"}]
    f = hits[0]
    assert f.severity is Severity.BLOCKING
    assert f.evidence["ibs0_per_het"] > IdentityConfig().parent_ibs0_max
    assert f.evidence["n_ibs0"] >= IdentityConfig().parent_ibs0_min_count
    # the true parent of C2 (M) is not accused
    assert all(set(x.subjects) != {"M", "C2"} for x in hits)
