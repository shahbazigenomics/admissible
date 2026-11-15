"""The parent-offspring IBS0 rule counts biallelic SNVs only and only ever warns.

Measured on 1000 Genomes: genotype noise (worst on indels) puts many true parent-child
pairs above the old 0.005 limit, and the one clear sibling pair sits inside that noise.
"""

from __future__ import annotations

import random

from admissible.checks.identity import IdentityConfig, check_identity
from admissible.model import Severity
from admissible.ped import read_ped
from admissible.vcfio import load_cohort

HEADER = (
    "##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tF\tM\tC\n"
)
PED = "FAM\tF\t0\t0\t1\t1\nFAM\tM\t0\t0\t2\t1\nFAM\tC\tF\tM\t1\t1\n"


def _gt(a, b):
    return "/".join(sorted((str(a), str(b))))


def _noisy(g, rng):
    return rng.choice(["0/0", "0/1", "1/1"]) if g else g


def _run(tmp_path, snv_noise=0.0, indel_noise=0.0, n_snv=8000, n_indel=3000, seed=3):
    rng = random.Random(seed)
    rows = [HEADER]
    pos = 1000
    for kind, n, noise in (("snv", n_snv, snv_noise), ("indel", n_indel, indel_noise)):
        ref, alt = ("A", "G") if kind == "snv" else ("AT", "A")
        for _ in range(n):
            p = rng.uniform(0.15, 0.5)
            f = [int(rng.random() < p) for _ in range(2)]
            m = [int(rng.random() < p) for _ in range(2)]
            c = [rng.choice(f), rng.choice(m)]
            gc = _gt(*c)
            if rng.random() < noise:
                gc = rng.choice(["0/0", "0/1", "1/1"])
            gts = [_gt(*f), _gt(*m), gc]
            if all(g == "0/0" for g in gts):
                gts[0] = "0/1"
            pos += 3
            rows.append(f"1\t{pos}\t.\t{ref}\t{alt}\t99\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")
    vcf = tmp_path / "t.vcf"
    vcf.write_text("".join(rows))
    ped = tmp_path / "t.ped"
    ped.write_text(PED)
    return check_identity(load_cohort([vcf]), ped=read_ped(ped), cfg=IdentityConfig())


def _hits(res):
    return [f for f in res.findings if f.code == "PARENT_OFFSPRING_IBS0_HIGH"]


def test_clean_trio_has_no_finding(tmp_path):
    assert not _hits(_run(tmp_path))


def test_indel_noise_does_not_trigger_the_rule(tmp_path):
    res = _run(tmp_path, snv_noise=0.0, indel_noise=0.15)
    assert not _hits(res)  # the same noise on SNVs would trip it, see below


def test_snv_noise_only_warns(tmp_path):
    res = _run(tmp_path, snv_noise=0.08)
    hits = _hits(res)
    assert hits, "a pair this noisy is above the sibling-like edge"
    assert all(f.severity is Severity.WARN for f in hits)
    assert not [f for f in res.findings if f.severity is Severity.BLOCKING]


def test_the_rule_reports_the_snv_figures(tmp_path):
    f = _hits(_run(tmp_path, snv_noise=0.08))[0]
    assert f.evidence["ibs0_per_het_snv"] > IdentityConfig().parent_ibs0_warn
    assert f.evidence["n_ibs0_snv"] >= IdentityConfig().parent_ibs0_min_count
    assert f.evidence["n_het_snv"] > 0
