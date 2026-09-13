#!/usr/bin/env python3
"""A second, small fixture: one family joint-called into a single VCF.

The cohort fixture is built from separate single-sample VCFs, which is the
realistic case and the case that motivated the tool - but in that layout
hom-reference genotypes are unobserved, so KING-robust kinship and IBS0 cannot
be computed.  This fixture exists to exercise the branch that can: genuine
Mendelian transmission, written as one multi-sample VCF with explicit ``0/0``
calls.

Expected, and asserted in the tests:

  parents (unrelated founders)  KING phi ~ 0.00, band "unrelated"
  parent-offspring              KING phi ~ 0.25, band "first-degree", IBS0 ~ 0
  full sibs                     KING phi ~ 0.25, band "first-degree", IBS0 > 0

The IBS0 contrast is the whole point: it separates the two first-degree
relationships, which a site-set overlap statistic cannot do at all.
"""

from __future__ import annotations

import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "joint"

SEED = 7
N_SITES = 20_000
SAMPLES = ["JF_FATHER", "JF_MOTHER", "JF_CHILD1", "JF_CHILD2", "JF_CHILD3"]
CONTIGS = [("1", 249250621), ("2", 243199373), ("X", 155270560), ("Y", 59373566)]


def build(outdir: Path | None = None) -> Path:
    global OUT
    OUT = Path(outdir) if outdir else OUT
    rng = random.Random(SEED)
    OUT.mkdir(parents=True, exist_ok=True)

    afs = [rng.uniform(0.05, 0.5) for _ in range(N_SITES)]
    haps: dict[str, list[tuple[int, int]]] = {}
    for founder in ("JF_FATHER", "JF_MOTHER"):
        haps[founder] = [
            (int(rng.random() < af), int(rng.random() < af)) for af in afs
        ]
    for child in ("JF_CHILD1", "JF_CHILD2", "JF_CHILD3"):
        haps[child] = [
            (rng.choice(haps["JF_FATHER"][i]), rng.choice(haps["JF_MOTHER"][i]))
            for i in range(N_SITES)
        ]

    positions = sorted(rng.sample(range(1_000_000, 240_000_000), N_SITES))
    lines = ["##fileformat=VCFv4.2", "##reference=GRCh37"]
    lines += [f"##contig=<ID={c},length={n}>" for c, n in CONTIGS]
    lines += [
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FILTER=<ID=LowQual,Description="Low quality">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(SAMPLES),
    ]
    for i, pos in enumerate(positions):
        cells = []
        keep = False
        for s in SAMPLES:
            a, b = haps[s][i]
            gt = f"{a}/{b}"
            keep = keep or (a or b)
            cells.append(f"{gt}:{rng.randint(20, 60)}")
        if not keep:
            continue  # a site nobody carries would not be emitted by a caller
        chrom = "1" if pos < 249_000_000 else "2"
        lines.append(
            f"{chrom}\t{pos}\t.\tA\tG\t900.0\tPASS\tAC=1\tGT:DP\t" + "\t".join(cells)
        )

    vcf = OUT / "joint_family.vcf"
    vcf.write_text("\n".join(lines) + "\n")

    ped = OUT / "joint_family.ped"
    ped.write_text(
        "\n".join(
            [
                "#FID\tIID\tPAT\tMAT\tSEX\tPHENO",
                "JF\tJF_FATHER\t0\t0\t1\t1",
                "JF\tJF_MOTHER\t0\t0\t2\t1",
                "JF\tJF_CHILD1\tJF_FATHER\tJF_MOTHER\t1\t2",
                "JF\tJF_CHILD2\tJF_FATHER\tJF_MOTHER\t2\t2",
                "JF\tJF_CHILD3\tJF_FATHER\tJF_MOTHER\t1\t2",
            ]
        )
        + "\n"
    )
    return vcf


if __name__ == "__main__":
    print(build())
