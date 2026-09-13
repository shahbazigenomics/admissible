#!/usr/bin/env python3
"""Generate the synthetic validation cohort.

**No real data.**  Every genotype here is simulated.  The cohort is built to
reproduce, in a file anyone can regenerate from this script, the class of
errors this tool exists to find:

  13 samples, 3 families (4 / 5 / 4), GRCh37, coding-only and PASS-only.

  * ``FAMB_04`` is declared a daughter in family B.  Its file is a copy of
    ``FAMC_04``'s file, identical except for a single byte (the sample name in
    the ``#CHROM`` line).  Jaccard 1.0000, agreement 1.0000.
  * ``FAMB_05`` is declared a daughter in family B.  It is the same individual
    as ``FAMC_03``, sequenced to a lower depth: Jaccard 0.7508, agreement
    0.9726, and a visibly smaller file.
  * Both of those samples are male by chrX heterozygosity (29.7% and 37.6%)
    while their pedigree slots declare them female.
  * Declared-unrelated pairs sit at Jaccard 0.34-0.39; true first-degree pairs
    at 0.48-0.53, which is a narrow separation on purpose - that thin margin
    is why the boundary is calibrated on the cohort instead of hard-coded.

Run ``python tests/fixtures/make_fixtures.py`` to (re)write the VCFs and the
PED, and print the achieved metrics.
"""

from __future__ import annotations

import random
import zlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
COHORT = HERE / "cohort"

SEED = 20260912

# GRCh37 lengths.  Only 1, X and Y are used to infer the build.
CONTIGS = [
    ("1", 249250621), ("2", 243199373), ("3", 198022430), ("4", 191154276),
    ("5", 180915260), ("6", 171115067), ("7", 159138663), ("8", 146364022),
    ("9", 141213431), ("10", 135534747), ("11", 135006516), ("12", 133851895),
    ("13", 115169878), ("14", 107349540), ("15", 102531392), ("16", 90354753),
    ("17", 81195210), ("18", 78077248), ("19", 59128983), ("20", 63025520),
    ("21", 48129895), ("22", 51304566), ("X", 155270560), ("Y", 59373566),
]
CONTIG_ORDER = {name: i for i, (name, _) in enumerate(CONTIGS)}
AUTOSOMES = [c for c, _ in CONTIGS if c.isdigit()]

# Region mix: 92% exonic / 3.6% intronic, i.e. a coding-only call set.
REGION_MIX = [
    ("exonic", 0.920), ("intronic", 0.036), ("splicing", 0.012),
    ("UTR3", 0.015), ("UTR5", 0.010), ("ncRNA_exonic", 0.007),
]

# Set sizes solved so that, at these sample sizes, cross-family Jaccard lands
# in 0.34-0.39 and within-family Jaccard in 0.48-0.53.  See module docstring.
N_COMMON = 8628          # sites every sample carries
N_FAMILY_EXTRA = 2400    # extra sites shared within one family
PRIVATE_RANGE = (4900, 5500)
N_PRIVATE_C04 = 6665     # makes FAMC_04 carry exactly 18,205 variants

HET, HOMALT = 1, 2
P_HET = 0.78

FAMILIES = {
    "FAM_A": ["FAMA_01", "FAMA_02", "FAMA_03", "FAMA_04"],
    "FAM_B": ["FAMB_01", "FAMB_02", "FAMB_03", "FAMB_04", "FAMB_05"],
    "FAM_C": ["FAMC_01", "FAMC_02", "FAMC_03", "FAMC_04"],
}
# iid -> (family, father, mother, declared_sex, affection)
PED_ROWS = {
    "FAMA_01": ("FAM_A", "0", "0", 1, 1),
    "FAMA_02": ("FAM_A", "0", "0", 2, 1),
    "FAMA_03": ("FAM_A", "FAMA_01", "FAMA_02", 1, 2),
    "FAMA_04": ("FAM_A", "FAMA_01", "FAMA_02", 2, 2),
    "FAMB_01": ("FAM_B", "0", "0", 1, 1),
    "FAMB_02": ("FAM_B", "0", "0", 2, 1),
    "FAMB_03": ("FAM_B", "FAMB_01", "FAMB_02", 1, 2),
    "FAMB_04": ("FAM_B", "FAMB_01", "FAMB_02", 2, 2),  # declared female: wrong
    "FAMB_05": ("FAM_B", "FAMB_01", "FAMB_02", 2, 2),  # declared female: wrong
    "FAMC_01": ("FAM_C", "0", "0", 1, 1),
    "FAMC_02": ("FAM_C", "0", "0", 2, 1),
    "FAMC_03": ("FAM_C", "FAMC_01", "FAMC_02", 1, 2),
    "FAMC_04": ("FAM_C", "FAMC_01", "FAMC_02", 1, 2),
}
# True sex of the genotypes, which is what the simulator builds.
TRUE_SEX = {
    "FAMA_01": "M", "FAMA_02": "F", "FAMA_03": "M", "FAMA_04": "F",
    "FAMB_01": "M", "FAMB_02": "F", "FAMB_03": "M", "FAMB_04": "M", "FAMB_05": "M",
    "FAMC_01": "M", "FAMC_02": "F", "FAMC_03": "M", "FAMC_04": "M",
}
PRIMARY = [s for s in TRUE_SEX if s not in ("FAMB_04", "FAMB_05")]

# chrX: (n_sites, n_het).  FAMC_03 and FAMC_04 carry the two profiles that end
# up, via the two swaps, in female-declared slots.
X_PROFILE = {
    "FAMC_03": (500, 188),   # 37.6%
    "FAMC_04": (505, 150),   # 29.7%
}
X_FEMALE = (500, 310)        # 62.0%
X_MALE = (500, 175)          # 35.0%
Y_COUNT = {"M": 6, "F": 2}
Y_OVERRIDE = {"FAMC_03": 5, "FAMC_04": 7}

B05_RETENTION = 0.7508
B05_DISCORDANCE = 0.0274


def make_positions(n: int, chroms: list[str], rng: random.Random) -> list[tuple[str, int]]:
    """Unique positions spread over the real span of each contig.

    Spread matters: a fixture whose chrX positions all sit in the first 3 Mb
    would land inside PAR1 and be excluded from the sex call, which would make
    the fixture test the exclusion logic instead of the sex logic.
    """
    lengths = dict(CONTIGS)
    out: list[tuple[str, int]] = []
    per = n // len(chroms) + 1
    for chrom in chroms:
        span = lengths[chrom]
        for pos in rng.sample(range(1_000_000, span - 100_000), per):
            out.append((chrom, pos))
    rng.shuffle(out)
    return out[:n]


def region_for(rng: random.Random) -> str:
    x = rng.random()
    acc = 0.0
    for name, w in REGION_MIX:
        acc += w
        if x <= acc:
            return name
    return "exonic"


def build(outdir: Path | None = None) -> dict:
    global COHORT
    COHORT = Path(outdir) if outdir else COHORT
    rng = random.Random(SEED)
    COHORT.mkdir(parents=True, exist_ok=True)

    # ---- allocate site pools -------------------------------------------
    priv_size = {
        s: (N_PRIVATE_C04 if s == "FAMC_04" else rng.randint(*PRIVATE_RANGE))
        for s in PRIMARY
    }
    n_auto = N_COMMON + N_FAMILY_EXTRA * 3 + sum(priv_size.values())
    auto_pos = make_positions(n_auto, AUTOSOMES, rng)
    x_pos = make_positions(1200, ["X"], rng)
    y_pos = make_positions(120, ["Y"], rng)

    cursor = 0

    def take(k: int) -> list[int]:
        nonlocal cursor
        ids = list(range(cursor, cursor + k))
        cursor += k
        return ids

    common = take(N_COMMON)
    fam_extra = {fam: take(N_FAMILY_EXTRA) for fam in FAMILIES}
    private = {s: take(priv_size[s]) for s in PRIMARY}

    all_pos = auto_pos + x_pos + y_pos
    x_offset, y_offset = len(auto_pos), len(auto_pos) + len(x_pos)
    regions = [region_for(rng) for _ in all_pos]
    alleles = [
        tuple(rng.sample("ACGT", 2)) for _ in all_pos
    ]  # (ref, alt)

    # ---- build primary samples -----------------------------------------
    calls: dict[str, dict[int, int]] = {}
    for s in PRIMARY:
        fam = PED_ROWS[s][0]
        sites = common + fam_extra[fam] + private[s]
        calls[s] = {i: (HET if rng.random() < P_HET else HOMALT) for i in sites}

        n_x, n_het_x = X_PROFILE.get(s, X_FEMALE if TRUE_SEX[s] == "F" else X_MALE)
        chosen_x = rng.sample(range(len(x_pos)), n_x)
        for rank, xi in enumerate(chosen_x):
            calls[s][x_offset + xi] = HET if rank < n_het_x else HOMALT

        n_y = Y_OVERRIDE.get(s, Y_COUNT[TRUE_SEX[s]])
        for yi in rng.sample(range(len(y_pos)), n_y):
            calls[s][y_offset + yi] = HOMALT

    # ---- swap 1: FAMB_04 is a byte-level copy of FAMC_04 ----------------
    calls["FAMB_04"] = dict(calls["FAMC_04"])

    # ---- swap 2: FAMB_05 is FAMC_03 at lower depth ----------------------
    src = calls["FAMC_03"]
    auto_src = [i for i in src if i < x_offset]
    x_src_het = [i for i in src if x_offset <= i < y_offset and src[i] == HET]
    x_src_hom = [i for i in src if x_offset <= i < y_offset and src[i] == HOMALT]
    y_src = [i for i in src if i >= y_offset]

    keep_total = round(B05_RETENTION * len(src))
    keep_x_het = round(B05_RETENTION * len(x_src_het))
    keep_x_hom = round(B05_RETENTION * len(x_src_hom))
    keep_y = round(B05_RETENTION * len(y_src))
    keep_auto = keep_total - keep_x_het - keep_x_hom - keep_y

    kept = (
        rng.sample(auto_src, keep_auto)
        + rng.sample(x_src_het, keep_x_het)
        + rng.sample(x_src_hom, keep_x_hom)
        + rng.sample(y_src, keep_y)
    )
    b05 = {i: src[i] for i in kept}
    # Genotype discordance, autosomes only, so the chrX profile is preserved.
    auto_kept = [i for i in kept if i < x_offset]
    for i in rng.sample(auto_kept, round(B05_DISCORDANCE * len(kept))):
        b05[i] = HOMALT if b05[i] == HET else HET
    calls["FAMB_05"] = b05

    # ---- write ----------------------------------------------------------
    order = sorted(range(len(all_pos)), key=lambda i: (CONTIG_ORDER[all_pos[i][0]], all_pos[i][1]))
    for s in TRUE_SEX:
        # FAMB_04's file must be FAMC_04's file with exactly one byte changed,
        # so it is rendered from the same records with the same depth stream and
        # only the sample name in the #CHROM line swapped.
        seed_name = "FAMC_04" if s == "FAMB_04" else s
        body = _render(s, seed_name, calls[s], order, all_pos, alleles, regions)
        (COHORT / f"{s}.vcf").write_text(body)

    _write_ped(COHORT / "cohort.ped")
    return _verify(calls, x_offset, y_offset)


def _render(sample, seed_name, calls, order, all_pos, alleles, regions) -> str:
    lines = ["##fileformat=VCFv4.2", "##reference=GRCh37"]
    lines += [f"##contig=<ID={c},length={n}>" for c, n in CONTIGS]
    lines += [
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">',
        '##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Max-likelihood allele count">',
        '##INFO=<ID=Func.refGene,Number=.,Type=String,Description="Region">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}",
    ]
    local = random.Random(zlib.crc32(seed_name.encode()))  # stable across runs
    for i in order:
        code = calls.get(i)
        if code is None:
            continue
        chrom, pos = all_pos[i]
        ref, alt = alleles[i]
        dp = local.randint(20, 60)
        if code == HET:
            gt, ac, ad = "0/1", 1, f"{dp // 2},{dp - dp // 2}"
        else:
            gt, ac, ad = "1/1", 2, f"0,{dp}"
        info = f"AC={ac};MLEAC={ac};DP={dp};Func.refGene={regions[i]}"
        lines.append(
            f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t500.0\tPASS\t{info}"
            f"\tGT:AD:DP:GQ\t{gt}:{ad}:{dp}:99"
        )
    return "\n".join(lines) + "\n"


def _write_ped(path: Path) -> None:
    rows = ["#FID\tIID\tPAT\tMAT\tSEX\tPHENO"]
    for iid, (fam, pat, mat, sex, aff) in PED_ROWS.items():
        rows.append(f"{fam}\t{iid}\t{pat}\t{mat}\t{sex}\t{aff}")
    path.write_text("\n".join(rows) + "\n")


def _verify(calls, x_offset, y_offset) -> dict:
    def jac(a, b):
        sa, sb = set(calls[a]), set(calls[b])
        return len(sa & sb) / len(sa | sb)

    def agree(a, b):
        shared = set(calls[a]) & set(calls[b])
        return sum(calls[a][i] == calls[b][i] for i in shared) / len(shared)

    def xhet(s):
        xs = [i for i in calls[s] if x_offset <= i < y_offset]
        return sum(calls[s][i] == HET for i in xs) / len(xs), len(xs)

    unrelated, sibs = [], []
    for a in TRUE_SEX:
        for b in TRUE_SEX:
            if a >= b or a in ("FAMB_04", "FAMB_05") or b in ("FAMB_04", "FAMB_05"):
                continue
            (sibs if PED_ROWS[a][0] == PED_ROWS[b][0] else unrelated).append(jac(a, b))

    size_b04 = (COHORT / "FAMB_04.vcf").stat().st_size
    size_c04 = (COHORT / "FAMC_04.vcf").stat().st_size
    ba = (COHORT / "FAMB_04.vcf").read_bytes()
    bc = (COHORT / "FAMC_04.vcf").read_bytes()
    byte_diff = sum(1 for x, y in zip(ba, bc, strict=False) if x != y) + abs(len(ba) - len(bc))

    return {
        "unrelated_jaccard_range": (round(min(unrelated), 4), round(max(unrelated), 4)),
        "sib_jaccard_range": (round(min(sibs), 4), round(max(sibs), 4)),
        "FAMB_04_vs_FAMC_04": (round(jac("FAMB_04", "FAMC_04"), 4),
                               round(agree("FAMB_04", "FAMC_04"), 4),
                               len(set(calls["FAMB_04"]) & set(calls["FAMC_04"]))),
        "FAMB_05_vs_FAMC_03": (round(jac("FAMB_05", "FAMC_03"), 4),
                               round(agree("FAMB_05", "FAMC_03"), 4)),
        "chrX_het_FAMB_04": round(xhet("FAMB_04")[0], 4),
        "chrX_het_FAMB_05": round(xhet("FAMB_05")[0], 4),
        "variant_counts": {s: len(calls[s]) for s in TRUE_SEX},
        "file_bytes_FAMB_04": size_b04,
        "file_bytes_FAMC_04": size_c04,
        "byte_differences": byte_diff,
    }


if __name__ == "__main__":
    import json

    print(json.dumps(build(), indent=2, default=str))
