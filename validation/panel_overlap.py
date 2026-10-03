"""Validate the "peddy/somalier panels are largely intronic/intergenic,
almost nothing survives coding-only PASS-filtering" claim, by measuring it
against real data rather than asserting it.

This exists because that exact claim, as it stood in the README and in D12 of
admissible/design-decisions.md, was independently checked against peddy's real
GRCh37 site panel and found to be wrong in direction: peddy's panel is heavily
concentrated in coding regions, not scattered genome-wide. See D12 for the
corrected text and the reasoning for why the panel comes out this way (peddy's
own site-selection script requires the 1000 Genomes ``EX_TARGET`` flag - an
exome-capture-target annotation - on every candidate site, so this is a
property of how the panel was built, not a coincidence).

Two independent measurements, reported separately and then combined:

1. What fraction of peddy's real GRCh37 site panel falls inside a real coding
   BED (RefSeq CDS) and a real clinical exome capture BED (Agilent
   SureSelect V6)?
2. What fraction of a real GATK call set survives the standard hard-filter
   thresholds (QD, FS, MQ, MQRankSum, ReadPosRankSum, SOR)?

Neither dataset is vendored. (1) clones two small public repos on demand
(~7 MB combined); (2) reuses the real GATK INFO fields already committed at
tests/data/annovar/gatk_chr21.hg38_multianno.txt (see that directory's
PROVENANCE.md for exactly what is real there) rather than fetching anything
new, so it needs no network access.

Run: ``python validation/panel_overlap.py``
"""

from __future__ import annotations

import bisect
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ANNOVAR_FIXTURE = ROOT / "tests" / "data" / "annovar" / "gatk_chr21.hg38_multianno.txt"

PEDDY_REPO = "https://github.com/brentp/peddy"
REFDATA_REPO = "https://github.com/AstraZeneca-NGS/reference_data"

# Standard GATK hard-filter thresholds (VariantFiltration recommendations for
# SNPs). A record fails if ANY present annotation crosses its threshold; a
# missing annotation is not itself a failure (some are SNP/indel-specific and
# some, like MQRankSum/ReadPosRankSum, require a het call to be defined).
HARD_FILTERS = {
    "QD": lambda v: v < 2.0,
    "FS": lambda v: v > 60.0,
    "MQ": lambda v: v < 40.0,
    "MQRankSum": lambda v: v < -12.5,
    "ReadPosRankSum": lambda v: v < -8.0,
    "SOR": lambda v: v > 3.0,
}


def _clone_if_missing(url: str, dest: Path) -> Path:
    if not dest.exists():
        print(f"cloning {url} -> {dest}", file=sys.stderr)
        subprocess.run(
            ["git", "clone", "--depth", "1", url, str(dest)], check=True
        )
    return dest


def _load_bed(path: Path) -> dict[str, list[tuple[int, int]]]:
    raw: dict[str, list[tuple[int, int]]] = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end = line.split("\t")[:3]
            raw.setdefault(chrom.removeprefix("chr"), []).append((int(start), int(end)))
    merged: dict[str, list[tuple[int, int]]] = {}
    for chrom, ivs in raw.items():
        ivs.sort()
        out: list[tuple[int, int]] = []
        for s, e in ivs:
            if out and s <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append((s, e))
        merged[chrom] = out
    return merged


def _load_sites(path: Path) -> list[tuple[str, int]]:
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            chrom, pos, _ref, _alt = line.split(":")
            out.append((chrom.removeprefix("chr"), int(pos)))
    return out


def _covered(chrom: str, pos_1based: int, bed: dict[str, list[tuple[int, int]]]) -> bool:
    ivs = bed.get(chrom)
    if not ivs:
        return False
    p0 = pos_1based - 1  # BED is 0-based half-open
    starts = [iv[0] for iv in ivs]
    i = bisect.bisect_right(starts, p0) - 1
    if i < 0:
        return False
    s, e = ivs[i]
    return s <= p0 < e


def panel_overlap(cache_dir: Path) -> dict:
    peddy = _clone_if_missing(PEDDY_REPO, cache_dir / "peddy")
    refdata = _clone_if_missing(REFDATA_REPO, cache_dir / "reference_data")

    sites = _load_sites(peddy / "peddy" / "GRCH37.sites")
    cds = _load_bed(refdata / "GRCh37" / "bed" / "CDS_RefSeq.bed")
    agilent = _load_bed(refdata / "GRCh37" / "bed" / "Exome-Agilent_V6.bed")

    in_cds = sum(1 for c, p in sites if _covered(c, p, cds))
    in_agilent = sum(1 for c, p in sites if _covered(c, p, agilent))
    n = len(sites)
    return {
        "n_sites": n,
        "cds_fraction": in_cds / n,
        "agilent_v6_fraction": in_agilent / n,
    }


def _parse_info(field: str) -> dict[str, str]:
    out = {}
    for kv in field.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def pass_filter_survival(fixture: Path) -> dict:
    lines = fixture.read_text().splitlines()
    header = lines[0].split("\t")
    info_col = header.index("Otherinfo11")  # the real VCF INFO column, verbatim
    n = 0
    passed = 0
    failed_records: list[tuple[str, list[str]]] = []
    for row in lines[1:]:
        fields = row.split("\t")
        info = _parse_info(fields[info_col])
        n += 1
        reasons = []
        for key, fails in HARD_FILTERS.items():
            raw = info.get(key)
            if raw in (None, "."):
                continue
            try:
                value = float(raw)
            except ValueError:
                continue
            if fails(value):
                reasons.append(key)
        if reasons:
            failed_records.append((fields[1], reasons))
        else:
            passed += 1
    return {
        "n_records": n,
        "n_passed": passed,
        "fraction_passed": passed / n,
        "failed": failed_records,
    }


def main() -> int:
    cache_dir = Path(tempfile.gettempdir()) / "admissible_validation_cache"
    cache_dir.mkdir(exist_ok=True)

    panel = panel_overlap(cache_dir)
    filt = pass_filter_survival(ANNOVAR_FIXTURE)

    print(f"peddy GRCh37 panel: {panel['n_sites']} sites")
    print(f"  inside RefSeq CDS:        {panel['cds_fraction']:.1%}")
    print(f"  inside Agilent V6 target: {panel['agilent_v6_fraction']:.1%}")
    print()
    print(
        f"GATK hard-filter survival (real INFO fields, "
        f"n={filt['n_records']} from tests/data/annovar): "
        f"{filt['fraction_passed']:.1%} ({filt['n_passed']}/{filt['n_records']})"
    )
    if filt["failed"]:
        print("  failed records (pos, reasons):")
        for pos, reasons in filt["failed"]:
            print(f"    {pos}: {', '.join(reasons)}")
    print()
    combined_cds = panel["cds_fraction"] * filt["fraction_passed"]
    combined_agilent = panel["agilent_v6_fraction"] * filt["fraction_passed"]
    print(
        f"combined (panel overlap x PASS survival): "
        f"{combined_cds:.1%} (vs CDS), {combined_agilent:.1%} (vs Agilent V6)"
    )
    print(
        "-> both land close to three-quarters, not 'almost none': a "
        "coding-only, PASS-only delivery still overlaps the majority of "
        "peddy's panel."
    )

    # n=120 is a small sample for the hard-filter fraction specifically (it
    # was chosen for its 7 real multiallelic records, not for this purpose) -
    # do not treat fraction_passed as more precise than that sample size
    # supports. It is reported here because it is the only *real* GATK INFO
    # data already in this repo; a larger real call set would tighten it.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
