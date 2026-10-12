#!/usr/bin/env python3
"""Compare admissible's callability JSON with the bedtools baseline for ONE sample.

    compare.py BASELINE_DIR ADMISSIBLE.json [--with-exclusion] [--out RESULT.tsv]

Exit 0: every compared number is identical.  1: at least one differs.  2: the inputs
could not be read or do not have the expected shape (reported, never raised).

Everything compared is an integer base count, so the expectation is exact equality,
not a tolerance: a difference of one base is a real disagreement about half-open
interval arithmetic and should be looked at, not rounded away.

Standard library only.  Not run against real data by its author; the comparison
logic itself is unit-tested on constructed inputs in tests/test_validation_compare.py.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def read_tsv(path: Path) -> list[list[str]]:
    return [ln.rstrip("\n").split("\t") for ln in path.read_text().splitlines() if ln.strip()]


def load_baseline(base: Path) -> tuple[dict[str, int], dict[str, tuple[int, int, int]]]:
    metrics = {k: int(v) for k, v in read_tsv(base / "baseline.tsv")}
    regions: dict[str, tuple[int, int, int]] = {}
    path = base / "regions.tsv"
    if path.exists():
        for name, bp, cbp, xbp in read_tsv(path):
            regions[name] = (int(bp), int(cbp), int(xbp))
    return metrics, regions


def compare(
    metrics: dict[str, int],
    regions: dict[str, tuple[int, int, int]],
    report: dict,
    with_exclusion: bool,
) -> tuple[list[tuple[str, str, object, object]], list[str]]:
    """Return (rows, problems). A row is (what, verdict, baseline, admissible)."""
    rows: list[tuple[str, str, object, object]] = []
    problems: list[str] = []

    def add(what: str, base: object, got: object) -> None:
        rows.append((what, "same" if base == got else "DIFFERENT", base, got))

    checks = report.get("checks") if isinstance(report, dict) else None
    if not isinstance(checks, list):
        return rows, ["the JSON is not an admissible report (no list of checks)"]
    cal = next(
        (c for c in checks if isinstance(c, dict) and c.get("check") == "callability"), None
    )
    if cal is None:
        return rows, ["no callability check in the JSON"]
    m = cal.get("metrics") or {}
    if "models" not in m or "target_bp" not in m:
        return rows, [f"callability did not compute fractions (status {cal.get('status')}): "
                      f"{cal.get('summary')}"]

    models = m["models"]
    model = models.get("recessive")
    if not model or model.get("callable_bp") is None:
        return rows, ["the recessive model was not computed; cannot compare"]

    add("target_bp", metrics["target_bp"], m["target_bp"])
    if with_exclusion:
        raw = (model.get("without_exclusion") or {}).get("callable_bp")
        add("callable_in_target_bp (no exclusion)", metrics["callable_in_target_bp"], raw)
        add("callable_in_target_minus_exclude_bp",
            metrics["callable_in_target_minus_exclude_bp"], model["callable_bp"])
        key = 2
    else:
        add("callable_in_target_bp", metrics["callable_in_target_bp"], model["callable_bp"])
        add("joint_callable_bp", metrics["callable_in_target_bp"], m.get("joint_callable_bp"))
        key = 1

    br = m.get("by_region")
    if not br:
        problems.append("no by_region block; per-region numbers not compared")
        return rows, problems
    units = br.get("units") or {}
    if len(units) != 1:
        problems.append(f"expected exactly one unit for a single sample, found {len(units)}")
        return rows, problems
    unit = next(iter(units.values()))
    listed = {r["region"]: r for r in unit.get("regions", [])}
    for name, tup in regions.items():
        row = listed.get(name)
        if row is None:
            if not unit.get("truncated"):
                rows.append((f"region {name}", "MISSING", tup[0], None))
            continue
        add(f"region {name} bp", tup[0], row["bp"])
        add(f"region {name} callable_bp", tup[key], row["callable_bp"])
    extra = sorted(set(listed) - set(regions))
    for name in extra:
        rows.append((f"region {name}", "ONLY_IN_ADMISSIBLE", None, listed[name]["bp"]))
    if unit.get("truncated"):
        problems.append(
            f"admissible listed {len(unit.get('regions', []))} of {unit.get('n_regions')} "
            f"regions; the rest were not compared"
        )
    return rows, problems


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("baseline_dir", type=Path)
    ap.add_argument("report", type=Path)
    ap.add_argument("--with-exclusion", action="store_true")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args(argv)
    try:
        metrics, regions = load_baseline(args.baseline_dir)
        report = json.loads(args.report.read_text())
    except (OSError, ValueError, KeyError) as exc:
        print(f"compare.py: cannot read inputs: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2
    try:
        rows, problems = compare(metrics, regions, report, args.with_exclusion)
    except (KeyError, TypeError, ValueError, AttributeError, IndexError) as exc:
        print(f"compare.py: unexpected report shape: {type(exc).__name__}: {exc}",
              file=sys.stderr)
        return 2
    lines = ["what\tverdict\tbaseline\tadmissible"] + [
        "\t".join(str(x) for x in r) for r in rows
    ]
    text = "\n".join(lines) + "\n"
    if args.out:
        args.out.write_text(text)
    bad = [r for r in rows if r[1] != "same"]
    print(f"{len(rows) - len(bad)} of {len(rows)} compared numbers identical")
    for r in bad[:20]:
        print(f"  {r[0]}: {r[1]} (baseline {r[2]}, admissible {r[3]})")
    for p in problems:
        print(f"  note: {p}")
    if not rows:
        return 2
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
