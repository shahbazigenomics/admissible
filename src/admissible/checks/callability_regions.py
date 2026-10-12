"""Per-region callability: one row per named target region, worst region per gene.

One number per run averages a dead exon away: a 10-exon gene with one exon wholly
uncovered is 90% callable, and "90% searched" is the wrong thing to say about it.
This module answers the question a single-gene analysis actually asks - was *this
exon* searched - for a target BED whose column 4 names its regions.

It is opt-in (``--by-region``) rather than detected, because capture-kit BEDs often
carry column-4 names already and switching the output on from their presence would
change what existing users get.

Two states are kept apart on purpose.  A region the samples cover with no callable
base is *uncallable*: it was looked at and could not be searched.  A region on a
contig the coverage file never mentions is *NOT_REPORTED*: nothing says whether it
was searched, so it gets no fraction, not a zero.  Reporting it as zero would claim
a measurement nobody made; reporting it as covered would claim the opposite.

Nothing here raises: a bad gene pattern, an unreadable BED or an empty target comes
back as a note and an empty table.
"""

from __future__ import annotations

import re
from collections.abc import Callable
from typing import Any

from ..intervals import (
    IntervalIndex,
    IntervalSet,
    NamedBed,
    NamedRegion,
    total_bp,
)

NOT_REPORTED = "NOT_REPORTED"

# Worst to best.  UNKNOWN is below EXPLORATORY: nothing at all can be said.
_ORDER = ["UNKNOWN", "EXPLORATORY", "SEARCHED_PARTLY", "SEARCHED_MOST"]


def worse_rung(a: str, b: str) -> str:
    """The lower of two rungs.  A rung this module does not know sorts as UNKNOWN."""
    ia = _ORDER.index(a) if a in _ORDER else 0
    ib = _ORDER.index(b) if b in _ORDER else 0
    return _ORDER[min(ia, ib)]


def compile_gene_pattern(pattern: str) -> tuple[re.Pattern[str] | None, str | None]:
    """Compile the configured gene pattern.  A bad one is a note, not an exception."""
    try:
        return re.compile(pattern, re.IGNORECASE), None
    except (re.error, TypeError) as exc:
        return None, (
            f"gene pattern {pattern!r} is not a valid regular expression ({exc}); "
            f"every region name was used as its own gene"
        )


def gene_of(name: str, pattern: re.Pattern[str] | None) -> str:
    """Region name with the trailing exon marker removed once; the name if nothing
    would be left."""
    if pattern is None:
        return name
    stripped = pattern.sub("", name, count=1)
    return stripped if stripped else name


def _interval_strings(iset: IntervalSet, limit: int) -> tuple[list[str], int]:
    flat = [(c, s, e) for c in sorted(iset) for s, e in iset[c]]
    return [f"{c}:{s}-{e}" for c, s, e in flat[:limit]], len(flat)


def _region_row(
    region: NamedRegion,
    gene: str,
    idx: IntervalIndex,
    reported_in: Callable[[str], bool],
    claim_rung_for: Callable[[float], str],
    claim_for: Callable[[float], str],
) -> dict[str, Any]:
    bp = region.bp
    contigs = sorted(region.intervals)
    absent = [c for c in contigs if not reported_in(c)]
    row: dict[str, Any] = {
        "region": region.name,
        "gene": gene,
        "named": region.named,
        "contigs": contigs,
        "bp": bp,
        "callable_bp": None,
        "fraction": None,
        "status": "ok",
        "fully_uncallable": False,
        "claim_rung": "UNKNOWN",
        "permitted_claim": None,
    }
    if absent:
        row["status"] = NOT_REPORTED
        row["not_reported_contigs"] = absent
        row["permitted_claim"] = (
            "cannot be assessed: the coverage file does not report "
            f"{', '.join(absent)}, so this region is neither searched nor shown "
            "to be unsearchable"
        )
        return row
    callable_bp = sum(
        idx.overlap_bp(c, s, e) for c, v in region.intervals.items() for s, e in v
    )
    row["callable_bp"] = callable_bp
    if bp <= 0:
        row["status"] = "empty"
        return row
    frac = callable_bp / bp
    row["fraction"] = frac
    row["fully_uncallable"] = callable_bp == 0
    row["claim_rung"] = claim_rung_for(frac)
    row["permitted_claim"] = claim_for(frac)
    return row


def analyse_regions(
    named: NamedBed,
    callable_set: IntervalSet,
    reported_in: Callable[[str], bool],
    *,
    gene_regex: str,
    claim_rung_for: Callable[[float], str],
    claim_for: Callable[[float], str],
    blind_cap: str,
    max_rows: int,
    max_listed: int,
) -> dict[str, Any]:
    """Region rows and per-gene roll-up for one callable set (one unit, one model).

    Rows are listed worst first and capped at ``max_rows``; the cap is reported, so
    a truncated list is never mistaken for a complete one.
    """
    pattern, pattern_note = compile_gene_pattern(gene_regex)
    idx = IntervalIndex(callable_set)

    rows: list[dict[str, Any]] = []
    for region in named.regions:
        rows.append(
            _region_row(
                region, gene_of(region.name, pattern), idx, reported_in,
                claim_rung_for, claim_for,
            )
        )

    by_gene: dict[str, list[dict[str, Any]]] = {}
    for r in rows:
        by_gene.setdefault(r["gene"], []).append(r)

    genes: list[dict[str, Any]] = []
    for gene, members in by_gene.items():
        ok = [m for m in members if m["status"] == "ok"]
        unreported = [m for m in members if m["status"] == NOT_REPORTED]
        blind = [m for m in ok if m["fully_uncallable"]]
        ok_bp = sum(m["bp"] for m in ok)
        ok_cbp = sum(m["callable_bp"] for m in ok)
        frac = (ok_cbp / ok_bp) if ok_bp else None
        worst = min(ok, key=lambda m: m["fraction"]) if ok else None
        if frac is None:
            rung, claim, capped = "UNKNOWN", None, False
            claim = (
                "cannot be assessed: no region of this gene is reported by the "
                "coverage file"
            )
        else:
            rung = claim_rung_for(frac)
            claim = claim_for(frac)
            capped = False
            if blind or unreported:
                new = worse_rung(rung, blind_cap)
                capped = new != rung
                rung = new
                bits = []
                if blind:
                    bits.append(
                        f"{len(blind)} region(s) had no callable base "
                        f"({', '.join(m['region'] for m in blind[:max_listed])})"
                    )
                if unreported:
                    bits.append(
                        f"{len(unreported)} region(s) are not reported by the "
                        f"coverage file "
                        f"({', '.join(m['region'] for m in unreported[:max_listed])})"
                    )
                claim = (
                    f"{claim}; the gene-level share hides that " + " and ".join(bits)
                    + (", so the claim is capped" if capped else "")
                )
        genes.append(
            {
                "gene": gene,
                "n_regions": len(members),
                "bp": sum(m["bp"] for m in members),
                "callable_bp": ok_cbp if ok else None,
                "fraction": frac,
                "worst_region": worst["region"] if worst else None,
                "worst_fraction": worst["fraction"] if worst else None,
                "n_fully_uncallable": len(blind),
                "fully_uncallable_regions": [m["region"] for m in blind[:max_listed]],
                "n_not_reported": len(unreported),
                "claim_rung": rung,
                "capped": capped,
                "permitted_claim": claim,
            }
        )

    def region_key(r: dict[str, Any]) -> tuple[int, float, str]:
        if r["status"] == NOT_REPORTED:
            return (0, 0.0, r["region"])
        return (1, r["fraction"] if r["fraction"] is not None else 2.0, r["region"])

    def gene_key(g: dict[str, Any]) -> tuple[int, float, str]:
        if g["fraction"] is None:
            return (0, 0.0, g["gene"])
        return (1, g["worst_fraction"], g["gene"])

    rows.sort(key=region_key)
    genes.sort(key=gene_key)
    n_rows, n_genes = len(rows), len(genes)
    kept_rows = rows[:max_rows]
    kept_genes = genes[:max_rows]

    # The uncallable part of a region is only worth computing for the rows we list.
    by_name = {r.name: r for r in named.regions}
    for r in kept_rows:
        frac = r["fraction"]
        if r["status"] == "ok" and frac is not None and frac < 1.0:
            gaps: IntervalSet = {}
            for c, ivs in by_name[r["region"]].intervals.items():
                found = [g for s, e in ivs for g in idx.gaps(c, s, e)]
                if found:
                    gaps[c] = found
            listed, n = _interval_strings(gaps, max_listed)
            r["uncallable_intervals"] = listed
            r["n_uncallable_intervals"] = n
            r["uncallable_bp"] = total_bp(gaps)

    notes = [pattern_note] if pattern_note else []
    return {
        "regions": kept_rows,
        "genes": kept_genes,
        "n_regions": n_rows,
        "n_genes": n_genes,
        "truncated": n_rows > len(kept_rows) or n_genes > len(kept_genes),
        "max_rows": max_rows,
        "n_fully_uncallable_regions": sum(1 for r in rows if r["fully_uncallable"]),
        "n_not_reported_regions": sum(1 for r in rows if r["status"] == NOT_REPORTED),
        "n_genes_capped": sum(1 for g in genes if g["capped"]),
        "n_genes_with_unsearched_region": sum(
            1 for g in genes if g["n_fully_uncallable"] or g["n_not_reported"]
        ),
        "notes": notes,
    }
