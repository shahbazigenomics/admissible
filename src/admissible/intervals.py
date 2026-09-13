"""Interval arithmetic for callability.

The only interesting operation here is ``coverage_at_least``: the set of
positions callable in at least *k* of *n* samples.  It exists because the joint
callable fraction is an intersection, never a product of per-sample fractions -
coverage failures are strongly correlated across samples (the same GC-rich first
exons, segmental duplications and low-efficiency probes fail in everyone), so
multiplying understates the truth, often by a factor of two or more.

It is parameterised by *k* rather than hard-wired to *n* because callability is a
property of a genetic model, not of a cohort: a fully penetrant model needs every
affected sample callable, while a phenocopy(n-1) model needs all but one, and one
global number cannot serve both.

Intervals are half-open ``[start, end)``, matching BED.
"""

from __future__ import annotations

import gzip
import os
from dataclasses import dataclass, field

from .contigs import normalize_contig

Interval = tuple[int, int]
IntervalSet = dict[str, list[Interval]]


@dataclass
class BedFile:
    path: str
    intervals: IntervalSet = field(default_factory=dict)
    n_lines: int = 0
    n_malformed: int = 0
    problems: list[str] = field(default_factory=list)
    # Column-4 values and the base pairs each covers.  Present so that a caller
    # can see which depth bins a mosdepth --quantize file actually contains,
    # rather than silently getting nothing when it asked for a bin that is not
    # there.
    labels: dict[str, int] = field(default_factory=dict)

    @property
    def total_bp(self) -> int:
        return total_bp(self.intervals)


def read_bed(
    path: str | os.PathLike[str], keep_label: str | None = None
) -> BedFile:
    """Parse a BED file.  Never raises; problems are recorded.

    ``keep_label`` filters on column 4, which is how a mosdepth ``--quantize``
    output is narrowed to one depth bin (e.g. ``"10:inf"``).
    """
    bed = BedFile(path=str(path))
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    try:
        fh = opener(path, "rt", errors="replace")  # type: ignore[operator]
    except OSError as exc:
        bed.problems.append(f"cannot open: {exc}")
        return bed
    raw: IntervalSet = {}
    try:
        with fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                bed.n_lines += 1
                f = line.split("\t")
                if len(f) < 3:
                    bed.n_malformed += 1
                    continue
                try:
                    start, end = int(f[1]), int(f[2])
                except ValueError:
                    bed.n_malformed += 1
                    continue
                if end <= start:
                    bed.n_malformed += 1
                    continue
                if len(f) >= 4 and f[3].strip():
                    lab = f[3].strip()
                    bed.labels[lab] = bed.labels.get(lab, 0) + (end - start)
                if keep_label is not None and (len(f) < 4 or f[3].strip() != keep_label):
                    continue
                raw.setdefault(normalize_contig(f[0]), []).append((start, end))
    except (OSError, EOFError, gzip.BadGzipFile) as exc:
        bed.problems.append(f"read stopped early: {type(exc).__name__}: {exc}")

    bed.intervals = {c: merge(v) for c, v in raw.items()}
    if bed.n_malformed:
        bed.problems.append(f"{bed.n_malformed} unparseable line(s) skipped")
    return bed


def merge(intervals: list[Interval]) -> list[Interval]:
    """Sort and coalesce touching or overlapping intervals."""
    if not intervals:
        return []
    out: list[Interval] = []
    for start, end in sorted(intervals):
        if out and start <= out[-1][1]:
            if end > out[-1][1]:
                out[-1] = (out[-1][0], end)
        else:
            out.append((start, end))
    return out


def total_bp(iset: IntervalSet) -> int:
    return sum(e - s for v in iset.values() for s, e in v)


def intersect(a: IntervalSet, b: IntervalSet) -> IntervalSet:
    """Positions present in both sets."""
    out: IntervalSet = {}
    for chrom, av in a.items():
        bv = b.get(chrom)
        if not bv:
            continue
        merged: list[Interval] = []
        i = j = 0
        while i < len(av) and j < len(bv):
            lo = max(av[i][0], bv[j][0])
            hi = min(av[i][1], bv[j][1])
            if lo < hi:
                merged.append((lo, hi))
            if av[i][1] < bv[j][1]:
                i += 1
            else:
                j += 1
        if merged:
            out[chrom] = merged
    return out


def coverage_at_least(sets: list[IntervalSet], k: int) -> IntervalSet:
    """Positions covered by at least ``k`` of the supplied sample interval sets.

    A sweep line over start/end events.  ``k == len(sets)`` gives the strict
    intersection; ``k == len(sets) - 1`` gives the phenocopy(n-1) denominator.
    """
    if k <= 0 or not sets:
        return {}
    k = min(k, len(sets))
    chroms = {c for s in sets for c in s}
    out: IntervalSet = {}
    for chrom in chroms:
        events: list[tuple[int, int]] = []
        for s in sets:
            for start, end in s.get(chrom, []):
                events.append((start, 1))
                events.append((end, -1))
        if not events:
            continue
        events.sort()
        spans: list[Interval] = []
        depth = 0
        run_start: int | None = None
        i = 0
        while i < len(events):
            pos = events[i][0]
            while i < len(events) and events[i][0] == pos:
                depth += events[i][1]
                i += 1
            if depth >= k and run_start is None:
                run_start = pos
            elif depth < k and run_start is not None:
                if pos > run_start:
                    spans.append((run_start, pos))
                run_start = None
        if spans:
            out[chrom] = merge(spans)
    return out
