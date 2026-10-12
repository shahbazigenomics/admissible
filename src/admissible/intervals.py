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
import zlib
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
    except (OSError, EOFError, gzip.BadGzipFile, zlib.error) as exc:
        # zlib.error is what a gzip with a valid header and a corrupt body raises
        # (it is not an OSError), so a bit-flipped download used to escape here.
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


def subtract(a: IntervalSet, b: IntervalSet) -> IntervalSet:
    """Positions in ``a`` that are not in ``b``.

    Used to take an exclusion set (low-mappability, repeats) out of a callable
    set, and to list the uncallable part of a region.  Subtraction distributes
    over both intersection and k-of-n coverage, so subtracting once from a joint
    set equals subtracting from every sample first.

    Neither argument is modified, and neither needs to be merged: both are
    normalised here, so a caller cannot get a wrong answer from unsorted or
    overlapping input.
    """
    out: IntervalSet = {}
    for chrom, av_raw in a.items():
        av = merge(list(av_raw))
        bv = merge(list(b.get(chrom, [])))
        if not bv:
            if av:
                out[chrom] = av
            continue
        res: list[Interval] = []
        j = 0
        for start, end in av:
            # b-intervals that end before this a-interval starts can never matter
            # again, because a is sorted.
            while j < len(bv) and bv[j][1] <= start:
                j += 1
            cur = start
            k = j
            while k < len(bv) and bv[k][0] < end:
                b_start, b_end = bv[k]
                if b_start > cur:
                    res.append((cur, b_start))
                if b_end > cur:
                    cur = b_end
                if cur >= end:
                    break
                k += 1
            if cur < end:
                res.append((cur, end))
        if res:
            out[chrom] = res
    return out


@dataclass
class NamedRegion:
    """One named target region: every BED line that shares a column-4 name.

    A name can sit on several lines (one exon listed under more than one
    transcript, or a region split across lines), so a region is a *set* of
    intervals, not a single one.  Two different names that overlap stay two
    regions: they are different questions, and merging them would hide the worse
    of the two.
    """

    name: str
    intervals: IntervalSet = field(default_factory=dict)
    # False when the name was made up from the coordinates because column 4 was
    # absent or empty on every line of this region.
    named: bool = True

    @property
    def bp(self) -> int:
        return total_bp(self.intervals)


@dataclass
class NamedBed:
    path: str
    regions: list[NamedRegion] = field(default_factory=list)
    n_lines: int = 0
    n_malformed: int = 0
    n_unnamed: int = 0
    problems: list[str] = field(default_factory=list)

    @property
    def has_names(self) -> bool:
        """True when every region carries a real column-4 name."""
        return bool(self.regions) and all(r.named for r in self.regions)

    @property
    def total_bp(self) -> int:
        """Bases in the union of all regions, so overlapping regions count once."""
        union: IntervalSet = {}
        for r in self.regions:
            for chrom, ivs in r.intervals.items():
                union.setdefault(chrom, []).extend(ivs)
        return total_bp({c: merge(v) for c, v in union.items()})


def read_named_bed(path: str | os.PathLike[str]) -> NamedBed:
    """Parse a BED file into named regions.  Never raises; problems are recorded.

    Regions keep the order in which their name first appears.  A line with no
    column-4 name still becomes a region, named ``contig:start-end`` and counted
    in ``n_unnamed``, so a caller that wants per-region output from an unnamed
    BED gets one row per line rather than nothing.  Negative starts are treated as
    malformed here (``read_bed`` keeps them, and is left as it was).
    """
    bed = NamedBed(path=str(path))
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    try:
        fh = opener(path, "rt", errors="replace")  # type: ignore[operator]
    except OSError as exc:
        bed.problems.append(f"cannot open: {exc}")
        return bed

    raw: dict[str, tuple[bool, IntervalSet]] = {}
    n_repeat_lines = 0
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
                if start < 0 or end <= start:
                    bed.n_malformed += 1
                    continue
                contig = normalize_contig(f[0])
                name = f[3].strip() if len(f) >= 4 else ""
                named = bool(name)
                if not named:
                    name = f"{contig}:{start}-{end}"
                    bed.n_unnamed += 1
                if name in raw:
                    n_repeat_lines += 1
                else:
                    raw[name] = (named, {})
                raw[name][1].setdefault(contig, []).append((start, end))
    except (OSError, EOFError, gzip.BadGzipFile, zlib.error) as exc:
        bed.problems.append(f"read stopped early: {type(exc).__name__}: {exc}")

    for name, (named, ivs) in raw.items():
        bed.regions.append(
            NamedRegion(
                name=name,
                intervals={c: merge(v) for c, v in ivs.items()},
                named=named,
            )
        )
    if bed.n_malformed:
        bed.problems.append(f"{bed.n_malformed} unparseable line(s) skipped")
    if n_repeat_lines:
        bed.problems.append(
            f"{n_repeat_lines} line(s) repeat an earlier name; their intervals were "
            f"combined into that region"
        )
    if bed.n_unnamed:
        bed.problems.append(
            f"{bed.n_unnamed} line(s) have no name in column 4 and were named by "
            f"their coordinates"
        )
    return bed
