"""Contig naming, genome-build inference, and pseudoautosomal intervals.

Build coordinates are constants, not a bundled database, so this file carries no
licence weight.  They are the standard PAR definitions for GRCh37 and GRCh38.
"""

from __future__ import annotations

from dataclasses import dataclass

# Canonical chromosome length -> build, used to infer the build from the header.
_CHR_LENGTHS = {
    "1": {249250621: "GRCh37", 248956422: "GRCh38"},
    "X": {155270560: "GRCh37", 156040895: "GRCh38"},
    "Y": {59373566: "GRCh37", 57227415: "GRCh38"},
}

# (start, end) inclusive, 1-based.
PAR = {
    "GRCh37": {
        "X": [(60001, 2699520), (154931044, 155260560)],
        "Y": [(10001, 2649520), (59034050, 59363566)],
    },
    "GRCh38": {
        "X": [(10001, 2781479), (155701383, 156030895)],
        "Y": [(10001, 2781479), (56887903, 57217415)],
    },
}

# X-transposed region: 99% identical to Yq, a well-known source of apparent
# heterozygosity in males.  Excluded by default from the sex call.
XTR = {
    "GRCh37": {"X": [(88400000, 92000000)]},
    "GRCh38": {"X": [(89140000, 92750000)]},
}

AUTOSOMES = tuple(str(i) for i in range(1, 23))


def normalize_contig(name: str) -> str:
    """'chr7' -> '7', 'chrM'/'MT' -> 'MT', leaves anything unrecognised alone."""
    n = name.strip()
    if n.lower().startswith("chr"):
        n = n[3:]
    if n.upper() in ("M", "MT"):
        return "MT"
    if n.upper() in ("X", "Y"):
        return n.upper()
    return n


def is_primary(contig: str) -> bool:
    """True for 1-22, X, Y, MT.  Alts, decoys and unplaced scaffolds are not."""
    c = normalize_contig(contig)
    return c in AUTOSOMES or c in ("X", "Y", "MT")


def infer_build(contig_lengths: dict[str, int]) -> str | None:
    """Infer GRCh37/GRCh38 from header contig lengths.  None if ambiguous."""
    votes: dict[str, int] = {}
    for name, length in contig_lengths.items():
        table = _CHR_LENGTHS.get(normalize_contig(name))
        if table and length in table:
            votes[table[length]] = votes.get(table[length], 0) + 1
    if not votes:
        return None
    ranked = sorted(votes.items(), key=lambda kv: -kv[1])
    if len(ranked) > 1 and ranked[0][1] == ranked[1][1]:
        return None  # contradictory header, refuse to guess
    return ranked[0][0]


@dataclass(frozen=True)
class SexRegions:
    """Which X/Y positions are usable for a sex call under a given build."""

    build: str | None
    exclude_par: bool = True
    exclude_xtr: bool = True

    def _blocked(self, contig: str) -> list[tuple[int, int]]:
        if self.build is None:
            return []
        out: list[tuple[int, int]] = []
        if self.exclude_par:
            out += PAR.get(self.build, {}).get(contig, [])
        if self.exclude_xtr:
            out += XTR.get(self.build, {}).get(contig, [])
        return out

    def usable(self, contig: str, pos: int) -> bool:
        c = normalize_contig(contig)
        if c not in ("X", "Y"):
            return False
        return not any(start <= pos <= end for start, end in self._blocked(c))
