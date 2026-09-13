"""PED parsing and the pairwise relationships a pedigree *claims*.

Nothing here touches genotypes.  This module produces the claim; the identity
check produces the observation; the reconciliation is what the tool reports.
"""

from __future__ import annotations

import gzip
import os
from dataclasses import dataclass, field
from enum import Enum


class Sex(str, Enum):
    MALE = "male"
    FEMALE = "female"
    UNKNOWN = "unknown"

    @classmethod
    def from_ped(cls, token: str) -> Sex:
        return {"1": cls.MALE, "2": cls.FEMALE}.get(token.strip(), cls.UNKNOWN)


class Affection(str, Enum):
    UNAFFECTED = "unaffected"
    AFFECTED = "affected"
    UNKNOWN = "unknown"

    @classmethod
    def from_ped(cls, token: str) -> Affection:
        return {"1": cls.UNAFFECTED, "2": cls.AFFECTED}.get(token.strip(), cls.UNKNOWN)


class Relationship(str, Enum):
    SELF = "self"
    PARENT_OFFSPRING = "parent-offspring"
    FULL_SIB = "full-sib"
    HALF_SIB = "half-sib"
    RELATED_UNSPECIFIED = "related-unspecified"  # same family, no stated link
    UNRELATED = "unrelated"

    @property
    def degree(self) -> float | None:
        """Expected kinship coefficient, or None where the PED does not say."""
        return {
            "self": 0.5,
            "parent-offspring": 0.25,
            "full-sib": 0.25,
            "half-sib": 0.125,
            "unrelated": 0.0,
        }.get(self.value)


@dataclass
class Individual:
    family_id: str
    iid: str
    father: str
    mother: str
    sex: Sex
    affection: Affection

    @property
    def parents(self) -> tuple[str, str]:
        return (self.father, self.mother)


@dataclass
class Pedigree:
    individuals: dict[str, Individual] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    _kin_cache: dict[tuple[str, str], float] = field(default_factory=dict, repr=False)

    @property
    def sample_ids(self) -> list[str]:
        return list(self.individuals)

    @property
    def families(self) -> dict[str, list[str]]:
        out: dict[str, list[str]] = {}
        for ind in self.individuals.values():
            out.setdefault(ind.family_id, []).append(ind.iid)
        return out

    def affected(self, family_id: str | None = None) -> list[str]:
        return [
            i.iid
            for i in self.individuals.values()
            if i.affection is Affection.AFFECTED
            and (family_id is None or i.family_id == family_id)
        ]

    def _depth(self, iid: str, seen: frozenset[str] = frozenset()) -> int:
        ind = self.individuals.get(iid)
        if ind is None or iid in seen:
            return 0
        parents = [p for p in ind.parents if p != "0" and p in self.individuals]
        if not parents:
            return 0
        return 1 + max(self._depth(p, seen | {iid}) for p in parents)

    def kinship(self, a: str, b: str) -> float | None:
        """Expected kinship coefficient from the pedigree alone.

        The recursive definition, so it is correct for relationships the coarse
        :class:`Relationship` enum cannot name: grandparent-grandchild (0.125),
        avuncular (0.125), spouses (0.0), half-sibs (0.125), and anything else
        the pedigree implies.  ``Relationship`` collapses all of those into
        "related-unspecified", which means a mislabelled spouse or grandparent
        could never be caught.

        Returns None when either individual is absent from the pedigree.
        """
        if a not in self.individuals or b not in self.individuals:
            return None
        return self._kinship(a, b, 0)

    def _kinship(self, a: str, b: str, depth_guard: int) -> float:
        if depth_guard > 64:  # pathological pedigree loop
            return 0.0
        key = (a, b) if a <= b else (b, a)
        cached = self._kin_cache.get(key)
        if cached is not None:
            return cached

        ia, ib = self.individuals.get(a), self.individuals.get(b)
        if ia is None or ib is None:
            return 0.0

        if a == b:
            fa, ma = (p if p in self.individuals else "0" for p in ia.parents)
            inbreeding = self._kinship(fa, ma, depth_guard + 1) if "0" not in (fa, ma) else 0.0
            value = 0.5 * (1.0 + inbreeding)
        else:
            # Recurse on the younger individual, so we always walk upwards.
            if self._depth(a) > self._depth(b):
                a, b = b, a
                ia, ib = ib, ia
            fb, mb = (p if p in self.individuals else "0" for p in ib.parents)
            if fb == "0" and mb == "0":
                value = 0.0  # b is a founder and is not a
            else:
                total = 0.0
                for parent in (fb, mb):
                    if parent != "0":
                        total += self._kinship(a, parent, depth_guard + 1)
                value = 0.5 * total
        self._kin_cache[key] = value
        return value

    def relationship(self, a: str, b: str) -> Relationship:
        if a == b:
            return Relationship.SELF
        ia, ib = self.individuals.get(a), self.individuals.get(b)
        if ia is None or ib is None:
            return Relationship.UNRELATED
        if ia.family_id != ib.family_id:
            return Relationship.UNRELATED
        if a in ib.parents or b in ia.parents:
            return Relationship.PARENT_OFFSPRING
        pa = {p for p in ia.parents if p != "0"}
        pb = {p for p in ib.parents if p != "0"}
        shared = pa & pb
        if len(shared) == 2:
            return Relationship.FULL_SIB
        if len(shared) == 1:
            return Relationship.HALF_SIB
        return Relationship.RELATED_UNSPECIFIED


def read_ped(path: str | os.PathLike[str]) -> Pedigree:
    """Parse a PED file.  Malformed lines are recorded as warnings, not raised."""
    ped = Pedigree()
    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        with opener(path, "rt", errors="replace") as fh:  # type: ignore[operator]
            for lineno, raw in enumerate(fh, 1):
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split()
                if len(fields) < 6:
                    ped.warnings.append(f"line {lineno}: expected >=6 columns, got {len(fields)}")
                    continue
                fam, iid, fa, mo, sex, aff = fields[:6]
                if iid in ped.individuals:
                    ped.warnings.append(f"line {lineno}: duplicate individual id {iid!r}")
                    continue
                ped.individuals[iid] = Individual(
                    family_id=fam,
                    iid=iid,
                    father=fa,
                    mother=mo,
                    sex=Sex.from_ped(sex),
                    affection=Affection.from_ped(aff),
                )
    except OSError as exc:
        ped.warnings.append(f"could not read {path}: {exc}")
        return ped

    # A parent named in the PED but never defined is a common, silent error.
    for ind in list(ped.individuals.values()):
        # Their own parent. Arises when a relabelling renames a row but not the
        # references to it, and it is not harmless: the kinship recursion takes
        # it at face value and returns 0.5 for a pair that is nothing of the
        # kind, so every reconciliation downstream is measured against a number
        # the pedigree never meant.
        if ind.iid in ind.parents:
            ped.warnings.append(
                f"{ind.iid} is listed as their own parent; the relationships "
                f"computed through them cannot be trusted"
            )
    for ind in list(ped.individuals.values()):
        for parent in ind.parents:
            if parent != "0" and parent not in ped.individuals:
                ped.warnings.append(
                    f"{ind.iid} names parent {parent!r} which has no row of its own"
                )
    return ped
