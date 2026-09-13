"""Check 5 - inheritance model sweep.

Runs every model and reports the **count profile**, because the shape of the
profile is itself diagnostic. A cohort with zero candidates under every model
says something quite different from one with zero under recessive and forty
under dominant-with-reduced-penetrance.

Three rules govern how the counts are reported, and each exists because the
obvious alternative is misleading:

**A count is never reported without its own callable fraction.** Zero candidates
over 20% of the target and zero over 95% are not the same finding, and printing
them identically is how "no monogenic cause" gets written down. Callability is a
property of the model, not of the cohort - a fully penetrant model needs every
affected sample callable, phenocopy(n-1) needs all but one - so each model
carries its own denominator, supplied by check 4.

**Counts are not rarity-filtered unless you supply the frequency.** Every model
below is a segregation filter, and segregation alone leaves thousands of common
variants standing. Rarity needs a population allele frequency, which this tool
does not bundle; point ``--af-field`` at a frequency already present in your own
annotation and the counts become meaningful. Without it the sweep still runs and
says plainly that the numbers are segregation-only.

**An exclusion clause is only satisfied by an observed genotype.** Every model
here says some version of "and the unaffected samples do not carry it". A sample
with no genotype at that site does not satisfy that clause - it says nothing at
all. Treating a missing genotype as a non-carrier is the same error as treating
zero candidates as evidence of absence, one level down, and it silently inflates
every count on any call set with a meaningful missing rate. Candidates whose
exclusion could not be verified are counted and reported separately, never merged
into the headline number.

**Models that cannot be evaluated report UNKNOWN, not zero.** De novo needs both
parents genotyped; compound heterozygous needs gene assignment; X-linked needs
chrX and known sex. Where the input cannot support a model, saying "0 candidates"
would be a lie of the most dangerous kind - it looks like evidence of absence.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..contigs import normalize_contig
from ..model import CheckResult, Finding, Severity, Status, never_raises, unknown
from ..ped import Affection, Pedigree, Sex
from ..vcfio import HET, HOMALT, HOMREF, MISSING, GenotypeMatrix

CHECK = "models"


@dataclass
class ModelConfig:
    max_af: float | None = None  # requires an af field in the user's own annotation
    max_examples: int = 20
    # A model needs at least this many affected samples to mean anything.
    min_affected: int = 1


@dataclass
class ModelResult:
    name: str
    status: str  # "computed" | "not-applicable"
    n_candidates: int = 0
    reason: str = ""
    examples: list[str] = field(default_factory=list)
    samples_required: int | None = None
    # Sites that fit the model only because a sample whose non-carriage the model
    # depends on had no genotype there. Never folded into n_candidates.
    n_exclusion_unverified: int = 0

    def to_dict(self) -> dict:
        d = {
            "model": self.name,
            "status": self.status,
            "n_candidates": self.n_candidates if self.status == "computed" else None,
            "n_exclusion_unverified": (
                self.n_exclusion_unverified if self.status == "computed" else None
            ),
            "samples_required": self.samples_required,
        }
        if self.reason:
            d["reason"] = self.reason
        if self.examples:
            d["examples"] = self.examples
        return d


def _carrier(code: int) -> bool:
    return code in (HET, HOMALT)


def _locus(matrix: GenotypeMatrix, site_id: int) -> str:
    chrom, pos, ref, alt = matrix.index.keys[site_id]
    return f"{chrom}:{pos}{ref}>{alt}"


@never_raises(CHECK)
def check_models(
    matrix: GenotypeMatrix | None = None,
    ped: Pedigree | None = None,
    callability: dict[str, dict] | None = None,
    cfg: ModelConfig | None = None,
    family_id: str | None = None,
) -> CheckResult:
    cfg = cfg or ModelConfig()
    if matrix is None or not matrix.samples:
        return unknown(CHECK, "no genotypes were supplied")
    if ped is None or not ped.individuals:
        return unknown(
            CHECK,
            "no pedigree supplied; every model here is a segregation filter and "
            "segregation cannot be evaluated without affected/unaffected status",
        )

    present = set(matrix.samples)
    affected = [s for s in ped.affected(family_id) if s in present]
    unaffected = [
        i.iid
        for i in ped.individuals.values()
        if i.affection is Affection.UNAFFECTED
        and i.iid in present
        and (family_id is None or i.family_id == family_id)
    ]
    if len(affected) < cfg.min_affected:
        return unknown(
            CHECK,
            f"only {len(affected)} affected sample(s) have genotypes; a segregation "
            f"sweep needs at least {cfg.min_affected}",
            n_affected=len(affected),
            n_unaffected=len(unaffected),
        )

    code = matrix.codes
    n_sites = matrix.n_sites

    def gt(sample: str, site: int) -> int:
        vec = code.get(sample)
        return vec[site] if vec is not None else MISSING

    # --- rarity ---------------------------------------------------------
    af = matrix.site_af
    rarity_applied = cfg.max_af is not None and bool(af)

    def rare(site: int) -> bool:
        if not rarity_applied:
            return True
        value = af.get(site)
        return value is None or value <= cfg.max_af

    # --- per-site context ------------------------------------------------
    chrom_of = [normalize_contig(k[0]) for k in matrix.index.keys]
    sex = {i.iid: i.sex for i in ped.individuals.values()}

    results: list[ModelResult] = []

    def _exclusion_verified(site: int, who: list[str]) -> bool:
        """True when every sample the model relies on NOT carrying was observed."""
        return all(gt(u, site) != MISSING for u in who)

    def record(name, predicate, required=None, applicable=True, reason="", excludes=None):
        if not applicable:
            results.append(
                ModelResult(name, "not-applicable", reason=reason, samples_required=required)
            )
            return
        who = unaffected if excludes is None else excludes
        hits: list[int] = []
        unverified = 0
        for site in range(n_sites):
            if not rare(site):
                continue
            if not predicate(site):
                continue
            if who and not _exclusion_verified(site, who):
                unverified += 1
            else:
                hits.append(site)
        results.append(
            ModelResult(
                name,
                "computed",
                n_candidates=len(hits),
                examples=[_locus(matrix, s) for s in hits[: cfg.max_examples]],
                samples_required=required,
                n_exclusion_unverified=unverified,
            )
        )

    autosomal = lambda s: chrom_of[s] not in ("X", "Y", "MT")  # noqa: E731

    def recessive(site: int) -> bool:
        if not autosomal(site):
            return False
        if not all(gt(a, site) == HOMALT for a in affected):
            return False
        return all(gt(u, site) != HOMALT for u in unaffected)

    def dominant(site: int) -> bool:
        if not autosomal(site):
            return False
        if not all(_carrier(gt(a, site)) for a in affected):
            return False
        return all(not _carrier(gt(u, site)) for u in unaffected)

    def dominant_reduced(site: int) -> bool:
        """Unaffected carriers are permitted; the affected must all carry."""
        return autosomal(site) and all(_carrier(gt(a, site)) for a in affected)

    def phenocopy(site: int) -> bool:
        """All but one affected sample fits a recessive model."""
        if not autosomal(site) or len(affected) < 2:
            return False
        n_hom = sum(gt(a, site) == HOMALT for a in affected)
        if n_hom != len(affected) - 1:
            return False
        return all(gt(u, site) != HOMALT for u in unaffected)

    def x_recessive(site: int) -> bool:
        if chrom_of[site] != "X":
            return False
        for a in affected:
            g = gt(a, site)
            if sex.get(a) is Sex.MALE:
                if g not in (HOMALT, HET):  # hemizygous callers emit either
                    return False
            elif g != HOMALT:
                return False
        for u in unaffected:
            if sex.get(u) is Sex.MALE and _carrier(gt(u, site)):
                return False
            if sex.get(u) is not Sex.MALE and gt(u, site) == HOMALT:
                return False
        return True

    def x_dominant(site: int) -> bool:
        if chrom_of[site] != "X":
            return False
        return all(_carrier(gt(a, site)) for a in affected) and all(
            not _carrier(gt(u, site)) for u in unaffected
        )

    def mitochondrial(site: int) -> bool:
        return chrom_of[site] == "MT" and all(_carrier(gt(a, site)) for a in affected)

    # de novo needs both parents genotyped for at least one affected sample
    trios = [
        (a, ped.individuals[a].father, ped.individuals[a].mother)
        for a in affected
        if a in ped.individuals
        and ped.individuals[a].father in present
        and ped.individuals[a].mother in present
    ]

    def de_novo(site: int) -> bool:
        # Both parents must be OBSERVED hom-reference. A missing parent genotype
        # is the single most common way a false de novo enters a candidate list.
        for child, father, mother in trios:
            if (
                _carrier(gt(child, site))
                and gt(father, site) == HOMREF
                and gt(mother, site) == HOMREF
            ):
                return True
        return False

    # compound het needs gene assignment
    genes = matrix.site_gene

    def compound_het_counts() -> tuple[int, list[str]]:
        by_gene: dict[str, list[int]] = {}
        for site in range(n_sites):
            g = genes.get(site)
            if not g or not autosomal(site) or not rare(site):
                continue
            if all(gt(a, site) == HET for a in affected) and all(
                gt(u, site) != HOMALT for u in unaffected
            ):
                by_gene.setdefault(g, []).append(site)
        hits = {g: v for g, v in by_gene.items() if len(v) >= 2}
        return len(hits), sorted(hits)[: cfg.max_examples]

    n_aff = len(affected)
    record("recessive", recessive, required=n_aff)
    record("dominant", dominant, required=n_aff)
    record("dominant-reduced-penetrance", dominant_reduced, required=n_aff, excludes=[])
    record(
        "x-linked-recessive",
        x_recessive,
        required=n_aff,
        applicable=any(c == "X" for c in chrom_of) and all(sex.get(a) for a in affected),
        reason="no chrX sites, or sex unknown for an affected sample",
    )
    record(
        "x-linked-dominant",
        x_dominant,
        required=n_aff,
        applicable=any(c == "X" for c in chrom_of),
        reason="no chrX sites in the call set",
    )
    record(
        "mitochondrial",
        mitochondrial,
        required=n_aff,
        applicable=any(c == "MT" for c in chrom_of),
        reason="no chrM/MT sites in the call set",
        excludes=[],
    )
    record(
        "de-novo",
        de_novo,
        required=1,
        applicable=bool(trios),
        reason="no affected sample has both parents genotyped",
        excludes=[],
    )
    record(
        "phenocopy(n-1)",
        phenocopy,
        required=max(n_aff - 1, 1),
        applicable=n_aff >= 2,
        reason="needs at least two affected samples",
    )

    if genes:
        n_ch, ch_examples = compound_het_counts()
        results.append(
            ModelResult(
                "compound-het",
                "computed",
                n_candidates=n_ch,
                examples=ch_examples,
                samples_required=n_aff,
            )
        )
    else:
        results.append(
            ModelResult(
                "compound-het",
                "not-applicable",
                reason=(
                    "no gene assignment in the input annotation, so variants cannot be "
                    "grouped by gene"
                ),
                samples_required=n_aff,
            )
        )
    results.append(
        ModelResult(
            "two-locus",
            "not-applicable",
            reason=(
                "not implemented: the pairwise search is combinatorial and needs an "
                "explicit multiple-testing treatment before any count is reportable"
            ),
        )
    )

    # --- attach each model's own callable fraction ------------------------
    profile: list[dict] = []
    for r in results:
        entry = r.to_dict()
        cov = (callability or {}).get(r.name)
        entry["callable_fraction"] = (cov or {}).get("fraction_of_target")
        entry["permitted_claim"] = (cov or {}).get("permitted_claim")
        profile.append(entry)

    computed = [r for r in results if r.status == "computed"]
    total = sum(r.n_candidates for r in computed)
    n_unknown = len(results) - len(computed)

    findings: list[Finding] = []
    notes: list[str] = []
    if not rarity_applied:
        notes.append(
            "counts are segregation-only and NOT rarity-filtered: no allele-frequency "
            "field was supplied (--af-field), so common variants that happen to "
            "segregate are still counted"
        )
    if callability is None:
        notes.append(
            "no callable fraction was available, so every count below is reported "
            "without its denominator; a count of zero means nothing on its own"
        )
    if n_unknown:
        notes.append(
            f"{n_unknown} model(s) could not be evaluated on these inputs and are "
            f"reported as not-applicable rather than as zero candidates"
        )

    if not unaffected:
        notes.append(
            "no unaffected sample is genotyped, so every 'and no unaffected carries "
            "it' clause is vacuous: the dominant and X-linked filters are doing no "
            "exclusion work and their counts are inflated accordingly"
        )
    by_name = {r.name: r for r in computed}
    dom, red = by_name.get("dominant"), by_name.get("dominant-reduced-penetrance")
    if dom and red and dom.n_candidates == red.n_candidates and dom.n_candidates:
        notes.append(
            f"dominant and dominant-reduced-penetrance returned the same count "
            f"({dom.n_candidates}); the two differ only in whether unaffected carriers "
            f"are tolerated, so an identical count means the cohort contains no "
            f"genotyped unaffected carrier to tell them apart"
        )

    unverified_total = sum(r.n_exclusion_unverified for r in computed)
    if unverified_total:
        worst = max(computed, key=lambda r: r.n_exclusion_unverified)
        notes.append(
            f"{unverified_total} site(s) across all models fit only because a sample "
            f"whose non-carriage the model depends on had NO genotype there - most in "
            f"'{worst.name}' ({worst.n_exclusion_unverified}). These are excluded from "
            f"the counts below. A missing genotype is not a negative result"
        )

    zero_models = [r.name for r in computed if r.n_candidates == 0]
    if zero_models:
        frac = None
        if callability:
            fracs = [
                v.get("fraction_of_target")
                for v in callability.values()
                if v.get("fraction_of_target") is not None
            ]
            frac = min(fracs) if fracs else None
        findings.append(
            Finding(
                code="NO_CANDIDATES",
                severity=Severity.WARN,
                message=(
                    f"{len(zero_models)} model(s) yielded no candidate: "
                    f"{', '.join(zero_models)}"
                    + (
                        f" - interpret against a callable fraction of {frac:.2f}"
                        if frac is not None
                        else " - with no callable fraction, this is not evidence of absence"
                    )
                ),
                evidence={
                    "models_with_zero": zero_models,
                    "callable_fraction": frac,
                    "rarity_filtered": rarity_applied,
                    "do_not_conclude": "no monogenic cause",
                    "next_steps": [
                        "report the callable fraction alongside every model that "
                        "yielded no candidate"
                    ],
                },
            )
        )

    status = Status.WARN if findings else Status.PASS
    summary = (
        f"{total} candidate(s) across {len(computed)} model(s)"
        + (f", {n_unknown} not evaluable" if n_unknown else "")
        + (f", {unverified_total} unverifiable" if unverified_total else "")
        + ("" if rarity_applied else "; not rarity-filtered")
    )
    return CheckResult(
        check=CHECK,
        status=status,
        summary=summary,
        findings=findings,
        metrics={
            "n_affected": n_aff,
            "n_unaffected": len(unaffected),
            "affected": sorted(affected),
            "rarity_filtered": rarity_applied,
            "n_exclusion_unverified_total": unverified_total,
            "max_af": cfg.max_af,
            "profile": profile,
        },
        notes=notes,
    )
