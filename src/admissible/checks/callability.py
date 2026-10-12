"""Check 4 - callability and permitted claims.

How much of the target did you actually search, and what may you therefore say?

Two things are settled here and must not be quietly undone:

**A VCF cannot answer this question.** A VCF says nothing about a position that
produced no call: zero coverage and reference-identical are the same absence of a
line.  Callability needs coverage evidence - ``mosdepth --quantize`` BEDs, or gVCF
reference blocks.  Given plain VCFs the honest output is UNKNOWN, not an estimate.

**The joint callable fraction is an intersection, never a product.** Multiplying
per-sample callable fractions assumes coverage failures are independent across
samples.  They are strongly positively correlated, so the product understates the
truth, often by a factor of two or more.  The report prints both numbers side by
side precisely so the size of that error is visible rather than argued about.

Callability is computed **per model**, because it is a property of the genetic
hypothesis and not of the cohort: a fully penetrant model requires every affected
sample to be callable, phenocopy(n-1) requires all but one.  A single global
number cannot serve both, and a model count reported without its own denominator
is not interpretable.

**A model's callable fraction is over every sample that model needs, not just the
affected ones.**  A recessive or de novo candidate is only assessable where the
parents were callable too: the inheritance sweep (check 5) drops any site at which
a relative it relies on was not observed.  A required sample with no coverage file
makes that model UNKNOWN; it is never dropped from the intersection, because that
would raise the fraction.

The claim a fraction permits says what was searched and what was not.  It never
says a cause was "excluded": the search is depth-only, cannot see heterozygous
deletions or other structural variants, and the unsearched share is unsearched.
The permitted-claim thresholds are conventions, not derived quantities.  They are
configurable and the report says so wherever it prints them.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field

from ..intervals import (
    IntervalSet,
    coverage_at_least,
    intersect,
    merge,
    read_bed,
    read_named_bed,
    subtract,
    total_bp,
)
from ..model import CheckResult, Finding, Severity, Status, never_raises, unknown
from ..ped import Pedigree
from ..stats import median
from .callability_regions import analyse_regions

CHECK = "callability"

# Machine-readable rungs.  A consumer should branch on these, not on the wording.
RUNG_MOST = "SEARCHED_MOST"
RUNG_PARTLY = "SEARCHED_PARTLY"
RUNG_EXPLORATORY = "EXPLORATORY"
RUNG_UNKNOWN = "UNKNOWN"

# Always true of this check, whatever the fraction.
LIMITS = (
    "callability here is depth-based: it cannot see heterozygous deletions or other "
    "structural variants, which leave depth above the floor",
)


@dataclass
class CallabilityConfig:
    # Depth at or above which a position counts as callable.  Bins are selected
    # by their own lower bound, not by matching one hard-coded label, because
    # `--quantize 0:10:` is only the most common of many binnings a user may
    # have run.
    min_depth: int = 10
    strong: float = 0.80
    weak: float = 0.50
    # Which samples each model needs callable.  These are the models whose check-5
    # predicate also requires the affected samples' parents to have been observed
    # (they exclude variants carried by unaffected relatives), the models that do
    # not, and the one that is a union over trios.
    with_parents_models: tuple[str, ...] = (
        "recessive",
        "compound-het",
        "dominant",
        "x-linked-recessive",
        "x-linked-dominant",
    )
    affected_only_models: tuple[str, ...] = (
        "dominant-reduced-penetrance",
        "mitochondrial",
    )
    trio_model: str = "de-novo"
    # Require the affected samples' parents to be callable too (the default, and
    # the safe reading for a family). False describes a proband-only or
    # affected-only search: the parents are not sequenced, or not part of the
    # analysis, and the fraction says so.
    require_relatives: bool = True
    # Treat every sample as an independent case: its own analysis unit, with no
    # relatives. For a cohort of unrelated people, or several singletons.
    independent: bool = False
    # All but one affected sample, plus the parents.
    phenocopy_model: str = "phenocopy(n-1)"
    # --- per-region output (opt-in) ---
    by_region: bool = False
    # Removed once from a region name to get its gene: a trailing separator and
    # digits, optionally "exon"/"ex"/"e" (BRCA1_exon3, BRCA1_3, BRCA1.3 -> BRCA1).
    # A name this does not match is its own gene.  Names are not parsed beyond this.
    gene_regex: str = r"[_.\-:](?:exon|ex|e)?\d+$"
    # A gene holding a region with no callable base, or one the coverage file does
    # not report, cannot be called better than this rung however high its fraction.
    blind_region_cap: str = RUNG_PARTLY
    # Worst-first rows kept per unit; the rest are counted, not listed.
    max_region_rows: int = 5000
    max_listed_intervals: int = 10
    # Target/kit sanity.  Only a target at least this large is treated as a
    # whole-exome-scale kit; a gene panel's target is legitimately a sliver of a
    # sample's callable bases.  Both numbers are UNCALIBRATED defaults, there so the
    # check can be switched on, not values derived from a benchmark.
    kit_min_bp: int = 5_000_000
    kit_share_warn: float = 0.50


@dataclass
class SampleCoverage:
    label: str
    path: str
    intervals: IntervalSet = field(default_factory=dict)
    problems: list[str] = field(default_factory=list)
    # How the depth bins in the file were treated, so the report can say it.
    bins_kept: list[str] = field(default_factory=list)
    bins_straddling: list[str] = field(default_factory=list)
    bins_dropped: list[str] = field(default_factory=list)
    unlabelled: bool = False
    # Contigs the file mentions at all (any label).  Absent from here means "not
    # reported", which is not the same as "reported and uncallable".
    contigs_reported: set[str] = field(default_factory=set)


def _bin_bounds(label: str) -> tuple[float, float] | None:
    """Parse a mosdepth quantize label such as ``10:inf`` or ``0:5``.

    Returns ``None`` for anything that is not a depth bin, which is how a
    plain callable-region BED with a name in column 4 is told apart from a
    quantize output.
    """
    lo_s, _, hi_s = label.partition(":")
    if not _:
        return None
    try:
        lo = float(lo_s)
        hi = float("inf") if hi_s in ("inf", "") else float(hi_s)
    except ValueError:
        return None
    return lo, hi


def load_coverage(
    spec: dict[str, str | os.PathLike[str]], cfg: CallabilityConfig | None = None
) -> list[SampleCoverage]:
    """Read one callable-region BED per sample.

    Accepts either a pre-extracted BED of callable regions, or a raw mosdepth
    ``--quantize`` output, in which case the configured depth bin is selected.
    """
    cfg = cfg or CallabilityConfig()
    out: list[SampleCoverage] = []
    for label, path in spec.items():
        bed = read_bed(path)
        problems = list(bed.problems)
        depth_bins = {k: _bin_bounds(k) for k in bed.labels}
        depth_bins = {k: v for k, v in depth_bins.items() if v is not None}

        if not depth_bins:
            # A plain callable-region BED: every interval in it is callable by
            # construction, and there is nothing to select.
            out.append(
                SampleCoverage(
                    label=label, path=str(path), intervals=bed.intervals,
                    problems=problems, unlabelled=True,
                    contigs_reported=set(bed.contigs_seen),
                )
            )
            continue

        kept, straddling, dropped = [], [], []
        for name, (lo, hi) in sorted(depth_bins.items(), key=lambda kv: kv[1]):
            if lo >= cfg.min_depth:
                kept.append(name)
            elif hi > cfg.min_depth:
                # The bin spans the threshold, so it cannot be resolved without
                # re-running mosdepth.  Excluding it understates callability,
                # which is the safe direction; the check says so.
                straddling.append(name)
            else:
                dropped.append(name)

        merged: IntervalSet = {}
        for name in kept:
            part = read_bed(path, keep_label=name)
            for contig, ivs in part.intervals.items():
                merged.setdefault(contig, []).extend(ivs)
        merged = {c: merge(v) for c, v in merged.items()}

        if not kept:
            problems.append(
                f"no depth bin at or above {cfg.min_depth}x: file has "
                f"{', '.join(sorted(depth_bins))}"
            )
        out.append(
            SampleCoverage(
                label=label, path=str(path), intervals=merged, problems=problems,
                bins_kept=kept, bins_straddling=straddling, bins_dropped=dropped,
                contigs_reported=set(bed.contigs_seen),
            )
        )
    return out


def claim_rung_for(fraction: float, cfg: CallabilityConfig) -> str:
    """Which rung a joint callable fraction reaches, using the configured thresholds."""
    if fraction >= cfg.strong:
        return RUNG_MOST
    if fraction >= cfg.weak:
        return RUNG_PARTLY
    return RUNG_EXPLORATORY


def claim_for(fraction: float, cfg: CallabilityConfig) -> str:
    """What the data permits you to say, stated as what was and was not searched."""
    rung = claim_rung_for(fraction, cfg)
    unsearched = f"{max(0.0, 1.0 - fraction):.1%} of the target was not searched"
    if rung == RUNG_MOST:
        return f"no candidate identified in the searched portion; {unsearched}"
    if rung == RUNG_PARTLY:
        return (
            f"no candidate identified in the searched portion; {unsearched}, "
            f"so power is limited"
        )
    return f"exploratory only; {unsearched}, and a causal variant may lie in it"


@dataclass
class ModelRequirement:
    """The samples one model needs callable, and how they combine.

    ``all_of`` must all be callable at a position.  ``k_of`` is ``(samples, k)``:
    at least ``k`` of those.  ``trios`` is a union: a position counts if at least
    one complete (child, father, mother) trio is callable there, because a de novo
    call needs only one child and both of that child's parents.
    """

    model: str
    required: list[str] = field(default_factory=list)
    all_of: list[str] = field(default_factory=list)
    k_of: tuple[list[str], int] | None = None
    trios: list[tuple[str, str, str]] | None = None
    min_samples: int = 0
    missing: list[str] = field(default_factory=list)
    # Why the model cannot be evaluated at all; None when it can.
    unresolved: str | None = None
    # True when the model could be evaluated, but over fewer samples or trios than
    # it would ideally use, so its fraction is a floor on what was searched.
    lower_bound: bool = False
    notes: list[str] = field(default_factory=list)

    @property
    def evaluable(self) -> bool:
        return self.unresolved is None


def affected_missing(affected: list[str], have: set[str]) -> str:
    return ", ".join(a for a in affected if a not in have)


def _dedupe(items: list[str]) -> list[str]:
    seen: set[str] = set()
    return [x for x in items if not (x in seen or seen.add(x))]


def resolve_requirements(
    ped: Pedigree | None,
    affected: list[str],
    have: set[str],
    cfg: CallabilityConfig,
) -> dict[str, ModelRequirement]:
    """Work out, per model, which samples must be callable and whether they are.

    ``affected`` is the pedigree's affected samples (or, with no pedigree, every
    sample that has coverage).  ``have`` is the set of labels with a coverage file.
    A parent who is not in the pedigree cannot be required, and is noted; a parent
    who is in the pedigree but has no coverage file makes the model UNKNOWN.
    """
    individuals = ped.individuals if ped else {}

    parents: list[str] = []
    no_parent_notes: list[str] = []
    for a in affected:
        ind = individuals.get(a)
        if ind is None:
            continue
        for who, pid in (("father", ind.father), ("mother", ind.mother)):
            if pid and pid != "0" and pid in individuals:
                parents.append(pid)
            else:
                no_parent_notes.append(f"{a} has no {who} in the pedigree")
    parents = [p for p in _dedupe(parents) if p not in affected]
    if not cfg.require_relatives:
        parents = []
        no_parent_notes = [
            "relatives are not required (affected-only mode): the fraction describes "
            "a search over the affected samples alone"
        ]

    out: dict[str, ModelRequirement] = {}

    def need(model: str, members: list[str]) -> ModelRequirement:
        req = ModelRequirement(model=model, required=list(members), all_of=list(members))
        req.min_samples = len(members)
        req.missing = [m for m in members if m not in have]
        if req.missing:
            req.unresolved = (
                f"no coverage file for {', '.join(req.missing)}, which {model} needs"
            )
        return req

    for name in cfg.with_parents_models:
        req = need(name, _dedupe(affected + parents))
        if no_parent_notes:
            req.notes.append(
                "; ".join(no_parent_notes)
                + ("" if not cfg.require_relatives else " (so they could not be required)")
            )
        out[name] = req
    for name in cfg.affected_only_models:
        out[name] = need(name, list(affected))

    # de novo: a union over trios
    trios = [
        (a, individuals[a].father, individuals[a].mother)
        for a in affected
        if a in individuals
        and individuals[a].father in individuals
        and individuals[a].mother in individuals
    ]
    dn = ModelRequirement(model=cfg.trio_model, min_samples=3)
    if not cfg.require_relatives:
        dn.unresolved = (
            "de novo needs a child and both parents, and affected-only mode does not "
            "use relatives"
        )
    elif not trios:
        dn.unresolved = (
            "no affected sample has both parents in the pedigree, so de novo cannot "
            "be assessed" + (" (no pedigree was supplied)" if not individuals else "")
        )
    else:
        usable = [t for t in trios if all(m in have for m in t)]
        dn.required = _dedupe([m for t in trios for m in t])
        dn.missing = _dedupe([m for t in trios for m in t if m not in have])
        if not usable:
            dn.unresolved = (
                f"no complete trio has coverage; no coverage file for "
                f"{', '.join(dn.missing)}, which {cfg.trio_model} needs"
            )
        else:
            dn.trios = usable
            if len(usable) < len(trios):
                dn.lower_bound = True
                dn.notes.append(
                    f"{len(trios) - len(usable)} of {len(trios)} trio(s) were skipped for "
                    f"want of coverage ({', '.join(dn.missing)}); the fraction is a floor"
                )
    out[cfg.trio_model] = dn

    # phenocopy: all but one affected sample, plus the parents
    if len(affected) >= 2:
        k = len(affected) - 1
        have_aff = [a for a in affected if a in have]
        req = ModelRequirement(
            model=cfg.phenocopy_model,
            required=_dedupe(affected + parents),
            all_of=list(parents),
            k_of=(have_aff, k),
            min_samples=k + len(parents),
        )
        req.missing = [m for m in req.required if m not in have]
        missing_parents = [m for m in parents if m not in have]
        if len(have_aff) < k:
            req.unresolved = (
                f"only {len(have_aff)} of {len(affected)} affected samples have coverage; "
                f"{cfg.phenocopy_model} needs at least {k}"
            )
        elif missing_parents:
            req.unresolved = (
                f"no coverage file for {', '.join(missing_parents)}, which "
                f"{cfg.phenocopy_model} needs"
            )
        elif len(have_aff) < len(affected):
            req.lower_bound = True
            req.notes.append(
                f"{affected_missing(affected, have)} has no coverage, which uses up the "
                f"one sample the model allows to be an exception; the fraction is a floor"
            )
        out[cfg.phenocopy_model] = req
    return out


def requirement_intervals(
    req: ModelRequirement, cov: dict[str, IntervalSet]
) -> IntervalSet:
    """The positions at which a model's required samples are callable."""
    if req.trios is not None:
        per_trio = [
            coverage_at_least([cov[c], cov[f], cov[m]], 3) for c, f, m in req.trios
        ]
        return coverage_at_least(per_trio, 1)
    parts: list[IntervalSet] = []
    if req.all_of:
        parts.append(coverage_at_least([cov[m] for m in req.all_of], len(req.all_of)))
    if req.k_of is not None:
        members, k = req.k_of
        parts.append(coverage_at_least([cov[m] for m in members], k))
    if not parts:
        return {}
    result = parts[0]
    for extra in parts[1:]:
        result = intersect(result, extra)
    return result


def _analysis_units(
    ped: Pedigree | None,
    affected_ped: list[str],
    labels: list[str],
    cfg: CallabilityConfig,
) -> list[tuple[str, str | None, list[str]]]:
    """The units a claim is made about: ``(unit id, family id, affected samples)``.

    A claim about "the family" and a claim about "this person" are different, and a
    cohort is many of the second, not one big first. So a unit is one PED family
    with at least one affected member, or, for independent cases, one sample.
    With no pedigree and no ``independent`` flag, every sample is one group, which
    is the safe reading for a real family whose PED was not supplied.
    """
    if cfg.independent:
        members = affected_ped or list(labels)
        return [(m, None, [m]) for m in members]
    if ped is not None and affected_ped:
        by_family: dict[str, list[str]] = {}
        for a in affected_ped:
            by_family.setdefault(ped.individuals[a].family_id, []).append(a)
        return [(fam, fam, members) for fam, members in by_family.items()]
    return [("all", None, list(labels))]


def _model_table(
    reqs: dict[str, ModelRequirement],
    cov_sets: dict[str, IntervalSet],
    fraction_of_target,
    cfg: CallabilityConfig,
    exclude: IntervalSet | None = None,
) -> dict[str, dict]:
    models: dict[str, dict] = {}
    for name, req in reqs.items():
        entry: dict = {
            "status": "computed" if req.evaluable else "unknown",
            "samples_required": req.min_samples,
            "required_samples": req.required,
            "missing_coverage": req.missing,
            "lower_bound": req.lower_bound,
            "relatives_required": cfg.require_relatives,
            "reason": req.unresolved,
            "notes": req.notes,
            "callable_bp": None,
            "fraction_of_target": None,
            "permitted_claim": None,
            "claim_rung": RUNG_UNKNOWN,
        }
        if req.evaluable:
            callable_set = requirement_intervals(req, cov_sets)
            m_bp, m_frac = fraction_of_target(callable_set)
            if exclude is not None:
                # Subtraction distributes over intersection and k-of-n, so removing
                # the exclusion once from the joint set equals removing it from every
                # sample first.  The claim follows the stricter, filtered number; the
                # raw one stays beside it so both can be reported.
                entry["without_exclusion"] = {
                    "callable_bp": m_bp, "fraction_of_target": m_frac,
                }
                m_bp, m_frac = fraction_of_target(subtract(callable_set, exclude))
            entry["callable_bp"] = m_bp
            entry["fraction_of_target"] = m_frac
            if m_frac is not None:
                entry["permitted_claim"] = claim_for(m_frac, cfg)
                entry["claim_rung"] = claim_rung_for(m_frac, cfg)
        models[name] = entry
    return models


def _unit_summary(
    unit_id: str, family: str | None, affected: list[str], models: dict[str, dict]
) -> dict:
    computed = {m: e for m, e in models.items() if e["fraction_of_target"] is not None}
    gating_model = gating_fraction = claim = None
    rung = RUNG_UNKNOWN
    if computed:
        gating_model = min(computed, key=lambda m: computed[m]["fraction_of_target"])
        gating_fraction = computed[gating_model]["fraction_of_target"]
        claim = computed[gating_model]["permitted_claim"]
        rung = computed[gating_model]["claim_rung"]
    return {
        "unit": unit_id,
        "family": family,
        "affected": affected,
        "status": "computed" if computed else "unknown",
        "gating_model": gating_model,
        "gating_fraction": gating_fraction,
        "claim_rung": rung,
        "permitted_claim": claim,
        "models": models,
    }


def _contig_reporter(req: ModelRequirement, contigs: dict[str, set[str]]):
    """Whether the samples a model needs *report* a contig at all.

    Follows the model's own structure: every ``all_of`` sample, at least ``k`` of a
    k-of-n set, or at least one complete trio.  A contig nobody reports is not
    "uncallable"; it is unmeasured, and the region table says so.
    """

    def reported(contig: str) -> bool:
        if req.trios is not None:
            return any(all(contig in contigs.get(x, ()) for x in t) for t in req.trios)
        ok = all(contig in contigs.get(m, ()) for m in req.all_of)
        if ok and req.k_of is not None:
            members, k = req.k_of
            ok = sum(contig in contigs.get(m, ()) for m in members) >= k
        return ok

    return reported


def _by_region_block(
    target_bed: str | os.PathLike[str],
    unit_rows: list[dict],
    unit_reqs: dict[str, dict[str, ModelRequirement]],
    cov_sets: dict[str, IntervalSet],
    contigs: dict[str, set[str]],
    exclude: IntervalSet | None,
    cfg: CallabilityConfig,
) -> tuple[dict, list[Finding], list[str]]:
    """Per-region rows for each unit's strictest assessable model."""
    named = read_named_bed(target_bed)
    notes = [f"target BED (by region): {p}" for p in named.problems]
    units: dict[str, dict] = {}
    for u in unit_rows:
        model = u["gating_model"]
        if model is None:
            units[u["unit"]] = {"model": None, "reason": "no model could be assessed"}
            continue
        req = unit_reqs[u["unit"]][model]
        callable_set = requirement_intervals(req, cov_sets)
        if exclude is not None:
            callable_set = subtract(callable_set, exclude)
        out = analyse_regions(
            named,
            callable_set,
            _contig_reporter(req, contigs),
            gene_regex=cfg.gene_regex,
            claim_rung_for=lambda f: claim_rung_for(f, cfg),
            claim_for=lambda f: claim_for(f, cfg),
            blind_cap=cfg.blind_region_cap,
            max_rows=cfg.max_region_rows,
            max_listed=cfg.max_listed_intervals,
        )
        notes += out.pop("notes")
        out["model"] = model
        units[u["unit"]] = out
    if not named.has_names and named.regions:
        notes.append(
            "the target BED has no column-4 names on some or all lines, so those "
            "regions are named by their coordinates and each is its own gene"
        )
    block = {
        "path": named.path,
        "n_lines": named.n_lines,
        "n_unnamed": named.n_unnamed,
        "has_names": named.has_names,
        "n_regions": len(named.regions),
        "gene_regex": cfg.gene_regex,
        "basis": "each unit's strictest assessable model",
        "units": units,
    }
    findings: list[Finding] = []
    flagged = {
        uid: v for uid, v in units.items()
        if v.get("n_fully_uncallable_regions") or v.get("n_not_reported_regions")
    }
    if flagged:
        ex = {}
        for uid, v in flagged.items():
            ex[uid] = {
                "model": v["model"],
                "n_fully_uncallable_regions": v["n_fully_uncallable_regions"],
                "n_not_reported_regions": v["n_not_reported_regions"],
                "n_genes_with_unsearched_region": v["n_genes_with_unsearched_region"],
                "n_genes_capped": v["n_genes_capped"],
                "genes": [
                    g["gene"] for g in v["genes"]
                    if g["n_fully_uncallable"] or g["n_not_reported"]
                ][:20],
            }
        findings.append(
            Finding(
                code="TARGET_REGIONS_NOT_SEARCHED",
                severity=Severity.WARN,
                message=(
                    f"in {len(flagged)} of {len(units)} unit(s), some target region has "
                    f"no callable base or is not reported by the coverage file; a gene "
                    f"holding one cannot be rated above {cfg.blind_region_cap}, "
                    f"whatever its overall fraction"
                ),
                evidence={"units": ex, "see": "metrics.by_region"},
            )
        )
    return block, findings, notes


@never_raises(CHECK)
def check_callability(
    coverage: list[SampleCoverage] | None = None,
    ped: Pedigree | None = None,
    target_bed: str | os.PathLike[str] | None = None,
    cfg: CallabilityConfig | None = None,
    family_id: str | None = None,
    exclude_bed: str | os.PathLike[str] | None = None,
) -> CheckResult:
    cfg = cfg or CallabilityConfig()

    if not coverage:
        return unknown(
            CHECK,
            "no coverage evidence supplied; a VCF cannot answer this question, so "
            "this check needs mosdepth --quantize BEDs or gVCFs",
        )

    affected_ped = list(ped.affected(family_id)) if ped else []
    affected_set = set(affected_ped)
    # With no pedigree (or none affected), every sample with coverage is treated as
    # affected, and nobody's parents are known.
    affected_all = affected_ped or [c.label for c in coverage]
    known = set(ped.sample_ids) if ped and affected_set else None
    # Coverage that names nobody in the pedigree is ignored (and reported below).
    cov_by_label: dict[str, SampleCoverage] = {
        c.label: c for c in coverage if known is None or c.label in known
    }
    have = set(cov_by_label)

    usable = [c for c in coverage if not affected_set or c.label in affected_set]
    if not usable:
        return unknown(
            CHECK,
            "none of the supplied coverage files correspond to an affected sample "
            "in the pedigree",
            coverage_labels=[c.label for c in coverage],
            affected=sorted(affected_set),
        )
    missing_cov = sorted(affected_set - {c.label for c in coverage}) if affected_set else []
    unmatched = (
        sorted({c.label for c in coverage} - set(ped.sample_ids if ped else []))
        if affected_set
        else []
    )

    # A file whose depth bins never reach the threshold carries no evidence of
    # callability at all.  Answering from it would be the exact error this check
    # exists to prevent, so refuse rather than report a fraction.
    empty = [c for c in cov_by_label.values() if not c.unlabelled and not c.bins_kept]
    if empty:
        return unknown(
            CHECK,
            f"no depth bin at or above {cfg.min_depth}x in "
            f"{', '.join(c.label for c in empty)}; re-run mosdepth with a "
            f"--quantize boundary at {cfg.min_depth}, e.g. "
            f"--quantize 0:{cfg.min_depth}:",
            bins_present={c.label: sorted(c.bins_dropped + c.bins_straddling)
                          for c in empty},
        )

    if not affected_set and not cfg.independent:
        notes_affected = (
            "no pedigree supplied, so every sample with coverage was treated as "
            "affected and as one group; supply a PED to compute this against the right "
            "samples, or pass --independent if they are unrelated cases"
        )
    else:
        notes_affected = ""

    target: IntervalSet | None = None
    target_problems: list[str] = []
    if target_bed is not None:
        tb = read_bed(target_bed)
        target = tb.intervals
        target_problems = tb.problems
    target_total = total_bp(target) if target else None

    pre_findings: list[Finding] = []
    pre_notes: list[str] = []
    exclude: IntervalSet | None = None
    exclusion_info: dict | None = None
    if exclude_bed is not None:
        xb = read_bed(exclude_bed)
        pre_notes += [f"exclusion BED: {p}" for p in xb.problems]
        if not xb.intervals:
            pre_findings.append(
                Finding(
                    code="EXCLUDE_BED_UNUSABLE",
                    severity=Severity.WARN,
                    message=(
                        f"the exclusion BED {exclude_bed} has no usable interval, so it "
                        f"was not applied and the fractions are unfiltered"
                    ),
                    evidence={"path": str(exclude_bed), "problems": xb.problems},
                )
            )
        else:
            exclude = xb.intervals
            if target:
                in_target = total_bp(intersect(exclude, target))
                exclusion_info = {
                    "path": str(exclude_bed),
                    "bp": total_bp(exclude),
                    "bp_in_target": in_target,
                    "fraction_of_target_excluded": (
                        in_target / target_total if target_total else None
                    ),
                }
                if in_target == 0:
                    pre_notes.append(
                        "the exclusion BED does not overlap the target, so it changed "
                        "no fraction; check that it is on the same build"
                    )
            else:
                pre_notes.append(
                    "an exclusion BED was supplied but there is no target BED, so it "
                    "could not be applied"
                )
                exclude = None

    def fraction_of_target(iset: IntervalSet) -> tuple[int, float | None]:
        bp = total_bp(intersect(iset, target)) if target else total_bp(iset)
        frac = (bp / target_total) if target_total else None
        return bp, frac

    role: dict[str, str] = {}
    for label in cov_by_label:
        role[label] = "affected" if label in affected_set or not affected_set else "relative"

    per_sample: dict[str, dict] = {}
    product = 1.0
    have_all_fractions = target_total is not None
    for c in cov_by_label.values():
        bp, frac = fraction_of_target(c.intervals)
        per_sample[c.label] = {
            "callable_bp": bp, "fraction_of_target": frac, "role": role[c.label],
        }
        if role[c.label] == "affected":
            if frac is None:
                have_all_fractions = False
            else:
                product *= frac

    # The affected-only intersection.  Kept because consumers read it, but it is
    # not what any model is held to: see ``models`` below.
    sets = [c.intervals for c in usable]
    n = len(sets)
    joint = coverage_at_least(sets, n)
    joint_bp, joint_frac = fraction_of_target(joint)

    # --- per model: the samples it needs, and what that intersection covers ----
    # ``models`` is pooled over every affected sample in the run. It is what check 5
    # reads, because that sweep also asks "the same variant in all of them".
    # The CLAIMS are made per analysis unit, below.
    ped_for_reqs = None if cfg.independent else ped
    cov_sets = {label: c.intervals for label, c in cov_by_label.items()}
    reqs = resolve_requirements(ped_for_reqs, affected_all, have, cfg)
    models = _model_table(reqs, cov_sets, fraction_of_target, cfg, exclude)

    unit_defs = _analysis_units(ped, affected_ped, [c.label for c in coverage], cfg)
    unit_reqs: dict[str, dict[str, ModelRequirement]] = {}
    if len(unit_defs) == 1:
        unit_rows = [_unit_summary(*unit_defs[0], models)]
        unit_reqs[unit_defs[0][0]] = reqs
    else:
        unit_rows = []
        for uid, fam, members in unit_defs:
            u_reqs = resolve_requirements(ped_for_reqs, members, have, cfg)
            unit_reqs[uid] = u_reqs
            unit_rows.append(
                _unit_summary(
                    uid, fam, members,
                    _model_table(u_reqs, cov_sets, fraction_of_target, cfg, exclude),
                )
            )
    multi = len(unit_rows) > 1
    for label in cov_by_label:
        if any(label in u["affected"] for u in unit_rows):
            role[label] = "affected"
            per_sample[label]["role"] = "affected"

    def _no_cov(table: dict[str, dict]) -> dict[str, dict]:
        return {
            m: e for m, e in table.items()
            if e["status"] == "unknown" and e["missing_coverage"]
        }

    def _undefined(table: dict[str, dict]) -> dict[str, dict]:
        return {
            m: e for m, e in table.items()
            if e["status"] == "unknown" and not e["missing_coverage"]
        }

    no_cov_by_unit = {u["unit"]: _no_cov(u["models"]) for u in unit_rows}
    undefined_by_unit = {u["unit"]: _undefined(u["models"]) for u in unit_rows}
    any_no_cov = any(no_cov_by_unit.values())

    findings: list[Finding] = list(pre_findings)
    notes: list[str] = list(pre_notes)
    if notes_affected:
        notes.append(notes_affected)
    if cfg.independent:
        notes.append(
            f"--independent: {len(unit_rows)} sample(s) were each treated as their own "
            f"case, with no relatives; a claim holds for one sample at a time, and "
            f"the pooled figure answers 'the same variant in all of them'"
        )
    if not cfg.require_relatives:
        notes.append(
            "affected-only mode: relatives were not required, so these fractions "
            "describe a search over the affected samples alone; a variant could not "
            "be checked against unaffected relatives at positions they did not cover"
        )
    notes += target_problems
    for c in cov_by_label.values():
        notes += [f"{c.label}: {p}" for p in c.problems]
    used = {m for r in reqs.values() for m in r.required}
    for label in sorted(have - used - set(affected_all)):
        notes.append(
            f"{label}: coverage was supplied but no model needs this sample, so it was "
            f"not used"
        )
    if missing_cov:
        findings.append(
            Finding(
                code="AFFECTED_WITHOUT_COVERAGE",
                severity=Severity.ERROR,
                message=(
                    f"{len(missing_cov)} affected sample(s) have no coverage file, so "
                    f"the affected-only fraction is the intersection over {n} samples, "
                    f"not over the {len(affected_set)} affected; every model that "
                    f"needs them is UNKNOWN"
                ),
                evidence={
                    "affected_without_coverage": missing_cov,
                    "affected_total": len(affected_set),
                    "samples_intersected": n,
                },
            )
        )
    if unmatched:
        findings.append(
            Finding(
                code="COVERAGE_LABEL_UNKNOWN",
                severity=Severity.WARN,
                message=(
                    f"coverage supplied under label(s) {', '.join(unmatched)}, which "
                    f"name nobody in the pedigree; those files were ignored"
                ),
                evidence={"unmatched_labels": unmatched,
                          "pedigree_samples": sorted(ped.sample_ids) if ped else []},
            )
        )

    if any_no_cov:
        # model -> samples missing, merged over units; per-unit detail alongside
        merged: dict[str, list[str]] = {}
        for table in no_cov_by_unit.values():
            for m, e in table.items():
                merged.setdefault(m, [])
                merged[m] = sorted(set(merged[m]) | set(e["missing_coverage"]))
        units_with = [u for u, t in no_cov_by_unit.items() if t]
        findings.append(
            Finding(
                code="MODEL_REQUIRED_SAMPLE_WITHOUT_COVERAGE",
                severity=Severity.WARN,
                message=(
                    f"{len(merged)} model(s) are UNKNOWN"
                    + (f" in {len(units_with)} of {len(unit_rows)} unit(s)" if multi else "")
                    + " because a sample they need has no coverage file: "
                    + "; ".join(f"{m} needs {', '.join(v)}" for m, v in merged.items())
                    + ". They are not computed over the samples that do have coverage, "
                    "because that would overstate what was searched"
                ),
                evidence={
                    "models": merged,
                    **({"units": {u: {m: e["missing_coverage"] for m, e in t.items()}
                                  for u, t in no_cov_by_unit.items() if t}} if multi else {}),
                    "next_steps": [
                        "supply --coverage LABEL=BED for each sample listed, or, if the "
                        "relatives are not sequenced, rerun with --affected-only to "
                        "describe a search over the affected samples alone"
                    ],
                },
            )
        )
    undefined_models: dict[str, str] = {}
    for table in undefined_by_unit.values():
        for m, e in table.items():
            undefined_models.setdefault(m, e["reason"])
    if undefined_models:
        findings.append(
            Finding(
                code="MODEL_REQUIREMENTS_UNDEFINED",
                severity=Severity.INFO,
                message=(
                    "the pedigree cannot define which samples these model(s) need, so they "
                    "are UNKNOWN: "
                    + "; ".join(f"{m} ({r})" for m, r in undefined_models.items())
                ),
                evidence={"models": undefined_models},
            )
        )

    for c in cov_by_label.values():
        if c.unlabelled:
            notes.append(
                f"{c.label}: the BED carries no depth bins, so every interval in "
                f"it was taken as callable on trust; admissible cannot verify a "
                f"{cfg.min_depth}x floor from it"
            )
        elif c.bins_straddling:
            notes.append(
                f"{c.label}: bin(s) {', '.join(c.bins_straddling)} span "
                f"{cfg.min_depth}x and were excluded, so this sample's callable "
                f"fraction is a lower bound; re-quantize with a boundary at "
                f"{cfg.min_depth} to resolve it"
            )
    for m, e in models.items():
        for note in e["notes"]:
            notes.append(f"{m}: {note}")
    notes.append(
        f"the {cfg.strong:.2f}/{cfg.weak:.2f} claim thresholds are conventions, not "
        f"derived quantities; they are configurable and should be stated as "
        f"conventions wherever this number is reported"
    )
    notes.extend(LIMITS)

    mode = "independent" if cfg.independent else ("family" if ped and affected_set else "pooled")
    analysis = {
        "mode": mode,
        "relatives_required": cfg.require_relatives,
        "n_units": len(unit_rows),
    }

    if target is None:
        notes.append(
            "no target BED supplied, so fractions could not be computed: the "
            "denominator would otherwise be the union of what was covered, which "
            "makes the number self-referential. Callable base counts are reported "
            "instead."
        )
        if cfg.by_region:
            notes.append("--by-region needs a target BED and was not applied")
        return CheckResult(
            check=CHECK,
            status=Status.UNKNOWN,
            summary=f"{joint_bp:,} bp callable in all {n} affected, no target to divide by",
            findings=findings
            + [
                Finding(
                    code="NO_TARGET_INTERVALS",
                    severity=Severity.WARN,
                    message=(
                        "callable base pairs were computed but no capture-kit target "
                        "BED was supplied, so 'fraction of target searched' is undefined"
                    ),
                    evidence={
                        "joint_callable_bp": joint_bp,
                        "per_sample": per_sample,
                        "next_steps": ["supply the capture kit's target BED"],
                    },
                )
            ],
            metrics={"per_sample": per_sample, "joint_callable_bp": joint_bp,
                     "models": models, "analysis": analysis},
            notes=notes,
        )

    if have_all_fractions:
        notes.append(
            f"the product of the per-sample fractions is {product:.3f}, while the true "
            f"intersection is {joint_frac:.3f}. The product assumes coverage failures "
            f"are independent across samples; they are not, which is why it must never "
            f"be used as the joint callable fraction."
        )

    # --- scope: contigs the coverage file does not mention --------------------
    contigs = {label: set(c.contigs_reported) for label, c in cov_by_label.items()}
    target_contigs = sorted(target)
    scope: dict[str, dict] = {}
    for label, seen in contigs.items():
        absent = [c for c in target_contigs if c not in seen]
        if absent:
            scope[label] = {
                "contigs": absent,
                "target_bp": sum(e - st for c in absent for st, e in target[c]),
            }
    if scope:
        listed = "; ".join(f"{k}: " + ", ".join(v["contigs"][:6]) for k, v in scope.items())
        findings.append(
            Finding(
                code="TARGET_CONTIG_NOT_IN_COVERAGE",
                severity=Severity.WARN,
                message=(
                    f"{len(scope)} coverage file(s) never mention some target contig(s) "
                    f"({listed}). "
                    f"The fractions count those bases as not searched, but the file "
                    f"does not say whether they were uncallable or simply not reported"
                ),
                evidence={
                    "per_sample": scope,
                    "next_steps": [
                        "re-run mosdepth over every contig in the target, or confirm "
                        "the coverage file is the one for this sample and build"
                    ],
                },
            )
        )

    # --- sanity: does this target look like the kit that produced the coverage? ---
    # Share of a sample's callable bases that fall inside the target.  Gated on
    # target size, because a gene panel is legitimately a sliver of any sample's
    # callable bases.  Thresholds are uncalibrated.
    kit: dict[str, float] = {}
    if target_total and target_total >= cfg.kit_min_bp:
        for label, c in cov_by_label.items():
            total = total_bp(c.intervals)
            if total > 0:
                share = per_sample[label]["callable_bp"] / total
                if share < cfg.kit_share_warn:
                    kit[label] = share
    if kit:
        findings.append(
            Finding(
                code="TARGET_KIT_MISMATCH_SUSPECTED",
                severity=Severity.WARN,
                message=(
                    f"only a small share of the callable bases of "
                    f"{', '.join(f'{k} ({v:.0%})' for k, v in kit.items())} fall inside "
                    f"the supplied target ({target_total:,} bp). That is what coverage "
                    f"from a different kit, a whole-genome run or a different build "
                    f"looks like. The build cannot be read from contig names, so this "
                    f"is a plausibility check, not a test; its thresholds are "
                    f"uncalibrated"
                ),
                evidence={
                    "share_callable_in_target": kit,
                    "threshold": cfg.kit_share_warn,
                    "applies_when_target_bp_at_least": cfg.kit_min_bp,
                    "calibrated": False,
                },
            )
        )

    by_region_block = None
    if cfg.by_region:
        by_region_block, region_findings, region_notes = _by_region_block(
            target_bed, unit_rows, unit_reqs, cov_sets, contigs, exclude, cfg
        )
        findings += region_findings
        notes += region_notes

    assessable = [u for u in unit_rows if u["gating_fraction"] is not None]
    by_rung: dict[str, int] = {}
    for u in unit_rows:
        by_rung[u["claim_rung"]] = by_rung.get(u["claim_rung"], 0) + 1
    base_metrics = {
        "n_affected_with_coverage": n,
        "target_bp": target_total,
        "joint_callable_bp": joint_bp,
        "joint_callable_fraction": joint_frac,
        "product_of_per_sample_fractions": product if have_all_fractions else None,
        "per_sample": per_sample,
        # pooled over every affected sample; check 5 reads this
        "models": models,
        "units": unit_rows,
        "analysis": analysis,
    }
    if by_region_block is not None:
        base_metrics["by_region"] = by_region_block
    if exclusion_info is not None:
        base_metrics["exclusion"] = exclusion_info
    if not assessable:
        return CheckResult(
            check=CHECK,
            status=Status.UNKNOWN,
            summary=(
                f"no model could be assessed: {len(models)} model(s) need a sample that has "
                f"no coverage or a relative the pedigree does not list"
            ),
            findings=findings,
            metrics={**base_metrics, "gating_model": None, "gating_fraction": None},
            notes=notes,
        )

    # The headline is the worst unit's strictest assessable model, because a single
    # number has to be safe to quote for every unit and every model at once.
    worst = min(assessable, key=lambda u: u["gating_fraction"])
    gating_model = worst["gating_model"]
    gating = worst["models"][gating_model]
    gating_frac = worst["gating_fraction"]
    claim = worst["permitted_claim"]
    rung = worst["claim_rung"]
    cohort = {
        "n_units": len(unit_rows),
        "units_by_rung": by_rung,
        "worst_unit": worst["unit"],
        "worst_fraction": gating_frac,
        "median_fraction": median([u["gating_fraction"] for u in assessable]),
        "n_units_unknown": len(unit_rows) - len(assessable),
    }

    if gating_frac < cfg.weak:
        severity, code = Severity.BLOCKING, "CALLABILITY_EXPLORATORY"
    elif gating_frac < cfg.strong:
        severity, code = Severity.WARN, "CALLABILITY_LIMITED"
    else:
        severity, code = Severity.INFO, "CALLABILITY_ADEQUATE"

    scope = (
        f"the worst of {len(unit_rows)} units ({worst['unit']}; "
        + ", ".join(f"{v} {k}" for k, v in sorted(by_rung.items()))
        + ")"
        if multi
        else "the target"
    )
    findings.append(
        Finding(
            code=code,
            severity=severity,
            message=(
                f"{gating_frac:.1%} of the target is callable at once in the samples the "
                f"strictest assessable model ({gating_model}) needs "
                f"({', '.join(gating['required_samples'])})"
                + (f", for {scope}" if multi else "")
                + f" - the data supports at most: \"{claim}\". {LIMITS[0].capitalize()}"
            ),
            subjects=list(gating["required_samples"]),
            evidence={
                "gating_unit": worst["unit"],
                "gating_model": gating_model,
                "gating_fraction": gating_frac,
                "claim_rung": rung,
                "joint_callable_fraction": joint_frac,
                "joint_callable_bp": joint_bp,
                "target_bp": target_total,
                "per_sample": per_sample,
                "product_of_per_sample_fractions": product if have_all_fractions else None,
                "permitted_claim": claim,
                "limits": list(LIMITS),
                "cohort": cohort,
                "models_unknown": sorted(
                    m for m, e in worst["models"].items() if e["status"] == "unknown"
                ),
                "thresholds_are_conventions": {"strong": cfg.strong, "weak": cfg.weak},
                **(
                    {
                        "do_not_conclude": "no monogenic cause",
                        "next_steps": [
                            "report the callable fraction alongside every negative result"
                        ],
                    }
                    if gating_frac < cfg.strong
                    else {}
                ),
            },
        )
    )

    status = (
        Status.FAIL
        if gating_frac < cfg.weak
        else Status.WARN
        if gating_frac < cfg.strong
        else Status.PASS
    )
    # A number that leaves a model unassessed for want of a coverage file the user
    # could have supplied is not a PASS whatever its value, and neither is a cohort
    # with a unit nothing could be said about. A model the pedigree cannot define
    # (de novo with no trio) is not that: check 5 calls it not-applicable, and
    # holding it against PASS would make PASS unreachable for any family without a
    # trio, or any single sample.
    if status is Status.PASS and (missing_cov or any_no_cov or cohort["n_units_unknown"]):
        status = Status.WARN
    return CheckResult(
        check=CHECK,
        status=status,
        # Short on purpose: the text report truncates a summary to one terminal line.
        summary=(
            f"{gating_frac:.2f} worst of {len(unit_rows)} units ({worst['unit']})"
            if multi
            else f"{gating_frac:.2f} of target searched ({gating_model})"
        ),
        findings=findings,
        metrics={
            **base_metrics,
            "gating_unit": worst["unit"],
            "gating_model": gating_model,
            "gating_fraction": gating_frac,
            "cohort": cohort,
        },
        notes=notes,
    )
