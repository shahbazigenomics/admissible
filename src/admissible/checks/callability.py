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
    total_bp,
)
from ..model import CheckResult, Finding, Severity, Status, never_raises, unknown
from ..ped import Pedigree

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
    # All but one affected sample, plus the parents.
    phenocopy_model: str = "phenocopy(n-1)"


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
            req.notes.append("; ".join(no_parent_notes) + " (so they could not be required)")
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
    if not trios:
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


@never_raises(CHECK)
def check_callability(
    coverage: list[SampleCoverage] | None = None,
    ped: Pedigree | None = None,
    target_bed: str | os.PathLike[str] | None = None,
    cfg: CallabilityConfig | None = None,
    family_id: str | None = None,
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

    if not affected_set:
        notes_affected = (
            "no pedigree supplied, so every sample with coverage was treated as "
            "affected; supply a PED to compute this against the right samples"
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
    reqs = resolve_requirements(ped, affected_all, have, cfg)
    cov_sets = {label: c.intervals for label, c in cov_by_label.items()}
    models: dict[str, dict] = {}
    for name, req in reqs.items():
        entry: dict = {
            "status": "computed" if req.evaluable else "unknown",
            "samples_required": req.min_samples,
            "required_samples": req.required,
            "missing_coverage": req.missing,
            "lower_bound": req.lower_bound,
            "reason": req.unresolved,
            "notes": req.notes,
            "callable_bp": None,
            "fraction_of_target": None,
            "permitted_claim": None,
            "claim_rung": RUNG_UNKNOWN,
        }
        if req.evaluable:
            m_bp, m_frac = fraction_of_target(requirement_intervals(req, cov_sets))
            entry["callable_bp"] = m_bp
            entry["fraction_of_target"] = m_frac
            if m_frac is not None:
                entry["permitted_claim"] = claim_for(m_frac, cfg)
                entry["claim_rung"] = claim_rung_for(m_frac, cfg)
        models[name] = entry

    findings: list[Finding] = []
    notes: list[str] = []
    if notes_affected:
        notes.append(notes_affected)
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

    no_cov_models = {
        m: e for m, e in models.items() if e["status"] == "unknown" and e["missing_coverage"]
    }
    undefined_models = {
        m: e for m, e in models.items() if e["status"] == "unknown" and not e["missing_coverage"]
    }
    if no_cov_models:
        findings.append(
            Finding(
                code="MODEL_REQUIRED_SAMPLE_WITHOUT_COVERAGE",
                severity=Severity.WARN,
                message=(
                    f"{len(no_cov_models)} model(s) are UNKNOWN because a sample they need "
                    f"has no coverage file: "
                    + "; ".join(
                        f"{m} needs {', '.join(e['missing_coverage'])}"
                        for m, e in no_cov_models.items()
                    )
                    + ". They are not computed over the samples that do have coverage, "
                    "because that would overstate what was searched"
                ),
                evidence={
                    "models": {m: e["missing_coverage"] for m, e in no_cov_models.items()},
                    "next_steps": [
                        "supply --coverage LABEL=BED for each sample listed, or accept "
                        "that these models were not assessed"
                    ],
                },
            )
        )
    if undefined_models:
        findings.append(
            Finding(
                code="MODEL_REQUIREMENTS_UNDEFINED",
                severity=Severity.INFO,
                message=(
                    "the pedigree cannot define which samples these model(s) need, so they "
                    "are UNKNOWN: "
                    + "; ".join(f"{m} ({e['reason']})" for m, e in undefined_models.items())
                ),
                evidence={"models": {m: e["reason"] for m, e in undefined_models.items()}},
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

    if target is None:
        notes.append(
            "no target BED supplied, so fractions could not be computed: the "
            "denominator would otherwise be the union of what was covered, which "
            "makes the number self-referential. Callable base counts are reported "
            "instead."
        )
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
            metrics={"per_sample": per_sample, "joint_callable_bp": joint_bp, "models": models},
            notes=notes,
        )

    if have_all_fractions:
        notes.append(
            f"the product of the per-sample fractions is {product:.3f}, while the true "
            f"intersection is {joint_frac:.3f}. The product assumes coverage failures "
            f"are independent across samples; they are not, which is why it must never "
            f"be used as the joint callable fraction."
        )

    computed = {m: e for m, e in models.items() if e["fraction_of_target"] is not None}
    if not computed:
        return CheckResult(
            check=CHECK,
            status=Status.UNKNOWN,
            summary=(
                f"no model could be assessed: {len(models)} model(s) need a sample that has "
                f"no coverage or a relative the pedigree does not list"
            ),
            findings=findings,
            metrics={
                "n_affected_with_coverage": n,
                "target_bp": target_total,
                "joint_callable_bp": joint_bp,
                "joint_callable_fraction": joint_frac,
                "product_of_per_sample_fractions": product if have_all_fractions else None,
                "per_sample": per_sample,
                "models": models,
                "gating_model": None,
                "gating_fraction": None,
            },
            notes=notes,
        )

    # The headline is the strictest model that could be assessed, because a single
    # number has to be safe to quote for every model at once.
    gating_model = min(computed, key=lambda m: computed[m]["fraction_of_target"])
    gating = computed[gating_model]
    gating_frac = gating["fraction_of_target"]
    claim = claim_for(gating_frac, cfg)
    rung = claim_rung_for(gating_frac, cfg)

    if gating_frac < cfg.weak:
        severity, code = Severity.BLOCKING, "CALLABILITY_EXPLORATORY"
    elif gating_frac < cfg.strong:
        severity, code = Severity.WARN, "CALLABILITY_LIMITED"
    else:
        severity, code = Severity.INFO, "CALLABILITY_ADEQUATE"

    findings.append(
        Finding(
            code=code,
            severity=severity,
            message=(
                f"{gating_frac:.1%} of the target is callable at once in the samples the "
                f"strictest assessable model ({gating_model}) needs "
                f"({', '.join(gating['required_samples'])}) - the data supports at most: "
                f"\"{claim}\". {LIMITS[0].capitalize()}"
            ),
            subjects=list(gating["required_samples"]),
            evidence={
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
                "models_unknown": sorted(
                    m for m, e in models.items() if e["status"] == "unknown"
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
    # could have supplied is not a PASS whatever its value. A model the pedigree
    # cannot define (de novo with no trio) is not that: check 5 calls it
    # not-applicable, and holding it against PASS would make PASS unreachable for
    # any family without a trio.
    if status is Status.PASS and (missing_cov or no_cov_models):
        status = Status.WARN
    return CheckResult(
        check=CHECK,
        status=status,
        # Short on purpose: the text report truncates a summary to one terminal line.
        summary=f"{gating_frac:.2f} of target searched ({gating_model})",
        findings=findings,
        metrics={
            "n_affected_with_coverage": n,
            "target_bp": target_total,
            "joint_callable_bp": joint_bp,
            "joint_callable_fraction": joint_frac,
            "product_of_per_sample_fractions": product if have_all_fractions else None,
            "per_sample": per_sample,
            "models": models,
            "gating_model": gating_model,
            "gating_fraction": gating_frac,
        },
        notes=notes,
    )
