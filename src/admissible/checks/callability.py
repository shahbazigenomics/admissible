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

# (minimum joint callable fraction, what the data permits you to say)
CLAIM_LADDER = [
    (0.80, "monogenic cause excluded for the models tested"),
    (0.50, "no candidate identified; limited power"),
    (0.00, "exploratory; a causal variant cannot be excluded"),
]


@dataclass
class CallabilityConfig:
    # Depth at or above which a position counts as callable.  Bins are selected
    # by their own lower bound, not by matching one hard-coded label, because
    # `--quantize 0:10:` is only the most common of many binnings a user may
    # have run.
    min_depth: int = 10
    strong: float = 0.80
    weak: float = 0.50
    # Models whose denominator is "every affected sample", vs the phenocopy
    # model, which needs all but one.
    full_penetrance_models: tuple[str, ...] = (
        "recessive",
        "compound-het",
        "dominant",
        "x-linked-recessive",
        "x-linked-dominant",
        "de-novo",
        "mitochondrial",
    )
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


def claim_for(fraction: float, cfg: CallabilityConfig) -> str:
    for floor, claim in CLAIM_LADDER:
        if fraction >= floor:
            return claim
    return CLAIM_LADDER[-1][1]


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

    affected = set(ped.affected(family_id)) if ped else set()
    usable = [c for c in coverage if not affected or c.label in affected]
    if not usable:
        return unknown(
            CHECK,
            "none of the supplied coverage files correspond to an affected sample "
            "in the pedigree",
            coverage_labels=[c.label for c in coverage],
            affected=sorted(affected),
        )
    # An affected sample with no coverage file silently shrinks the intersection
    # to the samples that do have one, which *raises* the reported fraction: the
    # number then answers a smaller question than the one being asked.
    missing_cov = sorted(affected - {c.label for c in coverage}) if affected else []
    unmatched = (
        sorted({c.label for c in coverage} - set(ped.sample_ids if ped else []))
        if affected
        else []
    )

    # A file whose depth bins never reach the threshold carries no evidence of
    # callability at all.  Answering from it would be the exact error this check
    # exists to prevent, so refuse rather than report a fraction.
    empty = [c for c in usable if not c.unlabelled and not c.bins_kept]
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

    if not affected:
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

    sets = [c.intervals for c in usable]
    n = len(sets)

    def fraction_of_target(iset: IntervalSet) -> tuple[int, float | None]:
        bp = total_bp(intersect(iset, target)) if target else total_bp(iset)
        frac = (bp / target_total) if target_total else None
        return bp, frac

    per_sample: dict[str, dict] = {}
    product = 1.0
    have_all_fractions = target_total is not None
    for c in usable:
        bp, frac = fraction_of_target(c.intervals)
        per_sample[c.label] = {"callable_bp": bp, "fraction_of_target": frac}
        if frac is None:
            have_all_fractions = False
        else:
            product *= frac

    joint = coverage_at_least(sets, n)
    joint_bp, joint_frac = fraction_of_target(joint)

    models: dict[str, dict] = {}
    for name in cfg.full_penetrance_models:
        models[name] = {
            "samples_required": n,
            "callable_bp": joint_bp,
            "fraction_of_target": joint_frac,
            "permitted_claim": claim_for(joint_frac, cfg) if joint_frac is not None else None,
        }
    if n >= 2:
        pheno = coverage_at_least(sets, n - 1)
        p_bp, p_frac = fraction_of_target(pheno)
        models[cfg.phenocopy_model] = {
            "samples_required": n - 1,
            "callable_bp": p_bp,
            "fraction_of_target": p_frac,
            "permitted_claim": claim_for(p_frac, cfg) if p_frac is not None else None,
        }

    findings: list[Finding] = []
    notes: list[str] = []
    if notes_affected:
        notes.append(notes_affected)
    notes += target_problems
    for c in usable:
        notes += [f"{c.label}: {p}" for p in c.problems]
    if missing_cov:
        findings.append(
            Finding(
                code="AFFECTED_WITHOUT_COVERAGE",
                severity=Severity.ERROR,
                message=(
                    f"{len(missing_cov)} affected sample(s) have no coverage file, so "
                    f"the fraction below is the intersection over {n} samples, not "
                    f"over the {len(affected)} the models actually require - it "
                    f"overstates what was searched"
                ),
                evidence={
                    "affected_without_coverage": missing_cov,
                    "affected_total": len(affected),
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

    for c in usable:
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
    notes.append(
        f"the {cfg.strong:.2f}/{cfg.weak:.2f} claim thresholds are conventions, not "
        f"derived quantities; they are configurable and should be stated as "
        f"conventions wherever this number is reported"
    )

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
            findings=[
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

    claim = claim_for(joint_frac, cfg)
    if have_all_fractions:
        notes.append(
            f"the product of the per-sample fractions is {product:.3f}, while the true "
            f"intersection is {joint_frac:.3f}. The product assumes coverage failures "
            f"are independent across samples; they are not, which is why it must never "
            f"be used as the joint callable fraction."
        )

    if joint_frac < cfg.weak:
        severity, code = Severity.BLOCKING, "CALLABILITY_EXPLORATORY"
    elif joint_frac < cfg.strong:
        severity, code = Severity.WARN, "CALLABILITY_LIMITED"
    else:
        severity, code = Severity.INFO, "CALLABILITY_ADEQUATE"

    findings.append(
        Finding(
            code=code,
            severity=severity,
            message=(
                f"{joint_frac:.1%} of the target is callable at once across all {n} "
                f"affected samples - the data supports at most: \"{claim}\""
            ),
            subjects=[c.label for c in usable],
            evidence={
                "joint_callable_fraction": joint_frac,
                "joint_callable_bp": joint_bp,
                "target_bp": target_total,
                "per_sample": per_sample,
                "product_of_per_sample_fractions": product if have_all_fractions else None,
                "permitted_claim": claim,
                "thresholds_are_conventions": {"strong": cfg.strong, "weak": cfg.weak},
                **(
                    {
                        "do_not_conclude": "no monogenic cause",
                        "next_steps": [
                            "report the callable fraction alongside every negative result"
                        ],
                    }
                    if joint_frac < cfg.strong
                    else {}
                ),
            },
        )
    )

    status = (
        Status.FAIL
        if joint_frac < cfg.weak
        else Status.WARN
        if joint_frac < cfg.strong
        else Status.PASS
    )
    # A fraction computed over fewer samples than the models require is not a
    # PASS whatever its value, because it is not the number it claims to be.
    if missing_cov and status is Status.PASS:
        status = Status.WARN
    return CheckResult(
        check=CHECK,
        status=status,
        summary=f"{joint_frac:.2f} of target searched across {n} affected samples",
        findings=findings,
        metrics={
            "n_affected_with_coverage": n,
            "target_bp": target_total,
            "joint_callable_bp": joint_bp,
            "joint_callable_fraction": joint_frac,
            "product_of_per_sample_fractions": product if have_all_fractions else None,
            "per_sample": per_sample,
            "models": models,
        },
        notes=notes,
    )
