"""Check 1 - sample identity.

Answers three questions, in this order, because each one poisons the next:

1. Is each sample the sex the pedigree says it is?
2. Is any sample a duplicate of any other sample **anywhere in the cohort**?
3. Do the observed genotype relationships match the declared pedigree?

Question 2 is cohort-wide on purpose.  A within-family duplicate scan cannot
find a swap that lives *between* two families, and that is the failure mode
that motivated this tool.

Two relatedness engines are used, and which one applies is a property of the
input, not a preference:

* ``jaccard`` / ``agreement`` over non-reference sites.  Always available.
  Depth- and pipeline-dependent, so the boundary is calibrated on the cohort
  rather than hard-coded.  It can separate "related" from "unrelated"; it
  cannot resolve relationship *degree*.
* ``KING-robust`` kinship and ``IBS0``.  Allele-frequency-free, so no external
  database is needed, and robust to depth.  Requires hom-ref genotypes to be
  observed, which is only true within a single multi-sample or joint-called
  VCF - in separate single-sample VCFs, "absent" and "hom-ref" are the same
  string of bytes.  Where it applies it resolves degree, and IBS0 separates
  parent-offspring from full sibs, which Jaccard cannot do at all.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..contigs import SexRegions, normalize_contig
from ..model import CheckResult, Finding, Severity, Status, never_raises
from ..ped import Pedigree, Relationship, Sex
from ..stats import mad, median, wilson_interval
from ..vcfio import HET, HOMALT, HOMREF, MISSING, GenotypeMatrix

CHECK = "identity"

# KING-robust degree bands (Manichaikul et al. 2010).  Allele-frequency-free.
KING_BANDS: list[tuple[float, float, str]] = [
    (0.354, 1.0, "duplicate/MZ"),
    (0.177, 0.354, "first-degree"),
    (0.0884, 0.177, "second-degree"),
    (0.0442, 0.0884, "third-degree"),
    (-1.0, 0.0442, "unrelated"),
]

FIRST_DEGREE = (Relationship.PARENT_OFFSPRING, Relationship.FULL_SIB)
AMBIGUOUS_DECLARATIONS = (Relationship.HALF_SIB, Relationship.RELATED_UNSPECIFIED)


@dataclass
class IdentityConfig:
    min_x_sites: int = 50
    min_shared_sites: int = 500
    min_informative_sites: int = 2000
    # A pair is a duplicate candidate at or above both of these.  Agreement is
    # the discriminating term: no true relative reaches ~0.95 concordance at
    # shared non-reference sites, because sibs and parent-offspring pairs
    # routinely differ in zygosity at sites they share.  Jaccard is deliberately
    # permissive, because the same individual sequenced twice at different depth
    # has a genuinely low site-set overlap.
    dup_agreement: float = 0.95
    dup_jaccard: float = 0.60
    # Cohort calibration refuses to reconcile below this separation.
    min_band_separation: float = 2.0
    # KING-robust degree boundaries (Manichaikul et al. 2010), used to decide
    # when an observed kinship is too far from the pedigree's to be noise.
    first_degree_min: float = 0.177
    second_degree_min: float = 0.0884
    unrelated_max: float = 0.0442
    # Fallback sex boundary if the cohort gives no clean gap to calibrate on.
    default_sex_boundary: float = 0.45
    exclude_par: bool = True
    exclude_xtr: bool = True


@dataclass
class SampleSets:
    assessed: set[int] = field(default_factory=set)
    het: set[int] = field(default_factory=set)
    homalt: set[int] = field(default_factory=set)
    homref: set[int] = field(default_factory=set)

    @property
    def nonref(self) -> set[int]:
        return self.het | self.homalt


@dataclass
class SexEvidence:
    sample: str
    n_x_sites: int
    n_x_het: int
    x_het_frac: float
    ci_low: float
    ci_high: float
    n_y_sites: int
    inferred: Sex
    reason: str

    def to_dict(self) -> dict:
        return {
            "sample": self.sample,
            "chrX_sites": self.n_x_sites,
            "chrX_het": self.n_x_het,
            "chrX_het_frac": self.x_het_frac,
            "chrX_het_ci95": [self.ci_low, self.ci_high],
            "chrY_sites": self.n_y_sites,
            "inferred_sex": self.inferred.value,
            "basis": self.reason,
        }


@dataclass
class PairEvidence:
    a: str
    b: str
    dense: bool
    n_shared: int
    n_union: int
    jaccard: float
    agreement: float
    n_informative: int = 0
    king_phi: float | None = None
    ibs0_rate: float | None = None
    n_nonref_a: int = 0
    n_nonref_b: int = 0

    @property
    def king_band(self) -> str | None:
        if self.king_phi is None:
            return None
        for lo, hi, label in KING_BANDS:
            if lo < self.king_phi <= hi:
                return label
        return "unrelated"

    def to_dict(self) -> dict:
        d = {
            "a": self.a,
            "b": self.b,
            "mode": "dense" if self.dense else "sparse-nonref",
            "n_shared_nonref": self.n_shared,
            "n_union_nonref": self.n_union,
            "n_nonref_a": self.n_nonref_a,
            "n_nonref_b": self.n_nonref_b,
            "jaccard": self.jaccard,
            "agreement": self.agreement,
        }
        if self.dense:
            d |= {
                "n_informative": self.n_informative,
                "king_robust_phi": self.king_phi,
                "king_band": self.king_band,
                "ibs0_rate": self.ibs0_rate,
            }
        return d


def build_sets(matrix: GenotypeMatrix) -> dict[str, SampleSets]:
    out: dict[str, SampleSets] = {}
    for sample in matrix.samples:
        vec = matrix.codes.get(sample)
        s = SampleSets()
        if vec is None:
            out[sample] = s
            continue
        for site_id, code in enumerate(vec):
            if code == MISSING:
                continue
            s.assessed.add(site_id)
            if code == HET:
                s.het.add(site_id)
            elif code == HOMALT:
                s.homalt.add(site_id)
            elif code == HOMREF:
                s.homref.add(site_id)
        out[sample] = s
    return out


def sex_evidence(
    matrix: GenotypeMatrix,
    sets: dict[str, SampleSets],
    cfg: IdentityConfig,
    build: str | None,
) -> tuple[dict[str, SexEvidence], float, list[str]]:
    """Per-sample X heterozygosity, PAR- and XTR-excluded, with a cohort boundary."""
    notes: list[str] = []
    regions = SexRegions(build=build, exclude_par=cfg.exclude_par, exclude_xtr=cfg.exclude_xtr)
    if build is None:
        notes.append(
            "genome build could not be inferred from the VCF headers; PAR and XTR "
            "were NOT excluded, which inflates male chrX heterozygosity"
        )

    x_sites: set[int] = set()
    y_sites: set[int] = set()
    for site_id, (chrom, pos, _ref, _alt) in enumerate(matrix.index.keys):
        c = normalize_contig(chrom)
        if c == "X" and (build is None or regions.usable(c, pos)):
            x_sites.add(site_id)
        elif c == "Y" and (build is None or regions.usable(c, pos)):
            y_sites.add(site_id)

    raw: dict[str, tuple[int, int, int]] = {}
    for sample in matrix.samples:
        s = sets[sample]
        nonref = s.nonref
        x_nonref = nonref & x_sites
        n_x = len(x_nonref)
        n_het = len(s.het & x_sites)
        n_y = len(nonref & y_sites)
        raw[sample] = (n_x, n_het, n_y)

    fracs = [h / x for x, h, _ in raw.values() if x >= cfg.min_x_sites]
    boundary = _calibrate_sex_boundary(fracs, cfg, notes)

    evidence: dict[str, SexEvidence] = {}
    for sample, (n_x, n_het, n_y) in raw.items():
        if n_x < cfg.min_x_sites:
            evidence[sample] = SexEvidence(
                sample, n_x, n_het, float("nan"), float("nan"), float("nan"), n_y,
                Sex.UNKNOWN,
                f"only {n_x} usable chrX variant sites (need {cfg.min_x_sites})",
            )
            continue
        frac = n_het / n_x
        lo, hi = wilson_interval(n_het, n_x)
        if lo > boundary:
            inferred, why = Sex.FEMALE, "chrX heterozygosity above the cohort boundary"
        elif hi < boundary:
            inferred, why = Sex.MALE, "chrX heterozygosity below the cohort boundary"
        else:
            inferred, why = (
                Sex.UNKNOWN,
                "95% interval straddles the cohort boundary; sex is not determined",
            )
        evidence[sample] = SexEvidence(sample, n_x, n_het, frac, lo, hi, n_y, inferred, why)

    if any(e.n_y_sites for e in evidence.values()):
        notes.append(
            "chrY call counts are reported as context only. Females routinely carry a "
            "handful of chrY calls from X-Y homologous mismapping, so small counts "
            "carry no evidential weight and are not used in the sex call."
        )
    return evidence, boundary, notes


def _calibrate_sex_boundary(fracs: list[float], cfg: IdentityConfig, notes: list[str]) -> float:
    """Prefer the cohort's own bimodal gap; fall back to a documented constant."""
    from ..stats import largest_gap

    gap = largest_gap(fracs, 0.40, 0.58)
    if gap and (gap[1] - gap[0]) >= 0.08:
        boundary = (gap[0] + gap[1]) / 2.0
        notes.append(
            f"sex boundary calibrated on this cohort at chrX het {boundary:.3f} "
            f"(clean gap between {gap[0]:.3f} and {gap[1]:.3f})"
        )
        return boundary
    notes.append(
        f"no clean bimodal gap in cohort chrX heterozygosity; using the default "
        f"boundary of {cfg.default_sex_boundary:.2f}. Consanguinity and long runs of "
        f"homozygosity depress female chrX heterozygosity, so verify borderline calls."
    )
    return cfg.default_sex_boundary


def pair_evidence(matrix: GenotypeMatrix, sets: dict[str, SampleSets]) -> list[PairEvidence]:
    out: list[PairEvidence] = []
    names = matrix.samples
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            sa, sb = sets[a], sets[b]
            nr_a, nr_b = sa.nonref, sb.nonref
            shared = nr_a & nr_b
            union = len(nr_a) + len(nr_b) - len(shared)
            n_shared = len(shared)
            agree = len(sa.het & sb.het) + len(sa.homalt & sb.homalt)
            ev = PairEvidence(
                a=a,
                b=b,
                dense=matrix.pair_is_dense(a, b),
                n_shared=n_shared,
                n_union=union,
                jaccard=(n_shared / union) if union else float("nan"),
                agreement=(agree / n_shared) if n_shared else float("nan"),
                n_nonref_a=len(nr_a),
                n_nonref_b=len(nr_b),
            )
            if ev.dense:
                inform = sa.assessed & sb.assessed
                ev.n_informative = len(inform)
                n_het_i = len(sa.het & sb.assessed)
                n_het_j = len(sb.het & sa.assessed)
                n_hethet = len(sa.het & sb.het)
                n_ibs0 = len(sa.homalt & sb.homref) + len(sb.homalt & sa.homref)
                denom = n_het_i + n_het_j
                ev.king_phi = ((n_hethet - 2 * n_ibs0) / denom) if denom else None
                ev.ibs0_rate = (n_ibs0 / ev.n_informative) if ev.n_informative else None
            out.append(ev)
    return out


def _classify_duplicate(ev: PairEvidence, matrix: GenotypeMatrix) -> tuple[str, str, dict]:
    """Describe a duplicate pair by *evidence*, not by a guessed cause.

    The tool reports what is observable - byte identity, site-set identity,
    call-count ratio - and lists the mechanisms compatible with it.  It does
    not assert a wet-lab versus file-handling cause, because the same
    observation is produced by re-running one FASTQ through a different
    pipeline version as by preparing a second library.
    """
    scan_a = matrix.scans.get(ev.a)
    scan_b = matrix.scans.get(ev.b)
    same_file = matrix.source.get(ev.a) == matrix.source.get(ev.b)
    detail: dict = {
        "identical_site_sets": ev.jaccard >= 0.9999,
        "same_source_file": same_file,
    }
    if scan_a and scan_b and not same_file:
        detail |= {
            "file_bytes_a": scan_a.file_size,
            "file_bytes_b": scan_b.file_size,
            "file_sha256_equal": bool(
                scan_a.file_sha256 and scan_a.file_sha256 == scan_b.file_sha256
            ),
            "file_byte_size_delta": abs(scan_a.file_size - scan_b.file_size),
        }

    if ev.jaccard >= 0.9999 and ev.n_nonref_a == ev.n_nonref_b:
        if detail.get("file_sha256_equal"):
            return (
                "DUPLICATE_BYTE_IDENTICAL",
                "the two files are byte-for-byte identical",
                detail,
            )
        same_size = detail.get("file_byte_size_delta") == 0 and not same_file
        return (
            "DUPLICATE_CONTENT_IDENTICAL",
            "identical variant content"
            + (
                " in two files of exactly the same byte length but different content "
                "hashes, which is what a relabelled copy of one file looks like"
                if same_size
                else " in files that are not byte-identical"
            )
            + " (compatible with: one dataset delivered twice under two names, or the "
            "same input re-run through the same pipeline)",
            detail,
        )
    return (
        "DUPLICATE_SAME_INDIVIDUAL",
        "same individual, different variant content (compatible with: two sequencing "
        "runs of one DNA, or one dataset re-processed with different filters or a "
        "different pipeline version)",
        detail,
    )


def _calibrate_relatedness(
    pairs: list[PairEvidence],
    ped: Pedigree,
    duplicate_pairs: set[tuple[str, str]],
    cfg: IdentityConfig,
) -> tuple[float | None, dict, list[str]]:
    """Derive the related/unrelated Jaccard boundary from the cohort's own declarations.

    Medians and MAD, not means and SD: the declarations are exactly what may be
    wrong, so the calibration must tolerate a minority of mislabelled pairs.

    Every pair involving a sample caught in a duplicate is excluded, not just
    the duplicate pair itself.  When two samples are the same individual we do
    not know which of the two labels is wrong, so neither sample's declared
    relationships can be trusted to calibrate anything.
    """
    notes: list[str] = []
    tainted = {s for pair in duplicate_pairs for s in pair}
    unrel, first = [], []
    for ev in pairs:
        if (ev.a, ev.b) in duplicate_pairs or ev.a in tainted or ev.b in tainted:
            continue
        if ev.n_shared < cfg.min_shared_sites:
            continue
        rel = ped.relationship(ev.a, ev.b)
        if rel is Relationship.UNRELATED:
            unrel.append(ev.jaccard)
        elif rel in FIRST_DEGREE:
            first.append(ev.jaccard)

    stats = {
        "n_declared_unrelated_pairs": len(unrel),
        "n_declared_first_degree_pairs": len(first),
        "unrelated_median": median(unrel) if unrel else None,
        "unrelated_mad": mad(unrel) if len(unrel) > 1 else None,
        "first_degree_median": median(first) if first else None,
        "first_degree_mad": mad(first) if len(first) > 1 else None,
    }
    if len(unrel) < 3 or len(first) < 3:
        notes.append(
            "too few declared pairs in one or both classes to calibrate a relatedness "
            "boundary on this cohort; pedigree reconciliation reported as UNKNOWN"
        )
        return None, stats, notes

    m_u, m_f = median(unrel), median(first)
    spread = (mad(unrel) or 0.0) + (mad(first) or 0.0)
    separation = (m_f - m_u) / spread if spread > 0 else float("inf")
    stats["separation"] = separation
    if m_f <= m_u or separation < cfg.min_band_separation:
        notes.append(
            f"declared unrelated (median J={m_u:.3f}) and first-degree (median "
            f"J={m_f:.3f}) pairs are not separable on this cohort "
            f"(separation {separation:.1f} < {cfg.min_band_separation}); pedigree "
            f"reconciliation reported as UNKNOWN rather than guessed"
        )
        return None, stats, notes

    boundary = (m_u + m_f) / 2.0
    stats["boundary"] = boundary
    notes.append(
        f"relatedness boundary calibrated on this cohort at Jaccard {boundary:.3f} "
        f"(declared-unrelated median {m_u:.3f}, declared-first-degree median {m_f:.3f})"
    )
    return boundary, stats, notes


@never_raises(CHECK)
def check_identity(
    matrix: GenotypeMatrix,
    ped: Pedigree,
    cfg: IdentityConfig | None = None,
    build: str | None = None,
) -> CheckResult:
    cfg = cfg or IdentityConfig()
    findings: list[Finding] = []
    notes: list[str] = []

    if not matrix.samples:
        from ..model import unknown

        return unknown(CHECK, "no samples were read from the supplied VCFs")

    if build is None:
        builds = {s.build for s in matrix.scans.values() if s.build}
        if len(builds) == 1:
            build = builds.pop()
        elif len(builds) > 1:
            findings.append(
                Finding(
                    code="BUILD_CONFLICT",
                    severity=Severity.ERROR,
                    message=f"VCFs declare more than one genome build: {sorted(builds)}",
                    evidence={
                        "builds": sorted(builds),
                        "next_steps": ["re-align or lift over so all samples share one build"],
                    },
                )
            )

    sets = build_sets(matrix)
    sex_ev, sex_boundary, sex_notes = sex_evidence(matrix, sets, cfg, build)
    notes += sex_notes
    pairs = pair_evidence(matrix, sets)

    # --- sample / pedigree membership -------------------------------------
    in_ped = set(ped.sample_ids)
    in_vcf = set(matrix.samples)
    if missing := sorted(in_vcf - in_ped):
        findings.append(
            Finding(
                code="SAMPLE_NOT_IN_PED",
                severity=Severity.WARN,
                message=f"{len(missing)} sample(s) present in the VCFs but absent from the PED",
                subjects=missing,
                evidence={"samples": missing},
            )
        )
    if absent := sorted(in_ped - in_vcf):
        findings.append(
            Finding(
                code="PED_SAMPLE_NOT_IN_VCF",
                severity=Severity.WARN,
                message=f"{len(absent)} pedigree member(s) have no genotypes",
                subjects=absent,
                evidence={"samples": absent},
            )
        )
    for w in ped.warnings:
        notes.append(f"PED: {w}")

    # --- sex ---------------------------------------------------------------
    n_sex_mismatch = 0
    for sample, ev in sex_ev.items():
        declared = ped.individuals[sample].sex if sample in ped.individuals else Sex.UNKNOWN
        if ev.inferred is Sex.UNKNOWN or declared is Sex.UNKNOWN:
            if ev.inferred is Sex.UNKNOWN and ev.n_x_sites < cfg.min_x_sites:
                findings.append(
                    Finding(
                        code="SEX_NOT_DETERMINED",
                        severity=Severity.INFO,
                        message=f"{sample}: {ev.reason}",
                        subjects=[sample],
                        evidence=ev.to_dict(),
                    )
                )
            continue
        if declared is not ev.inferred:
            n_sex_mismatch += 1
            findings.append(
                Finding(
                    code="SEX_MISMATCH",
                    severity=Severity.BLOCKING,
                    message=(
                        f"{sample}: pedigree says {declared.value}, genotypes say "
                        f"{ev.inferred.value} (chrX het {ev.x_het_frac:.1%} over "
                        f"{ev.n_x_sites} sites, 95% CI "
                        f"{ev.ci_low:.1%}-{ev.ci_high:.1%}; boundary {sex_boundary:.1%})"
                    ),
                    subjects=[sample],
                    evidence=ev.to_dict()
                    | {
                        "declared_sex": declared.value,
                        "boundary": sex_boundary,
                        "do_not_conclude": (
                            "any result that assumes the declared sex of this sample"
                        ),
                        "next_steps": ["resolve sample identity (blocking)"],
                    },
                )
            )

    # --- duplicates, cohort-wide ------------------------------------------
    duplicate_pairs: set[tuple[str, str]] = set()
    for ev in pairs:
        if ev.n_shared < cfg.min_shared_sites:
            continue
        if ev.agreement >= cfg.dup_agreement and ev.jaccard >= cfg.dup_jaccard:
            duplicate_pairs.add((ev.a, ev.b))
            code, explanation, detail = _classify_duplicate(ev, matrix)
            declared = ped.relationship(ev.a, ev.b)
            cross_family = declared is Relationship.UNRELATED
            findings.append(
                Finding(
                    code=code,
                    severity=Severity.BLOCKING,
                    message=(
                        f"{ev.a} and {ev.b} are the same individual "
                        f"(Jaccard {ev.jaccard:.4f}, agreement {ev.agreement:.4f} over "
                        f"{ev.n_shared} shared sites) but the pedigree declares them "
                        f"{declared.value}"
                        + (" ACROSS FAMILIES" if cross_family else "")
                        + f" - {explanation}"
                    ),
                    subjects=[ev.a, ev.b],
                    evidence=ev.to_dict()
                    | detail
                    | {
                        "declared_relationship": declared.value,
                        "cross_family": cross_family,
                        "do_not_conclude": "anything that treats these as two independent samples",
                        "next_steps": [
                            "resolve sample identity (blocking)",
                            "re-derive the affected/unaffected counts after the identity "
                            "of every sample is settled",
                        ],
                    },
                )
            )

    # --- pedigree reconciliation ------------------------------------------
    dense_pairs = [ev for ev in pairs if ev.dense and ev.king_phi is not None]
    boundary, calib, calib_notes = _calibrate_relatedness(pairs, ped, duplicate_pairs, cfg)
    notes += calib_notes
    n_ped_mismatch = 0

    if dense_pairs:
        median_sites = sorted(e.n_informative for e in dense_pairs)[len(dense_pairs) // 2]
        notes.append(
            f"{len(dense_pairs)} of {len(pairs)} pairs were jointly genotyped, so "
            f"KING-robust kinship and IBS0 were computed for them "
            f"(median {median_sites:,} informative sites per pair)"
        )
        if median_sites < 50_000:
            notes.append(
                f"at {median_sites:,} informative sites per pair the kinship estimator "
                f"cannot reliably separate first- from second-degree relatives for an "
                f"individual pair, so only gross discrepancies are reported as errors; "
                f"KING is usually run on >=100k common variants"
            )
        # Compare the observed kinship to the kinship the PEDIGREE implies, not
        # to a band label derived from a coarse relationship enum.  Two reasons.
        # First, the enum cannot name a grandparent or a spouse, so it could
        # never catch either being mislabelled.  Second, comparing band LABELS
        # is brittle at the edges: validated against the public CEPH 1463
        # pedigree, true full sibs came out at phi 0.153 and 0.176 - just under
        # the 0.177 first-degree boundary - and a label comparison called both
        # of them pedigree errors.  Realised IBD genuinely varies between full
        # sibs, and at ~20k sites the estimator is not sharp enough to resolve
        # first from second degree for an individual pair.
        #
        # So the rule is deliberately conservative: flag only discrepancies too
        # large to be noise, and stay silent on first-vs-second degree.
        for ev in dense_pairs:
            if (ev.a, ev.b) in duplicate_pairs:
                continue
            if ev.n_informative < cfg.min_informative_sites:
                continue
            expected = ped.kinship(ev.a, ev.b)
            if expected is None:
                continue
            observed = ev.king_phi
            declared = ped.relationship(ev.a, ev.b)

            close_but_looks_distant = (
                expected >= cfg.first_degree_min and observed < cfg.second_degree_min
            )
            distant_but_looks_close = (
                expected <= cfg.unrelated_max and observed >= cfg.first_degree_min
            )
            if not (close_but_looks_distant or distant_but_looks_close):
                continue

            n_ped_mismatch += 1
            findings.append(
                Finding(
                    code="PEDIGREE_MISMATCH",
                    severity=Severity.BLOCKING,
                    message=(
                        f"{ev.a}/{ev.b}: the pedigree implies kinship {expected:.4f} "
                        f"({declared.value}), genotypes give KING phi={observed:.4f} "
                        f"-> {ev.king_band}"
                    ),
                    subjects=[ev.a, ev.b],
                    evidence=ev.to_dict()
                    | {
                        "declared_relationship": declared.value,
                        "expected_kinship": expected,
                        "observed_kinship": observed,
                        "do_not_conclude": "any segregation result computed on this pedigree",
                        "next_steps": ["resolve sample identity (blocking)"],
                    },
                )
            )
    elif boundary is not None:
        for ev in pairs:
            if (ev.a, ev.b) in duplicate_pairs or ev.n_shared < cfg.min_shared_sites:
                continue
            declared = ped.relationship(ev.a, ev.b)
            if declared in AMBIGUOUS_DECLARATIONS:
                continue
            observed_related = ev.jaccard >= boundary
            if declared is Relationship.UNRELATED and observed_related:
                n_ped_mismatch += 1
                findings.append(
                    Finding(
                        code="PEDIGREE_MISMATCH",
                        severity=Severity.BLOCKING,
                        message=(
                            f"{ev.a}/{ev.b}: declared unrelated (different families) but "
                            f"Jaccard {ev.jaccard:.4f} sits above the cohort boundary "
                            f"{boundary:.4f}"
                        ),
                        subjects=[ev.a, ev.b],
                        evidence=ev.to_dict()
                        | {
                            "declared_relationship": declared.value,
                            "boundary": boundary,
                            "do_not_conclude": "any segregation result computed on this pedigree",
                            "next_steps": ["resolve sample identity (blocking)"],
                        },
                    )
                )
            elif declared in FIRST_DEGREE and not observed_related:
                n_ped_mismatch += 1
                findings.append(
                    Finding(
                        code="PEDIGREE_MISMATCH",
                        severity=Severity.BLOCKING,
                        message=(
                            f"{ev.a}/{ev.b}: declared {declared.value} but Jaccard "
                            f"{ev.jaccard:.4f} is below the cohort boundary {boundary:.4f}"
                        ),
                        subjects=[ev.a, ev.b],
                        evidence=ev.to_dict()
                        | {
                            "declared_relationship": declared.value,
                            "boundary": boundary,
                            "do_not_conclude": "any segregation result computed on this pedigree",
                            "next_steps": ["resolve sample identity (blocking)"],
                        },
                    )
                )

    if not dense_pairs:
        shared_file = len({matrix.source.get(s) for s in matrix.samples}) == 1
        if shared_file and len(matrix.samples) > 1:
            # The trap: a file that looks jointly called but carries no hom-ref.
            why = (
                "all samples are in ONE multi-sample file, but it contains no "
                "hom-reference genotypes - every non-carrier is './.'. So the file "
                "looks jointly called and is not: hom-reference is unobserved, and "
                "IBS0 depends on telling it apart from missing data"
            )
        else:
            why = (
                "samples came from separate single-sample VCFs, so hom-reference "
                "genotypes are unobserved"
            )
        notes.append(
            f"{why}. KING-robust kinship and IBS0 therefore cannot be computed, and "
            f"relationship DEGREE is not resolved - only related vs unrelated. "
            f"Re-merge from gVCFs retaining hom-reference calls to get degree."
        )

    if findings and any(f.severity is Severity.BLOCKING for f in findings):
        status = Status.FAIL
        parts = []
        if duplicate_pairs:
            parts.append(f"{len(duplicate_pairs)} duplicate pair(s)")
        if n_sex_mismatch:
            parts.append(f"{n_sex_mismatch} sex mismatch(es)")
        if n_ped_mismatch:
            parts.append(f"{n_ped_mismatch} pedigree inconsistencies")
        summary = ", ".join(parts) or "samples do not match the pedigree"
    elif any(f.severity in (Severity.WARN, Severity.ERROR) for f in findings):
        status = Status.WARN
        summary = "pedigree broadly consistent, with caveats"
    elif boundary is None and not dense_pairs:
        status = Status.UNKNOWN
        summary = "could not calibrate a relatedness boundary on this cohort"
    else:
        status = Status.PASS
        summary = "genotypes are consistent with the declared pedigree"

    return CheckResult(
        check=CHECK,
        status=status,
        summary=summary,
        findings=findings,
        metrics={
            "n_samples": len(matrix.samples),
            "n_sites": matrix.n_sites,
            "genome_build": build,
            "n_sex_mismatch": n_sex_mismatch,
            "n_pedigree_mismatch": n_ped_mismatch,
            "n_duplicate_pairs": len(duplicate_pairs),
            "sex_boundary": sex_boundary,
            "relatedness_calibration": calib,
            "sex": [e.to_dict() for e in sex_ev.values()],
            "pairs": [e.to_dict() for e in pairs],
        },
        notes=notes,
    )
