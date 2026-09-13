"""Check 3 - content provenance.  Is this file actually a whole exome?

Everything here is derived from the user's own file.  Region classes are read
from annotation that is *already in the VCF* (ANNOVAR's ``Func.refGene`` today);
nothing is fetched, and no annotation database is bundled.  If the VCF carries
no region annotation the coding-only question is reported UNKNOWN rather than
guessed, because there is no annotation-free way to distinguish "whole exome"
from "exome filtered to coding" - both are sparse, clustered and gene-shaped.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..intervals import IntervalSet, read_bed
from ..model import CheckResult, Finding, Severity, Status, never_raises
from ..vcfio import VcfScan

CHECK = "provenance"

CODING_CLASSES = ("exonic", "splicing", "exonic;splicing")
INTRONIC_CLASSES = ("intronic", "ncRNA_intronic")


@dataclass
class ProvenanceConfig:
    min_wes_variants: int = 40_000
    coding_fraction_alarm: float = 0.85
    intronic_fraction_alarm: float = 0.10
    # Below this on-target fraction, the caller was not restricted to the kit.
    on_target_alarm: float = 0.75


@never_raises(CHECK)
def check_provenance(
    scans: list[VcfScan],
    cfg: ProvenanceConfig | None = None,
    target_bed: str | None = None,
) -> CheckResult:
    cfg = cfg or ProvenanceConfig()
    findings: list[Finding] = []
    notes: list[str] = []
    verdicts: set[str] = set()
    per_file: list[dict] = []

    if not scans:
        from ..model import unknown

        return unknown(CHECK, "no VCFs were supplied")

    target: IntervalSet | None = None
    if target_bed:
        target = read_bed(target_bed).intervals or None
    on_target = _on_target_fractions(scans, target) if target else {}

    for scan in scans:
        st = scan.stats
        total = st.n_records
        called = st.n_pass + st.n_nonpass
        nonpass_frac = (st.n_nonpass / called) if called else None
        region_total = sum(st.region_counts.values())
        coding = sum(st.region_counts.get(c, 0) for c in CODING_CLASSES)
        intronic = sum(st.region_counts.get(c, 0) for c in INTRONIC_CLASSES)
        coding_frac = (coding / region_total) if region_total else None
        intronic_frac = (intronic / region_total) if region_total else None

        file_verdicts: list[str] = []
        if total and total < cfg.min_wes_variants:
            file_verdicts.append("SUBSET")
        if nonpass_frac == 0.0 and called:
            file_verdicts.append("PASS_FILTERED")
        if region_total and coding_frac is not None and intronic_frac is not None:
            if (
                coding_frac > cfg.coding_fraction_alarm
                and intronic_frac < cfg.intronic_fraction_alarm
            ):
                file_verdicts.append("CODING_ONLY")
        if not st.has_format_column or not scan.header.samples:
            file_verdicts.append("METRICS_STRIPPED")
        elif not (st.format_keys_seen - {"GT"}):
            file_verdicts.append("METRICS_STRIPPED")
        multisample = len(scan.header.samples) > 1
        if multisample and not any(scan.has_explicit_homref.values()):
            file_verdicts.append("HOMREF_STRIPPED")

        verdict = "WHOLE_EXOME" if not file_verdicts else "+".join(sorted(set(file_verdicts)))
        verdicts.add(verdict)
        per_file.append(
            {
                "path": scan.path,
                "verdict": verdict,
                "n_records": total,
                "n_reference_blocks": st.n_reference_blocks,
                "nonpass_fraction": nonpass_frac,
                "coding_fraction": coding_frac,
                "intronic_fraction": intronic_frac,
                "annotation_source": st.annotation_source,
                "format_keys": sorted(st.format_keys_seen),
                "n_malformed_lines": st.n_malformed,
            }
        )

        frac_on = on_target.get(scan.path)
        per_file[-1]["on_target_fraction"] = frac_on
        if frac_on is not None and frac_on < cfg.on_target_alarm:
            file_verdicts.append("NOT_INTERVAL_RESTRICTED")
            per_file[-1]["verdict"] = verdict = "+".join(sorted(set(file_verdicts)))
            verdicts.discard(
                "+".join(sorted({v for v in file_verdicts if v != "NOT_INTERVAL_RESTRICTED"}))
            )
            verdicts.add(verdict)

        name = scan.path
        if "NOT_INTERVAL_RESTRICTED" in file_verdicts:
            findings.append(
                Finding(
                    code="NOT_INTERVAL_RESTRICTED",
                    severity=Severity.WARN,
                    message=(
                        f"{name}: only {frac_on:.1%} of calls fall inside the capture "
                        f"target, so the caller was not restricted to the kit. The "
                        f"off-target majority is shallow by construction and will "
                        f"dominate any quality statistic computed over the whole file"
                    ),
                    evidence={
                        "on_target_fraction": frac_on,
                        "next_steps": [
                            "restrict quality statistics to the capture target, or "
                            "re-call with -L targets.bed"
                        ],
                    },
                )
            )
        if "SUBSET" in file_verdicts:
            findings.append(
                Finding(
                    code="SUBSET_SUSPECTED",
                    severity=Severity.WARN,
                    message=(
                        f"{name}: {total} variants, well under the ~{cfg.min_wes_variants} "
                        f"expected for a whole exome - this file was filtered or subset "
                        f"before you received it"
                    ),
                    evidence={
                        "n_records": total,
                        "expected_min": cfg.min_wes_variants,
                        "do_not_conclude": "no monogenic cause",
                        "next_steps": ["obtain or re-call the unfiltered call set"],
                    },
                )
            )
        if "CODING_ONLY" in file_verdicts:
            findings.append(
                Finding(
                    code="CODING_ONLY",
                    severity=Severity.WARN,
                    message=(
                        f"{name}: {coding_frac:.1%} coding and only {intronic_frac:.1%} "
                        f"intronic - this is a coding-only call set, not a whole exome"
                    ),
                    evidence={
                        "coding_fraction": coding_frac,
                        "intronic_fraction": intronic_frac,
                        "do_not_conclude": "no monogenic cause",
                        "next_steps": ["re-call without the coding-only interval file"],
                    },
                )
            )
        if "PASS_FILTERED" in file_verdicts:
            findings.append(
                Finding(
                    code="PASS_FILTERED",
                    severity=Severity.WARN,
                    message=(
                        f"{name}: 0.0% of {called} variants are non-PASS - the filter "
                        f"column carries no information, so the file was pre-filtered "
                        f"and you cannot assess what was discarded"
                    ),
                    evidence={
                        "n_pass": st.n_pass,
                        "n_nonpass": st.n_nonpass,
                        "next_steps": ["obtain the pre-filter call set to see what was removed"],
                    },
                )
            )
        if "HOMREF_STRIPPED" in file_verdicts:
            findings.append(
                Finding(
                    code="HOMREF_STRIPPED",
                    severity=Severity.WARN,
                    message=(
                        f"{name}: {len(scan.header.samples)} samples in one file, but not "
                        f"a single hom-reference genotype anywhere - every non-carrier is "
                        f"written as './.'. The file looks jointly called and is not: "
                        f"'./.' now means both 'hom-reference' and 'never covered', and "
                        f"nothing can tell them apart"
                    ),
                    evidence={
                        "n_samples": len(scan.header.samples),
                        "do_not_conclude": (
                            "that a sample lacking a variant here was tested and found "
                            "negative"
                        ),
                        "next_steps": [
                            "re-merge from the gVCFs retaining hom-reference calls, so "
                            "absence of a variant can be distinguished from absence of data"
                        ],
                    },
                )
            )
        if "METRICS_STRIPPED" in file_verdicts:
            findings.append(
                Finding(
                    code="METRICS_STRIPPED",
                    severity=Severity.ERROR,
                    message=(
                        f"{name}: no per-genotype metrics (FORMAT keys: "
                        f"{sorted(st.format_keys_seen) or 'none'}) - genotype "
                        f"confidence cannot be assessed at all"
                    ),
                    evidence={
                        "format_keys": sorted(st.format_keys_seen),
                        "do_not_conclude": "that any individual genotype is reliable",
                        "next_steps": ["obtain a VCF retaining GT:AD:DP:GQ:PL"],
                    },
                )
            )
        if st.n_malformed:
            notes.append(f"{name}: {st.n_malformed} unparseable line(s) skipped")
        for p in scan.problems[:5]:
            notes.append(f"{name}: {p}")
        if not st.region_counts:
            notes.append(
                f"{name}: no region annotation in the VCF, so the coding-only question "
                f"is UNKNOWN for this file"
            )

    status = Status.PASS if verdicts == {"WHOLE_EXOME"} else Status.WARN
    summary = (
        "whole exome"
        if status is Status.PASS
        else ", ".join(sorted(v for v in verdicts if v != "WHOLE_EXOME")).replace("+", " + ")
    )
    return CheckResult(
        check=CHECK,
        status=status,
        summary=summary,
        findings=findings,
        metrics={"files": per_file, "verdicts": sorted(verdicts)},
        notes=notes,
    )


def _on_target_fractions(
    scans: list[VcfScan], target: IntervalSet
) -> dict[str, float]:
    """Fraction of each file's records that fall inside the capture target.

    A low value means the caller was never restricted to the kit. That matters
    far more than it sounds: off-target calls are shallow by construction, so
    they dominate any depth or genotype-quality statistic computed over the
    whole file and make a perfectly adequate exome look like a failure: a depth
    distribution computed over an unrestricted call set describes the off-target
    reads, not the experiment. Restricting the same statistic to the target can
    change the verdict entirely, which is why this is reported rather than left
    for the reader to infer.
    """
    from ..annovar import iter_multianno_records, looks_like_multianno
    from ..vcfio import iter_records

    def contains(chrom: str, pos: int) -> bool:
        spans = target.get(chrom)
        if not spans:
            return False
        lo, hi = 0, len(spans) - 1
        while lo <= hi:
            mid = (lo + hi) // 2
            start, end = spans[mid]
            if pos < start:
                hi = mid - 1
            elif pos >= end:
                lo = mid + 1
            else:
                return True
        return False

    out: dict[str, float] = {}
    for scan in scans:
        if scan.path in out:
            continue
        reader = (
            iter_multianno_records(scan.path)
            if looks_like_multianno(scan.path)
            else iter_records(scan.path)
        )
        n = k = 0
        for rec in reader:
            n += 1
            if contains(rec.chrom, rec.pos - 1):  # BED is 0-based half-open
                k += 1
        if n:
            out[scan.path] = k / n
    return out
