"""Check 2 - genotype confidence.

Which individual genotypes are not supported by their own evidence?

The headline is the false-homozygote rule, because a ``1/1`` call is what puts a
variant into a recessive model, and a ``1/1`` call from two reads is the cheapest
way to lose a year.  Three things are asked of every homozygous-alternate call,
and they are kept apart because they have different causes and different fixes:

``LOW_DEPTH_HOM``
    There were not enough reads for "homozygous" to mean anything.  At DP=2 a
    true heterozygote yields two alt reads one time in four; the call is a
    coin-flip dressed as a genotype.

``HOM_CONTRADICTED_BY_LIKELIHOOD``
    The caller's own likelihoods do not prefer hom-alt over het by any margin
    worth having.  Read from ``PL``/``GL``, so it works for GATK, DeepVariant,
    bcftools, DRAGEN and Strelka alike.

``ALLELE_IMBALANCE_HOM``
    Depth was adequate and the call still is not clean - reference reads are
    present at a locus called homozygous alternate.  That is contamination, a
    paralogous mapping, or a copy-number event, and re-sequencing deeper will
    not fix it.  Merging this with LOW_DEPTH_HOM would hide a different problem
    behind the same label.

``MLEAC_CONTRADICTS_GT`` is a **fallback**, not a corroborator, and the reason is
structural rather than empirical: for a biallelic hom-alt call ``PL[het]`` *is*
GQ, by the definition of GQ as the difference between the best and second-best
genotype likelihoods.  So "the likelihoods contradict the genotype" and "GQ is
low" are the same statement, and MLEAC is a noisy third way of saying it.  On
GATK call sets it fires on a large fraction of all homozygous calls, and a flag
that fires on half the data is not a flag.

MLEAC is therefore only allowed to condemn a genotype when neither GQ nor PL is
available - which is exactly the case for a delivered analysis whose FORMAT block
has been stripped, and there it is the only surviving witness to genotype
uncertainty.  ``MLEAC_DISCORDANT`` is recorded everywhere regardless, so
``rule_overlap`` in the metrics re-measures the redundancy on whatever dataset is
in front of it instead of asking you to trust the paragraph above.

Allele balance is tested with an exact binomial against p=0.5 rather than a fixed
0.30-0.70 window, because a fixed window mis-scales with depth: 3/10 is
unremarkable for a true heterozygote and 60/200 is not, and both sit at 0.30.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field

from ..contigs import is_primary
from ..model import CheckResult, Finding, Severity, Status, never_raises
from ..stats import binom_two_sided_p
from ..vcfio import iter_records, unique_label

CHECK = "genotype"


@dataclass
class GenotypeConfig:
    low_depth: int = 10
    uncertain_zygosity: int = 20
    hom_min_depth: int = 8
    hom_min_allele_balance: float = 0.90
    min_pl_het_margin: int = 20
    min_gq: int = 20
    ab_alpha: float = 0.01
    min_qd: float = 2.0
    min_mq: float = 40.0
    max_fs: float = 60.0
    max_sor: float = 3.0
    cluster_window: int = 50
    cluster_min: int = 2
    max_examples: int = 40


@dataclass
class SampleTally:
    n_genotypes: int = 0
    n_het: int = 0
    n_hom_alt: int = 0
    flags: dict[str, int] = field(default_factory=dict)

    def bump(self, flag: str) -> None:
        self.flags[flag] = self.flags.get(flag, 0) + 1


def _num(value: str | None) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _allele_counts(
    cell: dict[str, str], info: dict[str, str]
) -> tuple[list[float] | None, str | None]:
    """Per-allele read counts, from whichever field the caller wrote them in.

    Returned **per allele** - index 0 is the reference, index *i* the *i*-th ALT -
    rather than as a ref/alt pair, because at a multiallelic site the pair is not
    a well-defined thing.  Summing every ALT and calling the result "alt support"
    says a ``1/1`` call backed by 5 reads out of 30 is in perfect balance, when
    the other 25 reads support a different allele entirely.

    ``AD`` is GATK's spelling and is far from universal: freebayes writes
    ``RO``/``AO``, and older samtools/bcftools pipelines write ``DP4``.  Reading
    only ``AD`` meant the allele-balance arm of this check silently did nothing
    on freebayes output - no ``AB_SKEW``, no ``ALLELE_IMBALANCE_HOM`` - while the
    report still looked complete.  Found by running the check on the public CEPH
    1463 call set.

    VarScan's single-valued ``AD`` (alt count only, with ``RD`` for reference) is
    deliberately not read as GATK's: the second element of the return value names
    the field actually used, so the source is auditable rather than guessed at.
    ``DP4`` carries no per-allele breakdown at all, so it yields two entries and
    is therefore only usable at biallelic sites.
    """
    ad = cell.get("AD", "")
    if ad:
        parts = [_num(p) for p in ad.split(",")]
        if len(parts) >= 2 and all(p is not None for p in parts):
            return [float(p) for p in parts], "AD"  # type: ignore[arg-type]

    ro, ao = _num(cell.get("RO")), cell.get("AO", "")
    if ro is not None and ao:
        alts = [_num(p) for p in ao.split(",")]
        if all(a is not None for a in alts):
            return [ro] + [float(a) for a in alts], "RO/AO"  # type: ignore[arg-type]

    dp4 = cell.get("DP4") or info.get("DP4", "")
    if dp4:
        parts = [_num(p) for p in dp4.split(",")]
        if len(parts) == 4 and all(p is not None for p in parts):
            return [parts[0] + parts[1], parts[2] + parts[3]], "DP4"  # type: ignore[operator]

    return None, None


def _pl_list(cell: dict[str, str]) -> list[int] | None:
    raw = cell.get("PL")
    scale = 1.0
    if not raw:
        # GL is log10-likelihood (negative); PL is -10*log10, normalised to 0.
        raw, scale = cell.get("GL"), 10.0
    if not raw:
        return None
    out: list[float] = []
    for part in raw.split(","):
        v = _num(part)
        if v is None:
            return None
        out.append(abs(v) * scale)
    if not out:
        return None
    floor = min(out)
    return [int(round(v - floor)) for v in out]


def evaluate_genotype(
    gt: str,
    cell: dict[str, str],
    info: dict[str, str],
    cfg: GenotypeConfig,
    single_sample: bool,
) -> tuple[list[str], dict]:
    """Flags raised by one genotype, plus the evidence behind them."""
    flags: list[str] = []
    ev: dict = {}

    alleles = [a for a in gt.replace("|", "/").split("/") if a.isdigit()]
    if not alleles:
        return ["NO_CALL"], ev
    is_het = len(set(alleles)) > 1
    is_hom_alt = not is_het and alleles[0] != "0"

    # INFO/DP is the depth summed over every sample at the site.  Using it as a
    # per-sample depth in a multi-sample VCF multiplies each sample's apparent
    # depth by the cohort size, which turns shallow homozygotes into
    # "adequate depth" ones - the precise failure this check exists to catch.
    # It is only a legitimate proxy when there is exactly one sample.
    dp = _num(cell.get("DP"))
    if dp is None and single_sample:
        dp = _num(info.get("DP"))
    gq = _num(cell.get("GQ"))
    counts, ad_source = _allele_counts(cell, info if single_sample else {})
    idx = [int(a) for a in alleles]
    have_counts = counts is not None and all(i < len(counts) for i in idx)
    # Balance is computed against the allele that was actually called, so a
    # multiallelic hom-alt backed by a minority of reads cannot read as clean.
    total = sum(counts) if have_counts else None  # type: ignore[arg-type]
    ab = (counts[idx[0]] / total) if (have_counts and total) else None  # type: ignore[index]
    pl = _pl_list(cell)
    ev |= {
        "DP": dp,
        "GQ": gq,
        "AD": cell.get("AD") or None,
        "allele_depth_source": ad_source,
        "allele_balance": ab,
    }

    if dp is not None:
        if dp < cfg.low_depth:
            flags.append("LOW_DEPTH")
        if dp < cfg.uncertain_zygosity:
            flags.append("UNCERTAIN_ZYGOSITY")
    if gq is not None and gq < cfg.min_gq:
        flags.append("LOW_GQ")

    if is_het and have_counts and total and total >= cfg.low_depth:
        # Only the two called alleles enter the test; reads supporting a third
        # allele are evidence about the site, not about this genotype's balance.
        a, b = counts[idx[0]], counts[idx[1]]  # type: ignore[index]
        pair = a + b
        p = binom_two_sided_p(int(b), int(pair)) if pair >= cfg.low_depth else 1.0
        ev["allele_balance_p"] = p
        if p < cfg.ab_alpha:
            flags.append("AB_SKEW")

    if is_hom_alt:
        # Silence is not evidence of adequacy.  A homozygous call with neither a
        # depth nor a likelihood to examine has not been checked, and must not be
        # tallied alongside the ones that passed a check.
        if dp is None and not _pl_list(cell) and ab is None:
            flags.append("HOM_NOT_ASSESSABLE")
        pl_het = None
        if pl and len(pl) >= 2:
            pl_het = pl[1]
            ev["PL_het"] = pl_het
        if dp is not None and dp < cfg.hom_min_depth:
            flags.append("LOW_DEPTH_HOM")
        if pl_het is not None and pl_het < cfg.min_pl_het_margin:
            flags.append("HOM_CONTRADICTED_BY_LIKELIHOOD")
        if (
            dp is not None
            and dp >= cfg.hom_min_depth
            and ab is not None
            and ab < cfg.hom_min_allele_balance
        ):
            flags.append("ALLELE_IMBALANCE_HOM")
        if single_sample:
            mleac = _num((info.get("MLEAC") or "").split(",")[0])
            ac = _num((info.get("AC") or "").split(",")[0])
            ev |= {"MLEAC": mleac, "AC": ac}
            if mleac is not None and ac is not None and mleac != ac:
                # Always recorded, so redundancy can be measured on any dataset.
                flags.append("MLEAC_DISCORDANT")
                # Only allowed to condemn a genotype when there is no direct
                # evidence to condemn it with.  Where GQ or PL exist, MLEAC is
                # almost perfectly redundant with them and fires on a large
                # fraction of all homozygous calls, so acting on it there would
                # bury the real signal.
                if gq is None and pl is None:
                    flags.append("MLEAC_CONTRADICTS_GT")
        if {
            "LOW_DEPTH_HOM",
            "HOM_CONTRADICTED_BY_LIKELIHOOD",
            "ALLELE_IMBALANCE_HOM",
            "MLEAC_CONTRADICTS_GT",
        } & set(flags):
            flags.append("FALSE_HOM_SUSPECT")

    return flags, ev


def _site_flags(info: dict[str, str], chrom: str, cfg: GenotypeConfig) -> list[str]:
    flags = []
    for key, limit, worse in (
        ("QD", cfg.min_qd, "lt"),
        ("MQ", cfg.min_mq, "lt"),
        ("FS", cfg.max_fs, "gt"),
        ("SOR", cfg.max_sor, "gt"),
    ):
        v = _num(info.get(key))
        if v is None:
            continue
        if (worse == "lt" and v < limit) or (worse == "gt" and v > limit):
            flags.append(f"SITE_{key}")
    dup = info.get("genomicSuperDups")
    if dup and dup not in (".", ""):
        flags.append("SEGMENTAL_DUPLICATION")
    if not is_primary(chrom):
        flags.append("NONPRIMARY_CONTIG")
    return flags


# The three rules that stand on their own, against which MLEAC is measured.
INDEPENDENT_HOM_RULES = ("LOW_DEPTH_HOM", "HOM_CONTRADICTED_BY_LIKELIHOOD", "ALLELE_IMBALANCE_HOM")
CONVENTIONAL_QC = ("LOW_DEPTH", "LOW_GQ")


@never_raises(CHECK)
def check_genotype(
    paths: list[str | os.PathLike[str]], cfg: GenotypeConfig | None = None
) -> CheckResult:
    cfg = cfg or GenotypeConfig()
    if not paths:
        from ..model import unknown

        return unknown(CHECK, "no files were supplied")

    from ..annovar import iter_multianno_records, looks_like_multianno

    tallies: dict[str, SampleTally] = {}
    examples: list[dict] = []
    n_records = 0
    n_with_format = 0
    files_without_format: list[str] = []
    partial_files: list[str] = []
    # Filtering state: a call set whose FILTER column is empty everywhere has
    # never been variant-quality filtered at all, which changes what every other
    # number in this check means.
    filter_seen = {"absent": 0, "pass": 0, "nonpass": 0}
    # Homozygotes deep enough that depth is not the explanation.
    deep_hom = {"n": 0, "flagged": 0, "depth_unknown": 0}
    # MLEAC redundancy measurement
    overlap = {"mleac_flagged": 0, "mleac_also_caught_by_rules": 0, "mleac_also_caught_by_qc": 0}
    # clustering state, per sample
    recent: dict[str, list[tuple[str, int]]] = {}
    cluster_hits: dict[str, int] = {}

    labels: dict[tuple[str, str], str] = {}

    for path in paths:
        saw_format = False
        saw_metrics = False
        reader = (
            iter_multianno_records(path) if looks_like_multianno(path) else iter_records(path)
        )
        for rec in reader:
            n_records += 1
            if not rec.filters:
                filter_seen["absent"] += 1
            elif rec.filters == ["PASS"]:
                filter_seen["pass"] += 1
            else:
                filter_seen["nonpass"] += 1
            if not rec.samples:
                continue
            saw_format = True
            single_sample = len(rec.samples) == 1
            site = _site_flags(rec.info, rec.chrom, cfg)
            for name, cell in rec.samples.items():
                gt = cell.get("GT", "")
                # Pipelines that call every sample "patient" are common; the label
                # has to distinguish them or two samples silently become one tally.
                lk = (str(path), name)
                if lk not in labels:
                    labels[lk] = unique_label(name, path, set(tallies))
                label = labels[lk]
                tally = tallies.setdefault(label, SampleTally())
                if set(cell) - {"GT"}:
                    saw_metrics = True
                flags, ev = evaluate_genotype(gt, cell, rec.info, cfg, single_sample)
                if "NO_CALL" in flags:
                    continue
                tally.n_genotypes += 1
                alleles = [a for a in gt.replace("|", "/").split("/") if a.isdigit()]
                is_hom_alt = len(set(alleles)) == 1 and alleles and alleles[0] != "0"
                if len(set(alleles)) > 1:
                    tally.n_het += 1
                elif is_hom_alt:
                    tally.n_hom_alt += 1
                    if ev.get("DP") is None:
                        # Depth unknown: this homozygote is neither adequate nor
                        # inadequate, and counting it either way is a claim the
                        # data does not support.
                        deep_hom["depth_unknown"] += 1
                    elif ev["DP"] >= cfg.low_depth:
                        deep_hom["n"] += 1
                        if "FALSE_HOM_SUSPECT" in flags:
                            deep_hom["flagged"] += 1

                # variant clustering within one sample
                window = recent.setdefault(label, [])
                window.append((rec.chrom, rec.pos))
                while window and (
                    window[0][0] != rec.chrom or rec.pos - window[0][1] > cfg.cluster_window
                ):
                    window.pop(0)
                if len(window) > cfg.cluster_min:
                    cluster_hits[label] = cluster_hits.get(label, 0) + 1

                all_flags = flags + site
                for f in all_flags:
                    tally.bump(f)

                if "MLEAC_DISCORDANT" in flags:
                    overlap["mleac_flagged"] += 1
                    if set(INDEPENDENT_HOM_RULES) & set(flags):
                        overlap["mleac_also_caught_by_rules"] += 1
                    if set(CONVENTIONAL_QC) & set(flags):
                        overlap["mleac_also_caught_by_qc"] += 1

                # Prefer examples that depth alone does not explain: a list of the
                # first 40 shallow calls in chromosome order tells you nothing,
                # while a hom call contradicted at adequate depth is a real lead.
                deep_enough = (ev.get("DP") or 0) >= cfg.low_depth
                room = (
                    len(examples) < cfg.max_examples
                    if deep_enough
                    else len(examples) < cfg.max_examples // 2
                )
                if "FALSE_HOM_SUSPECT" in flags and room:
                    examples.append(
                        {
                            "sample": label,
                            "locus": f"{rec.chrom}:{rec.pos}",
                            "ref": rec.ref,
                            "alt": ",".join(rec.alts),
                            "gt": gt,
                            "flags": sorted(set(all_flags)),
                            **{k: v for k, v in ev.items() if v is not None},
                        }
                    )
        if not saw_format:
            files_without_format.append(str(path))
        elif not saw_metrics:
            partial_files.append(str(path))
        else:
            n_with_format += 1

    if not tallies:
        from ..model import unknown

        return unknown(
            CHECK,
            "none of the supplied files carry per-genotype FORMAT fields, so genotype "
            "confidence cannot be assessed at all",
            n_records=n_records,
            files_without_format=files_without_format,
        )

    findings: list[Finding] = []
    notes: list[str] = []
    total_false_hom = sum(t.flags.get("FALSE_HOM_SUSPECT", 0) for t in tallies.values())
    total_hom = sum(t.n_hom_alt for t in tallies.values())

    if total_false_hom:
        by_sample = {
            s: t.flags.get("FALSE_HOM_SUSPECT", 0)
            for s, t in tallies.items()
            if t.flags.get("FALSE_HOM_SUSPECT")
        }
        findings.append(
            Finding(
                code="FALSE_HOM_SUSPECT",
                severity=Severity.ERROR,
                message=(
                    f"{total_false_hom} homozygous-alternate calls out of {total_hom} "
                    f"({total_false_hom / total_hom:.1%}) are not supported by their own "
                    f"evidence and must be excluded from recessive models until confirmed "
                    f"orthogonally"
                ),
                subjects=sorted(by_sample),
                evidence={
                    "n_false_hom": total_false_hom,
                    "n_hom_alt": total_hom,
                    "by_sample": by_sample,
                    "breakdown": {
                        rule: sum(t.flags.get(rule, 0) for t in tallies.values())
                        for rule in INDEPENDENT_HOM_RULES
                    },
                    "examples": examples,
                    "n_hom_alt_at_adequate_depth": deep_hom["n"],
                    "n_flagged_at_adequate_depth": deep_hom["flagged"],
                    "do_not_conclude": (
                        "that any homozygous call on this list supports a recessive model"
                    ),
                    "next_steps": [
                        "confirm every candidate homozygote by an orthogonal method "
                        "before building a recessive model on it"
                    ],
                },
            )
        )

    unassessable = sum(t.flags.get("HOM_NOT_ASSESSABLE", 0) for t in tallies.values())
    total_hom = sum(t.n_hom_alt for t in tallies.values())
    if unassessable and total_hom and unassessable / total_hom > 0.01:
        findings.append(
            Finding(
                code="HOM_NOT_ASSESSABLE",
                severity=Severity.ERROR,
                message=(
                    f"{unassessable} of {total_hom} homozygous-alternate calls "
                    f"({unassessable / total_hom:.0%}) carry neither a depth, a "
                    f"likelihood nor an allele count, so this check could not "
                    f"examine them at all - their absence from the flagged list "
                    f"is not evidence that they are sound"
                ),
                evidence={
                    "n_hom_not_assessable": unassessable,
                    "n_hom_alt": total_hom,
                    "needed": "FORMAT DP, or PL/GL, or AD / RO+AO / DP4",
                },
            )
        )

    total_filter = sum(filter_seen.values())
    if total_filter and filter_seen["absent"] / total_filter > 0.99:
        findings.append(
            Finding(
                code="NO_VARIANT_FILTERING",
                severity=Severity.WARN,
                message=(
                    f"the FILTER column is empty on {filter_seen['absent']} of "
                    f"{total_filter} records, so no variant-quality filtering (VQSR or "
                    f"hard filters) was ever applied to this call set - every quality "
                    f"judgement below is therefore being made for the first time here"
                ),
                evidence={
                    "filter_distribution": dict(filter_seen),
                    "next_steps": [
                        "apply VQSR or GATK hard filters before interpreting this call set"
                    ],
                },
            )
        )

    for label, t in sorted(tallies.items()):
        low = t.flags.get("LOW_DEPTH", 0)
        if t.n_genotypes and low / t.n_genotypes > 0.25:
            findings.append(
                Finding(
                    code="WIDESPREAD_LOW_DEPTH",
                    severity=Severity.WARN,
                    message=(
                        f"{label}: {low / t.n_genotypes:.0%} of genotypes are below DP="
                        f"{cfg.low_depth}, so zygosity is unreliable across much of this sample"
                    ),
                    subjects=[label],
                    evidence={"n_low_depth": low, "n_genotypes": t.n_genotypes},
                )
            )

    if overlap["mleac_flagged"]:
        caught = overlap["mleac_also_caught_by_rules"]
        unique = overlap["mleac_flagged"] - caught
        notes.append(
            f"MLEAC disagreed with the emitted genotype {overlap['mleac_flagged']} time(s); "
            f"{caught} of those were already flagged by depth, likelihood or allele balance, "
            f"leaving {unique} that only MLEAC caught"
        )
    if files_without_format:
        notes.append(
            f"{len(files_without_format)} input(s) carry no FORMAT fields and were skipped by "
            f"this check"
        )
    for label, n in sorted(cluster_hits.items()):
        if n:
            notes.append(
                f"{label}: {n} variant(s) sit in a cluster of >{cfg.cluster_min} within "
                f"{cfg.cluster_window} bp, a common signature of local misalignment"
            )

    if partial_files:
        notes.append(
            f"{len(partial_files)} input(s) carry no FORMAT block, so only the depth-based "
            f"and site-level rules could run on them; allele balance, GQ and genotype "
            f"likelihoods are unrecoverable there"
        )

    if total_false_hom:
        status = Status.WARN
        rate = total_false_hom / total_hom if total_hom else 0.0
        summary = f"{total_false_hom} of {total_hom} homozygous calls unsupported ({rate:.0%})"
        if deep_hom["n"]:
            summary += f"; {deep_hom['flagged']} at DP>={cfg.low_depth}"
    elif findings:
        status = Status.WARN
        summary = "genotype quality problems found"
    else:
        status = Status.PASS
        summary = f"{sum(t.n_genotypes for t in tallies.values())} genotypes, none contradicted"

    return CheckResult(
        check=CHECK,
        status=status,
        summary=summary,
        findings=findings,
        metrics={
            "n_records": n_records,
            "n_files_with_format": n_with_format,
            "rule_overlap": overlap,
            "filter_distribution": dict(filter_seen),
            "hom_at_adequate_depth": dict(deep_hom),
            "files_partially_assessable": partial_files,
            "per_sample": {
                s: {
                    "n_genotypes": t.n_genotypes,
                    "n_het": t.n_het,
                    "n_hom_alt": t.n_hom_alt,
                    "flags": dict(sorted(t.flags.items())),
                }
                for s, t in sorted(tallies.items())
            },
        },
        notes=notes,
    )
