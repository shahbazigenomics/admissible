"""The one-page verdict.

Two rules govern this renderer.  First, a line never states a conclusion the
check did not reach: an UNKNOWN check prints UNKNOWN, it does not print a
reassuring blank.  Second, every threshold that is a convention rather than a
derived quantity says so where it is printed.
"""

from __future__ import annotations

import textwrap

from ..model import Report, Severity, Status

WIDTH = 78
ROWS = [
    ("identity", "Identity"),
    ("provenance", "Provenance"),
    ("callability", "Callable"),
    ("genotype", "Genotype QC"),
    ("models", "Models"),
]


def render_text(report: Report, verbose: bool = False) -> str:
    out: list[str] = []
    label = "COHORT" if "families" in report.family_id else "FAMILY"
    head = f"{label} {report.family_id}"
    tail = f"VERDICT: {report.verdict.value}"
    out.append(f"{head}{' ' * max(2, WIDTH - len(head) - len(tail))}{tail}")
    out.append("")

    for key, label in ROWS:
        res = report.get(key)
        dots = label + " " + "." * max(2, 18 - len(label))
        if res is None:
            out.append(f"{dots} {'SKIPPED':<8} not run")
            continue
        room = WIDTH - len(dots) - 10
        summary = res.summary if len(res.summary) <= room else res.summary[: room - 1] + "…"
        out.append(f"{dots} {res.status.value:<8} {summary}")

    if report.do_not_conclude:
        out.append("")
        for claim in report.do_not_conclude:
            out.append(f'DO NOT CONCLUDE: "{claim}"')

    steps = report.next_steps
    if steps:
        out.append("")
        for i, (step, blocking) in enumerate(steps, 1):
            prefix = "NEXT STEPS: " if i == 1 else " " * 12
            mark = " (blocking)" if blocking and "blocking" not in step else ""
            out.append(
                "\n".join(
                    textwrap.wrap(
                        f"{prefix}{i}. {step}{mark}",
                        width=WIDTH,
                        subsequent_indent=" " * (len(prefix) + 3),
                        break_long_words=False,
                    )
                )
            )

    if verbose:
        out.append("")
        out.append("-" * WIDTH)
        for key, label in ROWS:
            res = report.get(key)
            if res is None:
                continue
            pairwise = _pairwise_evidence_lines(res) if key == "identity" else []
            if not res.findings and not res.notes and not pairwise:
                continue
            out.append("")
            out.append(f"{label.upper()}  [{res.status.value}]")
            for f in res.findings:
                if f.severity is Severity.INFO and res.status is Status.UNKNOWN:
                    continue
                out.append(_wrap(f"  [{f.severity.value}] {f.code}: {f.message}", "      "))
            for n in res.notes:
                out.append(_wrap(f"  note: {n}", "        "))
            out.extend(pairwise)

    return "\n".join(out) + "\n"


def _wrap(text: str, indent: str) -> str:
    return "\n".join(
        textwrap.wrap(text, width=WIDTH, subsequent_indent=indent, break_long_words=False)
    )


# Printed even for a pair that didn't cross any threshold hard enough to
# become a Finding - the numbers a borderline call was (or wasn't) made from
# used to be visible only in --json. "Close to a threshold" is deliberately
# asymmetric: a pair whose agreement alone is near/above the duplicate bar is
# notable regardless of how far its jaccard is from its own threshold, since
# that mismatch (high agreement, low jaccard) is exactly the failure mode
# HIGH_AGREEMENT_LOW_JACCARD exists to catch - "close to the jaccard
# threshold" would miss it entirely (see identity.py's dup_completeness_ratio
# comment for why jaccard alone is not trustworthy here).
_AGREEMENT_MARGIN = 0.05
_BOUNDARY_MARGIN = 0.05
_MAX_PAIRWISE_ROWS = 25


def _pairwise_evidence_lines(res) -> list[str]:
    pairs = (res.metrics or {}).get("pairs")
    if not pairs:
        return []
    thresholds = (res.metrics or {}).get("dup_thresholds") or {}
    dup_agreement = thresholds.get("agreement")
    boundary = ((res.metrics or {}).get("relatedness_calibration") or {}).get("boundary")

    notable = []
    for p in pairs:
        agreement, jaccard = p.get("agreement"), p.get("jaccard")
        near_dup = (
            dup_agreement is not None
            and agreement is not None
            and agreement >= dup_agreement - _AGREEMENT_MARGIN
        )
        near_boundary = (
            boundary is not None
            and jaccard is not None
            and abs(jaccard - boundary) <= _BOUNDARY_MARGIN
        )
        if near_dup or near_boundary:
            notable.append(p)
    if not notable:
        return []

    lines = ["", "  PAIRWISE EVIDENCE (verbose)"]
    for p in notable[:_MAX_PAIRWISE_ROWS]:
        declared = p.get("declared_relationship")
        declared_part = f"  (declared: {declared})" if declared else ""
        lines.append(
            f"    {p['a']} / {p['b']}    agreement={p['agreement']:.4f}  "
            f"jaccard={p['jaccard']:.4f}  n_shared={p['n_shared_nonref']}"
            f"{declared_part}"
        )
    if len(notable) > _MAX_PAIRWISE_ROWS:
        lines.append(
            f"    ... {len(notable) - _MAX_PAIRWISE_ROWS} more pair(s) meeting the "
            f"same criteria; see --json for the full pairwise table"
        )
    return lines
