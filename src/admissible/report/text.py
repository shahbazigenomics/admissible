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
            if res is None or (not res.findings and not res.notes):
                continue
            out.append("")
            out.append(f"{label.upper()}  [{res.status.value}]")
            for f in res.findings:
                if f.severity is Severity.INFO and res.status is Status.UNKNOWN:
                    continue
                out.append(_wrap(f"  [{f.severity.value}] {f.code}: {f.message}", "      "))
            for n in res.notes:
                out.append(_wrap(f"  note: {n}", "        "))

    return "\n".join(out) + "\n"


def _wrap(text: str, indent: str) -> str:
    return "\n".join(
        textwrap.wrap(text, width=WIDTH, subsequent_indent=indent, break_long_words=False)
    )
