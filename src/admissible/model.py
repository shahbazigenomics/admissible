"""Core result types.

Design contract for this whole package: **no check ever raises.**  A check that
cannot answer its question returns a well-formed result with ``Status.UNKNOWN``
and a stated reason.  "I could not determine this" is a valid, machine-readable
answer; a traceback is not.
"""

from __future__ import annotations

import functools
import traceback
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

SCHEMA_VERSION = "admissible/report/1"


class Status(str, Enum):
    """Outcome of a single check."""

    PASS = "PASS"
    WARN = "WARN"
    FAIL = "FAIL"
    UNKNOWN = "UNKNOWN"  # inputs cannot answer the question
    SKIPPED = "SKIPPED"  # not applicable, or not requested

    @property
    def rank(self) -> int:
        return {"PASS": 0, "SKIPPED": 1, "UNKNOWN": 2, "WARN": 3, "FAIL": 4}[self.value]


class Severity(str, Enum):
    INFO = "INFO"
    WARN = "WARN"
    ERROR = "ERROR"
    BLOCKING = "BLOCKING"  # nothing downstream is interpretable until this is fixed


class Verdict(str, Enum):
    INTERPRETABLE = "INTERPRETABLE"
    INTERPRETABLE_WITH_CAVEATS = "INTERPRETABLE WITH CAVEATS"
    NOT_INTERPRETABLE = "NOT INTERPRETABLE"
    UNDETERMINED = "UNDETERMINED"


@dataclass
class Finding:
    """One specific, evidenced statement about the data.

    ``evidence`` holds the numbers behind the claim so a reader can disagree
    with the threshold without re-running anything.
    """

    code: str
    severity: Severity
    message: str
    subjects: list[str] = field(default_factory=list)
    evidence: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "code": self.code,
            "severity": self.severity.value,
            "message": self.message,
            "subjects": list(self.subjects),
            "evidence": _jsonable(self.evidence),
        }


@dataclass
class CheckResult:
    check: str
    status: Status
    summary: str = ""
    findings: list[Finding] = field(default_factory=list)
    metrics: dict[str, Any] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    @property
    def blocking(self) -> bool:
        return any(f.severity is Severity.BLOCKING for f in self.findings)

    def to_dict(self) -> dict[str, Any]:
        return {
            "check": self.check,
            "status": self.status.value,
            "summary": self.summary,
            "blocking": self.blocking,
            "findings": [f.to_dict() for f in self.findings],
            "metrics": _jsonable(self.metrics),
            "notes": list(self.notes),
        }


@dataclass
class Report:
    family_id: str
    checks: list[CheckResult] = field(default_factory=list)
    inputs: dict[str, Any] = field(default_factory=dict)
    tool_version: str = "0.1.0"

    @property
    def verdict(self) -> Verdict:
        if any(c.blocking for c in self.checks):
            return Verdict.NOT_INTERPRETABLE
        statuses = [c.status for c in self.checks]
        if Status.FAIL in statuses:
            return Verdict.NOT_INTERPRETABLE
        if Status.WARN in statuses:
            return Verdict.INTERPRETABLE_WITH_CAVEATS
        if all(s in (Status.UNKNOWN, Status.SKIPPED) for s in statuses):
            return Verdict.UNDETERMINED
        return Verdict.INTERPRETABLE

    @property
    def do_not_conclude(self) -> list[str]:
        out: list[str] = []
        for c in self.checks:
            for f in c.findings:
                claim = f.evidence.get("do_not_conclude")
                if claim and claim not in out:
                    out.append(claim)
        return out

    @property
    def next_steps(self) -> list[tuple[str, bool]]:
        """(step, blocking): blocking steps first, then in the order raised."""
        seen: dict[str, bool] = {}
        order: dict[str, int] = {}
        for c in self.checks:
            for f in c.findings:
                for step in f.evidence.get("next_steps", []) or []:
                    blocking = f.severity is Severity.BLOCKING
                    seen[step] = seen.get(step, False) or blocking
                    order.setdefault(step, len(order))
        return sorted(seen.items(), key=lambda kv: (not kv[1], order[kv[0]]))

    def get(self, check: str) -> CheckResult | None:
        for c in self.checks:
            if c.check == check:
                return c
        return None

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA_VERSION,
            "tool_version": self.tool_version,
            "family_id": self.family_id,
            "verdict": self.verdict.value,
            "do_not_conclude": self.do_not_conclude,
            "next_steps": [{"step": s, "blocking": b} for s, b in self.next_steps],
            "inputs": _jsonable(self.inputs),
            "checks": [c.to_dict() for c in self.checks],
        }


def unknown(check: str, reason: str, **metrics: Any) -> CheckResult:
    """The canonical 'I cannot answer this from what you gave me' result."""
    return CheckResult(
        check=check,
        status=Status.UNKNOWN,
        summary=reason,
        findings=[Finding(code="CHECK_NOT_COMPUTABLE", severity=Severity.INFO, message=reason)],
        metrics=dict(metrics),
    )


def never_raises(check_name: str):
    """Decorator turning any escaped exception into an UNKNOWN CheckResult.

    This is a backstop, not a licence to skip error handling: a check that
    lands here has a bug, and the traceback is preserved in the evidence so it
    shows up in the JSON report rather than on someone's terminal.
    """

    def deco(fn):
        @functools.wraps(fn)
        def wrapper(*args, **kwargs) -> CheckResult:
            try:
                return fn(*args, **kwargs)
            except Exception as exc:  # noqa: BLE001 - deliberate catch-all
                return CheckResult(
                    check=check_name,
                    status=Status.UNKNOWN,
                    summary=f"check aborted: {type(exc).__name__}: {exc}",
                    findings=[
                        Finding(
                            code="CHECK_CRASHED",
                            severity=Severity.ERROR,
                            message=f"{type(exc).__name__}: {exc}",
                            evidence={"traceback": traceback.format_exc(limit=12)},
                        )
                    ],
                )

        return wrapper

    return deco


def _jsonable(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, Enum):
        return obj.value
    if isinstance(obj, float):
        return round(obj, 6)
    if isinstance(obj, (str, int, bool)) or obj is None:
        return obj
    return str(obj)
