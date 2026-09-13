from __future__ import annotations

import json

from ..model import Report


def render_json(report: Report, indent: int | None = 2) -> str:
    return json.dumps(report.to_dict(), indent=indent, sort_keys=False)
