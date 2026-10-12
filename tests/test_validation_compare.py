"""validation/callability/compare.py: the comparison logic, on constructed inputs.

The bedtools baseline and the driver script are NOT exercised here (no bedtools in
this environment); these tests only pin what compare.py does with the two things it is
handed, so that a disagreement or a malformed report is reported rather than missed.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "validation_compare",
    Path(__file__).resolve().parents[1] / "validation" / "callability" / "compare.py",
)
cmp = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(cmp)


def _baseline(tmp_path, target=1000, callable_=800, minus=700, regions=None):
    (tmp_path / "baseline.tsv").write_text(
        f"target_bp\t{target}\ncallable_in_target_bp\t{callable_}\n"
        f"callable_in_target_minus_exclude_bp\t{minus}\n"
    )
    regions = regions if regions is not None else {"E1": (500, 400, 350), "E2": (500, 400, 350)}
    (tmp_path / "regions.tsv").write_text(
        "".join(f"{k}\t{a}\t{b}\t{c}\n" for k, (a, b, c) in regions.items())
    )
    return tmp_path


def _report(target=1000, callable_=800, raw=None, rows=None, truncated=False, units=1):
    model = {"callable_bp": callable_}
    if raw is not None:
        model["without_exclusion"] = {"callable_bp": raw}
    rows = rows if rows is not None else [
        {"region": "E1", "bp": 500, "callable_bp": 400},
        {"region": "E2", "bp": 500, "callable_bp": 400},
    ]
    unit = {"regions": rows, "truncated": truncated, "n_regions": len(rows)}
    return {"checks": [{"check": "callability", "status": "PASS", "summary": "",
                        "metrics": {"target_bp": target, "joint_callable_bp": callable_,
                                    "models": {"recessive": model},
                                    "by_region": {"units": {f"u{i}": unit for i in range(units)}},
                                    }}]}


def _run(tmp_path, report, *flags):
    j = tmp_path / "r.json"
    j.write_text(json.dumps(report))
    return cmp.main([str(tmp_path), str(j), *flags])


def test_identical_numbers_exit_zero(tmp_path):
    _baseline(tmp_path)
    assert _run(tmp_path, _report()) == 0


def test_one_base_of_difference_is_a_failure_not_rounded_away(tmp_path):
    _baseline(tmp_path)
    assert _run(tmp_path, _report(callable_=801)) == 1


def test_a_region_that_differs_is_found(tmp_path):
    _baseline(tmp_path)
    rows = [{"region": "E1", "bp": 500, "callable_bp": 399},
            {"region": "E2", "bp": 500, "callable_bp": 400}]
    assert _run(tmp_path, _report(rows=rows)) == 1


def test_a_region_missing_from_admissible_is_a_failure_unless_truncated(tmp_path, capsys):
    _baseline(tmp_path)
    rows = [{"region": "E1", "bp": 500, "callable_bp": 400}]
    assert _run(tmp_path, _report(rows=rows)) == 1
    assert _run(tmp_path, _report(rows=rows, truncated=True)) == 0
    assert "not compared" in capsys.readouterr().out


def test_with_exclusion_compares_the_filtered_and_the_raw_figure(tmp_path):
    _baseline(tmp_path)
    rows = [{"region": "E1", "bp": 500, "callable_bp": 350},
            {"region": "E2", "bp": 500, "callable_bp": 350}]
    assert _run(tmp_path, _report(callable_=700, raw=800, rows=rows), "--with-exclusion") == 0
    assert _run(tmp_path, _report(callable_=700, raw=799, rows=rows), "--with-exclusion") == 1


def test_a_report_that_computed_nothing_is_reported_not_raised(tmp_path):
    _baseline(tmp_path)
    unknown = {"checks": [{"check": "callability", "status": "UNKNOWN",
                           "summary": "no target", "metrics": {}}]}
    assert _run(tmp_path, unknown) == 2
    assert _run(tmp_path, {"checks": []}) == 2


@pytest.mark.parametrize("junk", ["", "not json", "[]", '{"checks": 3}', '{"checks": [3]}'])
def test_garbage_json_never_raises(tmp_path, junk):
    _baseline(tmp_path)
    (tmp_path / "r.json").write_text(junk)
    assert cmp.main([str(tmp_path), str(tmp_path / "r.json")]) == 2


def test_missing_or_garbled_baseline_never_raises(tmp_path):
    (tmp_path / "r.json").write_text(json.dumps(_report()))
    assert cmp.main([str(tmp_path), str(tmp_path / "r.json")]) == 2  # no baseline.tsv
    (tmp_path / "baseline.tsv").write_text("target_bp\tNaN-ish\n")
    assert cmp.main([str(tmp_path), str(tmp_path / "r.json")]) == 2


def test_more_than_one_unit_is_refused_for_a_single_sample(tmp_path):
    _baseline(tmp_path)
    assert _run(tmp_path, _report(units=2)) in (0, 1)  # base numbers still compared
    out = tmp_path / "res.tsv"
    j = tmp_path / "r.json"
    j.write_text(json.dumps(_report(units=2)))
    cmp.main([str(tmp_path), str(j), "--out", str(out)])
    assert "region" not in out.read_text()  # no per-region rows were compared
