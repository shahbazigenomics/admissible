from __future__ import annotations

import json

from admissible.cli import main
from admissible.model import SCHEMA_VERSION


def test_audit_exits_one_and_names_the_problem(cohort, capsys, tmp_path):
    out_json = tmp_path / "report.json"
    code = main(
        [
            "audit",
            *[str(p) for p in cohort["vcfs"]],
            "--ped", str(cohort["ped"]),
            "--json", str(out_json),
        ]
    )
    assert code == 1
    text = capsys.readouterr().out
    assert "VERDICT: NOT INTERPRETABLE" in text
    assert "Identity" in text and "FAIL" in text
    assert "DO NOT CONCLUDE" in text
    assert "NEXT STEPS" in text
    # The blocking step must come first, whatever order findings were raised in.
    steps = [ln for ln in text.splitlines() if ". resolve sample identity" in ln]
    assert steps and steps[0].startswith("NEXT STEPS: 1.")

    doc = json.loads(out_json.read_text())
    assert doc["schema"] == SCHEMA_VERSION
    assert doc["verdict"] == "NOT INTERPRETABLE"
    assert {c["check"] for c in doc["checks"]} == {
        "identity", "provenance", "callability", "genotype", "models"
    }


def test_unknown_means_missing_data_not_missing_code(cohort, capsys):
    """All five checks are implemented, so UNKNOWN must reflect the inputs.

    Callability is UNKNOWN here because no coverage was supplied - a VCF cannot
    answer that question - while every check the inputs CAN support reports a
    real status.
    """
    main(["audit", *[str(p) for p in cohort["vcfs"]], "--ped", str(cohort["ped"])])
    text = capsys.readouterr().out
    callable_line = next(ln for ln in text.splitlines() if ln.startswith("Callable"))
    assert "UNKNOWN" in callable_line
    for label in ("Identity", "Provenance", "Genotype QC", "Models"):
        line = next(ln for ln in text.splitlines() if ln.startswith(label))
        assert "UNKNOWN" not in line, label


def test_provenance_flags_the_coding_only_pass_only_subset(cohort, capsys):
    code = main(["provenance", *[str(p) for p in cohort["vcfs"]]])
    text = capsys.readouterr().out
    assert "CODING_ONLY" in text
    assert "PASS_FILTERED" in text
    assert "SUBSET" in text
    assert code in (0, 1)


def test_every_line_fits_the_page(cohort, capsys):
    main(["audit", *[str(p) for p in cohort["vcfs"]], "--ped", str(cohort["ped"])])
    for line in capsys.readouterr().out.splitlines():
        assert len(line) <= 100


def test_missing_input_exits_two(tmp_path, capsys):
    assert main(["audit", str(tmp_path / "nope.vcf"), "--ped", str(tmp_path / "x.ped")]) == 2


def test_clean_family_passes_and_exits_zero(joint, capsys):
    code = main(["identity", str(joint["vcf"]), "--ped", str(joint["ped"])])
    assert code == 0
    assert "VERDICT: INTERPRETABLE" in capsys.readouterr().out
