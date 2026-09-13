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


def test_callability_runs_without_any_vcf(tmp_path, capsys):
    """Check 4 asks a question a VCF cannot answer, so it must not demand one.

    The inputs are a target BED and one coverage BED per sample - nothing else.
    An earlier version required a positional VCF here, which made the only
    check that never reads a VCF impossible to run on its own.
    """
    ped = tmp_path / "f.ped"
    ped.write_text(
        "#FID\tIID\tPAT\tMAT\tSEX\tPHENO\n"
        "F\tDAD\t0\t0\t1\t1\nF\tMUM\t0\t0\t2\t1\n"
        "F\tKID1\tDAD\tMUM\t1\t2\nF\tKID2\tDAD\tMUM\t2\t2\n"
    )
    target = tmp_path / "target.bed"
    target.write_text("chr1\t0\t1000\n")
    args = ["callability", "--ped", str(ped), "--target", str(target)]
    for s, end in (("KID1", 900), ("KID2", 800)):
        # written in genuine mosdepth --quantize 0:10: shape: both bins present
        p = tmp_path / f"{s}.quantized.bed"
        p.write_text(f"chr1\t0\t{end}\t10:inf\nchr1\t{end}\t1000\t0:10\n")
        args += ["--coverage", f"{s}={p}"]

    assert main(args) == 0
    out = capsys.readouterr().out
    assert "0.80 of target searched" in out


def test_min_depth_is_settable_from_the_command_line(tmp_path, capsys):
    """The depth floor is a clinical choice, so it must not be compiled in."""
    ped = tmp_path / "f.ped"
    ped.write_text(
        "#FID\tIID\tPAT\tMAT\tSEX\tPHENO\n"
        "F\tDAD\t0\t0\t1\t1\nF\tMUM\t0\t0\t2\t1\nF\tKID1\tDAD\tMUM\t1\t2\n"
    )
    target = tmp_path / "target.bed"
    target.write_text("chr1\t0\t1000\n")
    q = tmp_path / "KID1.quantized.bed"
    q.write_text("chr1\t0\t200\t0:10\nchr1\t200\t600\t10:30\nchr1\t600\t1000\t30:inf\n")
    base = ["callability", "--ped", str(ped), "--target", str(target),
            "--coverage", f"KID1={q}"]

    main(base)
    assert "0.80 of target searched" in capsys.readouterr().out
    main(base + ["--min-depth", "30"])
    assert "0.40 of target searched" in capsys.readouterr().out


def test_identity_runs_without_a_pedigree(cohort, capsys):
    """Duplicate detection and sex inference need no pedigree at all.

    An unpedigreed cohort is exactly when "are any two of these the same
    person" matters most, and check 1 previously crashed on ``ped=None`` after
    doing all of the pairwise work - on a 2,504-sample public VCF, after six
    minutes.
    """
    code = main(["identity", *[str(p) for p in cohort["vcfs"]], "-v"])
    out = capsys.readouterr().out
    assert code in (0, 1)
    # Nothing is claimed about a pedigree that was never supplied.
    assert "SAMPLE_NOT_IN_PED" not in out
    assert "pedigree broadly consistent" not in out
    assert "no pedigree supplied" in out
    # The planted duplicates are still found without any pedigree.
    assert "DUPLICATE" in out


def test_the_header_names_only_the_families_actually_supplied(tmp_path, capsys):
    """A whole-cohort PED reused for one family's run must not make the report
    claim the other families were looked at.

    Found by running the tool the way a first-time user would: three samples of
    one family, plus the cohort PED that was already on disk. The header read
    "3 families: FAM_A,FAM_B,FAM_C" for a run that saw only FAM_A - an overclaim
    on the headline line of a report whose whole purpose is to stop overclaims.
    """
    head = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA1\n"
    )
    vcf = tmp_path / "a.vcf"
    vcf.write_text(head + "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP:GQ\t0/1:30:99\n")
    ped = tmp_path / "cohort.ped"
    ped.write_text(
        "FAM_A\tA1\t0\t0\t1\t2\n"
        "FAM_B\tB1\t0\t0\t1\t2\n"
        "FAM_C\tC1\t0\t0\t2\t1\n"
    )

    main(["identity", str(vcf), "--ped", str(ped)])
    out = capsys.readouterr().out
    assert "FAM_A" in out
    assert "FAM_B" not in out and "FAM_C" not in out
    assert "3 families" not in out
