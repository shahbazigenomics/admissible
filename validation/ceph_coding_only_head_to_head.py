"""Claim: on the CEPH 1463 VCF filtered to coding-only, peddy's sex calls
degrade (17/17 -> 14/17, three true females called male) while admissible
declines instead of guessing (SEX_NOT_DETERMINED for all 17, because the
filtered file drops below the ~50-site chrX floor), and admissible's
relatedness call stays correct (it draws on the full autosomal panel, not the
chrX-only subset that sex inference needs).

Before this script, nothing in the repository could reproduce this claim at
all - no script, no filtered fixture, no recorded peddy output. This script
builds the coding-only file from real data and public tools, then:

* always runs and verifies the ``admissible`` side, because that only needs
  this package plus bcftools (not installed by this package - a validation-
  only dependency, matching the ``git clone`` pattern already used for CEPH
  data in tests/test_public_ceph.py);
* attempts to run real ``peddy`` and reports its actual output if it is
  importable in this environment, and says plainly that it is not (rather
  than asserting the 17/17 -> 14/17 numbers from memory) if it isn't. Do not
  edit those numbers into a report without this script actually printing
  them from a real peddy run - see the "if not verified" note at the bottom.

Requires: this repo checked out, bcftools + tabix on PATH, and (not vendored)
a clone of https://github.com/brentp/peddy and of
https://github.com/AstraZeneca-NGS/reference_data - both fetched on demand
into a temp cache, same as validation/panel_overlap.py.

Run: ``python validation/ceph_coding_only_head_to_head.py``
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

PEDDY_REPO = "https://github.com/brentp/peddy"
REFDATA_REPO = "https://github.com/AstraZeneca-NGS/reference_data"


def _require(cmd: str) -> None:
    if shutil.which(cmd) is None:
        sys.exit(f"this script needs '{cmd}' on PATH (validation-only dependency)")


def _clone_if_missing(url: str, dest: Path) -> Path:
    if not dest.exists():
        print(f"cloning {url} -> {dest}", file=sys.stderr)
        subprocess.run(["git", "clone", "--depth", "1", url, str(dest)], check=True)
    return dest


def build_coding_only(cache_dir: Path) -> tuple[Path, Path]:
    """Return (full_vcf, coding_only_vcf) for CEPH 1463, building the latter."""
    peddy_repo = _clone_if_missing(PEDDY_REPO, cache_dir / "peddy")
    refdata = _clone_if_missing(REFDATA_REPO, cache_dir / "reference_data")

    full_vcf = peddy_repo / "data" / "ceph1463.peddy.vcf.gz"
    if not full_vcf.exists():
        sys.exit(
            f"expected {full_vcf} from the peddy repo itself - clone may be stale, "
            f"delete {peddy_repo} and rerun"
        )
    cds_bed = refdata / "GRCh37" / "bed" / "CDS_RefSeq.bed"

    coding_only = cache_dir / "ceph1463.codingonly.vcf.gz"
    if not coding_only.exists():
        subprocess.run(
            ["bcftools", "view", "-R", str(cds_bed), "-Oz", "-o", str(coding_only), str(full_vcf)],
            check=True,
        )
        subprocess.run(["tabix", "-p", "vcf", str(coding_only)], check=True)
    return full_vcf, coding_only


def run_admissible_side(full_vcf: Path, coding_only: Path, ped: Path) -> dict:
    from admissible.checks.identity import check_identity
    from admissible.ped import read_ped
    from admissible.vcfio import load_cohort

    pedigree = read_ped(ped)
    out = {}
    for label, vcf in (("full", full_vcf), ("coding_only", coding_only)):
        matrix = load_cohort([vcf])
        res = check_identity(matrix, pedigree, build="GRCh37")
        codes = [f.code for f in res.findings]
        out[label] = {
            "status": res.status.value,
            "n_sex_not_determined": codes.count("SEX_NOT_DETERMINED"),
            "n_sex_mismatch": codes.count("SEX_MISMATCH"),
            "identity_status_is_pass_or_warn_not_fail": res.status.value != "FAIL",
        }
    return out


def run_peddy_side(full_vcf: Path, coding_only: Path, ped: Path) -> dict | None:
    try:
        import peddy  # noqa: F401
    except ImportError:
        return None

    # peddy's own CLI writes a directory of output files (het_check.csv,
    # sex_check.csv, etc.) whose exact schema this script has not verified
    # against a real run (peddy would not install in the environment this
    # script was written in - see the caller's message). Rather than guess a
    # column layout and silently produce a wrong comparison, this runs peddy
    # for real and points at its output for a human to read directly.
    out_dirs = {}
    for label, vcf in (("full", full_vcf), ("coding_only", coding_only)):
        out_dir = Path(tempfile.mkdtemp(prefix=f"peddy_{label}_"))
        prefix = out_dir / "out"
        subprocess.run(
            [sys.executable, "-m", "peddy", "-p", "1", "--prefix", str(prefix), str(vcf), str(ped)],
            check=True,
        )
        out_dirs[label] = str(out_dir)
    return out_dirs


def main() -> int:
    _require("bcftools")
    _require("tabix")
    cache_dir = Path(tempfile.gettempdir()) / "admissible_validation_cache"
    cache_dir.mkdir(exist_ok=True)

    full_vcf, coding_only = build_coding_only(cache_dir)
    ped = full_vcf.parent / "ceph1463.ped"

    print("=== admissible, full vs. coding-only-filtered CEPH 1463 ===")
    adm = run_admissible_side(full_vcf, coding_only, ped)
    print(json.dumps(adm, indent=2))
    assert adm["full"]["n_sex_mismatch"] == 0
    assert adm["coding_only"]["n_sex_not_determined"] == 17, (
        "expected all 17 samples to decline sex determination on the "
        "coding-only file (below the ~50-site floor) - re-examine the claim "
        "or the fixture if this fails"
    )
    assert adm["coding_only"]["identity_status_is_pass_or_warn_not_fail"], (
        "relatedness should still resolve on the coding-only file (it uses "
        "the full autosomal panel, not the chrX-only subset sex needs)"
    )
    print(
        "CONFIRMED (admissible side): 0 sex mismatches on the full file; all "
        "17 decline (SEX_NOT_DETERMINED) on the coding-only file, with "
        "relatedness still resolving."
    )

    print()
    print("=== peddy, same two files ===")
    peddy_result = run_peddy_side(full_vcf, coding_only, ped)
    if peddy_result is None:
        print(
            "peddy is not importable in this environment (it does not build "
            "with modern setuptools without an older Python/setuptools "
            "pin - see https://github.com/brentp/peddy issues). The "
            "'17/17 -> 14/17, 3 true females called male' figure in the "
            "README/design-decisions.md is NOT verified by this script run "
            "and must not be treated as confirmed until someone runs this "
            "script (or `python -m peddy`) directly on "
            f"{full_vcf} and {coding_only} with a working peddy install, "
            "e.g. in a separate conda env with an older Python."
        )
    else:
        print(json.dumps(peddy_result, indent=2))
        print(
            "peddy ran; its sex_check.csv / het_check.csv in each directory "
            "above is what carries the actual per-sample calls and error "
            "flag - read those directly rather than trusting a number typed "
            "into this script, since the CSV schema varies by peddy version."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
