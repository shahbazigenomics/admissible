from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tests" / "fixtures" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="session")
def cohort(tmp_path_factory) -> dict:
    """The 13-sample synthetic cohort, with the metrics the generator achieved."""
    outdir = tmp_path_factory.mktemp("cohort")
    metrics = _load("make_fixtures").build(outdir)
    return {"dir": outdir, "metrics": metrics, "ped": outdir / "cohort.ped",
            "vcfs": sorted(outdir.glob("FAM*.vcf"))}


@pytest.fixture(scope="session")
def joint(tmp_path_factory) -> dict:
    outdir = tmp_path_factory.mktemp("joint")
    vcf = _load("make_joint_family").build(outdir)
    return {"vcf": vcf, "ped": outdir / "joint_family.ped"}
