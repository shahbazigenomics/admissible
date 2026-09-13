"""Check 4 against real capture coverage, not invented coverage.

Every other test of check 4 uses BEDs written by hand or by a simulator, and a
simulator cannot settle the one question this check exists to answer: how far
does the product of per-sample callable fractions understate the joint one?  That
gap is a property of real probe-efficiency correlation.  Invented coverage
answers only with the correlation that was put into it.

The data is four unrelated 1000 Genomes exomes over chr20:1,400,000-1,500,000,
with the denominator taken independently from an Ensembl GTF.  See
``tests/data/1000g_exome_chr20/PROVENANCE.md`` for exactly what it is and what it
is not - in particular, 4,134 bp is a small measurement.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from admissible.checks.callability import check_callability, load_coverage
from admissible.model import Status
from admissible.ped import read_ped

DATA = Path(__file__).parent / "data" / "1000g_exome_chr20"
SAMPLES = ("HG00349", "HG00350", "HG00351", "HG00358")
TARGET = DATA / "chr20_coding_target.bed"


@pytest.fixture(scope="module")
def result(tmp_path_factory):
    ped = tmp_path_factory.mktemp("ped") / "cohort.ped"
    # Unrelated individuals; "affected" here only selects who is intersected.
    ped.write_text(
        "#FID\tIID\tPAT\tMAT\tSEX\tPHENO\n"
        + "".join(f"FIN\t{s}\t0\t0\t0\t2\n" for s in SAMPLES)
    )
    spec = {s: DATA / f"{s}.chr20_1400000-1500000.quantized.bed.gz" for s in SAMPLES}
    return check_callability(load_coverage(spec), read_ped(ped), target_bed=TARGET)


def test_real_mosdepth_quantize_output_is_read(result):
    """Genuine mosdepth output, read without a shim."""
    assert result.metrics["target_bp"] == 4134
    assert result.metrics["n_affected_with_coverage"] == 4
    per = result.metrics["per_sample"]
    assert per["HG00351"]["fraction_of_target"] == pytest.approx(0.4768, abs=5e-4)
    assert per["HG00358"]["fraction_of_target"] == pytest.approx(0.7975, abs=5e-4)


def test_the_product_understates_the_joint_fraction_on_real_capture(result):
    """The whole argument for check 4, measured rather than asserted.

    Product 0.2232 against a true intersection of 0.3815: the shortcut is 41%
    low, because the four samples lose coverage in the same places - the same
    GC-extreme exons and low-efficiency probes - and the product assumes they do
    not.
    """
    m = result.metrics
    assert m["joint_callable_bp"] == 1577
    assert m["joint_callable_fraction"] == pytest.approx(0.3815, abs=5e-4)
    assert m["product_of_per_sample_fractions"] == pytest.approx(0.2232, abs=5e-4)
    assert m["product_of_per_sample_fractions"] < m["joint_callable_fraction"]


def test_a_model_needing_one_fewer_sample_searches_far_more_of_the_target(result):
    """Callability is a property of the model, not of the cohort."""
    models = result.metrics["models"]
    assert models["recessive"]["fraction_of_target"] == pytest.approx(0.3815, abs=5e-4)
    assert models["phenocopy(n-1)"]["fraction_of_target"] == pytest.approx(
        0.6870, abs=5e-4
    )


def test_a_third_of_the_target_forbids_the_negative_claim(result):
    """38% searched cannot support "no monogenic cause", and says so."""
    assert result.status is Status.FAIL
    f = next(x for x in result.findings if x.code == "CALLABILITY_EXPLORATORY")
    assert f.evidence["do_not_conclude"] == "no monogenic cause"
    assert "cannot be excluded" in f.evidence["permitted_claim"]
