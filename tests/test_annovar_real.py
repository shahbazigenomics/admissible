"""The ANNOVAR adapter against genotypes it did not invent.

Every other test of this adapter feeds it a table written alongside the code,
which cannot catch the adapter misunderstanding what ANNOVAR actually writes.
The fixture here carries real GATK genotypes inside a reconstructed ANNOVAR
column layout - see ``tests/data/annovar/PROVENANCE.md`` for exactly which half
is which, and why no fully real file was available.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from admissible.annovar import iter_multianno_records, scan_multianno
from admissible.checks.genotype import check_genotype
from admissible.contigs import normalize_contig
from admissible.vcfio import HET, HOMALT, HOMREF, MISSING, SiteIndex

DATA = Path(__file__).parent / "data" / "annovar"
TABLE = DATA / "gatk_chr21.hg38_multianno.txt"
TRUTH = json.loads((DATA / "expected_genotypes.json").read_text())

_WANT = {"0/1": HET, "1/0": HET, "0|1": HET, "1|0": HET, "1/2": HET, "2/1": HET,
         "1/1": HOMALT, "1|1": HOMALT, "2/2": HOMALT, "0/0": HOMREF, "0|0": HOMREF}


@pytest.fixture(scope="module")
def scanned():
    index = SiteIndex()
    scan = scan_multianno(TABLE, index)
    return index, scan


def test_every_genotype_matches_the_one_in_the_file(scanned):
    """113 of 120 was not good enough, and the 7 were all the same kind."""
    index, scan = scanned
    codes = next(iter(scan.genotypes.values()))
    wrong = []
    for chrom, pos, ref, alt, gt in TRUTH:
        site = index.get((normalize_contig(chrom), pos, ref.upper(), alt.upper()))
        got = codes.get(site, MISSING) if site is not None else MISSING
        if got != _WANT[gt]:
            wrong.append((chrom, pos, gt, got))
    assert wrong == []


def test_a_multiallelic_heterozygote_is_not_a_homozygote(scanned):
    """AC=1,1 sums to two, so reconstructing from AC/AN turns 1/2 into 1/1.

    This is the failure the whole tool exists to prevent, and it was being
    committed by the reader before anything else got a chance to look.
    """
    index, scan = scanned
    codes = next(iter(scan.genotypes.values()))
    multiallelic = [t for t in TRUTH if "," in t[3]]
    assert len(multiallelic) == 7
    for chrom, pos, ref, alt, gt in multiallelic:
        assert gt == "1/2"
        site = index.get((normalize_contig(chrom), pos, ref.upper(), alt.upper()))
        assert codes[site] == HET


def test_the_preserved_format_block_is_reported_as_present(scanned):
    """`--vcfinput` keeps FORMAT and the sample column; saying otherwise made
    check 2 fall back to MLEAC while GQ, PL and AD sat unread in the file."""
    _, scan = scanned
    assert scan.stats.has_format_column is True
    assert {"GT", "AD", "DP", "GQ", "PL"} <= scan.stats.format_keys_seen


def test_check_2_can_examine_an_annovar_delivery(scanned):
    """The evidence is there, so the direct rules must fire, not the fallback."""
    res = check_genotype([TABLE])
    flags: dict[str, int] = {}
    for v in res.metrics["per_sample"].values():
        for k, n in v["flags"].items():
            flags[k] = flags.get(k, 0) + n
    # All three arms, none of which could fire when FORMAT was discarded.
    assert flags.get("HOM_CONTRADICTED_BY_LIKELIHOOD", 0) > 0   # needs PL
    assert flags.get("LOW_GQ", 0) > 0                           # needs GQ
    assert flags.get("AB_SKEW", 0) > 0                          # needs AD
    # And MLEAC stays out of it: it is only allowed to condemn a genotype when
    # there is no direct evidence to condemn it with.
    assert flags.get("MLEAC_CONTRADICTS_GT", 0) == 0


def test_records_stream_with_their_real_evidence(scanned):
    """iter_multianno_records is what check 2 reads; it must carry the cell."""
    rec = next(iter_multianno_records(TABLE))
    cell = next(iter(rec.samples.values()))
    assert cell["GT"]
    assert "DP" in cell and "PL" in cell


def test_a_table_without_vcfinput_still_falls_back_to_ac_an(tmp_path):
    """Tables produced without --vcfinput genuinely have no FORMAT block, and
    the AC/AN reconstruction remains the only thing available there."""
    p = tmp_path / "noformat.hg19_multianno.txt"
    p.write_text(
        "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\t"
        "Otherinfo1\tOtherinfo2\tOtherinfo3\tOtherinfo4\tOtherinfo5\t"
        "Otherinfo6\tOtherinfo7\tOtherinfo8\n"
        "chr1\t100\t100\tA\tG\texonic\tGENEX\t"
        "chr1\t100\t.\tA\tG\t50\tPASS\tAC=1;AN=2\n"
    )
    index = SiteIndex()
    scan = scan_multianno(p, index)
    assert scan.stats.has_format_column is False
    assert scan.stats.format_keys_seen == set()
    codes = next(iter(scan.genotypes.values()))
    assert codes[index.get(("1", 100, "A", "G"))] == HET
