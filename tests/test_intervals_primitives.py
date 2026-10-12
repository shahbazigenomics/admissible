"""0.2.0 commit 1: ``subtract`` and the named-BED reader.

These are primitives only. Nothing in check 4 uses them yet, so the first group of
tests pins the behaviour that must NOT change: ``read_bed`` still ignores names and
merges across them.
"""

from __future__ import annotations

import gzip
import random

import pytest

from admissible.intervals import (
    coverage_at_least,
    intersect,
    merge,
    read_bed,
    read_named_bed,
    subtract,
    total_bp,
)


def _bed(tmp_path, name, lines, gz=False):
    p = tmp_path / name
    text = "".join(lines)
    if gz:
        p.write_bytes(gzip.compress(text.encode()))
    else:
        p.write_text(text)
    return p


def _bases(iset):
    return {(c, i) for c, v in iset.items() for s, e in v for i in range(s, e)}


# --- existing behaviour that must not change -------------------------------


def test_read_bed_still_ignores_names_and_merges_across_them(tmp_path):
    p = _bed(tmp_path, "t.bed", ["chr1\t0\t100\tA_exon1\n", "chr1\t50\t150\tA_exon2\n"])
    assert read_bed(p).intervals == {"1": [(0, 150)]}


# --- subtract ---------------------------------------------------------------


def test_subtract_basic_cases():
    a = {"1": [(0, 100)]}
    assert subtract(a, {}) == {"1": [(0, 100)]}
    assert subtract(a, {"1": [(0, 100)]}) == {}
    assert subtract(a, {"1": [(0, 200)]}) == {}
    assert subtract(a, {"1": [(40, 60)]}) == {"1": [(0, 40), (60, 100)]}
    assert subtract(a, {"1": [(0, 10)]}) == {"1": [(10, 100)]}
    assert subtract(a, {"1": [(90, 100)]}) == {"1": [(0, 90)]}
    assert subtract(a, {"1": [(-50, 10), (95, 300)]}) == {"1": [(10, 95)]}


def test_subtract_is_half_open_so_touching_intervals_remove_nothing():
    assert subtract({"1": [(0, 10)]}, {"1": [(10, 20)]}) == {"1": [(0, 10)]}
    assert subtract({"1": [(10, 20)]}, {"1": [(0, 10)]}) == {"1": [(10, 20)]}


def test_subtract_only_touches_the_matching_contig():
    a = {"1": [(0, 10)], "2": [(0, 10)]}
    assert subtract(a, {"2": [(0, 10)], "3": [(0, 99)]}) == {"1": [(0, 10)]}


def test_subtract_does_not_modify_its_arguments_and_accepts_unmerged_input():
    a = {"1": [(50, 60), (0, 40), (30, 55)]}  # unsorted, overlapping
    b = {"1": [(10, 20), (15, 25)]}
    a0, b0 = {"1": list(a["1"])}, {"1": list(b["1"])}
    assert subtract(a, b) == {"1": [(0, 10), (25, 60)]}
    assert a == a0 and b == b0


def test_subtract_empty_inputs():
    assert subtract({}, {"1": [(0, 5)]}) == {}
    assert subtract({"1": []}, {"1": [(0, 5)]}) == {}


def _random_set(rng, contigs=("1", "2"), span=60, n=6):
    out = {}
    for c in contigs:
        ivs = []
        for _ in range(rng.randrange(0, n)):
            s = rng.randrange(0, span)
            ivs.append((s, s + rng.randrange(1, 15)))
        if ivs:
            out[c] = merge(ivs)
    return out


def test_subtract_matches_a_base_by_base_reference_on_random_input():
    rng = random.Random(20261012)
    for _ in range(400):
        a, b = _random_set(rng), _random_set(rng)
        got = subtract(a, b)
        assert _bases(got) == _bases(a) - _bases(b)
        # and the output is itself sorted, merged and non-empty per contig
        for v in got.values():
            assert v and v == merge(v)


def test_subtract_and_intersect_partition_a():
    rng = random.Random(7)
    for _ in range(200):
        a, b = _random_set(rng), _random_set(rng)
        assert total_bp(subtract(a, b)) + total_bp(intersect(a, b)) == total_bp(a)


def test_subtracting_once_from_the_joint_equals_subtracting_from_every_sample():
    """The property the exclusion BED relies on, for any k of n."""
    rng = random.Random(99)
    for _ in range(150):
        sets = [_random_set(rng) for _ in range(3)]
        x = _random_set(rng)
        for k in (1, 2, 3):
            once = subtract(coverage_at_least(sets, k), x)
            each = coverage_at_least([subtract(s, x) for s in sets], k)
            assert _bases(once) == _bases(each)


# --- named BED reader -------------------------------------------------------


def test_names_and_order_are_kept(tmp_path):
    p = _bed(tmp_path, "t.bed", [
        "chr7\t100\t200\tCFTR_exon2\n",
        "chr7\t300\t400\tCFTR_exon1\n",
        "chr1\t0\t10\tNOD2_exon1\n",
    ])
    bed = read_named_bed(p)
    assert [r.name for r in bed.regions] == ["CFTR_exon2", "CFTR_exon1", "NOD2_exon1"]
    assert bed.regions[0].intervals == {"7": [(100, 200)]}
    assert [r.bp for r in bed.regions] == [100, 100, 10]
    assert bed.has_names and bed.n_unnamed == 0 and bed.problems == []


def test_overlapping_different_names_stay_separate_but_total_counts_the_union(tmp_path):
    p = _bed(tmp_path, "t.bed", ["1\t0\t100\tA\n", "1\t50\t150\tB\n"])
    bed = read_named_bed(p)
    assert [r.bp for r in bed.regions] == [100, 100]
    assert bed.total_bp == 150  # not 200


def test_one_name_on_several_lines_is_one_region(tmp_path):
    p = _bed(tmp_path, "t.bed", [
        "1\t0\t10\tX\n", "1\t5\t20\tX\n", "2\t0\t5\tX\n", "1\t100\t110\tY\n",
    ])
    bed = read_named_bed(p)
    x = bed.regions[0]
    assert x.name == "X" and x.intervals == {"1": [(0, 20)], "2": [(0, 5)]} and x.bp == 25
    assert any("repeat an earlier name" in m for m in bed.problems)


def test_unnamed_lines_become_coordinate_named_regions_and_are_counted(tmp_path):
    p = _bed(tmp_path, "t.bed", ["chr1\t0\t10\n", "chr1\t20\t30\t\n", "chr1\t40\t50\tNAMED\n"])
    bed = read_named_bed(p)
    assert [r.name for r in bed.regions] == ["1:0-10", "1:20-30", "NAMED"]
    assert [r.named for r in bed.regions] == [False, False, True]
    assert bed.n_unnamed == 2
    assert not bed.has_names
    assert any("no name in column 4" in m for m in bed.problems)


def test_malformed_lines_are_skipped_and_counted(tmp_path):
    p = _bed(tmp_path, "t.bed", [
        "# comment\n", "track name=x\n", "\n",
        "chr1\t0\n",            # too few columns
        "chr1\ta\tb\tN\n",      # not integers
        "chr1\t10\t10\tN\n",    # empty interval
        "chr1\t20\t10\tN\n",    # reversed
        "chr1\t-5\t10\tN\n",    # negative start
        "chr1\t0\t10\tGOOD\n",
    ])
    bed = read_named_bed(p)
    assert [r.name for r in bed.regions] == ["GOOD"]
    assert bed.n_malformed == 5
    assert bed.n_lines == 6


def test_crlf_gzip_and_unicode_names(tmp_path):
    p = _bed(tmp_path, "t.bed.gz", ["chr1\t0\t10\tGÈNE_exon1\r\n", "chr1\t20\t30\tB\r\n"], gz=True)
    bed = read_named_bed(p)
    assert [r.name for r in bed.regions] == ["GÈNE_exon1", "B"]
    assert [r.bp for r in bed.regions] == [10, 10]


def test_spaces_inside_a_name_are_kept(tmp_path):
    p = _bed(tmp_path, "t.bed", ["1\t0\t10\texon 4 of GENE\n"])
    assert read_named_bed(p).regions[0].name == "exon 4 of GENE"


@pytest.mark.parametrize("reader", [read_bed, read_named_bed])
def test_missing_file_and_directory_are_problems_not_exceptions(tmp_path, reader):
    assert reader(tmp_path / "nope.bed").problems
    d = tmp_path / "dir.bed"
    d.mkdir()
    assert reader(d).problems


@pytest.mark.parametrize("reader", [read_bed, read_named_bed])
def test_empty_file_is_empty_not_an_error(tmp_path, reader):
    p = tmp_path / "e.bed"
    p.write_text("")
    r = reader(p)
    assert r.n_lines == 0 and r.problems == []


@pytest.mark.parametrize("reader", [read_bed, read_named_bed])
def test_corrupt_gzip_never_raises(tmp_path, reader):
    """A valid gzip header with a corrupt body raises zlib.error, which is not an
    OSError. read_bed used to let it through; found while writing this module."""
    good = gzip.compress(("chr1\t1\t100\tA\n" * 5000).encode())
    flipped = good[:30] + bytes([good[30] ^ 0xFF]) + good[31:]
    header_then_junk = b"\x1f\x8b\x08\x00" + b"\x07" * 400
    for i, blob in enumerate((flipped, good[: len(good) // 2], header_then_junk)):
        p = tmp_path / f"c{i}.bed.gz"
        p.write_bytes(blob)
        assert reader(p).problems


@pytest.mark.parametrize("reader", [read_bed, read_named_bed])
def test_byte_soup_and_huge_numbers_never_raise(tmp_path, reader):
    rng = random.Random(1)
    for i in range(40):
        blob = bytes(rng.randrange(256) for _ in range(rng.randrange(0, 3000)))
        for suffix in (".bed", ".bed.gz"):
            p = tmp_path / f"s{i}{suffix}"
            p.write_bytes(blob)
            reader(p)
    p = tmp_path / "huge.bed"
    p.write_text(
        "1\t0\t" + "9" * 40 + "\tBIG\n"       # larger than any genome, still an int
        "1\t0\t" + "9" * 6000 + "\tTOOBIG\n"  # beyond Python's int-from-str limit
        "1\t0\t10\t" + "N" * 100_000 + "\n"   # absurdly long name
        "1\t\x00\x00\t5\tNUL\n"
    )
    out = reader(p)
    assert out.n_lines == 4


def test_a_random_named_bed_agrees_with_read_bed_on_the_union(tmp_path):
    """Names must not change which bases are in the target."""
    rng = random.Random(5)
    for t in range(60):
        lines = []
        for _ in range(rng.randrange(1, 12)):
            s = rng.randrange(0, 500)
            lines.append(
                f"{rng.choice(['1', 'chr2'])}\t{s}\t{s + rng.randrange(1, 40)}\t"
                f"N{rng.randrange(0, 5)}\n"
            )
        p = _bed(tmp_path, f"r{t}.bed", lines)
        assert read_named_bed(p).total_bp == total_bp(read_bed(p).intervals)
