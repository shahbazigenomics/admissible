"""Giving the tool a pedigree without writing a PED by hand.

PED is positional and numeric: ``2`` means affected, ``1`` means unaffected, and
getting that backwards produces no error at all - check 5 just finds nothing to
segregate. Column order matters too. So the pedigree can also be supplied as the
table a lab already keeps, read by column name and in words, and the format is
decided by content rather than by file extension.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from admissible.cli import main, read_pedigree
from admissible.ped import Affection, Sex
from admissible.pedtable import read_ped_table


def write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


def test_a_spreadsheet_in_words_is_read(tmp_path):
    p = write(tmp_path, "family.csv", (
        "Family,Sample ID,Father,Mother,Sex,Affected\n"
        "IPC,IPC-1,,,Male,no\n"
        "IPC,IPC-2,,,Female,No\n"
        "IPC,IPC-3,IPC-1,IPC-2,M,YES\n"
        "IPC,IPC-4,IPC-1,IPC-2,f,affected\n"
    ))
    ped = read_ped_table(p)
    assert ped.warnings == []
    assert list(ped.individuals) == ["IPC-1", "IPC-2", "IPC-3", "IPC-4"]
    assert ped.individuals["IPC-3"].sex is Sex.MALE
    assert ped.individuals["IPC-4"].sex is Sex.FEMALE
    assert ped.individuals["IPC-4"].affection is Affection.AFFECTED
    assert ped.individuals["IPC-1"].affection is Affection.UNAFFECTED
    assert ped.individuals["IPC-1"].family_id == "IPC"
    # The relationships have to survive the translation, not just the fields.
    assert ped.kinship("IPC-3", "IPC-1") == pytest.approx(0.25)
    assert ped.kinship("IPC-3", "IPC-4") == pytest.approx(0.25)
    assert ped.kinship("IPC-1", "IPC-2") == pytest.approx(0.0)


def test_other_column_spellings_and_tabs(tmp_path):
    p = write(tmp_path, "f.tsv", (
        "sample\tdad\tmum\tgender\tstatus\n"
        "A\t.\t.\t1\tcontrol\n"
        "B\t.\t.\t2\tcontrol\n"
        "C\tA\tB\tM\tcase\n"
    ))
    ped = read_ped_table(p)
    assert ped.warnings == []
    assert ped.individuals["C"].affection is Affection.AFFECTED
    assert ped.kinship("C", "A") == pytest.approx(0.25)


def test_an_unreadable_value_becomes_unknown_and_says_so(tmp_path):
    """Never round a value to the plausible one: that is how 1/2 goes wrong."""
    p = write(tmp_path, "f.csv", (
        "sample,sex,affected\n"
        "A,intersex,probably\n"
    ))
    ped = read_ped_table(p)
    assert ped.individuals["A"].sex is Sex.UNKNOWN
    assert ped.individuals["A"].affection is Affection.UNKNOWN
    assert any("intersex" in w for w in ped.warnings)
    assert any("probably" in w for w in ped.warnings)


def test_a_founder_is_not_a_missing_parent(tmp_path):
    """A blank parent is "not in the study", not "a person whose row is lost"."""
    p = write(tmp_path, "f.csv", "sample,father,mother,sex,affected\nA,,,M,no\n")
    assert read_ped_table(p).warnings == []


def test_a_parent_with_no_row_is_reported(tmp_path):
    p = write(tmp_path, "f.csv", "sample,father,mother\nX,GHOST,\n")
    assert any("GHOST" in w for w in read_ped_table(p).warnings)


def test_a_missing_id_column_is_refused_not_guessed(tmp_path):
    p = write(tmp_path, "f.csv", "who,dad,mum\nA,,\n")
    ped = read_ped_table(p)
    assert ped.individuals == {}
    assert any("no sample-id column" in w for w in ped.warnings)


# --- format is decided by content, not by extension -----------------------


def test_a_headerless_ped_still_reads_as_a_ped(tmp_path):
    p = write(tmp_path, "h.ped", "FAM\tA\t0\t0\t1\t1\nFAM\tB\t0\t0\t2\t2\n")
    ped = read_pedigree(Path(p))
    assert list(ped.individuals) == ["A", "B"]
    assert ped.individuals["B"].affection is Affection.AFFECTED


def test_a_spreadsheet_saved_as_ped_is_still_read(tmp_path):
    """People name files whatever they like; the content is what decides."""
    p = write(tmp_path, "sheet.ped", "sample,father,mother,sex,affected\nA,,,M,yes\n")
    ped = read_pedigree(Path(p))
    assert ped.individuals["A"].affection is Affection.AFFECTED


def test_a_ped_saved_as_txt_is_still_read(tmp_path):
    p = write(tmp_path, "fam.txt", "FAM\tA\t0\t0\t1\t2\n")
    ped = read_pedigree(Path(p))
    assert ped.individuals["A"].affection is Affection.AFFECTED


# --- the template ---------------------------------------------------------


def test_the_template_carries_the_real_sample_names(tmp_path, capsys):
    """The ids have to match the VCF sample names, and nobody can see those
    without looking - which is why the commonest first-run failure is a PED
    whose ids match the lab's names for people instead."""
    head = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tpatient_03\tmum\n"
    )
    v = write(tmp_path, "fam.vcf", head + "chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\n")

    assert main(["ped-template", str(v)]) == 0
    out = capsys.readouterr().out
    assert "patient_03" in out and "mum" in out
    assert out.startswith("#FID\tIID\tPAT\tMAT\tSEX\tPHENO\n")

    # And the filled-in template must round-trip back in.
    filled = write(tmp_path, "filled.ped", out.replace(
        "FAM1\tpatient_03\t0\t0\t0\t0", "FAM1\tpatient_03\tdad\tmum\t1\t2"
    ))
    ped = read_pedigree(Path(filled))
    assert ped.individuals["patient_03"].affection is Affection.AFFECTED
    assert any("dad" in w for w in ped.warnings)      # named but not given a row


def test_the_csv_template_is_the_spreadsheet_shape(tmp_path, capsys):
    head = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    )
    v = write(tmp_path, "a.vcf", head + "chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n")
    assert main(["ped-template", str(v), "--csv"]) == 0
    out = capsys.readouterr().out
    assert out.splitlines()[0] == "sample,father,mother,sex,affected"
    assert out.splitlines()[1].startswith("S1,")
    # What it prints is what the reader accepts.
    p = write(tmp_path, "t.csv", out.replace("S1,,,,", "S1,,,M,yes"))
    assert read_ped_table(p).individuals["S1"].affection is Affection.AFFECTED


def test_a_person_listed_as_their_own_parent_is_reported(tmp_path):
    """Found while demonstrating the tool on CEPH 1463.

    Renaming a row without renaming the references to it leaves someone as
    their own mother. The kinship recursion takes that at face value and
    returns 0.5 for a pair that is nothing of the kind, so every reconciliation
    downstream is measured against a number the pedigree never meant - and
    nothing said so.
    """
    p = write(tmp_path, "self.csv", (
        "sample,father,mother,sex,affected\n"
        "MUM,,,F,no\n"
        "KID,DAD,KID,F,yes\n"
    ))
    ped = read_ped_table(p)
    assert any("own parent" in w for w in ped.warnings)


def test_the_ped_reader_flags_it_too(tmp_path):
    from admissible.ped import read_ped

    p = write(tmp_path, "self.ped", "F\tA\tA\t0\t1\t2\n")
    assert any("own parent" in w for w in read_ped(p).warnings)
