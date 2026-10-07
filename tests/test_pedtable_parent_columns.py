"""0.1.1: parent-column spellings in a headered pedigree table."""

from __future__ import annotations

from admissible.cli import read_pedigree

BAD_PED_HEADER = (
    "#family_id\tsample_id\tdad_id\tmom_id\tsex\tphenotype\told\n"
    "F1\tNA12880\tNA12889\tNA12890\t1\t-9\twas-NA12877\n"
    "F1\tNA12878\tNA12891\tNA12892\t2\t-9\tx\n"
    "F1\tNA12879\tNA12880\tNA12878\t2\t-9\tx\n"
)


def test_dad_id_and_mom_id_columns_keep_their_parent_links(tmp_path):
    f = tmp_path / "ped.txt"
    f.write_text(BAD_PED_HEADER)
    ped = read_pedigree(f)
    assert ped.individuals["NA12879"].father == "NA12880"
    assert ped.individuals["NA12879"].mother == "NA12878"
    assert ped.individuals["NA12880"].parents == ("NA12889", "NA12890")


def test_other_id_suffixed_spellings_are_recognised(tmp_path):
    f = tmp_path / "ped.csv"
    f.write_text("fid,iid,sire_id,dam_id,sex,pheno\nF,kid,dad1,mum1,1,2\n")
    ped = read_pedigree(f)
    assert ped.individuals["kid"].parents == ("dad1", "mum1")


def test_unrecognised_parent_columns_warn_instead_of_silently_losing_links(tmp_path):
    f = tmp_path / "ped.csv"
    f.write_text("fid,iid,progenitor_a,progenitor_b,sex,pheno\nF,kid,x,y,1,2\n")
    ped = read_pedigree(f)
    assert ped.individuals["kid"].parents == ("0", "0")
    text = " ".join(ped.warnings)
    assert "no father column" in text and "no mother column" in text
    assert "parent links are lost" in text
