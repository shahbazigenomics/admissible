"""Read a pedigree from the spreadsheet it is actually kept in.

PED is a positional format with numeric codes, and a positional format with
numeric codes is a format people get wrong.  ``2`` meaning *affected* and ``1``
meaning *unaffected* is the single most common mistake, and it fails silently:
check 5 simply finds nothing to segregate.  Column order matters too, so a
column inserted in the middle corrupts every row without any parse error.

So this module accepts the thing a lab already has - a table with a header row,
written in words.  Columns are located by name and values are read as a person
would write them (``M``, ``female``, ``yes``, ``affected``, a blank for a founder).
Nothing is guessed: a value that is not recognised becomes *unknown* and is
listed in ``Pedigree.warnings`` rather than being rounded to something plausible.

CSV, TSV and Excel are all read; Excel only if ``openpyxl`` happens to be
installed, since this package has no runtime dependencies and will not grow one
for a convenience.
"""

from __future__ import annotations

import csv
import os

from .ped import Affection, Individual, Pedigree, Sex

# Header spellings seen in real pedigree spreadsheets, lower-cased and stripped
# of spaces, underscores and dots before matching.
_ID = ("id", "iid", "sample", "sampleid", "samplename", "individual",
       "individualid", "person", "subject", "name")
_FATHER = ("father", "pat", "paternalid", "fatherid", "dad", "sire", "pathid")
_MOTHER = ("mother", "mat", "maternalid", "motherid", "mum", "mom", "dam")
_SEX = ("sex", "gender")
_AFFECTED = ("affected", "affection", "pheno", "phenotype", "status",
             "affectedstatus", "disease")
_FAMILY = ("family", "fid", "familyid", "fam", "pedigree", "kindred")

_MALE = {"1", "m", "male", "man", "boy", "father", "son", "xy"}
_FEMALE = {"2", "f", "female", "woman", "girl", "mother", "daughter", "xx"}
_YES = {"2", "y", "yes", "affected", "case", "true", "t", "patient", "ill"}
_NO = {"1", "n", "no", "unaffected", "control", "false", "f", "healthy", "normal"}
_BLANK = {"", ".", "0", "-", "na", "n/a", "nan", "none", "null", "unknown", "?"}


def _norm(text: str) -> str:
    return "".join(ch for ch in text.strip().lower() if ch.isalnum())


def _find(header: list[str], names: tuple[str, ...]) -> int | None:
    normed = [_norm(h) for h in header]
    for want in names:
        if want in normed:
            return normed.index(want)
    return None


def looks_like_ped_table(path: str | os.PathLike[str]) -> bool:
    """True for a table this module should read rather than ``read_ped``.

    Decided by content, not by file extension: a PED saved as ``family.txt`` and
    a spreadsheet exported as ``family.ped`` are both things people do.  The test
    is whether the first meaningful line names a sample-id column - which is what
    a header row has and a PED data row does not.
    """
    name = str(path).lower()
    if name.endswith((".xlsx", ".xls")):
        return True
    try:
        with open(path, errors="replace") as fh:
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                line = line.lstrip("#").strip()
                if not line:
                    continue
                tokens = [t for part in line.split() for t in part.split(",")]
                tokens = [t for tok in tokens for t in tok.split(";")]
                return any(_norm(t) in _ID for t in tokens)
    except OSError:
        return False
    return False


def _rows(path: str | os.PathLike[str]) -> list[list[str]]:
    name = str(path).lower()
    if name.endswith((".xlsx", ".xls")):
        try:
            from openpyxl import load_workbook  # noqa: PLC0415
        except ImportError as exc:  # pragma: no cover - depends on environment
            raise OSError(
                "reading .xlsx needs openpyxl (pip install openpyxl); "
                "or save the sheet as CSV, which needs nothing"
            ) from exc
        wb = load_workbook(path, read_only=True, data_only=True)
        return [
            ["" if c is None else str(c) for c in row]
            for row in wb[wb.sheetnames[0]].iter_rows(values_only=True)
        ]
    with open(path, newline="", errors="replace") as fh:
        sample = fh.read(8192)
        fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;|")
        except csv.Error:
            dialect = csv.excel_tab if "\t" in sample else csv.excel
        return [list(r) for r in csv.reader(fh, dialect)]


def read_ped_table(path: str | os.PathLike[str]) -> Pedigree:
    """Read a headered pedigree table.  Never raises; problems become warnings."""
    ped = Pedigree()
    try:
        rows = _rows(path)
    except OSError as exc:
        ped.warnings.append(f"could not read {path}: {exc}")
        return ped

    rows = [r for r in rows if any(str(c).strip() for c in r)]
    if not rows:
        ped.warnings.append(f"{path}: no rows")
        return ped

    header = [str(c) for c in rows[0]]
    i_id = _find(header, _ID)
    if i_id is None:
        ped.warnings.append(
            f"{path}: no sample-id column found; expected a header naming one of "
            f"{', '.join(_ID[:5])}"
        )
        return ped
    i_fa, i_mo = _find(header, _FATHER), _find(header, _MOTHER)
    i_sex, i_aff = _find(header, _SEX), _find(header, _AFFECTED)
    i_fam = _find(header, _FAMILY)
    for label, idx in (("sex", i_sex), ("affected", i_aff)):
        if idx is None:
            ped.warnings.append(
                f"{path}: no {label} column; every sample will be {label}-unknown"
            )

    def cell(row: list[str], idx: int | None) -> str:
        return str(row[idx]).strip() if idx is not None and idx < len(row) else ""

    def blank(value: str) -> bool:
        return _norm(value) in _BLANK or value.strip() == ""

    for lineno, row in enumerate(rows[1:], 2):
        iid = cell(row, i_id)
        if not iid:
            continue
        if iid in ped.individuals:
            ped.warnings.append(f"row {lineno}: duplicate sample id {iid!r}")
            continue

        raw_sex = cell(row, i_sex)
        sex = Sex.UNKNOWN
        if not blank(raw_sex):
            key = _norm(raw_sex)
            if key in _MALE:
                sex = Sex.MALE
            elif key in _FEMALE:
                sex = Sex.FEMALE
            else:
                ped.warnings.append(
                    f"row {lineno}: sex {raw_sex!r} not understood for {iid}; "
                    f"recorded as unknown rather than guessed"
                )

        raw_aff = cell(row, i_aff)
        affection = Affection.UNKNOWN
        if not blank(raw_aff):
            key = _norm(raw_aff)
            if key in _YES:
                affection = Affection.AFFECTED
            elif key in _NO:
                affection = Affection.UNAFFECTED
            else:
                ped.warnings.append(
                    f"row {lineno}: affected status {raw_aff!r} not understood for "
                    f"{iid}; recorded as unknown rather than guessed"
                )

        father = cell(row, i_fa)
        mother = cell(row, i_mo)
        family = cell(row, i_fam)
        ped.individuals[iid] = Individual(
            family_id=family if not blank(family) else "FAM",
            iid=iid,
            father="0" if blank(father) else father,
            mother="0" if blank(mother) else mother,
            sex=sex,
            affection=affection,
        )

    for ind in list(ped.individuals.values()):
        # Their own parent. Arises when a relabelling renames a row but not the
        # references to it, and it is not harmless: the kinship recursion takes
        # it at face value and returns 0.5 for a pair that is nothing of the
        # kind, so every reconciliation downstream is measured against a number
        # the pedigree never meant.
        if ind.iid in ind.parents:
            ped.warnings.append(
                f"{ind.iid} is listed as their own parent; the relationships "
                f"computed through them cannot be trusted"
            )
    for ind in list(ped.individuals.values()):
        for parent in ind.parents:
            # "0" is the absent-parent sentinel, not a person who is missing.
            if parent != "0" and parent not in ped.individuals:
                ped.warnings.append(
                    f"{ind.iid} names parent {parent!r}, who has no row of their own"
                )
    if not ped.individuals:
        ped.warnings.append(f"{path}: a header was found but no sample rows")
    return ped


def write_ped(ped: Pedigree, out) -> None:
    """Write a Pedigree back out as a plain PED, for tools that want one."""
    sex_code = {Sex.MALE: "1", Sex.FEMALE: "2", Sex.UNKNOWN: "0"}
    aff_code = {Affection.UNAFFECTED: "1", Affection.AFFECTED: "2", Affection.UNKNOWN: "0"}
    out.write("#FID\tIID\tPAT\tMAT\tSEX\tPHENO\n")
    for ind in ped.individuals.values():
        out.write(
            f"{ind.family_id}\t{ind.iid}\t{ind.father}\t{ind.mother}\t"
            f"{sex_code[ind.sex]}\t{aff_code[ind.affection]}\n"
        )
