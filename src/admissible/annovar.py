"""Read ANNOVAR ``*_multianno.txt`` tables as if they were single-sample VCFs.

Why this exists: in practice the VCFs are often the thing that got lost.  A
delivered analysis frequently consists of annotated tables, with the call sets
they came from long gone - and an identity audit that cannot read the only
files you still have is useless precisely when it is needed most.

ANNOVAR keeps the original VCF columns in its ``Otherinfo`` block, so CHROM,
POS, REF, ALT, FILTER and INFO all survive.  What does *not* survive is the
FORMAT/sample block, so the genotype has to be reconstructed - and it can be,
exactly, for a single-sample call set: with ``AN=2``, ``AC=1`` is a heterozygote
and ``AC=2`` a homozygous alternate.

That reconstruction is only valid for a single-sample VCF.  In a multi-sample
call set AC is a cohort allele count and says nothing about any one individual,
so this reader will not assume: it checks ``AN`` on every record and marks any
record with ``AN != 2`` as unassessable.

Tables stripped down to ``AC`` alone, with no ``AN`` at all, do occur.  There the
assumption is tested rather than made: AC is used only if it never exceeds 2
anywhere in the file, which a cohort allele count in a family study would.  If it
does exceed 2, every genotype derived that way is discarded.
"""

from __future__ import annotations

import os

from .contigs import is_primary, normalize_contig
from .vcfio import (
    HET,
    HOMALT,
    HOMREF,
    MISSING,
    SiteIndex,
    VcfHeader,
    VcfScan,
    VcfStats,
    _sha256,
    encode_gt,
)

REGION_COLUMNS = ("Func.refGene", "Func.refGeneWithVer", "Func.knownGene")


def looks_like_multianno(path: str | os.PathLike[str]) -> bool:
    name = os.path.basename(str(path)).lower()
    return name.endswith(".txt") and "multianno" in name


def default_label(path: str | os.PathLike[str]) -> str:
    """``.../FAM01/03/x_multianno.txt`` -> ``FAM01-03``.

    These tables carry no sample name, so the directory layout is the only
    identifier available.  It is also usually the one the analyst actually uses.
    """
    p = os.path.normpath(str(path))
    parent = os.path.basename(os.path.dirname(p))
    grandparent = os.path.basename(os.path.dirname(os.path.dirname(p)))
    if parent and grandparent:
        return f"{grandparent}-{parent}"
    return parent or os.path.basename(p)


def _build_from_filename(path: str) -> str | None:
    low = os.path.basename(path).lower()
    if "hg19" in low or "grch37" in low:
        return "GRCh37"
    if "hg38" in low or "grch38" in low:
        return "GRCh38"
    return None


def _locate_columns(header: list[str], first_row: list[str]) -> dict[str, int] | None:
    """Find the original VCF columns inside the Otherinfo block.

    Located by content rather than by a fixed offset, because the number of
    Otherinfo columns varies with the ANNOVAR invocation.
    """
    other = [i for i, h in enumerate(header) if h.startswith("Otherinfo")]
    if not other:
        return None
    last = max(other)
    # Validate a candidate layout against the data rather than trusting a fixed
    # offset: an INFO field stripped down to a single key (it happens) has no
    # semicolon, so the anchor has to be the CHROM/POS pair, not the INFO shape.
    for i in reversed(other):
        if i >= len(first_row) or "=" not in first_row[i] or i - 7 < 0:
            continue
        cols = {
            "chrom": i - 7, "pos": i - 6, "ref": i - 4,
            "alt": i - 3, "filter": i - 1, "info": i,
        }
        if max(cols.values()) >= len(first_row):
            continue
        if not (first_row[cols["pos"]].isdigit() and first_row[cols["chrom"]]):
            continue
        # `table_annovar.pl --vcfinput` - the standard invocation - keeps the
        # original FORMAT and sample columns immediately after INFO.  They are
        # the genotype, verbatim; reconstructing one from AC/AN while the real
        # thing sits two columns to the right is how a 1/2 call becomes a
        # homozygote.
        if i + 2 <= last and i + 2 < len(first_row):
            fmt = first_row[i + 1].strip()
            keys = fmt.split(":")
            if keys and keys[0] == "GT" and _looks_like_genotype(first_row[i + 2]):
                cols["format"] = i + 1
                cols["sample"] = i + 2
        return cols
    return None


def _looks_like_genotype(cell: str) -> bool:
    gt = cell.split(":")[0].strip()
    if not gt or len(gt) > 7:
        return False
    sep = "/" if "/" in gt else "|" if "|" in gt else None
    parts = gt.split(sep) if sep else [gt]
    return all(p == "." or p.isdigit() for p in parts) and bool(parts)


def _ac_total(ac: str) -> int | None:
    total = 0
    for part in ac.split(","):
        part = part.strip()
        if not part.isdigit():
            return None
        total += int(part)
    return total


def genotype_from_info(info: dict[str, str]) -> tuple[int, int | None]:
    """Reconstruct a single-sample genotype from AC/AN, or AF as a fallback.

    Returns ``(code, provisional_ac)``.  A non-None ``provisional_ac`` means the
    genotype was derived from AC with no AN to confirm the call set is
    single-sample; the caller must validate that assumption across the whole
    file before trusting it, because in a multi-sample VCF AC is a cohort allele
    count and carries no information about any one individual.
    """
    an = info.get("AN")
    if an is not None:
        if an.strip() != "2":
            return MISSING, None  # not a diploid single-sample record
        total = _ac_total(info.get("AC", ""))
        if total is None:
            return MISSING, None
        return {0: HOMREF, 1: HET, 2: HOMALT}.get(total, MISSING), None

    if "AC" in info:
        total = _ac_total(info["AC"])
        if total is not None:
            return {0: HOMREF, 1: HET, 2: HOMALT}.get(total, MISSING), total

    af = info.get("AF", "").split(",")[0].strip()
    try:
        value = float(af)
    except ValueError:
        return MISSING, None
    if abs(value - 0.5) < 0.1:
        return HET, None
    if value >= 0.9:
        return HOMALT, None
    return MISSING, None


GENE_COLUMNS = ("Gene.refGene", "Gene.refGeneWithVer", "Gene.knownGene")


def scan_multianno(
    path: str | os.PathLike[str],
    index: SiteIndex,
    label: str | None = None,
    af_key: str | None = None,
) -> VcfScan:
    """Read one multianno table into ``index``.  Never raises."""
    sample = label or default_label(path)
    header_obj = VcfHeader(samples=[sample])
    stats = VcfStats(has_format_column=False)
    size, digest = _sha256(path)
    scan = VcfScan(
        path=str(path), header=header_obj, stats=stats,
        build=_build_from_filename(str(path)), file_size=size, file_sha256=digest,
    )
    scan.genotypes[sample] = {}
    scan.has_explicit_homref[sample] = False

    try:
        fh = open(path, errors="replace")
    except OSError as exc:
        scan.problems.append(f"cannot open: {exc}")
        return scan

    with fh:
        try:
            head = fh.readline().rstrip("\n").split("\t")
        except OSError as exc:
            scan.problems.append(f"read stopped early: {exc}")
            return scan
        region_col = next((head.index(c) for c in REGION_COLUMNS if c in head), None)
        if region_col is not None:
            stats.annotation_source = f"ANNOVAR:{head[region_col]}"
        gene_col = next((head.index(c) for c in GENE_COLUMNS if c in head), None)
        af_col = head.index(af_key) if af_key and af_key in head else None

        cols: dict[str, int] | None = None
        n_unassessable = 0
        provisional: list[tuple[int, int]] = []
        max_ac = 0
        try:
            for lineno, raw in enumerate(fh, 2):
                row = raw.rstrip("\n").split("\t")
                if len(row) < 6:
                    stats.n_malformed += 1
                    continue
                if cols is None:
                    cols = _locate_columns(head, row)
                    if cols is None:
                        scan.problems.append(
                            "could not locate the original VCF columns in the Otherinfo "
                            "block; genotypes cannot be reconstructed from this table"
                        )
                        return scan
                try:
                    chrom_raw = row[cols["chrom"]]
                    pos = int(row[cols["pos"]])
                    ref = row[cols["ref"]].upper()
                    alt = row[cols["alt"]].upper()
                    filt = row[cols["filter"]].strip()
                    info_raw = row[cols["info"]]
                except (IndexError, ValueError):
                    stats.n_malformed += 1
                    if len(scan.problems) < 20:
                        scan.problems.append(f"line {lineno}: unparseable record")
                    continue

                stats.n_records += 1
                chrom = normalize_contig(chrom_raw)
                stats.contig_counts[chrom] = stats.contig_counts.get(chrom, 0) + 1
                if not is_primary(chrom_raw):
                    stats.n_nonprimary_contig += 1
                if filt in ("", "."):
                    stats.n_filter_absent += 1
                elif filt == "PASS":
                    stats.n_pass += 1
                else:
                    stats.n_nonpass += 1
                if "," in alt:
                    stats.n_multiallelic += 1

                info: dict[str, str] = {}
                for kv in info_raw.split(";"):
                    k, _, v = kv.partition("=")
                    stats.info_keys_seen.add(k)
                    if k in ("AC", "AN", "AF", "MLEAC", "MLEAF", "DP"):
                        info[k] = v
                    elif af_col is None and af_key and k == af_key:
                        info["__af"] = v

                if region_col is not None and region_col < len(row):
                    cls = row[region_col].strip()
                    if cls:
                        stats.region_counts[cls] = stats.region_counts.get(cls, 0) + 1

                cell: dict[str, str] = {}
                if "sample" in cols and cols["sample"] < len(row):
                    keys = row[cols["format"]].split(":")
                    vals = row[cols["sample"]].split(":")
                    cell = dict(zip(keys, vals, strict=False))
                    stats.format_keys_seen.update(keys)
                    stats.has_format_column = True

                if cell.get("GT"):
                    code, provisional_ac = encode_gt(cell["GT"]), None
                    if code == HOMREF:
                        scan.has_explicit_homref[sample] = True
                else:
                    code, provisional_ac = genotype_from_info(info)
                if code == MISSING:
                    n_unassessable += 1
                    continue
                site_id = index.id_for((chrom, pos, ref, alt))
                if gene_col is not None and gene_col < len(row):
                    g = row[gene_col].strip()
                    if g and g != ".":
                        scan.site_gene[site_id] = g.split(",")[0].split(";")[0]
                raw_af = (
                    row[af_col] if af_col is not None and af_col < len(row) else info.get("__af")
                )
                if raw_af:
                    try:
                        scan.site_af[site_id] = float(raw_af.split(",")[0])
                    except ValueError:
                        pass
                if provisional_ac is not None:
                    max_ac = max(max_ac, provisional_ac)
                    provisional.append((site_id, code))
                    continue
                if code == HOMREF:
                    scan.has_explicit_homref[sample] = True
                scan.genotypes[sample][site_id] = code
        except OSError as exc:
            scan.problems.append(f"read stopped early: {type(exc).__name__}: {exc}")

    if provisional:
        # No AN anywhere.  If AC never exceeds 2 across the whole table, the call
        # set is diploid single-sample and AC does determine the genotype.  If it
        # ever exceeds 2 it is a cohort count and every such genotype is void.
        if max_ac <= 2:
            for site_id, code in provisional:
                if code == HOMREF:
                    scan.has_explicit_homref[sample] = True
                scan.genotypes[sample][site_id] = code
            scan.problems.append(
                f"no AN field in this table; genotypes for {len(provisional)} record(s) "
                f"were reconstructed from AC alone, which is valid only because AC never "
                f"exceeds 2 here (max observed {max_ac}), i.e. the call set is "
                f"single-sample diploid"
            )
        else:
            n_unassessable += len(provisional)
            scan.problems.append(
                f"no AN field and AC reaches {max_ac}, so this is a cohort allele count, "
                f"not a per-sample genotype; {len(provisional)} record(s) are unassessable"
            )

    if n_unassessable:
        scan.problems.append(
            f"{n_unassessable} record(s) had no usable AN=2/AC, so no genotype could be "
            f"reconstructed for them"
        )
    # A table produced without --vcfinput genuinely has no FORMAT block, and
    # there depth, GQ and allele balance are not recoverable - a real limitation
    # of the file, not a defect. Where ANNOVAR did keep the sample column, those
    # fields are present and are reported as present.
    if not scan.stats.has_format_column:
        scan.stats.format_keys_seen = set()
    return scan


def iter_multianno_records(path: str | os.PathLike[str], label: str | None = None):
    """Yield :class:`~admissible.vcfio.Record` objects from a multianno table.

    The FORMAT block is gone in these files, so allele balance, GQ and PL are
    simply not recoverable - but ``DP``, ``MLEAC``, ``QD``, ``MQ``, ``FS`` and
    ``SOR`` all survive in INFO, and the genotype is reconstructable.  That is
    enough to run the depth-based half of the genotype check on a call set whose
    VCFs no longer exist, which is a common and otherwise hopeless situation.
    """
    from .vcfio import HET, HOMALT, HOMREF, MISSING, Record

    sample = label or default_label(path)
    gt_text = {HOMREF: "0/0", HET: "0/1", HOMALT: "1/1"}
    try:
        fh = open(path, errors="replace")
    except OSError:
        return
    with fh:
        head = fh.readline().rstrip("\n").split("\t")
        cols = None
        for lineno, raw in enumerate(fh, 2):
            row = raw.rstrip("\n").split("\t")
            if len(row) < 6:
                continue
            if cols is None:
                cols = _locate_columns(head, row)
                if cols is None:
                    return
            try:
                chrom_raw = row[cols["chrom"]]
                pos = int(row[cols["pos"]])
                ref = row[cols["ref"]].upper()
                alt = row[cols["alt"]].upper()
                filt = row[cols["filter"]].strip()
                info_raw = row[cols["info"]]
            except (IndexError, ValueError):
                continue
            info = {}
            for kv in info_raw.split(";"):
                k, _, v = kv.partition("=")
                info[k] = v

            cell: dict[str, str] = {}
            if "sample" in cols and cols["sample"] < len(row):
                cell = dict(
                    zip(
                        row[cols["format"]].split(":"),
                        row[cols["sample"]].split(":"),
                        strict=False,
                    )
                )
            if cell.get("GT"):
                # The real genotype, with the real DP/GQ/AD/PL beside it, so
                # check 2 can examine an ANNOVAR-delivered analysis instead of
                # declaring it unassessable.
                code = encode_gt(cell["GT"])
            else:
                code, _provisional = genotype_from_info(info)
                cell = {"GT": gt_text[code]}
            if code == MISSING:
                continue
            yield Record(
                chrom=normalize_contig(chrom_raw),
                pos=pos,
                ref=ref,
                alts=alt.split(","),
                qual=None,
                filters=[x for x in filt.split(";") if x and x != "."],
                info=info,
                samples={sample: cell},
                lineno=lineno,
            )
