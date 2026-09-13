"""A deliberately paranoid, dependency-free VCF reader.

Why not cyvcf2/pysam?  Because this tool's whole contract is "never raise on
malformed input", and the fast C readers do the opposite: they abort (or worse)
on exactly the broken files we exist to describe.  A tool that dies on a
truncated VCF cannot report that the VCF is truncated.  VCF is line-oriented,
so a tolerant pure-Python parser costs us speed and buys us the contract plus a
zero-dependency ``pip install``.  A cyvcf2 fast path can be added later behind
an extra, once the semantics are pinned down by tests.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import re
from array import array
from dataclasses import dataclass, field

from .contigs import infer_build, is_primary, normalize_contig

MISSING = 3  # genotype code for no-call / unassessed
HOMREF, HET, HOMALT = 0, 1, 2

SiteKey = tuple[str, int, str, str]

_CONTIG_RE = re.compile(r"ID=([^,>]+).*?length=(\d+)", re.IGNORECASE)
_ID_RE = re.compile(r"ID=([^,>]+)")

# Gene-symbol keys written by the common annotators, in preference order.
GENE_INFO_KEYS = ("Gene.refGene", "Gene.refGeneWithVer", "Gene.knownGene", "GENE", "Gene")

# Sub-field names that carry a gene symbol inside a VEP ``CSQ`` or SnpEff ``ANN``
# block, in preference order: a symbol groups compound heterozygotes the way a
# clinician reads them, an Ensembl id does so correctly but unreadably.
_CSQ_GENE_FIELDS = ("SYMBOL", "Gene_Name", "Gene", "Gene_ID", "GeneID")
_CSQ_CONSEQUENCE_FIELDS = ("Consequence", "Annotation")
_CSQ_FORMAT_RE = re.compile(
    r"(?:Format|Functional annotations):\s*'?\s*(.+?)\s*'?\s*(?:\"|>)", re.IGNORECASE
)

# Sequence Ontology terms mapped onto the coarse region classes check 3 reasons
# about.  Only the two classes that decide "is this a whole exome or a
# coding-only extract" are mapped; everything else is counted under its own term
# and contributes to the denominator, so an unmapped term can never be silently
# read as coding.
_SO_CODING = frozenset({
    "missense_variant", "synonymous_variant", "stop_gained", "stop_lost",
    "start_lost", "frameshift_variant", "inframe_insertion", "inframe_deletion",
    "protein_altering_variant", "coding_sequence_variant", "incomplete_terminal_codon_variant",
    "stop_retained_variant", "start_retained_variant",
    "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant",
    "splice_donor_5th_base_variant", "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
})
_SO_INTRONIC = frozenset({"intron_variant", "non_coding_transcript_variant"})


def csq_field_index(header_line: str, wanted: tuple[str, ...]) -> tuple[str, int] | None:
    """Locate one named sub-field inside a VEP ``CSQ`` / SnpEff ``ANN`` block."""
    m = _ID_RE.search(header_line)
    if not m or m.group(1) not in ("CSQ", "ANN"):
        return None
    fmt = _CSQ_FORMAT_RE.search(header_line)
    if not fmt:
        return None
    cols = [c.strip() for c in fmt.group(1).split("|")]
    for want in wanted:
        for i, c in enumerate(cols):
            if c.lower() == want.lower():
                return m.group(1), i
    return None


def csq_gene_index(header_line: str) -> tuple[str, int] | None:
    """Locate the gene column inside a VEP/SnpEff annotation block.

    Both annotators declare their pipe-delimited layout in the INFO header line
    rather than fixing it, so the position has to be read from the file.  Only
    ANNOVAR-style flat keys were understood before, which left compound-het
    unevaluable on VEP- or SnpEff-annotated VCFs - that is, on most of them.
    """
    return csq_field_index(header_line, _CSQ_GENE_FIELDS)


def encode_gt(gt: str) -> int:
    """Map a GT string to 0 hom-ref / 1 het / 2 hom-alt / 3 missing.

    ``1/2`` encodes as het: for kinship arithmetic what matters is that the two
    alleles differ, not that neither is the reference.  Haploid calls (``1`` on
    male X/Y from a ploidy-aware caller) encode as hom-alt, which keeps them out
    of the heterozygous counts where they belong.
    """
    if not gt:
        return MISSING
    sep = "/" if "/" in gt else ("|" if "|" in gt else None)
    parts = gt.split(sep) if sep else [gt]
    vals: list[int] = []
    for p in parts:
        p = p.strip()
        if not p or not p.isdigit():
            return MISSING
        vals.append(int(p))
    if not vals:
        return MISSING
    if len(set(vals)) > 1:
        return HET
    return HOMREF if vals[0] == 0 else HOMALT


@dataclass
class VcfHeader:
    samples: list[str] = field(default_factory=list)
    contig_lengths: dict[str, int] = field(default_factory=dict)
    info_keys: set[str] = field(default_factory=set)
    format_keys: set[str] = field(default_factory=set)
    filter_ids: set[str] = field(default_factory=set)
    reference: str | None = None
    n_header_lines: int = 0
    # (INFO key, index of the gene sub-field) for a VEP CSQ or SnpEff ANN block.
    csq_gene: tuple[str, int] | None = None
    # The same, for the consequence sub-field, which is what check 3 needs to
    # answer "is this a whole exome or a coding-only extract".
    csq_consequence: tuple[str, int] | None = None


@dataclass
class VcfStats:
    """Everything checks 2-4 will want, collected in the single read pass."""

    n_records: int = 0
    n_malformed: int = 0
    n_reference_blocks: int = 0  # gVCF <NON_REF>-only rows
    n_pass: int = 0
    n_nonpass: int = 0
    n_filter_absent: int = 0
    n_multiallelic: int = 0
    n_nonprimary_contig: int = 0
    contig_counts: dict[str, int] = field(default_factory=dict)
    # Region classes taken from whatever annotation is already in the user's own
    # file (ANNOVAR Func.refGene today).  Nothing is fetched or bundled.
    region_counts: dict[str, int] = field(default_factory=dict)
    annotation_source: str | None = None
    info_keys_seen: set[str] = field(default_factory=set)
    format_keys_seen: set[str] = field(default_factory=set)
    has_format_column: bool = False


@dataclass
class VcfScan:
    path: str
    header: VcfHeader
    stats: VcfStats
    build: str | None
    file_size: int
    file_sha256: str
    problems: list[str] = field(default_factory=list)
    # per-sample, keyed by site id assigned by the SiteIndex that drove the scan
    genotypes: dict[str, dict[int, int]] = field(default_factory=dict)
    has_explicit_homref: dict[str, bool] = field(default_factory=dict)
    # Per-site annotation taken from the user's OWN file. Nothing is fetched and
    # no database is bundled; if the annotation is absent the dependent checks
    # report UNKNOWN rather than guessing.
    site_gene: dict[int, str] = field(default_factory=dict)
    site_af: dict[int, float] = field(default_factory=dict)

    @property
    def is_multisample(self) -> bool:
        return len(self.header.samples) > 1


class SiteIndex:
    """Interns (chrom, pos, ref, alt) so genotype vectors can be plain arrays."""

    def __init__(self) -> None:
        self._ids: dict[SiteKey, int] = {}
        self.keys: list[SiteKey] = []

    def __len__(self) -> int:
        return len(self.keys)

    def id_for(self, key: SiteKey) -> int:
        got = self._ids.get(key)
        if got is None:
            got = len(self.keys)
            self._ids[key] = got
            self.keys.append(key)
        return got

    def get(self, key: SiteKey) -> int | None:
        return self._ids.get(key)


def _sha256(path: str | os.PathLike[str]) -> tuple[int, str]:
    h = hashlib.sha256()
    size = 0
    try:
        with open(path, "rb") as fh:
            while chunk := fh.read(1 << 20):
                size += len(chunk)
                h.update(chunk)
    except OSError:
        return 0, ""
    return size, h.hexdigest()


def _open_text(path: str | os.PathLike[str]):
    if str(path).endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", errors="replace")
    return open(path, errors="replace")


def scan_vcf(
    path: str | os.PathLike[str], index: SiteIndex, af_key: str | None = None
) -> VcfScan:
    """Read one VCF into ``index``.  Never raises; problems land in ``problems``."""
    header = VcfHeader()
    stats = VcfStats()
    size, digest = _sha256(path)
    scan = VcfScan(
        path=str(path), header=header, stats=stats, build=None,
        file_size=size, file_sha256=digest,
    )

    try:
        fh = _open_text(path)
    except OSError as exc:
        scan.problems.append(f"cannot open: {exc}")
        return scan

    gt_idx_cache: dict[str, int | None] = {}
    try:
        with fh:
            for lineno, raw in enumerate(fh, 1):
                line = raw.rstrip("\n")
                if not line:
                    continue
                if line.startswith("##"):
                    header.n_header_lines += 1
                    _parse_meta(line, header)
                    continue
                if line.startswith("#CHROM"):
                    header.n_header_lines += 1
                    cols = line.split("\t")
                    if len(cols) > 9:
                        header.samples = cols[9:]
                    stats.has_format_column = len(cols) > 8
                    for s in header.samples:
                        scan.genotypes.setdefault(s, {})
                        scan.has_explicit_homref.setdefault(s, False)
                    continue

                fields = line.split("\t")
                if len(fields) < 8:
                    stats.n_malformed += 1
                    if len(scan.problems) < 20:
                        scan.problems.append(f"line {lineno}: {len(fields)} columns, expected >=8")
                    continue
                ok = _consume_record(fields, header, stats, scan, index, gt_idx_cache, af_key)
                if not ok:
                    stats.n_malformed += 1
                    if len(scan.problems) < 20:
                        scan.problems.append(f"line {lineno}: unparseable record")
    except (OSError, EOFError, gzip.BadGzipFile) as exc:
        # A truncated bgzip file lands here.  That is a finding, not a crash.
        scan.problems.append(f"read stopped early: {type(exc).__name__}: {exc}")

    scan.build = infer_build(header.contig_lengths)
    return scan


def _parse_meta(line: str, header: VcfHeader) -> None:
    if line.startswith("##contig="):
        m = _CONTIG_RE.search(line)
        if m:
            try:
                header.contig_lengths[m.group(1)] = int(m.group(2))
            except ValueError:
                pass
    elif line.startswith("##INFO="):
        if m := _ID_RE.search(line):
            header.info_keys.add(m.group(1))
        if (found := csq_gene_index(line)) is not None:
            header.csq_gene = found
        if (found := csq_field_index(line, _CSQ_CONSEQUENCE_FIELDS)) is not None:
            header.csq_consequence = found
    elif line.startswith("##FORMAT="):
        if m := _ID_RE.search(line):
            header.format_keys.add(m.group(1))
    elif line.startswith("##FILTER="):
        if m := _ID_RE.search(line):
            header.filter_ids.add(m.group(1))
    elif line.startswith("##reference="):
        header.reference = line.split("=", 1)[1].strip()


def _consume_record(
    fields: list[str],
    header: VcfHeader,
    stats: VcfStats,
    scan: VcfScan,
    index: SiteIndex,
    gt_idx_cache: dict[str, int | None],
    af_key: str | None = None,
) -> bool:
    chrom_raw, pos_raw, _id, ref, alt, _qual, filt = fields[:7]
    try:
        pos = int(pos_raw)
    except ValueError:
        return False

    alts = alt.split(",")
    # gVCF reference block: nothing to genotype, but worth counting.
    if alts == ["<NON_REF>"] or alt == ".":
        stats.n_reference_blocks += 1
        return True

    stats.n_records += 1
    chrom = normalize_contig(chrom_raw)
    stats.contig_counts[chrom] = stats.contig_counts.get(chrom, 0) + 1
    if not is_primary(chrom_raw):
        stats.n_nonprimary_contig += 1
    if len([a for a in alts if a != "<NON_REF>"]) > 1:
        stats.n_multiallelic += 1

    f = filt.strip()
    if f in ("", "."):
        stats.n_filter_absent += 1
    elif f == "PASS":
        stats.n_pass += 1
    else:
        stats.n_nonpass += 1

    gene_value: str | None = None
    af_value: float | None = None
    for kv in fields[7].split(";"):
        k, _, v = kv.partition("=")
        stats.info_keys_seen.add(k)
        if gene_value is None and k in GENE_INFO_KEYS and v and v != ".":
            gene_value = v.split(",")[0].split("\\x3b")[0].strip()
        if (
            gene_value is None
            and header.csq_gene is not None
            and k == header.csq_gene[0]
            and v
        ):
            # A CSQ/ANN value is one block per transcript, comma-separated; the
            # blocks at one site usually name the same gene, and where they do
            # not, the first is the annotator's own primary choice.  Grouping
            # compound heterozygotes needs one gene per site, so take it.
            idx = header.csq_gene[1]
            for block in v.split(","):
                parts = block.split("|")
                if idx < len(parts) and parts[idx].strip():
                    gene_value = parts[idx].strip()
                    break
        if af_key and k == af_key and v:
            try:
                af_value = float(v.split(",")[0])
            except ValueError:
                af_value = None
        if k in ("Func.refGene", "Func.refgene", "Func_refGene") and v:
            stats.annotation_source = "ANNOVAR:Func.refGene"
            for cls in v.split("\\x3b"):  # ANNOVAR escapes ';' inside INFO values
                cls = cls.strip()
                if cls:
                    stats.region_counts[cls] = stats.region_counts.get(cls, 0) + 1
        elif (
            header.csq_consequence is not None
            and k == header.csq_consequence[0]
            and v
            # ANNOVAR's own class, where present, wins: mixing two annotators'
            # answers into one tally would make the fractions meaningless.
            and not (stats.annotation_source or "").startswith("ANNOVAR")
        ):
            # VEP and SnpEff answer the same question as Func.refGene, in
            # Sequence Ontology terms and one block per transcript.  A site is
            # classed by the most consequential term any transcript gives it,
            # which is how a clinician reads it and how ANNOVAR's single value
            # behaves - taking the first transcript instead would class a coding
            # variant as intronic whenever VEP happened to list an intron-
            # containing transcript first.
            stats.annotation_source = f"VEP/SnpEff:{header.csq_consequence[0]}"
            idx = header.csq_consequence[1]
            terms: set[str] = set()
            for block in v.split(","):
                parts = block.split("|")
                if idx < len(parts) and parts[idx].strip():
                    terms.update(t.strip() for t in parts[idx].split("&") if t.strip())
            if terms & _SO_CODING:
                cls = "exonic"
            elif terms & _SO_INTRONIC:
                cls = "intronic"
            else:
                cls = "other"
            stats.region_counts[cls] = stats.region_counts.get(cls, 0) + 1

    if len(fields) < 10 or not header.samples:
        return True  # sites-only VCF; check 3 reports FORMAT stripped

    fmt = fields[8]
    if fmt not in gt_idx_cache:
        parts = fmt.split(":")
        for k in parts:
            stats.format_keys_seen.add(k)
        gt_idx_cache[fmt] = parts.index("GT") if "GT" in parts else None
    gt_idx = gt_idx_cache[fmt]
    if gt_idx is None:
        return True

    # ``<NON_REF>`` is gVCF bookkeeping - the symbolic "any other allele" - not an
    # observed alternate.  Leaving it in the site key gives the same variant two
    # different identities depending on whether it arrived in a gVCF or a plain
    # VCF, so a gVCF and a VCF of the SAME person share no sites at all and the
    # pair is reported unrelated.  Failing towards "unrelated" is the one
    # direction an identity check must never fail in.
    real_alts = [a for a in alts if a.upper() != "<NON_REF>"] or alts
    key: SiteKey = (chrom, pos, ref.upper(), ",".join(real_alts).upper())
    site_id = index.id_for(key)
    if gene_value:
        scan.site_gene[site_id] = gene_value
    if af_value is not None:
        scan.site_af[site_id] = af_value
    for sample, cell in zip(header.samples, fields[9:], strict=False):
        sub = cell.split(":")
        if gt_idx >= len(sub):
            continue
        code = encode_gt(sub[gt_idx])
        if code == MISSING:
            continue
        if code == HOMREF:
            scan.has_explicit_homref[sample] = True
        scan.genotypes[sample][site_id] = code
    return True


@dataclass
class GenotypeMatrix:
    """Cohort genotypes as fixed-length byte vectors over a shared site index."""

    index: SiteIndex
    samples: list[str] = field(default_factory=list)
    codes: dict[str, array] = field(default_factory=dict)
    source: dict[str, str] = field(default_factory=dict)
    dense: dict[str, bool] = field(default_factory=dict)
    scans: dict[str, VcfScan] = field(default_factory=dict)
    # Every file that was read, in the order given, whether or not it contained
    # any samples.  ``scans`` is keyed by sample and therefore cannot represent a
    # sites-only VCF at all - and a sites-only VCF is exactly what check 3 exists
    # to name, so it must not disappear on the way there.
    files: list[VcfScan] = field(default_factory=list)
    site_gene: dict[int, str] = field(default_factory=dict)
    site_af: dict[int, float] = field(default_factory=dict)

    @property
    def n_sites(self) -> int:
        return len(self.index)

    def pair_is_dense(self, a: str, b: str) -> bool:
        """Only true when both samples were jointly assessed at the same sites.

        In separate single-sample VCFs a site absent from a sample's file is
        indistinguishable between "hom-ref" and "never callable there".  That
        distinction is exactly what IBS0 and KING-robust depend on, so those
        estimators are only valid within one multi-sample (or joint) VCF.
        """
        return (
            self.source.get(a) == self.source.get(b)
            and self.dense.get(a, False)
            and self.dense.get(b, False)
        )


def unique_label(sample: str, path: str | os.PathLike[str], taken: set[str]) -> str:
    """A stable, human-readable label when several files share a sample name.

    Pipelines that write every sample to ``Patient_GATK.vcf.gz`` with the
    internal sample name ``patient`` are common, and basenames collide there
    too.  Walk outwards through the path until something distinguishes them,
    so the label reads like ``patient@5`` rather than an opaque index.
    """
    p = os.path.normpath(str(path))
    parent = os.path.basename(os.path.dirname(p))
    grandparent = os.path.basename(os.path.dirname(os.path.dirname(p)))
    for cand in (
        sample,
        f"{sample}@{parent}",
        f"{sample}@{grandparent}/{parent}",
        f"{sample}@{grandparent}/{parent}/{os.path.basename(p)}",
        f"{sample}@{p}",
    ):
        if cand and cand not in taken:
            return cand
    i = 2
    while f"{sample}#{i}" in taken:
        i += 1
    return f"{sample}#{i}"


def load_cohort(
    paths: list[str | os.PathLike[str]], af_key: str | None = None
) -> GenotypeMatrix:
    """Scan every VCF into one matrix.  Colliding sample names are disambiguated."""
    index = SiteIndex()
    raw: dict[str, dict[int, int]] = {}
    matrix = GenotypeMatrix(index=index)

    from .annovar import looks_like_multianno, scan_multianno  # avoids a cycle

    for path in paths:
        scan = (
            scan_multianno(path, index, af_key=af_key)
            if looks_like_multianno(path)
            else scan_vcf(path, index, af_key=af_key)
        )
        matrix.files.append(scan)
        matrix.site_gene.update(scan.site_gene)
        matrix.site_af.update(scan.site_af)
        for sample in scan.header.samples:
            name = unique_label(sample, path, set(raw))
            if name != sample:
                scan.problems.append(f"sample name {sample!r} already seen; stored as {name!r}")
            raw[name] = scan.genotypes.get(sample, {})
            matrix.samples.append(name)
            matrix.source[name] = str(path)
            matrix.dense[name] = scan.has_explicit_homref.get(sample, False)
            matrix.scans[name] = scan

    n = len(index)
    for name, calls in raw.items():
        vec = array("b", bytes([MISSING])) * n if n else array("b")
        for site_id, code in calls.items():
            vec[site_id] = code
        matrix.codes[name] = vec
    return matrix


# --------------------------------------------------------------------------
# Streaming record access, for checks that need per-genotype FORMAT fields.
#
# The cohort matrix deliberately keeps only genotype codes, because identity
# work needs nothing else and holding full records for a cohort is wasteful.
# Genotype confidence needs AD, DP, GQ and PL, so it streams the files again
# rather than inflating the matrix for every caller.
# --------------------------------------------------------------------------


@dataclass
class Record:
    chrom: str
    pos: int
    ref: str
    alts: list[str]
    qual: float | None
    filters: list[str]
    info: dict[str, str]
    samples: dict[str, dict[str, str]]
    lineno: int

    @property
    def is_biallelic(self) -> bool:
        return len([a for a in self.alts if a != "<NON_REF>"]) == 1


def iter_records(path: str | os.PathLike[str]):
    """Yield parsed records.  Malformed lines are skipped, never raised."""
    samples: list[str] = []
    try:
        fh = _open_text(path)
    except OSError:
        return
    try:
        with fh:
            for lineno, raw in enumerate(fh, 1):
                line = raw.rstrip("\n")
                if not line:
                    continue
                if line.startswith("#CHROM"):
                    cols = line.split("\t")
                    samples = cols[9:] if len(cols) > 9 else []
                    continue
                if line.startswith("#"):
                    continue
                f = line.split("\t")
                if len(f) < 8:
                    continue
                try:
                    pos = int(f[1])
                except ValueError:
                    continue
                alts = f[4].split(",")
                if alts == ["<NON_REF>"] or f[4] == ".":
                    continue
                try:
                    qual = float(f[5])
                except ValueError:
                    qual = None
                info: dict[str, str] = {}
                for kv in f[7].split(";"):
                    k, _, v = kv.partition("=")
                    info[k] = v
                per_sample: dict[str, dict[str, str]] = {}
                if len(f) > 9 and samples:
                    keys = f[8].split(":")
                    for name, cell in zip(samples, f[9:], strict=False):
                        per_sample[name] = dict(zip(keys, cell.split(":"), strict=False))
                yield Record(
                    chrom=normalize_contig(f[0]),
                    pos=pos,
                    ref=f[3].upper(),
                    alts=[a.upper() for a in alts],
                    qual=qual,
                    # "." means "no filtering applied", which is not the same as
                    # "failed a filter" - conflating them hides an unfiltered
                    # call set behind a wall of false non-PASS records.
                    filters=[x for x in f[6].split(";") if x and x != "."],
                    info=info,
                    samples=per_sample,
                    lineno=lineno,
                )
    except (OSError, EOFError, gzip.BadGzipFile):
        return
