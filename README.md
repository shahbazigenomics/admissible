# admissible

**Audit whether exome data can support an interpretation — before you interpret it.**

`admissible` does not tell you whether a variant is pathogenic. It tells you whether
your evidence is admissible: whether the samples are who the pedigree says they are,
whether the file you were given is the file you think you were given, and whether you
actually searched enough of the exome to be allowed to report a negative.

```
$ admissible audit family/*.vcf --ped family.ped

COHORT FAM_B                                        VERDICT: NOT INTERPRETABLE

Identity .......... FAIL     2 duplicate pair(s), 2 sex mismatch(es)
Provenance ........ WARN     CODING_ONLY + PASS_FILTERED + SUBSET
Callable .......... UNKNOWN  no coverage supplied; a VCF cannot answer this
Genotype QC ....... WARN     318 of 2104 homozygous calls unsupported (15%)
Models ............ WARN     0 candidates under 4 models; not rarity-filtered

DO NOT CONCLUDE: "no monogenic cause"

NEXT STEPS: 1. resolve sample identity (blocking)
            2. re-call without the coding-only interval file
```

---

## Why this exists

A candidate variant held up for a year turned out to be a false homozygous call
supported by two reads. Chasing that led to a second discovery: two exomes filed
under one family belonged to different people. Chasing *that* led to a third: the
VCFs had been pre-filtered to coding-only and PASS-only before delivery, so only a
fraction of the exome was ever searchable in the first place.

None of the tools in that pipeline — a commercial interpretation platform, a
standard annotator, a by-the-book GATK workflow — raised any of it. They all answer
the same question: *is this variant pathogenic?* None of them answers the question
that had to be answered first: *can I trust this data, and did I actually search
enough?*

That gap is this tool.

## What it does

**Input:** a family's VCFs + a PED file (+ optionally mosdepth coverage summaries)
**Output:** a one-page verdict, in human-readable text and machine-readable JSON

It never interprets a variant. It audits the evidence base that an interpretation
would rest on, across five checks:

| # | Check | Question | 0.1.0 |
|---|-------|----------|-------|
| 1 | Identity | Is each sample who the pedigree says it is? | implemented |
| 2 | Genotype | Which genotypes are not supported by their own evidence? | implemented |
| 3 | Provenance | Is this file actually a whole exome? | implemented |
| 4 | Callability | How much of the target did you really search? | implemented |
| 5 | Models | What does the inheritance-model count profile look like? | implemented |

`UNKNOWN` in a report therefore always means *your inputs cannot answer this*, never
*this is not written yet*. Callability reports `UNKNOWN` given only VCFs, because a
VCF genuinely cannot answer it; the two-locus model reports `not-applicable` because
the pairwise search needs an explicit multiple-testing treatment first.

## Install

```bash
git clone https://github.com/shahbazigenomics/admissible
cd admissible && pip install -e ".[dev]"
```

Not on PyPI yet, so there is deliberately no `pip install admissible` line here
to fail on.

Zero required runtime dependencies, Python 3.10+.

## Design decisions worth arguing with

These are the choices most likely to be wrong. They are stated here so they can be
challenged rather than discovered.

**Duplicate detection is cohort-wide, never within-family.** A within-family scan
cannot find a swap that lives *between* two families, which is the failure mode that
motivated the tool.

**Two relatedness engines, and the input decides which applies.** Jaccard and
genotype agreement over non-reference sites always work, but they are depth- and
pipeline-dependent — so the boundary is calibrated on your own cohort rather than
hard-coded, using medians and MAD so that a minority of mislabelled pairs cannot
drag it. They can separate related from unrelated; they cannot resolve degree.

Where samples were *jointly* genotyped, the tool switches to **KING-robust kinship
and IBS0**. Both are allele-frequency-free, so they need no external database, and
both are far more robust to depth. IBS0 also separates parent–offspring from full
sibs, a distinction site-set overlap cannot make at all — though see the public
validation below for how far that holds. On the bundled joint fixture Jaccard puts
parent–offspring at 0.570 and full sibs at 0.592 while IBS0 puts them at 0.000 and
0.025; on real data, where genotype error puts a floor under parent–offspring IBS0,
the separation is clear on average but the two distributions overlap at the edges.

The catch is that this only works within one multi-sample or joint-called VCF. In
separate single-sample VCFs, "hom-reference" and "never callable here" are the same
absence of a line, and IBS0 depends on telling them apart. So with single-sample
VCFs the tool reports degree as unresolved instead of guessing.

**Sex is called from chrX heterozygosity with PAR and XTR excluded, reported as an
interval.** chrY call counts are shown as context and take no part in the call:
females routinely carry a handful of chrY calls from X–Y homologous mismapping, so
single-digit counts carry no evidential weight. Consanguinity and long runs of
homozygosity depress female chrX heterozygosity, so the boundary is calibrated on
the cohort's own bimodality where one exists, and borderline samples are returned
as *not determined* rather than forced into a call.

**Duplicate pairs are reported by evidence, not by cause.** The tool states what is
observable — byte identity, site-set identity, call-count ratio, shared-site
agreement — and lists the mechanisms compatible with it. It will not assert that a
duplicate came from the wet lab rather than from file handling: re-running one FASTQ
through a different pipeline version produces the same observation as preparing a
second library.

**No bundled annotation database.** The tool reads only your own files. Region
classes come from annotation already present in your VCF; nothing is fetched and
nothing is shipped. Where a check genuinely needs a population allele frequency, it
takes a user-supplied source and reports `UNKNOWN` without one, rather than quietly
substituting a worse definition of "rare".

**A count is never reported without its denominator.** Zero candidates over 20% of
the target and zero over 95% are not the same finding, and printing them identically
is how "no monogenic cause" gets written down. Every model in the sweep carries its
own callable fraction, and a model the inputs cannot support reports
`not-applicable` rather than zero — because zero looks like evidence of absence.

**Nothing raises.** A tool whose purpose is to describe broken files must not die on
one. Truncated downloads, double-gzipped files, sites-only VCFs, headerless VCFs and
files that are not VCFs at all all produce a typed failure with a stated reason. The
test suite fuzzes for this.

## Giving it the pedigree

A PED file works, but PED is positional and numeric — `2` means affected, `1`
means unaffected — and getting that backwards produces no error at all: check 5
simply finds nothing to segregate. Two easier routes:

**Start from your VCFs.** The ids in a pedigree have to match the sample names
*inside* the VCFs, which nobody can see without looking, and a mismatch is the
commonest way a first run goes quiet:

```bash
admissible ped-template family/*.vcf.gz > family.ped        # names already filled in
admissible ped-template family/*.vcf.gz --csv > family.csv  # to fill in a spreadsheet
```

**Or hand it the spreadsheet you already keep.** Columns are found by name and
values read in words; `--ped` takes it directly:

```csv
Family,Sample ID,Father,Mother,Sex,Affected
IPC,IPC-1,,,Male,no
IPC,IPC-2,,,Female,no
IPC,IPC-3,IPC-1,IPC-2,M,yes
IPC,IPC-4,IPC-1,IPC-2,F,yes
```

`sample`/`id`/`individual`, `father`/`pat`/`dad`, `mother`/`mat`/`mum`,
`sex`/`gender`, `affected`/`status`/`phenotype`, `family`/`fid` are all
recognised; so are `M`/`F`, `male`/`female`, `1`/`2`, `yes`/`no`,
`case`/`control`. CSV, TSV and (with openpyxl) xlsx. The format is decided by
what is in the file, not by its extension.

A value that is *not* recognised becomes **unknown** and is reported — never
rounded to the plausible one, because rounding is how `1`/`2` goes wrong
silently in the first place. A blank parent means "not in the study"; a parent
named but given no row of their own is flagged, since that is a real and easy
mistake.

Leave `--ped` off entirely and `admissible identity` still does sex inference
and cohort-wide duplicate detection — which is the right first move on a cohort
you have not pedigreed yet.

## Running the checks separately

```bash
admissible identity   family/*.vcf --ped family.ped -v
admissible provenance family/*.vcf
admissible audit      family/*.vcf --ped family.ped --json report.json
```

Exit codes: `0` nothing blocking, `1` a check failed or raised a blocking finding,
`2` the inputs could not answer the question.

## Test fixtures

All fixtures are **synthetic**. No patient data is in this repository, and none ever
will be. The generators rebuild them deterministically:

```bash
python tests/fixtures/make_fixtures.py       # 13 samples, 3 families, two planted swaps
python tests/fixtures/make_joint_family.py   # one jointly called family, real transmission
pytest
```

The cohort fixture plants, and the test suite asserts recovery of: a pair of samples
in different families that are the same individual (Jaccard 1.0000, agreement 1.0000,
files of identical byte length differing at exactly one byte); a second cross-family
pair that is one individual sequenced twice at different depth (Jaccard 0.7508,
agreement 0.9726); two male samples sitting in pedigree slots declared female; and a
coding-only, PASS-only, subset call set. Declared-unrelated pairs sit at Jaccard
0.34–0.37 and true first-degree pairs at 0.48–0.52 — a deliberately thin margin,
because that thinness is the argument for calibrating the boundary instead of
hard-coding it.

## Validation on public data

Synthetic fixtures written by the author of the code are weak evidence, so
checks 1, 2 and 5 are also run against **CEPH pedigree 1463** — the Utah
three-generation family, as distributed with
[peddy](https://github.com/brentp/peddy): a real joint-called VCF (17 samples,
freebayes, `GL` rather than `PL`, no `##contig` headers), a real published
pedigree, and a deliberately corrupted copy of that pedigree that peddy ships as
its own canonical failure case.

| property | result |
|---|---|
| sex calls vs the published pedigree | 16 called, **0 wrong**; 1 returned *undetermined* at the boundary rather than guessed |
| truly unrelated pairs called related | **0 of 11** |
| first-degree pairs recognised as related | **all**, minimum φ 0.153 |
| false duplicate pairs across 136 comparisons | **0** |
| the correct pedigree | accepted, no findings |
| peddy's corrupted pedigree (father/daughter swapped) | caught, both swapped samples named |

Two things this exposed that the synthetic fixtures could not, both now
reflected in the code:

**Band-label comparison is brittle.** True full sibs in this family come out at
φ 0.153 and 0.176, just under the 0.177 first-degree boundary — realised IBD
genuinely varies between sibs. An earlier version compared band *labels* and
reported both as pedigree errors. Reconciliation now compares the observed
kinship to the kinship the pedigree *implies*, computed recursively, and flags
only discrepancies too large to be noise. It deliberately says nothing about
first-vs-second degree.

**Second-degree relatedness is not reliable per pair at this site count.**
Grandparent–grandchild pairs average the textbook 0.125, but the weakest falls
below the unrelated band floor. The tool reports the informative-site count and
says so, rather than implying a precision it does not have.

Running checks 2 and 5 on the same file found two defects that the synthetic
fixtures could not, because a fixture written alongside the code inherits the
code's assumptions about what a VCF looks like:

**Allele depth is not always spelled `AD`.** freebayes writes `RO`/`AO`, older
samtools pipelines write `DP4`. Reading only `AD` left the entire allele-balance
arm of check 2 silently inert on freebayes output — no `AB_SKEW`, no
`ALLELE_IMBALANCE_HOM` — while the report still looked complete. On CEPH those
two now fire 2,598 and 63 times respectively, and allele balance finds 20
false-homozygote candidates that the likelihood arm alone does not.

**Gene assignment is not always an ANNOVAR key.** VEP writes `CSQ` and SnpEff
writes `ANN`, both pipe-delimited with the layout declared in the header rather
than fixed. Reading only flat keys made compound-heterozygous report
*not-applicable* on most annotated VCFs, CEPH included; it now computes.

**A file with no samples is still a file.** The cohort loader keyed every scan
by sample name, so a sites-only VCF — no FORMAT column, no genotypes — vanished
before reaching check 3, and the CLI answered *"no VCFs were supplied"* about a
file it had just read. A sites-only export is the clearest possible case of
`METRICS_STRIPPED`, so the one case the verdict exists for was the one case it
could never report.

A third thing this exposed is not a bug but a reporting duty. Three declared-
affected siblings in CEPH yield 29 apparent de novo variants over ~20,000 sites.
At a germline rate near 1.3 × 10⁻⁸ per base per generation a whole exome expects
well under one true de novo per proband, so essentially all 29 are genotyping
error. The count is what the segregation filter finds and is not wrong — but
printed bare it reads as a mutation count, so it now travels with that caveat.

Check 3 is validated differently, because provenance has no ground truth in a
file that was never tampered with: one known transformation at a time is applied
to the real CEPH file — a coding-only extract, a PASS-only delivery, a merge that
rewrites `0/0` as `./.`, a sites-only export — and the tool must name that
transformation and no other. It does, and the unmodified file is called a subset
and nothing else.

Reproduce with `git clone --depth 1 https://github.com/brentp/peddy /tmp/peddy`
then `pytest tests/test_public_ceph.py` (the tests skip if the data is absent).

### Real capture coverage

Check 4's central claim — that the joint callable fraction is an intersection and
never a product — is measured rather than asserted, on four unrelated 1000
Genomes exomes (Baylor, SRP000808) over chr20:1,400,000–1,500,000, against a
denominator taken from an Ensembl GTF rather than from the coverage itself:

| | fraction of target |
|---|---|
| HG00349 / HG00350 / HG00351 / HG00358 | 0.7383 / 0.7951 / 0.4768 / 0.7975 |
| joint, all four callable at once | **0.3815** |
| product of the per-sample fractions | **0.2232** |
| three of four (phenocopy model) | 0.6870 |

The product is 41% low, because the four lose coverage in the same places. At
0.38 of the target searched, the strongest permitted claim is *exploratory; a
causal variant cannot be excluded* — which is what the tool reports.

4,134 bp over four samples is a small measurement and carries real sampling
noise. It demonstrates the effect on real capture; it does not estimate its size
in general.

### The ANNOVAR reader

A delivered analysis is often nothing but annotated tables, so the adapter that
reads them matters. Run against real GATK genotypes wrapped in ANNOVAR's
documented `--vcfinput` layout, it reconstructed 113 of 120 genotypes correctly
and got 7 wrong — **every multiallelic `1/2` call became a homozygous
alternate**, because `AC=1,1` sums to two. A heterozygote read as a homozygote
is the one failure this tool exists to prevent, and the reader was committing it
before any check got to look.

The cause was the same one that ran through checks 2, 3 and 5: reading a derived
proxy while the real field sat unread. `--vcfinput` keeps the original FORMAT and
sample columns; the genotype is now read from them, which is 120 of 120, and
`DP`, `GQ`, `PL` and `AD` come with it — so check 2 can examine an
ANNOVAR-delivered analysis using likelihoods and allele balance instead of
falling back to `MLEAC`. Tables produced without `--vcfinput` have no FORMAT
block at all, and there the `AC`/`AN` reconstruction is still the only option.

## What this tool cannot see

Callability here is depth-based, which makes it blind to structural variation. A
heterozygous deletion leaves depth comfortably above any sensible floor, so the
region reads as callable; the surviving allele is then called homozygous with a
genuinely clean allele balance, so check 2 does not flag it either. Both checks
pass and the genotype is still wrong. Where a CNV is a plausible mechanism —
`IL10RA`/`IL10RB` in early-onset IBD, for one — a read-depth or split-read CNV
caller answers a question this tool does not ask.

The target BED is also taken on trust. Supply the wrong capture kit and every
fraction in check 4 is wrong; only a gross mismatch trips `NOT_INTERVAL_RESTRICTED`.

## Related tools

`admissible` is not a replacement for any of these, and where they apply you should
prefer them:

- [somalier](https://github.com/brentp/somalier) — identity, relatedness, sex
- [peddy](https://github.com/brentp/peddy) — the same, VCF-native
- [slivar](https://github.com/brentp/slivar) — inheritance models, in-house AF
- [NGSCheckMate](https://github.com/parklab/NGSCheckMate) — identity at the FASTQ stage
- [mosdepth](https://github.com/brentp/mosdepth) — coverage summaries

One caveat drove a design decision here. somalier and peddy score against a fixed
panel of common genome-wide sites, largely intronic and intergenic. A coding-only,
PASS-only call set overlaps almost none of them — so on exactly the data class this
tool exists to detect, they have little to work with. The native genotype engine is
therefore the primary path here, not a fallback, and adapters for these tools will
report their overlapping-site count and decline to run below a floor.

## Licence

MIT. See [LICENSE](LICENSE).
