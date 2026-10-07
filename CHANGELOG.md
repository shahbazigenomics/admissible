# Changelog

All notable changes to this project are documented here.
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.1] — unreleased

Fixes found after 0.1.0 by an AI adversarial review and by running the tool on a
real 13-sample, 3-family familial-IBD WES cohort. **Everyone installing from PyPI
before this release has 0.1.0, which lacks all of them.** (Set the release date
here when this is tagged.)

### Fixed

- **A real cross-format duplicate produced no finding.** The duplicate gate needed
  `agreement >= 0.95` *and* `jaccard >= 0.60`; Jaccard collapses when one file has
  far fewer calls than the other (a native VCF versus a table rebuilt from AC/AN).
  A real pair with agreement 1.000 and Jaccard 0.181 passed silently. New WARN
  `HIGH_AGREEMENT_LOW_JACCARD`, shown when agreement clears the bar, Jaccard does
  not, and the non-ref counts differ by at least `dup_completeness_ratio` (default
  2.0). It does not auto-classify a duplicate: a human verifies.
- **Pairwise evidence was only in the JSON.** `-v` now prints a `PAIRWISE EVIDENCE`
  table in the text report for pairs near the duplicate or relatedness thresholds.
- **Mixed genotype sources were pooled silently.** New WARN `MIXED_GENOTYPE_SOURCE`
  names which samples have reconstructed (AC/AN) genotypes when a cohort mixes
  them with native VCF genotypes. Per-source calibration is deliberately not done.
- A header-only VCF returned WARN instead of UNKNOWN from the provenance check,
  which let its audit exit 0 instead of 2 ("could not read the input").
- The ANNOVAR AC/AN fallback read a multiallelic `1/2` heterozygote (`AC=1,1`) as
  homozygous-alt, and accepted an internally inconsistent `AC=2,1` as a het; any AC
  whose parts do not sum to 2 now reads as missing.
- `compound_het_counts` could count a site whose unaffected genotype was missing,
  inflating the count.
- A two-node pedigree cycle (A's father is B, B's father is A) was accepted.
- Cohort-wide pairwise relatedness had no cap on sample count; a header declaring
  thousands of samples could run for a long time. Capped by
  `max_samples_for_pairwise` (2,000) with a reported finding.
- A cohort PED no longer names families that were not read, and the report header
  lists only the families actually supplied.

### Changed

- README and design notes: the claim that somalier/peddy score against a panel that
  coding-only call sets barely overlap was measured and found backwards (85.4% of
  peddy's GRCh37 panel is in RefSeq CDS); corrected, with the scripts to reproduce
  it in `validation/`.
- Python 3.14 added to CI; a Dockerfile and a project website added; PyPI install
  instructions (`pipx` recommended) added.

### Known limitations (found after 0.1.0, not fixed in this release)

- The duplicate warning needs a 2x difference in non-ref counts. A same-person pair
  with similar-sized call sets and Jaccard below 0.60 gets no finding (seen on
  synthetic data: agreement 1.000, Jaccard 0.54).
- Monozygotic twins reach duplicate-level agreement and are reported as
  `DUPLICATE_SAME_INDIVIDUAL` (blocking); there is no way to declare a twin pair
  (seen on synthetic data).
- Pairs sharing fewer than 500 non-ref sites are skipped for duplicate detection
  without a message.
- With no clean gap in the cohort's chrX heterozygosity the sex boundary falls back
  to 0.45; females with depressed chrX heterozygosity (consanguinity, long runs of
  homozygosity) can then be called male. Borderline calls should be verified.
- With no pedigree, a duplicate message still says "the pedigree declares them
  unrelated ACROSS FAMILIES".

## [0.1.0] — 2026-09-13

First release. All five checks implemented.

### Added

- **Check 1 — sample identity.** Runs with or without a pedigree: sex inference
  and cohort-wide duplicate detection need none, and an unpedigreed cohort is
  when "are any two of these the same person" matters most. Sex inference from chrX heterozygosity with PAR and
  XTR excluded, reported as a Wilson interval with an explicit *not determined* state.
  Cohort-wide duplicate detection (never within-family) with evidence-based typing.
  Relatedness by KING-robust kinship and IBS0 where samples were jointly genotyped, and
  by Jaccard/agreement with a cohort-calibrated boundary otherwise. Pedigree
  reconciliation against recursively computed expected kinship.
- **Check 2 — genotype confidence.** Allele depths read from `AD`, freebayes
  `RO`/`AO` or `DP4`, whichever the caller wrote, with the field used recorded in
  the evidence; VarScan's single-valued `AD` is deliberately not read as GATK's.
  False-homozygote detection separating
  `LOW_DEPTH_HOM`, `HOM_CONTRADICTED_BY_LIKELIHOOD` and `ALLELE_IMBALANCE_HOM`. Allele
  balance by a two-sided binomial test rather than a fixed window, computed in
  log space so that deep coverage (panel or amplicon depths in the thousands)
  cannot overflow it. Site-level quality flags,
  segmental-duplication and clustered-variant flags, and detection of call sets that
  were never variant-quality filtered.
- **Check 3 — content provenance.** `CODING_ONLY`, `PASS_FILTERED`, `SUBSET`,
  `METRICS_STRIPPED`, `HOMREF_STRIPPED` and `NOT_INTERVAL_RESTRICTED` verdicts,
  from annotation already present in the input. Region class is read from
  ANNOVAR keys or from VEP/SnpEff consequence terms, so the coding-only question
  is answerable on ordinarily annotated VCFs. A sites-only VCF reaches the check
  rather than disappearing for having no samples.
- **Check 4 — callability.** Depth bins are selected by their own lower bound, so
  any `--quantize` binning works and none of them can silently fall back to
  counting shallow bases as callable; a binning with no boundary at the depth
  floor reports UNKNOWN, and one whose bin straddles the floor is excluded and
  declared a lower bound. Joint callable fraction as a true interval intersection,
  per inheritance model, from `mosdepth --quantize` output. Prints the product of the
  per-sample fractions alongside it to show how far that common shortcut is wrong.
- **Check 5 — inheritance model sweep.** Gene assignment from ANNOVAR keys, VEP
  `CSQ` or SnpEff `ANN` (the pipe-delimited layout is read from the header rather
  than assumed), so compound-heterozygous is evaluable on ordinarily annotated
  VCFs. A de-novo count carries a caveat stating that an unfiltered count at
  exome scale is dominated by genotyping error. Nine models with per-model callable
  fractions; models the inputs cannot support report `not-applicable`, never zero.
- ANNOVAR `*_multianno.txt` reader, for analyses whose VCFs no longer exist.
  Where `table_annovar.pl --vcfinput` preserved the original FORMAT and sample
  columns, the genotype is read from them and `DP`/`GQ`/`PL`/`AD` are available
  to check 2; the `AC`/`AN` reconstruction remains the fallback for tables that
  genuinely carry no FORMAT block.
- **Pedigree input without hand-writing a PED.** `admissible ped-template
  *.vcf` prints a skeleton carrying the sample names as they appear inside the
  VCFs — the part that has to match — and `--csv` prints the spreadsheet shape
  instead. `--ped` also accepts a headered CSV/TSV/xlsx written in words
  (`M`/`female`, `yes`/`control`, a blank for a founder), with columns located
  by name. An unrecognised value becomes *unknown* and is reported rather than
  rounded to the plausible one. Format is decided by content, not by extension.
- Text and JSON reports (schema `admissible/report/1`), per-check CLI subcommands,
  pipeline exit codes.

### Validated against

- Public CEPH pedigree 1463 (17 samples, three generations, freebayes + VEP,
  GRCh37). Check 1: 0 sex errors, 0 false duplicates across 136 pairs, the correct
  pedigree accepted and a deliberately corrupted copy of it caught. Checks 2 and 5
  run on the same file, which is what exposed the `RO`/`AO` and `CSQ` gaps above.
  Check 3 is validated by applying one known transformation at a time to that
  same real file — coding-only extract, PASS-only delivery, a merge that drops
  hom-reference, a sites-only export — and requiring the tool to name that
  transformation and no other. Reproducible by anyone — see the validation
  section of the README.
- Real capture coverage: four unrelated 1000 Genomes exomes over
  chr20:1,400,000–1,500,000, with the denominator taken independently from an
  Ensembl GTF rather than from the coverage. Per-sample callable fractions
  0.4768–0.7975; joint fraction across all four 0.3815; product of the
  per-sample fractions 0.2232 — the shortcut is 41% low on real capture. 4,134 bp
  is a small measurement and carries real sampling noise; see
  `tests/data/1000g_exome_chr20/PROVENANCE.md`.
- Real GATK genotypes inside a reconstructed ANNOVAR column layout, which is
  what caught the multiallelic reconstruction error above. See
  `tests/data/annovar/PROVENANCE.md` for which half of that fixture is real.
- Synthetic fixtures with planted sample swaps, regenerated deterministically by
  scripts in `tests/fixtures/`.

### Known limitations

- Two samples homozygous for *different* alternate alleles at the same
  multiallelic site are not counted as IBS0, because the genotype encoding keeps
  only hom-ref / het / hom-alt. On CEPH 1463 this affects 0 pairs across 31
  multiallelic sites, which is why it is recorded rather than fixed; a call set
  rich in multiallelics would need a per-allele encoding.

- The two-locus model is not implemented: the pairwise search needs an explicit
  multiple-testing treatment before any count is reportable.
- Relationship degree beyond *related vs unrelated* is not adjudicated at typical
  exome site counts. See the validation section of the README.
- Tested on GRCh37 and GRCh38 contig naming; exercised in anger only on GATK,
  freebayes and ANNOVAR output.
