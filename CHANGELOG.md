# Changelog

All notable changes to this project are documented here.
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] — unreleased

First release. All five checks implemented.

### Added

- **Check 1 — sample identity.** Sex inference from chrX heterozygosity with PAR and
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
  balance by exact binomial test rather than a fixed window. Site-level quality flags,
  segmental-duplication and clustered-variant flags, and detection of call sets that
  were never variant-quality filtered.
- **Check 3 — content provenance.** `CODING_ONLY`, `PASS_FILTERED`, `SUBSET` and
  `METRICS_STRIPPED` verdicts, from annotation already present in the input.
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
- Text and JSON reports (schema `admissible/report/1`), per-check CLI subcommands,
  pipeline exit codes.

### Validated against

- Public CEPH pedigree 1463 (17 samples, three generations, freebayes + VEP,
  GRCh37). Check 1: 0 sex errors, 0 false duplicates across 136 pairs, the correct
  pedigree accepted and a deliberately corrupted copy of it caught. Checks 2 and 5
  run on the same file, which is what exposed the `RO`/`AO` and `CSQ` gaps above.
  Reproducible by anyone — see the validation section of the README.
- Real capture coverage: four unrelated 1000 Genomes exomes over
  chr20:1,400,000–1,500,000, with the denominator taken independently from an
  Ensembl GTF rather than from the coverage. Per-sample callable fractions
  0.4768–0.7975; joint fraction across all four 0.3815; product of the
  per-sample fractions 0.2232 — the shortcut is 41% low on real capture. 4,134 bp
  is a small measurement and carries real sampling noise; see
  `tests/data/1000g_exome_chr20/PROVENANCE.md`.
- Synthetic fixtures with planted sample swaps, regenerated deterministically by
  scripts in `tests/fixtures/`.

### Known limitations

- The two-locus model is not implemented: the pairwise search needs an explicit
  multiple-testing treatment before any count is reportable.
- Relationship degree beyond *related vs unrelated* is not adjudicated at typical
  exome site counts. See the validation section of the README.
- Tested on GRCh37 and GRCh38 contig naming; exercised in anger only on GATK,
  freebayes and ANNOVAR output.
