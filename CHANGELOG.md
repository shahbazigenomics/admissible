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
- **Check 2 — genotype confidence.** False-homozygote detection separating
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
- **Check 5 — inheritance model sweep.** Nine models with per-model callable
  fractions; models the inputs cannot support report `not-applicable`, never zero.
- ANNOVAR `*_multianno.txt` reader, for analyses whose VCFs no longer exist.
- Text and JSON reports (schema `admissible/report/1`), per-check CLI subcommands,
  pipeline exit codes.

### Validated against

- Public CEPH pedigree 1463 (17 samples, three generations, freebayes, GRCh37):
  0 sex errors, 0 false duplicates across 136 pairs, the correct pedigree accepted
  and a deliberately corrupted copy of it caught. Reproducible by anyone — see the
  validation section of the README.
- Synthetic fixtures with planted sample swaps, regenerated deterministically by
  scripts in `tests/fixtures/`.

### Known limitations

- Check 4 has been exercised end-to-end on genuine `mosdepth --quantize 0:10:`
  output, but the BAMs behind that output were simulated. Simulated coverage can
  confirm that the quantize format is read correctly and that the intersection
  arithmetic is right; it cannot establish how far the product understates the
  joint fraction on a real capture, because that gap is set by the real
  correlation structure of probe efficiency. Treat the printed product-vs-
  intersection gap as demonstrated in principle and unmeasured in practice until
  it has been run on a real exome.
- The two-locus model is not implemented: the pairwise search needs an explicit
  multiple-testing treatment before any count is reportable.
- Relationship degree beyond *related vs unrelated* is not adjudicated at typical
  exome site counts. See the validation section of the README.
- Tested on GRCh37 and GRCh38 contig naming; exercised in anger only on GATK,
  freebayes and ANNOVAR output.
