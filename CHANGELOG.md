# Changelog

All notable changes to this project are documented here.
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- **A declared female called male now comes with the evidence that separates a sample
  swap from a real female with a homozygous X.** chrX heterozygosity cannot tell
  them apart (consanguinity and long runs of homozygosity depress it), so the call
  and the BLOCKING severity are unchanged. The `SEX_MISMATCH` message and evidence
  now add her autosomal het fraction next to the cohort median (runs of homozygosity
  lower it genome-wide, a swap does not) and, when the cohort's males have any, her
  chrY call count. The autosomal het fraction is also reported for every sample.
  The comparison group is samples from other families when there are at least three,
  because consanguinity lowers autosomal het for a whole family and comparing a
  member with her own relatives would hide it; otherwise it is all other samples and
  the message says relatives weaken the comparison.
  Checked on a simulated cohort only; no consanguineous data was available, so no
  interpretation threshold is applied.
- An attempted fix that refused male calls above a fixed chrX het of 0.30 was
  dropped: male het is 0.12 on CEPH 1463 but 0.35 in the test fixture, so a fixed
  ceiling is wrong for some datasets.

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
  not. It does not auto-classify a duplicate: a human verifies. The first version
  also required the non-ref counts to differ by at least `dup_completeness_ratio`
  (2.0); that requirement is dropped, because the same person called by two
  pipelines gives similar-sized call sets with low Jaccard (synthetic: agreement
  1.000, Jaccard 0.54) and was still missed. The ratio now only chooses which
  explanation the message offers.
- **A PED header naming parents `dad_id` / `mom_id` silently lost every parent
  link.** Parent columns were matched against a fixed list that lacked those
  spellings, so every sample was read as a founder and relationship checks ran
  against a pedigree with no relationships in it. A trailing `id` is now accepted on
  any listed spelling, and a header with no recognisable father or mother column
  now warns that parent links are lost. (Found by running peddy's
  `ceph1463.bad.ped`, which uses those headers.)
- **A sibling declared as a parent was not caught.** Kinship is 0.25 for both, so
  the kinship comparison cannot tell them apart. New BLOCKING
  `PARENT_OFFSPRING_NOT_SUPPORTED` fires for a declared parent-offspring pair with
  at least 20 IBS0 sites and more than 0.005 IBS0 per heterozygous call (a true
  parent and child have almost none). On CEPH 1463, true parent-offspring pairs were
  at most 0.0008 and full sibs at least 0.0145. **This changes a 0.1.0 result:** on
  the full published CEPH 1463 pedigree with peddy's VCF it flags 7 declared
  parent-offspring pairs (NA12877 with NA12882, NA12883, NA12884, NA12886, NA12888,
  NA12893, and NA12889 with NA12877), which have sibling-like IBS0 in that file. The
  cause (pedigree or this extract) is not established. peddy's own
  `ceph1463.good.ped` keeps only two of those children and passes.
- **Pairs with under 500 shared sites were skipped silently.** New INFO
  `DUPLICATE_CHECK_SKIPPED` states how many pairs were not examined, so a clean
  result is not read as covering them.
- **No-pedigree duplicate message was wrong.** With no pedigree, or a sample missing
  from it, the message said "the pedigree declares them unrelated ACROSS FAMILIES".
  It now says the pedigree does not list both samples; the claim is made only when
  both are in the pedigree.
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

- A blocking duplicate between two samples that the pedigree lists as same-sex full
  siblings (or relatives of unspecified degree) now also says it is compatible with
  monozygotic twins, but only if the pair is already known to be identical twins;
  otherwise treat it as a swap or duplicate. Still blocking, and the evidence carries
  `possible_mz_twins`. Nothing is downgraded automatically.
- README and design notes: the claim that somalier/peddy score against a panel that
  coding-only call sets barely overlap was measured and found backwards (85.4% of
  peddy's GRCh37 panel is in RefSeq CDS); corrected, with the scripts to reproduce
  it in `validation/`.
- Python 3.14 added to CI; a Dockerfile and a project website added; PyPI install
  instructions (`pipx` recommended) added.

### Known limitations (found after 0.1.0, not fixed in this release)

- Monozygotic twins reach duplicate-level agreement and are still reported as
  `DUPLICATE_SAME_INDIVIDUAL` (blocking); there is no way to declare a twin pair
  yet. A `--mz-twins` declaration is planned for 0.2.0.
- With no clean gap in the cohort's chrX heterozygosity the sex boundary falls back
  to 0.45; females with depressed chrX heterozygosity (consanguinity, long runs of
  homozygosity) can then be called male. Borderline calls should be verified.
  Hardening this changes calls, so it is planned for 0.2.0.
- The parent-offspring IBS0 check only looks at pairs the pedigree declares as
  parent and child. It does not yet catch the reverse (siblings declared, a parent
  and child in fact) and has not been tested on consanguineous families.
- The duplicate thresholds (agreement 0.95, Jaccard 0.60, 500 shared sites) are
  conventions, not derived values. On real CEPH 1463 data the highest non-duplicate
  agreement was 0.826.

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
