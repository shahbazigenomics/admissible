# Check 4 (callability) in 0.2.0: status per problem

Branch `feat/callability-0.2.0`, off 4328c4c (v0.1.1). Identity, provenance,
genotype and models behave as in 0.1.1. Terms: **confirmed** = the problem existed
in 0.1.1, found by reading its code or running it, with a regression test added;
**fixed** = changed and covered by a test; **documented only** = written down, not
changed.

**Read this first: none of 0.2.0 has been validated against real capture data or
against an independent tool.** The tests are synthetic. The `validation/callability/`
scripts were written but not run (no bedtools, mosdepth or kit data in the
environment where they were written).

| # | problem | verdict | what was done | what is still open |
|---|---|---|---|---|
| 1 | one number per run | confirmed, fixed | `--by-region`: one row per named region, worst region per gene, gene-level claim. Unnamed output unchanged (test). | Opt-in, so a caller must ask. Region rows are for each unit's strictest model, not every model. Gene names are a regex on the region name, not annotation. |
| 2 | claim ladder overstates | confirmed, fixed | "excluded" removed; claims say what was searched. Thresholds were *ignored* by the claim (a bug found while fixing this) and now drive it. A gene with a wholly uncallable or unreported region cannot exceed `SEARCHED_PARTLY` (configurable). | The 0.80/0.50 thresholds remain conventions. The cap is on the gene row; the headline is deliberately not capped (a whole exome always has some blind exon). |
| 3 | joint fraction counted affected only | confirmed, fixed | Each model is held to the samples it needs (affected plus parents, trio union, phenocopy n-1), listed in `required_samples`; a missing required file makes the model `UNKNOWN`, not a higher number. Added in the same work: judging per family / per sample (`--independent`) and `--affected-only`, because pooling an unrelated cohort reported 0.38 for twelve cases that were each 0.93. | Unaffected relatives other than parents are not required, although check 5 drops sites where they were unobserved, so the share can exceed what check 5 could verify. |
| 4 | exclusion BED | implemented | `--exclude-bed`; excluded bases count as not searched, denominator stays the full target, claim follows the stricter figure, `without_exclusion` kept. | Counting excluded bases only in the numerator is a design choice, untested on real data. No exclusion track is shipped. |
| 5 | target/kit sanity | partly as specified | WARN findings, never a hard fail: target contig missing from the coverage file (and `NOT_REPORTED` kept apart from "uncallable"), and a kit-size-gated share-in-target check. | **Build mismatch is not detectable** from contig names, so it is not detected. The share check's thresholds (5 Mb gate, 0.50) are uncalibrated and unvalidated; it may false-alarm or miss. |
| 6 | hemizygous chrX | documented only | `docs/hemizygous-chrx-floor.md`: current bias (conservative, X-linked fractions in males biased low) and a sex-aware-floor proposal with what it needs. | Not implemented. The 0.5 depth factor is a prior, not a measurement. |
| 7 | validation scripts | written, **not run** | `validation/callability/`: bedtools baseline, comparator, driver, two-kit requirement. Comparator logic is unit-tested. | Everything end to end. Shell scripts passed `bash -n` and `shellcheck` only. |

## Known limit, applies to every claim above

Depth-based callability is blind to heterozygous deletions and other structural
variants: a deleted allele leaves depth above any floor, so the region reads as
callable. No fraction here speaks to them. This is stated beside the claim text in
the README and in every finding's `limits`.

## Not confirmed

- That the old claim text was *wrong* for any particular dataset: it overstated the
  evidence by construction, but no real case was shown to have been misled by it.
- That the new warnings fire usefully on real kits.
- That VarLens consumes the new JSON as intended; the field names are stable only in
  the sense that no existing name was removed or renamed except `CLAIM_LADDER` and
  `full_penetrance_models` (see CHANGELOG).

## Unrelated flaw noticed, not fixed

`vcfio.load_cohort` has the same corrupt-gzip `zlib.error` flaw that `read_bed` had.
Outside this brief, so left as it was.
