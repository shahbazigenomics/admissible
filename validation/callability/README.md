# Check 4 against an independent baseline (bedtools)

**Status: written, not run.** The environment these were written in had neither
`bedtools` nor `mosdepth`, and no capture-kit data. Nothing in admissible 0.2.0 has
been validated by these scripts, so no claim in the README or CHANGELOG rests on
them. The only part exercised so far is `compare.py`'s comparison logic, on
constructed inputs (`tests/test_validation_compare.py`). The shell scripts have
passed `bash -n` and `shellcheck` and nothing more; expect to fix small things on first use.

## What it checks

For each sample, whether admissible's callable-base counts equal the counts obtained
by plain interval arithmetic in `bedtools`/`awk`, which shares no code with it:

| compared | admissible | baseline |
|---|---|---|
| target size | `target_bp` | `bedtools merge` of the target |
| callable in target | `models.recessive.callable_bp`, `joint_callable_bp` | `bedtools intersect` |
| callable minus exclusion | `models.recessive.callable_bp` (with `--exclude-bed`) and `without_exclusion` | `bedtools subtract` |
| per region | `by_region` rows | `bedtools intersect -wao`, summed per name |

Every number is an integer base count, so the expectation is **exact equality**.
A one-base difference is a disagreement about half-open interval arithmetic, not
noise. Run on a single sample with a one-line PED (`FAM sample 0 0 0 2`), so the
multi-sample, per-unit and pedigree logic is *not* what this exercises; those are
covered by synthetic tests only.

## What it does not show

- It checks arithmetic against arithmetic. It does **not** show that depth ≥ 10x is
  the right definition of callable, that the 0.80/0.50 thresholds mean anything, or
  that the kit-mismatch thresholds (uncalibrated) are sensible.
- Both sides strip `chr` and map `M` to `MT`, so contig normalisation is a shared
  assumption.
- Regions with the same name are assumed not to overlap each other; if they do, the
  baseline double-counts and a difference is reported that is not an admissible bug.
- Two or more kits are required (`run_validation.sh` reports INCOMPLETE otherwise),
  because one kit cannot show the result is not an accident of that kit's BED.

## How to run

1. Pick samples and their kits' target BEDs, with column-4 region names. This
   directory deliberately names no dataset: none has been chosen and checked.
2. Produce each sample's quantized coverage (check the flags against your mosdepth
   version; these were not run here):

       mosdepth --quantize 0:10: PREFIX sample.bam     # -> PREFIX.quantized.bed.gz

3. Copy `kits.example.tsv` to `kits.tsv`, fill it in, then:

       ./run_validation.sh kits.tsv out/

Per sample it leaves `baseline/`, `admissible.json`, `result.tsv` (every compared
number with a verdict) and `compare.txt` in `out/KIT/SAMPLE/`, and `versions.txt` in
`out/`. Exit 0 means every number was identical on at least two kits; 1 a
disagreement; 2 INCOMPLETE; 3 a tool is missing.

Requirements: `bedtools`, `python3`, `admissible` on `PATH`. bash 3.2 (macOS) is
supported by intent but untested.
