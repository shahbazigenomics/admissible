# Recorded validation run: peddy on CEPH 1463, full vs. coding-only

Committed so the D12 sex-check claim has evidence in the repo, not just a
script someone has to rerun to check. Produced by
`validation/ceph_coding_only_head_to_head.py`'s peddy-side path.

**Run provenance:** peddy 0.4.8 (installed via `conda install -c bioconda
peddy`, since peddy's legacy `setup.py` does not build against modern
setuptools/pip), run 2026-09-28 by Amir on macOS against the real CEPH 1463
VCF from the `peddy` repo itself and a RefSeq-CDS-filtered copy of it built
with `bcftools view -R`.

`ceph_sex_check_full.csv` and `ceph_sex_check_coding_only.csv` below are
verbatim copies of peddy's own `out.sex_check.csv` output for each run.

## Result

Full VCF: all 17 samples' `predicted_sex` matches `ped_sex`, `error` is
`False` for all 17.

Coding-only VCF: 14 of 17 still match; the other 3 (`NA12881`, `NA12887`,
`NA12890`) are all pedigree-declared `female` predicted `male` by peddy
(`error=True`), with `het_ratio` dropping under 1.0 once the site count is
too small for the het/hom-alt ratio to still read as female. No sample went
the other direction (male predicted female).

This confirms, rather than merely repeats, the "17/17 -> 14/17, three true
females called male" claim in the README and D12 of
`admissible/design-decisions.md` — it was flagged as unverified as of
2026-09-28 morning and confirmed later the same day once a working peddy
install was available. `admissible`'s own behavior on the same two files
(0 sex mismatches full; clean `SEX_NOT_DETERMINED` decline, not a wrong
guess, on all 17 in the coding-only file) is unchanged and was already
confirmed by the same script without needing peddy.
