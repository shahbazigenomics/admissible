# Validation scripts

Not run in CI (they need network access and, for one of them, `bcftools`/
`tabix` on `PATH`) and not part of the installed package - `pip install
admissible` does not pull these in or run them. They exist to make two
specific claims in the README and in `admissible/design-decisions.md`
reproducible from committed code rather than asserted from memory.

## `panel_overlap.py`

Measures what fraction of peddy's real GRCh37 site panel (23,770 sites) falls
inside a real coding BED (RefSeq CDS) and a real clinical exome capture BED
(Agilent SureSelect V6), and what fraction of a real GATK call set survives
the standard hard-filter thresholds (QD/FS/MQ/MQRankSum/ReadPosRankSum/SOR).
Clones two small public repos on first run (cached under
`$TMPDIR/admissible_validation_cache`); the hard-filter measurement reuses the
real GATK INFO fields already committed at `tests/data/annovar/` and needs no
network access of its own.

This directly reproduces the correction recorded in D12: peddy's panel is
concentrated in coding regions (~85-87% overlap), not "largely intronic and
intergenic" as an earlier draft of this README claimed - see D12 for how that
was found to be wrong and why the panel comes out this way.

```
python validation/panel_overlap.py
```

## `ceph_coding_only_head_to_head.py`

Builds a coding-only-filtered copy of the real CEPH 1463 VCF and runs
`admissible`'s own identity check against both the full and filtered files,
confirming: 0 sex mismatches on the full file, and a clean decline
(`SEX_NOT_DETERMINED` on all 17 samples, not a wrong guess) on the
coding-only file, once chrX site count drops below the calibrated floor -
while relatedness still resolves, since it draws on the full autosomal panel.

The companion claim - that peddy's own sex calls degrade on the same
coding-only file (17/17 -> 14/17) rather than declining - needs a real peddy
install to verify, which does not build with modern setuptools in the
environment this script was written in. The script runs peddy automatically
if it is importable and tells you plainly if it is not, rather than asserting
a number nobody in this repo's CI has actually reproduced. **If you have a
working peddy install, please run this and report back what its
`sex_check.csv` actually says** so that specific comparison can be stated as
verified rather than reported (see the script's own output for exactly what
that involves).

```
python validation/ceph_coding_only_head_to_head.py
```
