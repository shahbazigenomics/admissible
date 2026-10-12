# Hemizygous chrX and the depth floor — documented, not implemented

**Status: proposal. Nothing here is in the code.**

## What happens today

Check 4 applies one depth floor (`--min-depth`, default 10x) to every position of
every sample. For a male, the non-pseudoautosomal part of chrX has one copy, so the
same library gives roughly half the depth it gives on an autosome (an expectation to
be measured, not a figure this project has measured). A 10x floor on a male X
therefore asks for something like 20x-equivalent library depth.

Consequence: callable fractions for X-linked models in male samples are biased **low**.
That is the conservative direction, so it does not cause an overclaim, but it can make
an adequately sequenced X-linked gene read as `EXPLORATORY`. It matters little for the
whole-target headline (X is a few percent of an exome) and a lot for a single X-linked
candidate gene, which is exactly what `--by-region` is for.

## Proposal

A **sex-aware floor**, off by default, behind an explicit flag (for example
`--sex-aware-floor`):

- For a sample recorded as male, positions in chrX outside the PAR use
  `min_depth_x_male` (prior: `ceil(min_depth / 2)`); everywhere else, and for every
  female and every sample of unknown sex, the ordinary floor applies.
- The reasoning is that a hemizygous call needs support for one allele, not for two
  (a heterozygote needs both), and the male library delivers half the reads per copy.
  Both halves of that are plausible and neither is measured here.

## What it needs before it could be enabled

1. **Sex source.** The PED sex column only. Check 4 should not infer sex from data:
   check 1 already does that and may disagree with the PED. If the PED sex is 0 or
   missing, keep the ordinary floor and say so. If check 1 disagrees with the PED for
   a sample, the lowered floor should not be applied to it.
2. **Build.** PAR coordinates differ by build and a BED does not reveal the build, so
   `--build` must be given; without it the floor is not changed.
3. **Quantization.** A mosdepth file only resolves depth at its own bin boundaries.
   Lowering the floor to 5x needs a boundary at 5 (`--quantize 0:5:10:`); otherwise
   the bin straddles the floor and, as today, is excluded and flagged as a lower
   bound.
4. **Per-position floors in the interval code.** The callable set of a sample becomes
   the union of "≥ ordinary floor outside male non-PAR X" and "≥ lowered floor inside
   it". That is a region split before the existing intersection, not a change to the
   intersection itself.
5. **Evidence, before any default.** On real male WES, (a) the ratio of X to
   autosomal depth, and (b) hemizygous genotype concordance against a truth set at
   depths 5–9. The factor 0.5 is a prior. If (b) does not support calling at 5x, the
   proposal should not ship.
6. **Disclosure.** When the flag is used, the claim text and the finding must say
   that a lowered floor was applied to which samples, because the result is then more
   permissive than the default, which is the direction this project is wary of.

## Not covered

Females with X-linked recessive disease need both alleles resolved and keep the
ordinary floor. Y and the mitochondrial genome are separate questions (the
mitochondrion has no copy-number-2 assumption at all, and heteroplasmy is not a depth
question).
