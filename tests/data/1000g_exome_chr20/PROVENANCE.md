# Real exome coverage — 1000 Genomes, chr20:1,400,000–1,500,000

These are the only files in this repository derived from real sequencing. They
exist because check 4 could otherwise only ever be tested on coverage invented
by the same person who wrote the code, and invented coverage cannot settle the
question check 4 exists to answer: *how far does the product of per-sample
callable fractions understate the joint one on a real capture?* That gap is set
by real probe-efficiency correlation, and nothing simulated can establish it.

## What these are

`*.quantized.bed.gz` — `mosdepth 0.3.8 --by <target> --quantize 0:10:` output for
four unrelated Finnish individuals from the 1000 Genomes exome project
(Baylor College of Medicine, study SRP000808, 101 bp Illumina reads), restricted
to chr20:1,400,000–1,500,000 (GRCh38):

    HG00349  HG00350  HG00351  HG00358

The alignments they were computed from are the public test subsets distributed by
[nf-core/test-datasets](https://github.com/nf-core/test-datasets) under
`data/genomics/homo_sapiens/illumina/bam/`. Only the coverage summaries are
vendored here, not the alignments.

`chr20_coding_target.bed` — the denominator, and deliberately **not** derived
from the coverage: protein-coding CDS intervals from Ensembl release 111
(`Homo_sapiens.GRCh38.111_chr20.gtf`, also via nf-core/test-datasets), merged and
padded 50 bp either side to approximate probe placement. 19 intervals, 4,134 bp.

## What they establish, and what they do not

Measured on this data: per-sample callable fractions 0.4768–0.7975, joint
callable fraction across all four **0.3815**, and the product of the per-sample
fractions **0.2232**. The product understates the truth by 41%, on real capture,
which is the claim check 4 is built around.

4,134 bp over four samples is a small measurement and the figure carries real
sampling noise; it is a demonstration on real data, not an estimate of the gap
in general. A larger and less noisy figure needs a full exome, which is not
something this repository can vendor.

The four individuals are unrelated. The pedigree used in the test declares them
affected purely so that check 4 has affected samples to intersect — it asserts no
relationship and changes no coverage arithmetic.
