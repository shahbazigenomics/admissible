# ANNOVAR multianno fixture — real genotypes, reconstructed column layout

`gatk_chr21.hg38_multianno.txt` is 120 records, and it exists because the
ANNOVAR adapter was the last reader in this repository never run against
anything but a table written by the same person who wrote the adapter.

## What is real and what is not

**Real:** every genotype, every FORMAT field (`GT:AD:DP:GQ:PL`), every INFO
field and every position. They come verbatim from a public GATK
HaplotypeCaller call set —
`data/genomics/homo_sapiens/illumina/gatk/haplotypecaller_calls/test2_haplotc.vcf.gz`
in [nf-core/test-datasets](https://github.com/nf-core/test-datasets), sample
`disease_103`, chr21, GRCh38. Seven of the 120 records are genuine `1/2`
multiallelic calls; they are why this fixture was selected rather than the
first 120 rows.

**Not real:** the ANNOVAR column wrapper around them — the `Func.refGene`,
`Gene.refGene` and `ExonicFunc.refGene` values are placeholders, and the
`Otherinfo1`–`Otherinfo13` layout was reconstructed from ANNOVAR's documented
`--vcfinput` output format. No genuine `*_multianno.txt` file is published
anywhere reachable from this environment, so the layout is faithful to the
specification rather than to a file someone actually produced. If you have a
real one, running it through `admissible genotype` is worth more than this
fixture.

`expected_genotypes.json` is the ground truth: `[chrom, pos, ref, alt, GT]` per
record, read straight out of the sample column.

## What this fixture caught

The adapter reconstructed every genotype from `AC`/`AN` even when ANNOVAR had
preserved the real `FORMAT` and sample columns two positions to the right. On
these 120 records that reconstruction got **113 right and 7 wrong** — every
multiallelic `1/2` call became a homozygous alternate, because `AC=1,1` sums to
two. A heterozygote read as a homozygote is the single failure this whole tool
was built to prevent, and it fed straight into the recessive model in check 5.

Reading the genotype that is actually in the file gives 120 of 120, and hands
check 2 the `DP`, `GQ`, `PL` and `AD` it previously did without.
