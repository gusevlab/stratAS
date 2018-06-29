# :first_quarter_moon: stratAS

Allele specific analyses across cell states and conditions

## Workflow

### Input data / dependencies

* A *.vcf file containing phased individual genotypes and an "AS" field listing the allelic counts.
* (Optional, for tumors) A list of structural variant boundaries for which to estimate local distribution parameters.
* R, and the `VGAM`, `optparse` libraries installed.

### Recommended pre-processing and allelic counts

* Phase your genotype data with EAGLE/HRC.
* Process all sequence data with the [WASP](https://github.com/bmvdgeijn/WASP) mapping pipeline.
* Use the GATK ASEReadCounter to count reads then convert into vcf (see `pipeline/` scripts).

### `params.R` : estimating prior parameters

`params.R` infers the read count distribution (beta binomial) parameters and needs to be run on a whole genome vcf one individual at a time. The input is an `--inp_counts` file which must contain the headers `CHR POS HAP REF.READS ALT.READS` where `HAP` is the `0|1` or `1|0` vcf haplotype code.

A vcf file containing one individual is converted to counts as follows:
```
zcat $VCF \
| grep -v '#' \
| cut -f 1,2,10 \
| tr ':' '\t' \
| awk '$3 != "1|1" && $3 != "0|0"' \
| cut -f 1-3,7 \
| tr ',' '\t' \
| sed 's/chr//' \
| awk 'BEGIN { print "CHR POS HAP REF.READS ALT.READS" } $4 + $5 > 0' \
 > $OUT.counts
```

For tumor data, local CNV-specific parameters are also estimated, and the `--inp_cnv` file must be provided with headers `CHR P0 P1` listing the boundaries of CNV regions.

inference is then performed by running:
```
Rscript params.R --min_reads 5 --inp_counts $OUT.counts --inp_cnv $CNV --out $OUT
```

A `$OUT.global.params` file is generated containing the parameters (with header `PHI MU N`) and (optionally) an `$OUT.local.params` file is generated containing the positions and parameters for each CNV (with header `CHR P0 P1 PHI MU N`).

### `stratas.R` : testing for differential AS

`stratas.R` computes the AS statistics from a VCF of all individuals and the prior parameters estimated above.

A vcf file containing all individuals is converted to counts and split into batches as follows:
```
zcat $VCF\
| grep -v '#' \
| cut -f 1-5,9- \
| tr ':' '\t' \
| awk '{ for(i=1;i<=5;i++) { printf "%s ",$i; } for(i=6;i<=10;i++) { if($i == "GT") gtnum=i-6; if($i=="AS") asnum=i-6; } for(i=11;i<=NF;i++) { if( (i-11) % 5 == asnum || (i-11) % 5 == gtnum ) printf " %s",$i; } print ""; }' \
| sed 's/[|,]/ /g' | split -d -l 15000 - $OUT.MAT.split.
```

Each batch is then processed as follows:
```
Rscript stratas.R \
--input $OUT.MAT.split.00 \
--samples SAMPLES.ID \
--peaks gencode.protein_coding.transcripts.bed \
--global_param GLOBAL.params \
--local_param LOCAL.params \
```

### Output data

The ASE test is printed to screen with each line containing the following entries:

| Column | Description |
| --- | --- |
| CHR | Chromosome |
| POS | Position of test SNP |
| RSID | ID of test SNP |
| P0 | Start of gene/peak  |
| P1 | End of gene/peak |
| NAME | Name of gene/peak |
| CENTER | Center position of peak (or TSS for gene) |
| N.HET | # of heterozygous individuals tested |
| N.READS | # of reads tested in total |
| BINOM.EST | Allelic fraction estimate from standard binomial test across both conditions [*] |
| BINOM.P | Binomial test for imbalance across both conditions [*] |
| BBINOM.EST | Allelic fraction estimate from beta binomial test across both conditions |
| BBINOM.P | Beta-binomial test for imbalance across both conditions  |
| FISHER.EST | Fisher's test odd's ratio for difference between conditions [*] |
| BINOM.C0.P | Binomial test for imbalance in condition 0 [*] |
| BBINOM.C0.P | Beta-binomial test for imbalance in condition 0 |
| BINOM.C1.P | Binomial test for imbalance in condition 1 [*] |
| BBINOM.C1.P | Beta-binomial test for imbalance in condition 1 |
| FISHER.DIFF.P | Fisher's test difference between conditions [*] |
| BBINOM.DIFF.P | Beta-binomial test for difference between conditions |

*\* these fields are only reported if `--binom` is set* 

## Example

An example locus with significant AS associations can be run by calling:

```
Rscript stratas.R \
--input example/ENSG00000075240.12.mat \
--samples example/KIRC.ALL.AS.PHE \
--peaks example/ENSG00000075240.12.bed \
--global_param example/KIRC.ALL.AS.CNV \
--local_param=example/KIRC.ALL.AS.CNVLOCAL \
< stratas.R
```

## Work in progress

* Permutation: The permutation test permutes read counts with respect to alleles. This is a null of counts randomly sampled from the observed count distribution. If a small number of individuals have unusually high read counts and imbalance they will dominate the observed imbalance and appear highly significant even in permutation. In these cases, failing the permutation test may rule out true biological signal due to individual outliers.

* Unlike an eQTL, we can test for AS signal within each individual separately but this is currently not being evaluated. It may be interesting to know what fraction of individual heterozygotes exhibit significant imbalance (for example, through the Storey & Tibshirani pi statistic). Alternatively, individual standard errors can be approximated from each individual test and heterogeneity assessed using standard meta-analysis statistics.

* For a test SNP that is not in the target feature, it may be of interest to report the allelic imbalance for homozygous individuals (that have heterozygous read-carrying variants). If the test SNP explains all of the imbalance at the locus we would expect these homozygous individuals to follow the null. This could also be a way to identify secondary associations.