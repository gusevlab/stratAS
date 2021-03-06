# :first_quarter_moon: stratAS

Allele specific analyses across cell states and conditions

## Workflow

### Input data / dependencies

* A `.vcf` file containing phased individual genotypes and an "AS" field listing the allelic counts. This gets turned into a flat file to be loaded in the `--input` flag (see below).
* A samples file with the columns `ID` and `CONDITION` that is in the **same order** as the vcf input. `CONDITION` is currently restricted to 0/1. This file got into the `--samples` flag.
* A `.bed` file listing the physical positions of the features to be tested for allelic imbalance. This goes into the `--peaks` flag.
* Global/local parameter files generated by `params.R` (see below). At minimum a global overdispersion parameter file with the columns `ID` and `PHI` with entries in any order goes into the `--global_params` flag.
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

For tumor data, local CNV-specific parameters are also estimated, and the `--inp_cnv` file must be provided with headers `CHR P0 P1` listing the boundaries of CNV regions. This file can additionally contain a `CNV` column listing the CNV estimate (centered to zero as in TCGA calls) for inclusion as a covariate in the final analysis.

inference is then performed by running:
```
Rscript params.R --min_cov 5 --inp_counts $OUT.counts --inp_cnv $CNV --out $OUT
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

By default, the ASE test is printed to screen with each line containing the following entries:

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
| ALL.AF | Allelic fraction estimate from beta binomial test across both conditions |
| ALL.BBINOM.P | Beta-binomial test for imbalance across both conditions  |
| C0.AF | Allelic fraction estimate from condition 0 |
| C0.BBINOM.P | Beta-binomial test for imbalance in condition 0 |
| C1.AF | Allelic fraction estimate from condition 1 |
| C1.BBINOM.P | Beta-binomial test for imbalance in condition 1 |
| DIFF.BBINOM.P | Beta-binomial test for difference between conditions |

Enabling the `--binom` flag additionally runs a standard binomial test, and produces the following columns:

| Column | Description |
| --- | --- |
| ALL.BINOM.P | Binomial test for imbalance across both conditions|
| ALL.C0.BINOM.P | Binomial test for imbalance in condition 0 |
| ALL.C1.BINOM.P | Binomial test for imbalance in condition 1 |
| FISHER.OR | Fisher's test odd's ratio for difference between conditions|
| FISHER.DIFF.P | Fisher's test difference between conditions |

Enabling the `--bbreg` flag additionally runs a beta binomial regression with CNV as covariate, and produces the following columns:

| Column | Description |
| --- | --- |
| ALL.BBREG.P | Beta binomial regression (with covariates) for imbalance across both conditions |
| DIFF.BBREG.P | Beta binomial regression (with covariates) for imbalance difference between conditions |
| CNV.BBREG.P | Beta binomial regression (with covariates) for imbalance along CNV covariate |

Enabling the `--indiv` flag additionally produces the following columns:

| Column | Description |
| --- | --- |
| IND.C0 | Number of each condition 0 individual included in this test (comma separated) |
| IND.C0.COUNT.REF | condition 0 REF allele counts of each individual included in this test (comma separated) |
| IND.C0.COUNT.ALT | condition 0 ALT allele counts of each individual included in this test (comma separated) |
| IND.C1 | Number of each condition 1 individual included in this test (comma separated) |
| IND.C1.COUNT.REF | condition 1 REF allele counts of each individual included in this test (comma separated) |
| IND.C1.COUNT.ALT | condition 1 ALT allele counts of each individual included in this test (comma separated) |

## Example

An example locus with significant AS associations can be run by calling:

```
Rscript stratas.R \
--input example/ENSG00000075240.12.mat \
--samples example/KIRC.ALL.AS.PHE \
--peaks example/ENSG00000075240.12.bed \
--global_param example/KIRC.ALL.AS.CNV \
--local_param=example/KIRC.ALL.AS.CNVLOCAL
```
## Detailed Parameters

| Parameter | Description |
| --- | --- |
| `--input` | Path to input file |
| `--samples` | Path to sample identifier file, must have ID and CONDITION columns |
| `--peaks` | Path to file containing peak/gene boundaries, must contain CHR P0 P1 NAME CENTER columns |
| `--global_param` | Path to global overdispersion parameter file |
| `--local_param` | Path to local overdispersion parameter file |
| `--out` | Path to output |
| `--window` | Window (in bp) for SNPs to test around the peak boundary |
| `--perm` | # of permutations to shuffle the allele labels (0=off) |
| `--perm_cond` | # of permutations to shuffle the condition labels (0=off) |
| `--min_cov` | Individuals must have at least this many reads (for both alleles) to be tested |
| `--min_maf` | Minimum minor allele frequency for test SNP |
| `--min_het` | Minimum minor heterozygous frequency for test SNP |
| `--max_rho` | Maximum local/global over-dispersion parameter for which to include individual in test |
| `--binom` | Also perform a standard binomial test |
| `--bbreg` | Also perform a beta binomial regression with local CNV status as a covariate (must also provide `--local_param` file) |
| `--indiv` | Also report the per-individual allele fractions (Warning, this can produce large files) |
| `--exclude` | The mimium distance between SNPs allowed in the haplotype (to exclude variants in the same read) |

## Notes:

* Permutation: The permutation test permutes read counts with respect to alleles. This is a null of counts randomly sampled from the observed count distribution. If a small number of individuals have unusually high read counts and imbalance they will dominate the observed imbalance and appear highly significant even in permutation. In these cases, failing the permutation test may rule out true biological signal due to individual outliers.

* Unlike an eQTL, we can test for AS signal within each individual separately but this is currently not being evaluated. It may be interesting to know what fraction of individual heterozygotes exhibit significant imbalance (for example, through the Storey & Tibshirani pi statistic). Alternatively, individual standard errors can be approximated from each individual test and heterogeneity assessed using standard meta-analysis statistics.

* For a test SNP that is not in the target feature, it may be of interest to report the allelic imbalance for homozygous individuals (that have heterozygous read-carrying variants). If the test SNP explains all of the imbalance at the locus we would expect these homozygous individuals to follow the null. This could also be a way to identify secondary associations.
