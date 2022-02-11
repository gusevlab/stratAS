# :first_quarter_moon: stratAS

A method for allele specific analyses across cell states and conditions. `stratAS` implements population-scale tests and predictive models for allele-specific activity (typically RNA or ATAC data).

## Workflow

### Input data / dependencies

* A `.vcf` file containing phased individual genotypes and an "AS" field listing the allelic counts. This gets turned into a flat file to be loaded in the `--input` flag (see below).
* A samples file with the columns `ID` and `CONDITION` that is in the **same order** as the vcf input. `CONDITION` is currently restricted to 0/1. This file goes into the `--samples` flag.
* A `.bed` file listing the physical positions of the features to be tested for allelic imbalance. This file goes into the `--peaks` flag.
* Global/local parameter files generated by `params.R` (see below). At minimum a global overdispersion parameter file with the columns `ID` and `PHI` with entries in any order goes into the `--global_params` flag.
* R, and the `VGAM`, `optparse` libraries installed.

### Recommended pre-processing and allelic counts

* Phase and impute your genotype data with EAGLE on the Haplotype Reference Consortium imputation server.
* Process all sequence data with the [WASP](https://github.com/bmvdgeijn/WASP) mapping pipeline.
* Use the GATK ASEReadCounter to count reads then convert into vcf (see `pipeline/` scripts).

### Estimating prior parameters

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

### Testing for population allele specificity

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

### Learning predictive models for TWAS/RWAS/CWAS analysis

`stratas.R` can compute multi-variate predictive models for downstream integration with GWAS data, using both allelic and total activity.

Inputs are similar to the tests for allele specificity described above and additionally require a total activity matrix, covariates, and the predict flags:

```
Rscript stratas.R \
--input $OUT.MAT.split.00 \
--samples SAMPLES.ID \
--peaks gencode.protein_coding.transcripts.bed \
--global_param GLOBAL.params \
--local_param LOCAL.params \
--total_matrix KIRC.gexp \
--covar KIRC.gexp.cov \
--predict_snps HM3.extract \
--predict --predict_only > $OUT.profile
```

This generates predictive model files for every feature in the TWAS/FUSION format, which can then be integrated with GWAS data using the [FUSION](http://gusevlab.org/projects/fusion/) software. The generated model types inside each model file are `"lasso", "lasso.as", "lasso.combined", "top1.as", "top1", "top1.combined"`; where `lasso` versus `top` corresponds to full locus penalized versus top SNP models; and `*`, `*.as`, `*.combined` indicates models using total activity, allele-specific activity, or both jointly. An `$OUT.profile` file is generated which lists model performance characteristics.

 **Important**: This analysis differs from the individual allele-specific tests in two ways: (1) the CONDITION value is not used and one model is trained on all samples; (2) models are still trained even if no allele-specific informatino is available. To restrict to models with both sources of signal (which are likely to be more accurate as a group), run the allele-specific tests described above and retain only peaks that have non-NA values.

Example total activity input files (KIRC.gexp, KIRC.gexp.cov, and HM3.extract) based on TCGA data are provided in the `example/` directory.

## Output data

By default, the outputs are printed to screen with each line containing the following entries:

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

Enabling the `--total_matrix` flag additionally runs a standard eQTL test (and interaction term), producing the following columns:

| Column | Description |
| --- | --- |
| ALL.TOTAL.Z | Test statisic from linear model across all samples |
| ALL.TOTAL.P | P-value from linear model across all samples |
| C0.TOTAL.Z | Test statisic from linear model across condition 0 |
| C0.TOTAL.P | P-value from linear model across condition 0 |
| C1.TOTAL.Z | Test statisic from linear model across condition 1 |
| C1.TOTAL.P | P-value from linear model across condition 1 |
| DIFF.TOTAL.Z | Test statistic for difference/interaction between conditions |
| DIFF.TOTAL.P | P-value for difference/interaction between conditions |

Enabling the `--combine` flag additionally estimates the (inverse-variance weighted) combined statistics from the allele-specific and total models, producing the following columns:

| Column | Description |
| --- | --- |
| ALL.COMBINE.Z | Test statisic from linear model across all samples |
| ALL.COMBINE.P | P-value from linear model across all samples |
| C0.COMBINE.Z | Test statisic from linear model across condition 0 |
| C0.COMBINE.P | P-value from linear model across condition 0 |
| C1.COMBINE.Z | Test statisic from linear model across condition 1 |
| C1.COMBINE.P | P-value from linear model across condition 1 |
| DIFF.COMBINE.Z | Test statistic for difference between conditions |
| DIFF.COMBINE.P | P-value for difference between conditions |

Enabling the `--binom` flag additionally runs a standard binomial test (and Fisher's test for interaction), producing the following columns:

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

All parameters can be displayed by running `--help`.

### Basic input/output

| Parameter | |
| --- | --- |
| `--input` | Path to input file |
| `--samples` | Path to sample identifier file, must have ID and CONDITION columns |
| `--peaks` | Path to file containing peak/gene boundaries, must contain CHR P0 P1 NAME CENTER columns |
| `--global_param` | Path to global overdispersion parameter file |
| `--local_param` | Path to local overdispersion parameter file |
| `--out` | Path to output |

### Inclusion/exclusion

| Parameter | |
| --- | --- |
| `--exclude` | The mimium distance between SNPs allowed in the haplotype (to exclude variants in the same read) |
| `--keep` | Path to file listing samples to retain for analyses |
| `--min_cov` | Individuals must have at least this many reads (for both alleles) to be tested |
| `--min_het` | Minimum minor heterozygous frequency for test SNP |
| `--min_maf` | Minimum minor allele frequency for test SNP |
| `--max_rho` | Maximum local/global over-dispersion parameter for which to include individual in test |
| `--window` | Window (in bp) for SNPs to test around the peak boundary |

### Total (non-allelic) input

| Parameter | |
| --- | --- |
| `--covar` | Path to covariates for total activity |
| `--total_matrix` | Path to matrix of total activity, enables the linear model. Note: if --local_param is on, then CNV is included as covariate |
| `--total_rn` | Rank normalize the total expression phenotype |

### Prediction

| Parameter | |
| --- | --- |
| `--predict` | Build predictive models for each feature (for downstream TWAS/RWAS/CWAS analysis) |
| `--predict_only` | Do not perform individual allele-specific tests |
| `--predict_snps` | Path to file listing SNPs for which predictor weights should be fit (typically common HapMap3 SNPs) |

### Testing

| Parameter | |
| --- | --- |
| `--bbreg` | Also perform a beta binomial regression with local CNV status as a covariate (must also provide `--local_param` file) |
| `--binom` | Also perform a standard binomial test |
| `--collapse_reads` | Merge all sites for each individual prior to analysis (not recommended) |
| `--combine` | Output combined BBINOM and TOTAL statistics by Stouffer's method, must have all relevant flags for basic and total input |
| `--mbased` | Also perform the MBASED test for differences (requires MBASED libraries) |

### Copy number

| Parameter | |
| --- | --- |
| `--cond_cnv_mean` | Use MU as the covariate for BBREG, otherwise uses CNV as the covariate (requires --bbreg) |
| `--fill_cnv` | Set individuals with missing CNV calls to diploid and rho=0.01 (assuming they are neutral) |
| `--mask_cnv_cutoff` | Mask out any sites that have an absolute CNV value above this cutoff (NA = no masking) |
| `--min_cnv_round` | Absolute CNV values below this cutoff get set to zero |

### Simulation

| Parameter | |
| --- | --- |
| `--sim` | Simulate imbalance based on the input data and test it. Allelic-fraction specified by --sim_af |
| `--sim_af1` | Specify the allelic fraction for CONDITION==1 |
| `--sim_af0` | Specify the allelic fraction for CONDITION==0 |
| `--sim_cnv` | Add local CNVs to simulations |
| `--sim_cnv_allelic` | CNVs always impact one allele |

### Extras

| Parameter | |
| --- | --- |
| `--indiv` | Also report the per-individual allele fractions (Warning, this can produce large files, typically used for visualization of one feature) |
| `--perm` | # of permutations to shuffle the allele labels (0=off) |
| `--perm_cond` | # of permutations to shuffle the condition labels (0=off) |
| `--seed` | Random seed (for simulations/permtuations) |
