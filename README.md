# stratAS :first_quarter_moon:

Allele specific analyses across cell states and conditions

## Workflow

### Input data / dependencies

* A *.vcf file containing phased individual genotypes and an "AS" field listing the allelic counts.
* (Optional, for tumors) A list of structural variant boundaries for which to estimate local distribution parameters.
* R, and the `VGAM`, `optparse` libraries installed.

### Pre-processing and allelic counts

For QC, we recommend analyzing all sequence data with the [WASP](https://github.com/bmvdgeijn/WASP) mapping pipeline. For allelic quantification, we recommend using the GATK ASEReadCounter.

### Estimating prior parameters with `params.R`

`params.R` infers the read count distribution (beta binomial) parameters and needs to be run on a whole genome vcf one individual at a time. For tumor data, local CNV-specific parameters are also estimated.

### Analysis with `stratas.R`

`stratas.R` computes the actual AS statistics from a VCF of all individuals and the prior parameters estimated above.

### Output data
