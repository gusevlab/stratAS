# stratAS :first_quarter_moon:

Allele specific analyses across cell states and conditions

## Workflow

### Input data / dependencies

* A *.vcf file containing phased individual genotypes and an "AS" field listing the allelic counts.
* (Optional, for tumors) A list of structural variant boundaries for which to estimate local distribution parameters.
* R, and the `VGAM`, `optparse` libraries installed.

### Pre-processing and allelic counts

### Estimating prior parameters with `params.R`

Run on each full vcf to get global and local distributional parameters.

### Analysis with `stratas.R`

Run on combined vcf to get AS statistics.

### Output data
