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

`params.R` infers the read count distribution (beta binomial) parameters and needs to be run on a whole genome vcf one individual at a time. The input is an `--inp_counts` file which must contain the headers `CHR POS HAP REF.READS ALT.READS` where `HAP` is the `0|1` or `1|0` vcf haplotype code.

A vcf file is converted to counts as follows:
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

A `$OUT.global.params` file is generated containing the parameters and (optionally) an `$OUT.local.params` file is generated containing the positions and parameters for each CNV.

### Analysis with `stratas.R`

`stratas.R` computes the actual AS statistics from a VCF of all individuals and the prior parameters estimated above.

### Output data
