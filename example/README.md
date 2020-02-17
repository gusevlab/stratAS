### `params.R` : estimating prior parameters

```
Rscript params.R \
--inp_counts TCGA-CZ-5468.TUMOR.AS.counts.gz \
--inp_cnv TCGA-CZ-5468.TUMOR.AS.counts.cnv \
--out TCGA-CZ-5468.TUMOR.AS.param \
--group 10 \
--id TCGA-CZ-5468.TUMOR
```
### `stratas.R` : testing for differential AS

```
Rscript stratas.R \
--input example/ENSG00000075240.12.mat \
--samples example/KIRC.ALL.AS.PHE \
--peaks example/ENSG00000075240.12.bed \
--global_param example/KIRC.ALL.AS.CNV \
--local_param=example/KIRC.ALL.AS.CNVLOCAL
```

