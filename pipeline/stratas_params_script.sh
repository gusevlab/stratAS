AHS="/bcb/agusevlab/ashetty"
SAMPLE=$1
echo started

zcat $AHS/stratas_prep_files/${SAMPLE}_ase.sorted.named.vcf.gz \
| grep -v '#' \
| cut -f 1,2,10 \
| tr ':' '\t' \
| awk '$3 != "1|1" && $3 != "0|0"' \
| cut -f 1-3,7 \
| tr ',' '\t' \
| sed 's/chr//' \
| awk 'BEGIN { print "CHR POS HAP REF.READS ALT.READS" } $4 + $5 > 0' \
 > $AHS/stratas_prep_files/${SAMPLE}.counts

/apps/R-3.4.1share/bin/Rscript /bcb/agusevlab/stratAS/params.R \
--min_cov 5 \
--inp_counts /bcb/agusevlab/ashetty/stratas_prep_files/${SAMPLE}.counts \
--out /bcb/agusevlab/ashetty/stratas_input_files/${SAMPLE}

set +o posix

paste <(printf "%s\n%s" ID "${SAMPLE}") <(cat /bcb/agusevlab/ashetty/stratas_input_files/"${SAMPLE}".global.params) > /bcb/agusevlab/ashetty/stratas_input_files/$SAMPLE.modified.global.params
