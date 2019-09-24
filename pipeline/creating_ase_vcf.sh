#!/bin/bash

AHS="/bcb/agusevlab/ashetty"
VCF=$1
SAMPLE=$2
# This code will add onto the haplotypes as the last column in the vcf a column (which is delimiated by a :) the read count of either halpotype

echo "a"

set +o posix

join -a1 -1 1 -2 1 \
<(cat "$AHS"/all_chrs_"$VCF".recode.vcf | awk '{ printf $1; printf $2; printf " "; print $0 }' | sed 's/chr//' | awk '{ printf "%015d", $1; print " ", $2, $3, $4, $5, $6, $7, $8, $9, $10, $11 }' | sort) \
<(tail -n+2 "$AHS"/allele_counts/"$SAMPLE".csv | awk '{ printf $1; printf $2; printf " "; printf "%d %d\n",$6,$7 }' | sed 's/chr//' | awk '{printf "%015d", $1; print " ", $2, $3}' | sort) \
| cut -d ' ' -f2- \
| awk '{ if(substr($1,1,1) != "#" ) { if (NF != 10) print $1,$2,$3,$4,$5,$6,$7,$8,$9":AS",$10":"$11","$12; else print $1,$2,$3,$4,$5,$6,$7,$8,$9":AS",$10":0,0"; } else { print $0 } }' \
| tr ' ' '\t' > "$AHS"/stratas_prep_files/"${SAMPLE}"_ase.vcf

echo b

cat \
<(echo '##fileformat=VCFv4.2') \
<(cat "$AHS"/stratas_prep_files/"$SAMPLE"_ase.vcf | grep "##" | grep -v "##bcftools") \
<(echo '##FORMAT=<ID=AS,Number=2,Type=Integer,Description="Read counts for the ref and alt alleles">') \
<(cat $AHS/stratas_prep_files/${SAMPLE}_ase.vcf | grep -v "##") \
| sed "/^#CHROM/s/${VCF}/${SAMPLE}/g" > $AHS/stratas_prep_files/${SAMPLE}_ase.named.vcf

echo c

/apps/bcftools-1.6/bin/bcftools sort $AHS/stratas_prep_files/${SAMPLE}_ase.named.vcf > $AHS/stratas_prep_files/${SAMPLE}_ase.sorted.named.vcf

echo d

/apps/htslib-1.6/bin/bgzip $AHS/stratas_prep_files/${SAMPLE}_ase.sorted.named.vcf
/apps/htslib-1.6/bin/tabix $AHS/stratas_prep_files/${SAMPLE}_ase.sorted.named.vcf.gz

echo e


