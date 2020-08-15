#!/bin/bash

## This code will add onto the haplotypes as the last column in the vcf a column (which is delimiated by a :) the read count of either halpotype
## ---

# Tool paths:
BCFTOOLS="/apps/bcftools-1.6/bin/bcftools"
BGZIP="/apps/htslib-1.6/bin/bgzip"
TABIX="/apps/htslib-1.6/bin/tabix"
# ---

# Input (vcf.gz) with genotypes for this individual
VCF=$1
# Input (.csv) with alleles for this individual
CSV=$2
# Output prefix
OUT=$3
# Sample identifier to put into the VCF header
SAMPLE=$4
# ---

set +o posix

# Construct the header
zcat $VCF | awk -vs=$SAMPLE 'BEGIN{ OFS="\t" } { if(substr($1,1,1) != "#") exit 1; if($1 == "#CHROM") { print "##FORMAT=<ID=AS,Number=2,Type=Integer,Description=\"Read counts for the ref and alt alleles\">"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,s; } else print $0; }' > $OUT.vcf

# Join the body
join -1 1 -2 1 \
<(zcat $VCF | awk 'substr($1,1,1) != "#" { printf $1":"$2" "; print $0 }' | sed 's/chr//' | sort -k1,1) \
<(cat $CSV | tail -n+2 | awk '{ printf $1":"$2" "; printf "%d %d\n",$6,$7 }' | sed 's/chr//' | sort -k1,1) \
| cut -d ' ' -f2- \
| awk 'BEGIN{ OFS="\t" } { if (NF != 10) print $1,$2,$3,$4,$5,$6,$7,$8,$9":AS",$10":"$11","$12; else print $1,$2,$3,$4,$5,$6,$7,$8,$9":AS",$10":0,0"; }' \
>> $OUT.vcf

# Sort and index
$BCFTOOLS sort $OUT.vcf | $BGZIP -c > $OUT.sorted.vcf.gz
$TABIX $OUT.sorted.vcf.gz
