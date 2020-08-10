#!/bin/sh

# ---
# Parameters for this script

# Sample identifier from the ASVCF file
# Example: TCGA-CZ-5468
ID=$1

# AS VCF prefix for each chromosome (assuming split by chromosome)
# {1-22}.vcf.gz will be appended to this to load the file
VCF=$2

# CNV file (in TCGA/FireCloud format)
# Example file: KIRC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.seg.txt (unzip if gz)
# Example download: http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/KIRC/20160128/gdac.broadinstitute.org_KIRC.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.Level_3.2016012800.0.0.tar.gz
CNV=$3

# Output file prefix
OUT=$4

# stratas directory
DIR="../"
# ---

# Get the counts from each chromosome
for CHR in `seq 1 22`; do
bcftools view --samples $ID --no-header ${VCF}${CHR}.vcf.gz \
| cut -f 1,2,10 | tr ':' '\t' | awk '{ print $1,$2,$3,$NF }' | tr ',' '\t' | tr ' ' '\t'
done | awk 'BEGIN{ print "CHR POS HAP REF.READS ALT.READS" } $3 != "0|0" && $3 != "1|1" && ($4+$5) > 0 { print $0 }' | tr ' ' '\t' > $OUT.counts

# Get the CNV segments
cat $CNV | grep $ID | awk 'BEGIN { print "ID CHR P0 P1 CNV" } { print $1,"chr"$2,$3,$4,$NF }' | tr ' ' '\t' > $OUT.cnv

# Compute CNVs for each region w/ more than 50 SNPs
Rscript $DIR/params.R --inp_counts $OUT.counts --inp_cnv $OUT.cnv --out $OUT.par_region --min_cov 5 --id $ID  --min_snps 50

# For any regions w/ fewer than 50 SNPs, rerun aggregating across regions
cat $OUT.par_region.local.params | awk 'NR == 1 || $6 == "NA"' | cut -f 1-5 > $OUT.rerun.cnv
Rscript $DIR/params.R --inp_counts $OUT.counts --inp_cnv $OUT.rerun.cnv --out $OUT.par_group --min_cov 5 --id $ID --group 10 --min_snps 50

# Merge the two outputs
mv $OUT.par_region.global.params $OUT.global.params
cat $OUT.par_region.local.params | grep -v NA \
| cat - $OUT.par_group.local.params | sort -k2,2 -k3,3n | uniq > $OUT.local.params

# Clean up
rm $OUT.par_region.* $OUT.par_group.* $OUT.rerun.cnv

