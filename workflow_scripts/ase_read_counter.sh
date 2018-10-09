# Remember that when running ReorderSAM that there whouls be a .fasta files and a corresponding .dict file, both with the same stem 
# Remember that when running samtools sort, there should be an -o tag for the output file, otherwise it will fail to be outputted
AHS="/bcb/agusevlab/ashetty"
VCF=$1
SAMPLE=$2

/bcb/agusevlab/miniconda2/bin/java -jar /apps/picard-tools-2.2.2/dist/picard.jar AddOrReplaceReadGroups \
I=${AHS}/step_7_output_files/${SAMPLE}.sorted.bam \
O=${AHS}/ase_read_counter_output/${SAMPLE}.red.sorted.bam \
RGLB=lib \
RGPL=illumina \
RGPU=run \
RGSM=${SAMPLE}

samtools index $AHS/ase_read_counter_output/${SAMPLE}.red.sorted.bam

/opt/jre1.7.0_60/bin/java -Djava.io.tmpdir=/tmp -jar /apps/gatk3.5/protected/gatk-package-distribution/target/gatk-package-distribution-3.5.jar \
   -R $AHS/new_reference_files/ucsc.hg19.fasta \
   -T ASEReadCounter \
   -o $AHS/allele_counts/${SAMPLE}.csv \
   -I $AHS/ase_read_counter_output/${SAMPLE}.red.sorted.bam \
   -sites $AHS/all_chrs_${VCF}.recode.vcf \
   -U ALLOW_N_CIGAR_READS \
   -minDepth 10

