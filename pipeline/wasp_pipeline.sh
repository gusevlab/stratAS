# This is the WASP pipeline script
# Python 2.7, pytables 2.x, pysam and numpy are all required depdendencies for this - install them!
# If my minoconda installation in the lab directory is safe, then these should all be installed already
# This currently requires a folder for the reference fasta - if this does not exist, then this will need to be downloaded
# step{2..7}_output_files/ folders also are needed - hopefully this can be implemented so that we will have only one folder for WASP  

#!/bin/sh

# request Bourne shell as shell for job
#$ -S /bin/sh

AHS="/bcb/agusevlab/ashetty"
WASP="/bcb/agusevlab/WASP"
DATA="/bcb/agusevlab/DATA/PRCA_CHIP"

FASTQ=$1
VCF=$2
SAMPLE=$3

echo step 2 has started

# STEP 2
/apps/bwa-0.7.13/bwa aln -t 8 \
$AHS/new_reference_files/ucsc.hg19.fasta \
$AHS/input_fastq_files/${FASTQ} > $AHS/step_2_output_files/${SAMPLE}.sai

/apps/bwa-0.7.13/bwa samse \
$AHS/new_reference_files/ucsc.hg19.fasta \
$AHS/step_2_output_files/${SAMPLE}.sai \
$AHS/input_fastq_files/${FASTQ} | samtools view -b -q 10 - > $AHS/step_2_output_files/${SAMPLE}.bam

samtools sort -T $AHS/step_2_output_files/sort_temp/${SAMPLE} -O bam $AHS/step_2_output_files/${SAMPLE}.bam > $AHS/step_2_output_files/${SAMPLE}.sorted.bam

samtools index $AHS/step_2_output_files/${SAMPLE}.sorted.bam 

echo step 3 has started
# STEP 3
echo ${VCF} >  $AHS/step_3_output_files/${VCF}.txt

/bcb/agusevlab/miniconda2/bin/python2.7 $WASP/mapping/find_intersecting_snps.py \
	--is_sorted \
	--output_dir $AHS/step_3_output_files \
	--snp_tab $AHS/hdf5_files/snp_tab.h5 \
	--snp_index $AHS/hdf5_files/snp_index.h5 \
	--haplotype $AHS/hdf5_files/haps.h5 \
	--samples $AHS/step_3_output_files/${VCF}.txt \
	$AHS/step_2_output_files/${SAMPLE}.sorted.bam

echo step 4 has started
# STEP 4
/apps/bwa-0.7.13/bwa aln -t 8 \
$AHS/new_reference_files/ucsc.hg19.fasta \
$AHS/step_3_output_files/${SAMPLE}.sorted.remap.fq.gz > $AHS/step_4_output_files/${SAMPLE}.sai

/apps/bwa-0.7.13/bwa samse \
$AHS/new_reference_files/ucsc.hg19.fasta \
$AHS/step_4_output_files/${SAMPLE}.sai \
$AHS/step_3_output_files/${SAMPLE}.sorted.remap.fq.gz | samtools view -b -q 10 - > ${AHS}/step_4_output_files/${SAMPLE}.bam

samtools sort -T $AHS/step_4_output_files/sort_temp/${SAMPLE} -O bam $AHS/step_4_output_files/${SAMPLE}.bam > $AHS/step_4_output_files/${SAMPLE}.sorted.bam

samtools index $AHS/step_4_output_files/${SAMPLE}.sorted.bam 

echo step 5 has started
# STEP 5
/bcb/agusevlab/miniconda2/bin/python2.7 $WASP/mapping/filter_remapped_reads.py \
	$AHS/step_3_output_files/${SAMPLE}.sorted.to.remap.bam \
	$AHS/step_4_output_files/${SAMPLE}.sorted.bam \
	$AHS/step_5_output_files/${SAMPLE}.keep.bam

echo step 6 has started
# STEP 6
samtools merge $AHS/step_6_output_files/${SAMPLE}.keep.merge.bam \
              $AHS/step_5_output_files/${SAMPLE}.keep.bam  \
              $AHS/step_3_output_files/${SAMPLE}.sorted.keep.bam

echo step 7 has started
# STEP 7
# Be aware the step_6 file is not truly sortted!
/bcb/agusevlab/miniconda2/bin/python2.7 $WASP/mapping/rmdup.py $AHS/step_6_output_files/${SAMPLE}.keep.merge.bam $AHS/step_7_output_files/${SAMPLE}.sorted.bam

echo completed!
