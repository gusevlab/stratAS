# Allele-specific Quantification Workflow

These scripts provide an example workflow to QC and quantify allele-specific reads. Note these do not work out of the box as there are hard coded references to inputs/outputs, but should provide a general idea of the analysis steps.

## The input is:
* a FASTQ file containing your unmapped reads (or skip to step 1, line 32 if you have already done mapping)
* a phased VCF file containing all genotypes for the individual you are quantifying
* a [SAMPLE] identifier which will be used for output/intermediaries

## The output is:
* A VCF file containing allele-specific read counts after the WASP pipeline has been run to adjust for mapping bias.

## This workflow is run for each sample:

1. `wasp_pipeline.sh [FASTQ] [VCF] [SAMPLE]` takes bam files through the WASP pipeline.
2. `ase_read_counter.sh [VCF] [SAMPLE]` takes the bam file and counts the allele specific reads.
3. `creating_ase_vcf.sh [VCF] [SAMPLE]` takes the read counts and merges them with a .vcf file for haplotype phasing. 
4. `stratas_params_script.sh` provides a convient script to run the params.R script and should produce a single file which can be directly loaded into stratAS

After completion, all of the VCFs can be merged into a single VCF, all of the params estimates can be concatenated, and the files run through stratAS for allele-specific testing. 
