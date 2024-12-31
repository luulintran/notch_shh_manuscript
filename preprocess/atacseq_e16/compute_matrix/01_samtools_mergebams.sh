#!/bin/bash
#BSUB -J samtools_merge
#BSUB -o logs/samtools_merge.%J.out
#BSUB -e logs/samtools_merge.%J.err
#BSUB -n 12
#BSUB -R rusage[mem=50]

#make direcotories if they don't exist
mkdir -p logs
mkdir -p merged_filtered_bams

. /usr/share/Modules/init/bash
module load modules modules-init
module load samtools

# change to filtered_bams directory. Let me know if this directory doesn't exist.
cd filtered_bams || { echo "Directory not found: results/filtered_bams"; exit 1; }

# make arrays with the file names for each set of replicates that you want to merge together
ctrl_replicates=("filtered_bams/CTRL_REP1.mLb.clN.sorted_filtered.bam" "filtered_bams/CTRL_REP2.mLb.clN.sorted_filtered.bam" "filtered_bams/CTRL_REP3.mLb.clN.sorted_filtered.bam")
nicd_replicates=("filtered_bams/NICD_REP1.mLb.clN.sorted_filtered.bam" "filtered_bams/NICD_REP2.mLb.clN.sorted_filtered.bam" "filtered_bams/NICD_REP3.mLb.clN.sorted_filtered.bam")

# merge BAM files for control replicates
samtools merge merged_filtered_bams/CTRL_filtered_merged.bam "${ctrl_replicates[@]}"

# merge BAM files for NICD replicates
samtools merge merged_filtered_bams/NICD_filtered_merged.bam "${nicd_replicates[@]}"