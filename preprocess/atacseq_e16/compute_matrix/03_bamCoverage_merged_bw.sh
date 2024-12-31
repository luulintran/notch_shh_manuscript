#!/bin/bash
#BSUB -J bamCoverage
#BSUB -o logs/bamCoverage_bw.%J.out
#BSUB -e logs/bamCoverage_bw.%J.err
#BSUB -n 12
#BSUB -R rusage[mem=50]

mkdir -p logs
mkdir -p merged_bigwig

. /usr/share/Modules/init/bash
module load modules modules-init
module load python

pip install deeptools

# generate BigWig files from merged BAM files
bamCoverage -b merged_filtered_bams/CTRL_filtered_merged.bam -o merged_bigwig/CTRL_merged.bw --binSize 10 --normalizeUsing RPKM
bamCoverage -b merged_filtered_bams/NICD_filtered_merged.bam -o merged_bigwig/NICD_merged.bw --binSize 10 --normalizeUsing RPKM