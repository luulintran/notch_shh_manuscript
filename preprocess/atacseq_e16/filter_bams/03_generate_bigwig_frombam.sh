#!/bin/bash
#BSUB -o logs/bamCoverage.%J.out
#BSUB -e logs/bamCoverage.%J.err
#BSUB -n 12
#BSUB -R rusage[mem=50]

mkdir -p logs

. /usr/share/Modules/init/bash
module load modules modules-init
module load python

# install deepTools
pip install deeptools

# make arrays with the file names for each replicate for each condition, ctrl and nicd. Use filtered bams
ctrl_replicates=("CTRL_REP1.mLb.clN.sorted_filtered.bam" "CTRL_REP2.mLb.clN.sorted_filtered.bam" "CTRL_REP3.mLb.clN.sorted_filtered.bam")
nicd_replicates=("NICD_REP1.mLb.clN.sorted_filtered.bam" "NICD_REP2.mLb.clN.sorted_filtered.bam" "NICD_REP3.mLb.clN.sorted_filtered.bam")

# define function for generating BigWig files with bamCoverage from filtered bam files
generate_bigwig() {
    local replicate_files=("${!1}")
    local extension=$2

# iterate through each bam file in replicate_files
    for bam_file in "${replicate_files[@]}"; do
        # output BigWig file name
        output_bigwig="${bam_file%.bam}${extension}"  # replace .bam from input BAM file and add extension

        # generate BigWig file
        bamCoverage -b "$bam_file" -o "$output_bigwig" --binSize 10 --normalizeUsing RPKM

        echo "Generated BigWig file: $output_bigwig"
    done #end loop
}

# generate BigWig files for control replicates. 1st arg= ctrl_replicates. 2nd arg = extension .bw
generate_bigwig ctrl_replicates[@] ".bw"

# generate BigWig files for NICD replicates
generate_bigwig nicd_replicates[@] ".bw"