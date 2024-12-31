#!/bin/bash
#BSUB -J macs2_peak
#BSUB -o logs/macs2_peakcalling.%J.out
#BSUB -e logs/macs2_peakcalling.%J.err
#BSUB -n 12
#BSUB -R rusage[mem=50]

mkdir -p logs
mkdir -p macs2

. /usr/share/Modules/init/bash
module load modules modules-init
module load python

# install macs2
pip install macs2

#directories
#made new directory in results/ called "filtered_bams", and moved filtered .bam and .bai files here
filtered_bams_dir="filtered_bams"
output_dir="macs2"

#make output directory if it doesn't exist
mkdir -p "$output_dir"

# go through each bam file in the filtered_bams directory
for bam_file in "$filtered_bams_dir"/*.bam; do
    base_name=$(basename "$bam_file" .bam)

    # Perform MACS2 peak calling
    macs2 callpeak \
        -t "$bam_file" \
        -f BAMPE \
        --keep-dup all \
        --call-summits \
        -g mm \
        -n "${base_name}" \
        --outdir "$output_dir"

    echo "MACS2 peak calling completed for $bam_file"
done