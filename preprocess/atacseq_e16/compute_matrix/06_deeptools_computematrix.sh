#!/bin/bash
#BSUB -J deeptools
#BSUB -o logs/computematrix_heatmap.%J.out
#BSUB -e logs/computematrix_heatmap.%J.err
#BSUB -n 12
#BSUB -R rusage[mem=50]

mkdir -p logs

#load modules and install deeptools
. /usr/share/Modules/init/bash
module load modules modules-init
module load python

# install macs2
pip install deeptools

#directories
filtered_bigwig_dir="merged_bigwig"
filtered_bed_dir="merged_filtered_bed"
output_dir1="deeptools_output"
output_dir2="$output_dir1/computematrix"

#make output directory if it doesn't exist
mkdir -p "$output_dir1"
mkdir -p "$output_dir2"

#array of bigWig and bed files
declare -A files
files=( ["CTRL"]="CTRL_merged.bw CTRL_filtered_merged_peaks.bed"
        ["NICD"]="NICD_merged.bw NICD_filtered_merged_peaks.bed" )

#process bigWig and BED files into matrix, heatmap, and plotprofile
for key in "${!files[@]}"; do
    IFS=' ' read -r bigwig_file bed_file <<< "${files[$key]}"
    bigwig_file="$filtered_bigwig_dir/$bigwig_file"
    bed_file="$filtered_bed_dir/$bed_file"
    base_name="$key"

    if [[ -f "$bigwig_file" && -f "$bed_file" ]]; then
        # Compute matrix
        computeMatrix reference-point \
            -S "$bigwig_file" \
            -R "$bed_file" \
            --referencePoint center \
            -b 1000 -a 1000 \
            -bs 1 \
            -p max \
            --outFileName "$output_dir2/${base_name}.matrix.gz"
    else
        echo "Files for $key not found, skipping..."
    fi
done

