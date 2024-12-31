#!/bin/bash
#BSUB -J merged_bed
#BSUB -o logs/narrowPeak_to_bed.%J.out
#BSUB -e logs/narrowPeak_to_bed.%J.err
#BSUB -n 12
#BSUB -R rusage[mem=50]

mkdir -p logs

#directories
merged_macs2_dir="merged_filtered_macs2"
output_dir="merged_filtered_bed"

#make output directory if it doesn't exist
mkdir -p "$output_dir"


# go through each peak file in the directory
for peak_file in "$merged_macs2_dir"/*.narrowPeak; do
    base_name=$(basename "$peak_file" .narrowPeak)

    # take out first three columns of narrowPeak file and make a new file
    #awk prints columns and separates them by tabs \t. Output Field Separator
    #print first three columns of current line and put into new file
    awk 'BEGIN{OFS="\t"}{print $1, $2, $3}' "$peak_file" > "$output_dir/${base_name}.bed"
done

