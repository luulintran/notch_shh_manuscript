#!/bin/bash
#BSUB -o logs/samtools.%J.out  
#BSUB -e logs/samtools.%J.err  
#BSUB -n 12
#BSUB -R rusage[mem=50]

mkdir -p logs
  
. /usr/share/Modules/init/bash
module load modules modules-init
module load samtools

# go through each filtered BAM file in the current directory
for file in *_filtered.bam; do
    # index BAM file
    samtools index "$file"
    echo "Indexed $file"
done