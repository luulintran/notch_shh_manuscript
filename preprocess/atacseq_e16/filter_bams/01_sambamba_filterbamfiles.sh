#!/bin/bash  
#BSUB -o logs/sambamba.%J.out  
#BSUB -e logs/sambamba.%J.err  
#BSUB -n 12
#BSUB -R rusage[mem=50]  

mkdir -p logs
  
. /usr/share/Modules/init/bash
module load modules modules-init
module load sambamba

#go through each .bam file in the current directory
for file in *.bam; do
    #make new output file name with _filtered.bam
    output_file="${file%.bam}_filtered.bam"
    
    #filter fragment lengths less than 100 bp
    sambamba view -f bam -F 'template_length < 100 and template_length > -100' -t 8 "$file" > "$output_file"
    
    echo "Filtered $file to $output_file"
done
