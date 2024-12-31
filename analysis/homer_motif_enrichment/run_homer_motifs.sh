#!/bin/bash  
#BSUB -J homer
#BSUB -o logs/homer.findMotifs.%J.out  
#BSUB -e logs/homer.findMotifs.%J.err  
#BSUB -n 12
#BSUB -R rusage[mem=50]  

mkdir -p logs
  
. /usr/share/Modules/init/bash
module load modules modules-init
module load homer

#Paths to BED files and genome FASTA file
BED1="diffbind_output/CTRL_e16_enriched.bed"
BED2="diffbind_output/NICD_e16_enriched.bed"
GENOME="../genome_file/mm10.fa"

#Output directory
OUTPUT_DIR="homer_output"

#make output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

#Run HOMER findMotifsGenome for the first and second BED file

findMotifsGenome.pl "$BED1" "$GENOME" "$OUTPUT_DIR/CTRL_motifs"

findMotifsGenome.pl "$BED2" "$GENOME" "$OUTPUT_DIR/NICD_motifs"
