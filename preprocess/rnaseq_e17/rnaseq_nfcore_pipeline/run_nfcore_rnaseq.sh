#!/bin/bash  
#BSUB -J nfcore_rnaseq
#BSUB -o logs/nfcore.%J.out  
#BSUB -e logs/nfcore.%J.err  
#BSUB -n 12
#BSUB -R rusage[mem=50]  

mkdir -p logs
  
. /usr/share/Modules/init/bash
module load modules modules-init
module load java/18
module load singularity/3.9.2

 
~/bin/nextflow run nf-core/rnaseq -r 3.14.0 -c lsf.config -profile singularity -params-file nf-params.json