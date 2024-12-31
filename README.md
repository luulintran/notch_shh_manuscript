## Introduction
This repository contains the code and analyses associated with the manuscript titled "Epigenetic priming of neural progenitors by Notch enhances Sonic hedgehog signaling and establishes gliogenic competence". It includes scripts for setting up the project environment, performing differential analysis, generating figures, and creating tables. The code is organized to ensure a clear workflow, where figures depend on both analysis and tables.

Raw and processed data are in the process of being deposited in the Gene Expression Omnibus. The accession number can be found here once we receive it. 

To reproduce the full analysis pipeline used in the manuscript or to reproduce any of the figures, follow the steps below.

---

## Overview
For our manuscript, we generated data using RNA-seq and ATAC-seq. To preprocess the data from each assay, we used the Nf-core RNA-seq (https://nf-co.re/rnaseq/3.14.0/) and ATAC-seq (https://nf-co.re/atacseq/2.1.2/) pipelines. Then we performed downstream analysis using the processed data. For ATAC-seq, we used Diffbind (v 3.12.0) (using DESeq2 for differential analysis) and ChIPSeeker (v 1.38.0) R packages. For RNA-seq, we used DESeq2 (v 1.42.1), pheatmap (v 1.012) and clusterProfiler (v 4.10.10) R packages. This repo is meant for reproducing the figures in the manuscript starting with processed data. Information and scripts related to preprocessing the raw data can be found below.

---

## Table of Contents

- [Setup](#setup)
- [Preprocessing](#preprocessing)
- [Data](#data)
- [Analysis](#analysis)
- [Figures](#figures)
- [Tables](#tables)


---

## Setup
The majority of this project was done using R (v 4.3.2) and some shell scripting. 

Before running any analysis or generating figures, you need to set up the environment using the `setup/` scripts. This will ensure all dependencies are installed.

### To run the setup:

1. Clone the repository:
   ```bash
   git clone https://github.com/luulintran/notch_shh_manuscript.git
   cd notch_shh_manuscript
   ```

2. Run the setup script to install the required dependencies:
    ```
    Rscript setup/set_up_environment.R
    ```

## Preprocessing

This repo is meant to regenerate figures starting with processed data. However, if you are interested in the preprocessing steps, these scripts can be found in the `preprocess/` directory. The `preprocess/` directory contains scripts to run the Nf-core pipelines for preprocessing RNA-seq and ATAC-seq data, as well as filtering fragments for ATAC-seq.

Note: You will need to include a config file to run the Nf-core pipeline scripts. 

## Data

The `data/` directory contains both raw and processed data, which can be downloaded once they have been deposited in GEO.

Note: After running the preprocessing steps and filtering ATAC-seq fragments, move the processed data to `data/processed_data/` for downstream analysis.

## Analysis

The `analysis/` directory contains scripts for performing downstream analysis using DESeq2 or Diffbind depending on the specific figure you plan to generate. You need to run the analysis scripts before generating any figures or tables, as the results from these analyses will be used in both.

If you are interested in regenerating all the figures and tables in the manuscript, run each script in the order they are numbered. If you are only interested in regenerating the figures and tables related to one assay (such as RNA-seq), you only need to run the corresponding `analysis/` script.

### HOMER motif enrichment analysis
Included in the `analysis/` directory are shell scripts for performing motif enrichment analysis with HOMER (v 4.9) (http://homer.ucsd.edu/homer/). These are separate from the rest of the workflow, but the .txt files resulting from HOMER `findmotifsGenome.pl` to locate motif instances are needed to run `tables/tableS4-5.R`.  To run findmotifsGenome.pl to search for motif instances, you need motif files for the transcription factor you are interested in. You can read more about this on HOMER's website. 

Some of the results of the motif enrichment analysis can be found in Figure 4. The rest of the results are not included in this repo or in the manuscript, but can be reproduced using the BED files generated from `analysis/04_make_bedfiles.R` and the HOMER script `analysis/homer_motif_enrichment/run_homer_motifs.sh`, which will need a genome file.

## Figures

The `figures/` directory contains scripts for generating the figure panels included in the manuscript. Some figure panels depend on tables being generated first. This requirement should be noted at the top of the script. Be sure to follow the necessary steps in the Tables section before running certain figure scripts.

The name of scripts will correspond to the figure panel. Each figure directory (such as `fig2/`) will have numbered scripts starting with '01_....R'. Run these scripts in order to regenerate plots for the entire figure. 

Note: To regenerate the TSS enrichment plot in Figure 3 (`04_fig3c_tssenrichmentplot.R`), you will need to generate matrix files. Shell scripts that were used to make matrix files can be found in `preprocess/compute_matrix` and are numbered in the order that they should be run. Briefly, these will merge the bam files for each condition (to get average across replicates), index them, make bigwig files, call peaks, and make bed files. The bigwig and bed files are input for deeptools `computeMatrix` command.

## Tables

The `tables/scripts/` directory contains scripts for generating the supplemental tables in the manuscript. These tables are based on the results of the differential analysis. Some figures rely on tables to be created first, so make sure the analysis is completed before running these scripts. Generated tables will be saved in `tables/`.

Note: Tables 4 and 5 require .txt files resulting from Homer `findMotifsGenome.pl` which was used to locate instances of motifs (scripts can be found in `analysis/homer/`). Make sure these .txt files are in `processed_data/atacseq_e16/homer_output/` to run the scripts.

Note: The csv files resulting from `tables/scripts/TableS4-5.R` were combined to make the xlxs files for Supplemental Tables 4 and 5. For example, the csv files for table_S4 were combined into one xlxs file, 'Table_S4_neurogenic_motifs.xlxs' with each csv file representing a tab, and the same for 'Table_S5_gliogenic_motifs.xlxs'.

## Notes
Ensure you have the appropriate data files in `data/raw_data/` and `data/processed_data/`.

The scripts in `analysis/`, `figures/`, and `tables/` depend on the output from previous steps, so follow the order of operations carefully and look at the comments at the top of the scripts or notes in the README.
