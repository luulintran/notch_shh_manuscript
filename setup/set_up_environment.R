# setup_environment.R

# Load renv for package management
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
library(renv)

# Restore the environment from the renv.lock file
if (file.exists("renv.lock")) {
  cat("\nRestoring environment from renv.lock...\n")
  renv::restore()
} else {
  stop("renv.lock file not found. Make sure it exists in the project directory.")
}

# Set up directory structure
required_directories <- c(
  "data/raw_data",
  "data/processed_data",
  "data/meta_data",
  "data/motif_files",
  "data/processed_data/atacseq_e16",
  "data/processed_data/rnaseq_e16",
  "data/processed_data/rnaseq_e17",
  "data/processed_data/rnaseq_ctrl_e13",
  "data/processed_data/atacseq_e16/deeptools_output",
  "data/processed_data/atacseq_e16/diffbind_results_bed",
  "data/processed_data/atacseq_e16/filtered_bams",
  "data/processed_data/atacseq_e16/filtered_macs2",
  "data/processed_data/atacseq_e16/r_objects",
  "data/processed_data/homer_output",
  "data/processed_data/rnaseq_e16/raw_counts",
  "data/processed_data/rnaseq_e16/r_objects",
  "data/processed_data/rnaseq_e17/raw_counts",
  "data/processed_data/rnaseq_e17/r_objects/",
  "data/processed_data/rnaseq_ctrl_e13/raw_counts",
  "data/processed_data/rnaseq_ctrl_e13/r_objects/",
  "analysis",
  "figures",
  "preprocess",
  "setup",
  "tables"
)

for (dir in required_directories) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat(paste0("Created directory: ", dir, "\n"))
  }
}

# Load libraries
library(tidyverse)
library(dplyr)
library(readr)
library(ggplot2)
library(extrafont)
library(pheatmap)
library(data.table)
library(reshape2)
library(car)

library(DESeq2)
library(DiffBind)
library(ChIPseeker)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(RColorBrewer)

# import Arial font
font_import(pattern = "Arial", prompt = FALSE)

cat("\nEnvironment setup complete.\n")