# SET UP
library(DESeq2)
library(tidyverse)
library(dplyr)
library(readr)


# READ IN DATA:------------------------------------------------------------------------------------------
# raw gene counts data from star_salmon; remove gene_name column
cts <- read.csv("data/processed_data/rnaseq_e17/gene_counts/rnaseq_e17_gene_counts.tsv", sep= "\t", row.names = 1) # makes the rownames the first col
cts <- as.data.frame(cts[, -c(1)]) # remove the gene_name col

# sample metadata
coldata <- read.csv("data/meta_data/rnaseq_e17/coldata.csv", row.names = 1)

# CREATE DESEQ2 DATASET:---------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = round(cts), 
                              colData = coldata, 
                              design = ~ condition)

dds$condition <- relevel(dds$condition, ref = "control")

# RUN DESEQ2 ANALYSIS: ----------------------------------------------------------------------------------
dds <- DESeq(dds)

# SAVE DESEQ2 OBJECT TO RDS FOR LATER: -------------------------------------------------------------------
saveRDS(dds, file = "data/processed_data/rnaseq_e17/r_objects/deseq2_dds_e17.rds")

print("DESeq2 analysis done. RDS file saved in r_objects/")