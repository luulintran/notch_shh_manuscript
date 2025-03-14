# Run this after running 'analysis/01_deseq2_e16.R'

# SET UP
library(DESeq2)
library(org.Mm.eg.db)
library(tidyverse)
library(readr)
library(dplyr)



# LOAD DESEQ2 DDS OBJECT FROM RDS FILE: ----------------------------------------
dds <- readRDS("data/processed_data/rnaseq_e16/r_objects/deseq2_dds_e16.rds")

# STORE DESEQ2 RESULTS: --------------------------------------------------------
res <- results(dds)

# ANNOTATE RESULTS WITH GENE SYMBOLS AND ENTREZ IDS: ---------------------------
ensembl_ids <- rownames(res)

# annotate with gene symbols using org.Mm.eg.db package
res$symbol <- mapIds(org.Mm.eg.db,
                     keys=ensembl_ids,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
# annotate with entrez ids
res$entrez <- mapIds(org.Mm.eg.db,
                     keys=ensembl_ids,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

# Re-order res based on padj
resOrdered <- res[order(res$padj),]

# Move symbol and entrez columns to the front
resOrdered <- resOrdered[, c("symbol", 
                             "entrez", 
                             setdiff(colnames(resOrdered), 
                                     c("symbol", "entrez")))]

# Rename resOrdered to a simpler name as dataframe
out <- as.data.frame(resOrdered)

# Export CSV file with gene symbols and entrez ids
write.csv(out, file = "tables/table_S1_rnaseq_e16_deseq2_results.csv")

print("Table_S1 was generated and saved in tables/")
