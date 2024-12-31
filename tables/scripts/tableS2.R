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

# Save as dataframe
out <- as.data.frame(resOrdered)

# MAKE LIST OF GENES RELATED TO SHH PATHWAY: -----------------------------------
Shh_genes <- c("Smo", "Boc", "Cdon", "Gas1", "Gli1", "Gli2", "Gli3", "Sufu", 
               "Disp1", "Iqce", "Efcab7", "Ptch1", "Ptch2", "Hhip", "Sufu", 
               "Gsk3b", "Ck1", "Pcaf", "Cul3", "Kif7")

# FILTER DESEQ2 RESULTS FOR SHH GENES: -----------------------------------------

shh_deseq2_results <- out %>%
  filter(symbol %in% Shh_genes)


# Make rownames a column and name ensembl because it contains ensembl id's
shh_deseq2_results <- rownames_to_column(shh_deseq2_results, var = "ensembl")

# SAVE TO FILE: ----------------------------------------------------------------
write.csv(
  shh_deseq2_results, 
  "tables/table_S2_e16_shh_genes.csv", 
  row.names = F)

