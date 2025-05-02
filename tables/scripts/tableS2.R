# Run 'analysis/01_deseq2_e16.R' if you haven't already
#source("analysis/01_deseq2_e16.R")

# DEFINE FILES AND PATHS: ------------------------------------------------------
# input deseq2 results rds file
rds_deseq2_results <- "data/processed_data/rnaseq_e16/r_objects/deseq2_dds_e16.rds"

# output directory for supplementary tables
output_dir_tables <- "tables"

# output filename for table
filename <- "table_S2_e16_shh_genes.csv"

# MAKE LIST OF GENES RELATED TO SHH PATHWAY: -----------------------------------
Shh_genes <- c("Smo", "Boc", "Cdon", "Gas1", "Gli1", "Gli2", "Gli3", "Sufu", 
               "Disp1", "Iqce", "Efcab7", "Ptch1", "Ptch2", "Hhip", "Sufu", 
               "Gsk3b", "Ck1", "Pcaf", "Cul3", "Kif7")

# LOAD DESEQ2 DDS OBJECT FROM RDS FILE: ----------------------------------------
dds <- readRDS(rds_deseq2_results)

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

# FILTER DESEQ2 RESULTS FOR SHH GENES: -----------------------------------------

shh_deseq2_results <- out %>%
  filter(symbol %in% Shh_genes)


# Make rownames a column and name ensembl because it contains ensembl id's
shh_deseq2_results <- rownames_to_column(shh_deseq2_results, var = "ensembl")

# SAVE TO FILE: ----------------------------------------------------------------
write.csv(
  shh_deseq2_results, 
  file = file.path(output_dir_tables, filename), 
  row.names = F)

print(paste0(filename, " was generated and saved in ", output_dir_tables))