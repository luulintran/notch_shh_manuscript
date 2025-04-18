# Run 'analysis/04_deseq2_e13-vs-e17.R' if you haven't already
#source("analysis/04_deseq2_e13-vs-e17.R")

# DEFINE FILES AND PATHS: ------------------------------------------------------
# input deseq2 results rds file
rds_deseq2_results <- "data/processed_data/rnaseq_ctrl_e13/r_objects/deseq2_dds_e13-vs-e17.rds"

# output directory
output_dir_tables <- "tables"

# output filename
filename <- "table_S8_rnaseq_e13-vs-e17_deseq2_results.csv"

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

# Rename resOrdered to a simpler name as dataframe
out <- as.data.frame(resOrdered)

# Export CSV file with gene symbols and entrez ids
write.csv(out, 
          file = file.path(
            output_dir_tables, 
            filename)
          )

print(paste0(filename, " was generated and saved in ", output_dir_tables))

