# RUN DESEQ2 ANALYSIS SCRIPT: --------------------------------------------------
# Run 'analysis/01_deseq2_e16.R' if you have not already
# source("analysis/01_deseq2_e16.R")

# DEFINE FILES AND PATHS: ------------------------------------------------------
# input deseq2 results rds file
rds_deseq2_results <- "data/processed_data/rnaseq_e16/r_objects/deseq2_dds_e16.rds"

# output directory for supplementary tables
output_dir_tables <- "tables"

# output filename for table
filename <- "table_S1_rnaseq_e16_deseq2_results.csv"

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
                                     c("symbol", "entrez")
                                     )
                             )
                         ]

# Rename resOrdered to a simpler name as dataframe
out <- as.data.frame(resOrdered)

# Export CSV file with gene symbols and entrez ids
write.csv(out, 
          file = file.path(
            output_dir_tables, 
            filename)
          )

print(paste0(filename, " was generated and saved in ", output_dir_tables))

      