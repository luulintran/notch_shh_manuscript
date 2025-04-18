# READ DDS OBJECT FOLLOWING DESEQ2 ANALYSIS FROM RDS FILE: ---------------------
dds <- readRDS(rds_deseq2_results)
# Store results in res
res <- results(dds)

# ANNOTATE DESEQ2 RESULTS WITH GENE SYMBOLS AND ENTREZ IDS: --------------------
ensembl_ids <- rownames(res)

# annotate with gene symobols using org.Mm.eg.db package
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
                         
# PREPARE DATA FOR HEATMAP SHOWING SPECIFIC GENE SET: --------------------------
# Remove NAs from padj before filtering for significant genes with padj < 0.05
sig_res <- resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.05, ]

# Extract rownames/ensembl id's of significant DEGs
sig_genes <- rownames(sig_res)

# Map ensembl id's to symbols
sig_symbols <- mapIds(org.Mm.eg.db,
                      keys = sig_genes,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")

# Filter sig_genes list for genes in combined_gene list by symbol
selected_genes <- sig_genes[sig_symbols %in% combined_gene_list]

                         

# EXTRACT TRANSFORMED VALUES: --------------------------------------------------
vsd <- vst(dds, blind=FALSE)

# Subset assay matrix for the selected genes
mat <- assay(vsd)[selected_genes, ]

# Update row names with gene symbols for the selected genes
rownames(mat) <- sig_symbols[sig_symbols %in% combined_gene_list]

# Center the data
mat <- mat - rowMeans(mat)

# Make Heatmap in specific order ***********************************************
# Only include genes present in the matrix
ordered_genes <- desired_order[desired_order %in% rownames(mat)]

# Subset the matrix to include only the ordered genes
mat_ordered <- mat[match(ordered_genes, rownames(mat)), ]

# MAKE HEATMAP WITH ORDERED GENES: ---------------------------------------------
ordered_combined_gene_list_heatmap <- pheatmap(
  mat_ordered,
  # Set to FALSE to use specified order
  cluster_rows = FALSE,  
  cluster_cols = TRUE, 
  # Show gene symbols
  show_rownames = TRUE,  
  annotation_col = as.data.frame(
    colData(vsd)[, "condition", drop=FALSE]),
  color = colorRampPalette(c(control_color, "white", mutant_color))(50))

# SAVE FIGURE: -----------------------------------------------------------------
pdf(file.path(output_dir_figures, filename_heatmap), width = 6, height = 8)
print(ordered_combined_gene_list_heatmap)

dev.off()

print(paste0(filename_heatmap, " generated and saved in ", output_dir_figures))