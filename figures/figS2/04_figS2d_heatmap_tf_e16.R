# Run after running 'analysis/01_deseq2_e16.R' 

library(pheatmap)
library(DESeq2)
library(org.Mm.eg.db)
library(tidyverse)
library(readr)
library(dplyr)

# READ DDS OBJECT FOLLOWING DESEQ2 ANALYSIS FROM RDS FILE: --------------------
dds <- readRDS("data/processed_data/rnaseq_e16/r_objects/deseq2_dds_e16.rds")
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
                                     c("symbol", "entrez")))
]

# MAKE LIST OF TRANSCRIPTION FACTORS: ------------------------------------------
# These TFs have motifs enriched in sites with decreased or 
# increased accessibility

# Heatmap with specific genes of TFs from atac-seq motif enrichment 
TFs <- c('Lhx1','Lhx2', 'Nfix', 'Nfia', 'Nfic', 'Sox2', 'Sox8', 'Sox9','Sox10', 
         'Neurog1', 'Neurog2', 'Neurod2', 'Mef2c', 'Mef2d', 'Eomes', 'Tbr1')

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
selected_genes <- sig_genes[sig_symbols %in% TFs]

# EXTRACT TRANSFORMED VALUES: **************************************************
vsd <- vst(dds, blind=FALSE)

# Subset assay matrix for the selected genes
mat <- assay(vsd)[selected_genes, ]

# Update row names with gene symbols for the selected genes
rownames(mat) <- sig_symbols[sig_symbols %in% TFs]

# Center the data
mat <- mat - rowMeans(mat)


# ******************************************************************************
desired_order <- c('Lhx1','Lhx2', 'Nfix', 'Nfia', 'Nfic', 'Sox2', 'Sox8', 
                   'Sox9','Sox10', 'Neurog1', 'Neurog2', 'Neurod2', 'Mef2c', 
                   'Mef2d', 'Eomes', 'Tbr1') 

# Ensure the order vector only includes genes present in the matrix
ordered_genes <- desired_order[desired_order %in% rownames(mat)]

# Subset the matrix to include only the ordered genes
mat_ordered <- mat[match(ordered_genes, rownames(mat)), ]

# Make heatmap with ordered genes
tf_heatmap <- pheatmap(
  mat_ordered,
  cluster_rows = FALSE, 
  cluster_cols = TRUE, 
  show_rownames = TRUE, 
  annotation_col = as.data.frame(colData(vsd)[, "condition", drop=FALSE]),
  color = colorRampPalette(c("#CDA2CC", "white", "#FF6B35"))(50))

filename = "figures/figS2/figS2d_heatmap_tf_e16.pdf"
pdf(filename, width = 4, height = 5)
print(tf_heatmap)

dev.off()

print("heatmap saved in figures/figS2")