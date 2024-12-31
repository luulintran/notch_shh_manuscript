
# Run after running 'analysis/03_deseq2_e17.R'
library(pheatmap)
library(DESeq2)
library(org.Mm.eg.db)
library(dplyr)
library(tidyverse)
library(readr)


# READ DDS OBJECT FOLLOWING DESEQ2 ANALYSIS FROM RDS FILE: ---------------------
dds <- readRDS("data/processed_data/rnaseq_e17/r_objects/deseq2_dds_e17.rds")
res <- results(dds)

# ANNOTATE RESULTS WITH GENE SYMBOLS AND ENTREZ IDS: ---------------------------
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

# MAKE LIST OF SELECTED SHH PATHWAY GENES: -------------------------------------
Shh_gene_list <- c("Smo", "Boc", "Cdon", 
             "Gas1", "Gli1", "Gli2", "Gli3", 
             "Sufu", "Disp1", "Iqce", "Efcab7", 
             "Ptch1", "Ptch2", "Hhip", "Sufu", 
             "Gsk3b", "Ck1", "Pcaf", "Cul3", "Kif7")

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

# Filter sig_genes list for genes in Shh_gene_list by symbol
selected_genes <- sig_genes[sig_symbols %in% Shh_gene_list]

# EXTRACT TRANSFORMED VALUES: **************************************************
vsd <- vst(dds, blind=FALSE)

# Subset assay matrix for the selected genes
mat <- assay(vsd)[selected_genes, ]

# Update row names with gene symbols for the selected genes
rownames(mat) <- sig_symbols[sig_symbols %in% Shh_gene_list]

# Center the data
mat <- mat - rowMeans(mat)


# PLOT HEATMAP: ----------------------------------------------------------------
heatmap_sig_shh <- pheatmap(
  mat,
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  show_rownames = TRUE, 
  annotation_col = as.data.frame(colData(vsd)[, "condition", drop=FALSE]),
  color = colorRampPalette(c("#CDA2CC", "white", "#FF6B35"))(50)) 

filename = "figures/fig5/fig5h_heatmap.pdf"
pdf(filename, width = 5, height = 6)
print(heatmap_sig_shh)

dev.off()

print("Heatmap saved in figures/fig5")
