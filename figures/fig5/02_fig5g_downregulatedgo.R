# READ DESEQ2 RESULTS RDS FILE: ----------------------------------------
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

deseq2_results <- as.data.frame(resOrdered)

# FILTER DESEQ2 RESULTS FOR SIGNIFICANT DOWNREGULATED GENES: ---------------------
genes_down_list <- deseq2_results %>%
  dplyr:: filter(log2FoldChange < 0 & padj < 0.05)

# This is a df of all the significant downregulated genes in order by padj
genes_down_list <- genes_down_list[order(genes_down_list$padj), ]

# Display number of downregulated genes
cat("Number of downregulated genes:", nrow(genes_down_list))

# GO ANALYSIS: -----------------------------------------------------------------
# Extract the gene symbols into a list
down_genes <- genes_down_list$symbol

GO_results_down <- enrichGO(gene = down_genes,
                          OrgDb = "org.Mm.eg.db", #annotation database
                          keyType = "SYMBOL", #gene id type
                          ont = "BP") #ontology: BP (biological process)

# Save GO results for downregulated genes as a dataframe
GO_down_genes_df <- as.data.frame(GO_results_down)

# PLOT GO RESULTS FOR GLIOGENIC TERMS: -----------------------------------------
# Filter the dataframe for the relevant GO terms
GO_down_genes_df_specific <- GO_down_genes_df %>% 
  dplyr::filter(Description %in% GO_down_specific_terms)

# Calculate -log10(p.adjust) values
GO_down_genes_df_specific$log_p.adjust <- -log10(GO_down_genes_df_specific$p.adjust)

# Create the bar plot
GO_down_barplot <- ggplot(GO_down_genes_df_specific, 
                             aes(x = reorder(Description, log_p.adjust), 
                                 y = log_p.adjust)) + 
  geom_bar(stat = "identity", fill = control_color) + 
  coord_flip() +  # make horizontal barplot
  labs(title = "Downregulated GO terms", 
       x = "GO Term", 
       y = "-log10(p.adjust)") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# SAVE PLOT: -------------------------------------------------------------------
pdf(file.path(output_dir_figures, filename_downregulatedgo), width = 5, height = 3)
print(GO_down_barplot)

dev.off()

print(paste0(filename_downregulatedgo, " generated and saved in ", output_dir_figures))
