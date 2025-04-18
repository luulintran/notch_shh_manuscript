# Run 'analysis/03_deseq2_e17.R' and 'tables/scripts/tableS6.R' if you haven't already

# READ DESEQ2 RESULTS RDS FILE: ------------------------------------------------
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

# FILTER DESEQ2 RESULTS FOR SIGNIFICANT UPREGULATED GENES: ---------------------
genes_up_list <- deseq2_results %>%
  dplyr:: filter(log2FoldChange > 0 & padj < 0.05)

# This is a df of all the significant upregulated genes in order by padj
genes_up_list <- genes_up_list[order(genes_up_list$padj), ]

# Display number of upregulated and downregulated genes
cat("Number of upregulated genes:", nrow(genes_up_list))

# GO ANALYSIS: -----------------------------------------------------------------
# Extract the gene symbols into a list
up_genes <- genes_up_list$symbol

GO_results_up <- enrichGO(gene = up_genes,
                          OrgDb = "org.Mm.eg.db", #annotation database
                          keyType = "SYMBOL", #gene id type
                          ont = "BP") #ontology: BP (biological process)

# Save GO results for upregulated genes as a dataframe
GO_up_genes_df <- as.data.frame(GO_results_up)

# PLOT GO RESULTS FOR SPECIFIC TERMS: ------------------------------------------

# Filter the dataframe for the relevant GO terms
GO_up_genes_df_specific <- GO_up_genes_df %>% 
  dplyr::filter(Description %in% GO_up_specific_terms)

# Calculate -log10(p.adjust) values
GO_up_genes_df_specific$log_p.adjust <- -log10(GO_up_genes_df_specific$p.adjust)

# Make the bar plot
GO_up_barplot <- ggplot(GO_up_genes_df_specific, 
                           aes(x = reorder(Description, log_p.adjust),
                               y = log_p.adjust)) + 
  geom_bar(stat = "identity", fill = mutant_color) +
  coord_flip() +  # Flip coordinates to make it horizontal
  labs(title = "Upregulated GO Terms", 
       x = "GO Term", 
       y = "-log10(p.adjust)") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# SAVE PLOT: -------------------------------------------------------------------
pdf(file.path(output_dir_figures, filename_upregulatedgo), width = 5, height = 3)
print(GO_up_barplot)

dev.off()

print(paste0(filename_upregulatedgo, " generated and saved in ", output_dir_figures))