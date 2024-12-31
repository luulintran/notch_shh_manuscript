# Run after running 'analysis/03_deseq2_e17.R' and 'tables/scripts/ableS6.R'

# SET UP
library(clusterProfiler)
library(org.Mm.eg.db)
library(DESeq2)
library(dplyr)
library(ggplot2)

# READ DESEQ2 RESULTS CSV FILE: ----------------------------------------
dds <- readRDS("data/processed_data/rnaseq_e17/r_objects/deseq2_dds_e17.rds")

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
genes_down_list <- deseq2_results %>%
  dplyr:: filter(log2FoldChange < 0 & padj < 0.05)

# This is a df of all the significant upregulated genes in order by padj
genes_down_list <- genes_down_list[order(genes_down_list$padj), ]

# Display number of upregulated and downregulated genes
cat("Number of downregulated genes:", nrow(genes_down_list))

# GO ANALYSIS: -----------------------------------------------------------------
# Extract the gene symbols into a list
down_genes <- genes_down_list$symbol

GO_results_down <- enrichGO(gene = down_genes,
                          OrgDb = "org.Mm.eg.db", #annotation database
                          keyType = "SYMBOL", #gene id type
                          ont = "BP") #ontology: BP (biological process)

# Save GO results for upregulated genes as a dataframe
GO_down_genes_df <- as.data.frame(GO_results_down)

# PLOT GO RESULTS FOR GLIOGENIC TERMS: -----------------------------------------
# Define the GO terms of interest
GO_neurogenesis <- c('regulation of neurogenesis', 
                     'dendrite development', 
                     'positive regulation of neurogenesis', 
                     'regulation of synapse structure or activity', 
                     'regulation of neuron differentiation')

# Filter the dataframe for the relevant GO terms
GO_down_genes_df_neuron <- GO_down_genes_df %>% 
  dplyr::filter(Description %in% GO_neurogenesis)

# Calculate -log10(p.adjust) values
GO_down_genes_df_neuron$log_p.adjust <- -log10(GO_down_genes_df_neuron$p.adjust)

# Create the bar plot
GO_neuron_barplot <- ggplot(GO_down_genes_df_neuron, 
                             aes(x = reorder(Description, log_p.adjust), 
                                 y = log_p.adjust)) + 
  geom_bar(stat = "identity", fill = "#ca8bc9") + 
  coord_flip() +  # make horizontal barplot
  labs(title = "Downregulated GO terms", 
       x = "GO Term", 
       y = "-log10(p.adjust)") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# SAVE PLOT: -------------------------------------------------------------------
filename = "figures/fig5/fig5g_downregulatedgo_barplot.pdf"
pdf(filename, width = 5, height = 3)
print(GO_neuron_barplot)

dev.off()

print("GO bar plot generated and saved in figures/fig5")





