# Run after running 'analysis/03_deseq2_e17.R' and 'tables/scripts/tableS6.R'

# SET UP
library(org.Mm.eg.db)
library(clusterProfiler)
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

# PLOT GO RESULTS FOR GLIOGENIC TERMS: -----------------------------------------

GO_up_smo_cilia <- c('smoothened signaling pathway', 
                     'cilium assembly', 'cilium organization')

# Filter the dataframe for the relevant GO terms
GO_up_genes_df_smo_cilia <- GO_up_genes_df %>% 
  dplyr::filter(Description %in% GO_up_smo_cilia)

# Calculate -log10(p.adjust) values
GO_up_genes_df_smo_cilia$log_p.adjust <- -log10(GO_up_genes_df_smo_cilia$p.adjust)

# Create horizontal  bar plot 
GO_smo_cilia_barplot <- ggplot(GO_up_genes_df_smo_cilia, 
                               aes(x = reorder(Description, log_p.adjust), 
                                   y = log_p.adjust)) + 
  geom_bar(stat = "identity", fill = "#f48c67") + 
  coord_flip() +  
  labs(title = "Upregulated GO Terms", 
       x = "GO Term", 
       y = "-log10(p.adjust)") + 
  theme_minimal() + 
  theme(axis.text.x = element_text( hjust = 1))


# SAVE PLOT: -------------------------------------------------------------------
filename = "figures/figS2/figS2c_upregulatedgo_barplot.pdf"
pdf(filename, width = 5, height = 3)
print(GO_smo_cilia_barplot)

dev.off()

print("GO bar plot generated and saved in figures/figS2")