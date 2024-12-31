# Run after running 'analysis/03_deseq2_e17.R' and 'tables/scripts/tableS6.R'

# SET UP
library(ggplot2)
library(DESeq2)
library(extrafont)
library(tidyverse)
library(readr)
library(dplyr)

# Load CSV file of DESeq2 results, ordered by padj and containing gene symbol 
# and entrez id's.
res_df <- read.csv("tables/table_S6_rnaseq_e17_deseq2_results.csv")

# MAKE VOLCANO PLOT WITH SPECIFIC GENES LABELED: -------------------------------

# Make a list of selected genes. Here, I want to show Notch genes, 
# a few neurogenic genes, and a few progenitor and gliogenic genes.
specific_genes <- c('Dll3', 'Dll1', 'Neurod4', 'Hey1', 
                    'Neurog2', 'Hey1', 'Hes5', 'Notch1', 
                    'Sox2', 'Sox9', 'Notch2')

# filter the deseq2 results (res_df) dataframe to include only the genes 
# in the specific_genes list based on symbol column
specific_genes_res <- res_df %>%
  dplyr:: filter(symbol %in% specific_genes)

# define significance threshold (padj 0.05)
alpha <- 0.05

# Use mutate() to make color_group column in res_df based on log2FoldChange and 
# padj values

# Here, I want to make the dots that are positive log2foldchange to be one 
# color, and the negative log2foldchange to be a different color.

# I also want all dots with padj > 0.05 to be grey (non-significant)
res_df <- res_df %>% 
  dplyr::mutate(
    color_group = 
      # if log2FoldChange < 0 and padj < alpha, then make blue
      ifelse(log2FoldChange < 0 & padj < alpha, "#CDA2CC", 
             # if log2FoldChange > 0 and padj < alpha, then make pink
             ifelse(log2FoldChange > 0 & padj < alpha, "#FF6B35", 
                    # else, make grey
                    "grey"))
  )


specific_genes_res <- specific_genes_res %>%
  dplyr:: mutate(
    color_group =
      ifelse(log2FoldChange < 0 & padj < alpha, "#CDA2CC",
             ifelse(log2FoldChange > 0 & padj < alpha, "#FF6B35", 
                    "grey"))
    )


# ******************************************************************************
# Plot volcano plot with different colors for log2foldchange negative and 
# positive values, label only the genes in the specific_genes_res dataframe,
# make non significant values with padj > 0.05 grey

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = color_group), alpha = 0.8, size = 0.5) + 
  scale_color_manual(values = c("#CDA2CC", "#FF6B35", "grey")) + 
  theme_minimal(base_family = "Arial", base_size = 8) + 
  labs(
    title = "NICD vs CTRL E17.5",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  ) + 
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  ) + 
  # Outline specific genes in black
  geom_point(data = specific_genes_res, aes(fill = color_group), 
             color = "black", size = 0.5, shape = 21, stroke = 0.5) + 
  # Set fill color for specific genes
  scale_fill_manual(values = c("#CDA2CC", "#FF6B35", "grey")) + 
  # Add gene labels with ggrepel
  ggrepel::geom_text_repel(
    data = specific_genes_res, 
    aes(label = paste0("italic('", symbol, "')")), 
    size = 2, 
    point.padding = 0.3, 
    max.overlaps = 20, 
    segment.color = "black", 
    segment.size = 0.5, 
    parse = TRUE
  ) + 
  theme(legend.position = "none") + 
  coord_cartesian(xlim = c(-10, 10), ylim = c(0, 100))


# Set filename
filename = "figures/figS2/figS2b_volcanoplot.pdf"
pdf(filename, width = 2.5, height = 2.5)
print(volcano_plot)

dev.off()

print("Volcano plot generated and saved in figures/figS2/")
