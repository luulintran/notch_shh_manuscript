# READ IN DATA: ----------------------------------------------------------------
# Load CSV file of DESeq2 results, ordered by padj and containing gene symbol 
# and entrez id's.
res_df <- read.csv(csv_deseq2_results)

# MAKE VOLCANO PLOT WITH SPECIFIC GENES LABELED: -------------------------------

# filter the deseq2 results (res_df) dataframe to include only the genes 
# in the specific_genes list based on symbol column
specific_genes_res <- res_df %>%
  dplyr:: filter(symbol %in% specific_genes)

# Use mutate() to make color_group column in res_df based on log2FoldChange and 
# padj values

# Here, I want to make the dots that are positive log2foldchange to be one 
# color, and the negative log2foldchange to be a different color.

# I also want all dots with padj > 0.05 to be grey (non-significant)
res_df <- res_df %>% 
  dplyr::mutate(
    color_group = 
      # if log2FoldChange < 0 and padj < alpha, then make blue
      ifelse(log2FoldChange < 0 & padj < alpha, downreg_color, 
             # if log2FoldChange > 0 and padj < alpha, then make pink
             ifelse(log2FoldChange > 0 & padj < alpha, upreg_color, 
                    # else, make grey
                    "grey"))
  )


specific_genes_res <- specific_genes_res %>%
  dplyr:: mutate(
    color_group =
      ifelse(log2FoldChange < 0 & padj < alpha, downreg_color,
             ifelse(log2FoldChange > 0 & padj < alpha, upreg_color, 
                    "grey"))
    )


# ******************************************************************************
# Plot volcano plot with different colors for log2foldchange negative and 
# positive values, label only the genes in the specific_genes_res dataframe,
# make non significant values with padj > 0.05 grey

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = color_group), alpha = 0.8, size = 0.5) + 
  scale_color_manual(values = c(downreg_color, upreg_color, "grey")) + 
  theme_minimal(base_family = "Arial", base_size = 8) + 
  labs(
    title = "WT E13.5 vs E17.5",
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
  scale_fill_manual(values = c(downreg_color, upreg_color, "grey")) + 
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


# SAVE
pdf(file.path(output_dir_figures, filename_volcanoplot), width = 2.5, height = 2.5)
print(volcano_plot)

dev.off()

print(paste0(filename_volcanoplot, " generated and saved in ", output_dir_figures))
