# READ DDS OBJECT FROM RDS FILE: -----------------------------------------------
dds <- readRDS(
  file = rds_deseq2_results)

# EXTRACT TRANSFORMED VALUES: --------------------------------------------------
vsd <- vst(dds, blind=FALSE)

# CUSTOMIZE A PCA PLOT: --------------------------------------------------------
# extract PCA data from vsd
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

# get percentage of variance for first two principal components
percentVar <- round(100 * attr(pcaData, "percentVar"))

# plot PCA
pcaplot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4) +
  labs(title = "PCA",
       x = paste0("PC1[", percentVar[1], "%]"),
       y = paste0("PC2[", percentVar[2], "%]")) +
  scale_color_manual(values = c("e13" = control_color, 
                                "e17" = mutant_color)) + 
  theme_classic(base_family = "Arial", base_size = 12) +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16)
  )

# SAVE 
pdf(file.path(output_dir_figures, filename_pca), width = 3, height = 3)
print(pcaplot)

dev.off()

print(paste0(filename_pca, " saved in ", output_dir_figures))