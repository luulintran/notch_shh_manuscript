# Run after running 'analysis/01_deseq2_e16.R'

# SET UP
library(ggplot2)
library(DESeq2)
library(extrafont)

# import Arial font
font_import(pattern = "Arial", prompt = FALSE)

# READ DDS OBJECT FROM RDS FILE: -----------------------------------------------
dds <- readRDS(
  file = "data/processed_data/rnaseq_e16/r_objects/deseq2_dds_e16.rds")

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
  scale_color_manual(values = c("control" = "#76BAE0", "NICD" = "#B8396B")) +
  theme_classic(base_family = "Arial", base_size = 12) +
  theme(
    axis.line = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16)
  )

# Set filename to save
filename <- "figures/fig2/fig2b_pcaplot.pdf"
pdf(filename, width = 3, height = 3)
print(pcaplot)

dev.off()

print("PCA plot was generated and saved in figures/fig2/")