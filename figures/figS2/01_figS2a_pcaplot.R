# Run after running 'analysis/03_deseq2_e17.R'

# SET UP
library(ggplot2)
library(DESeq2)
library(extrafont)

# import Arial font
font_import(pattern = "Arial", prompt = FALSE)

# Run after running 'analysis/03_deseq2_e17.R'

# READ DDS OBJECT FROM RDS FILE: -----------------------------------------------
dds <- readRDS(
  file = "data/processed_data/rnaseq_e17/r_objects/deseq2_dds_e17.rds")

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
  scale_color_manual(values = c("control" = "#008080", "NICD" = "#dc9b06")) +
  theme_classic(base_family = "Arial", base_size = 12) +
  theme(
    axis.line = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16)
  )

# Set filename to save
filename <- "figures/figS2/figS2a_pcaplot.pdf"
pdf(filename, width = 3, height = 3)
print(pcaplot)

dev.off()

print("PCA plot saved in figures/figS2/")