# Run after running scripts in 'preprocess/atacseq_e16_compute_matrix/'
# and put the matrix files in 'processed_data/atacseq_e16/deeptools_output'

# SET UP
library(data.table)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(readr)
library(dplyr)

# READ IN MATRIX FILES: -----------------------------------------------------------------------------------
# Unzip gz files
# Skip the JSON metadata line
ctrl_matrix <- fread("zcat < data/processed_data/atacseq_e16/deeptools_output/computematrix/CTRL.matrix.gz", skip = 1)
nicd_matrix <- fread("zcat < data/processed_data/atacseq_e16/deeptools_output/computematrix/NICD.matrix.gz", skip = 1) 

# PREPARE DATA: -----------------------------------------------------------------------------------
# Extract the numeric data for plotting starting at column 7
ctrl_data <- ctrl_matrix[, 7:ncol(ctrl_matrix), with = FALSE]
nicd_data <- nicd_matrix[, 7:ncol(nicd_matrix), with = FALSE]

# Calculate the mean across all regions for each position
ctrl_mean <- colMeans(ctrl_data, na.rm = TRUE)
nicd_mean <- colMeans(nicd_data, na.rm = TRUE)

# Create a dataframe for plotting
plot_data <- data.frame(
  Position = seq(-1000, 1000, length.out = length(ctrl_mean)),  # Set the x-axis from -1000 to 1000
  CTRL = ctrl_mean,
  NICD = nicd_mean
)

# Melt the dataframe for ggplot2
plot_data_melt <- melt(plot_data, id.vars = "Position", variable.name = "Condition", value.name = "Score")

# PLOT TSS ENRICHMENT: -----------------------------------------------------------------------------------
# Plot using ggplot2 and save as png

tss_enrichment_plot <- ggplot(plot_data_melt, aes(x = Position, 
                                                  y = Score, 
                                                  color = Condition)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "dashed") +  
  labs(title = "TSS Enrichment Plot",
       x = "Position relative to TSS",
       y = "Enrichment Score") +
  theme_classic() +
  scale_color_manual(values = c("#AFA2C2", "#833D66")) +  
  theme(legend.position = "none")


filename = "figures/fig3/fig3c_tssenrichmentplot.pdf"
pdf(filename, width = 3, height = 3)
print(tss_enrichment_plot)

dev.off()

print("TSS enrichment plot saved in figures/fig3/")