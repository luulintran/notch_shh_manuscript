# PREPARE DATA: ----------------------------------------------------------------
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
plot_data_melt <- melt(plot_data, 
                       id.vars = "Position", 
                       variable.name = "Condition", 
                       value.name = "Score"
                       )

# PLOT TSS ENRICHMENT: ---------------------------------------------------------
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
  scale_color_manual(values = c(control_color, mutant_color)) +  
  theme(legend.position = "none")


# SAVE PLOT: -------------------------------------------------------------------
pdf(file.path(output_dir_figures, filename_tssenrichment), width = 3, height = 3)
print(tss_enrichment_plot)

dev.off()

print(paste0(filename_tssenrichment, " saved in ", output_dir_figures))