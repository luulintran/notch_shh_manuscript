## Run after running 'analysis/02_diffbind_e16.R'

library(Diffbind)

# READ DBA OBJECT FOLLOWING DIFFBIND ANALYSIS FROM RDS FILE: -------------------
# Make sure it's the one after running diffbind analysis
dbObj.DB <- readRDS(
  "data/processed_data/atacseq_e16/r_objects/diffbind_dbObj.rds")

# MAKE HEATMAP: ----------------------------------------------------------
# Set file name and pdf settings for saving
filename <- "figures/fig3/fig3d_heatmap.pdf"
pdf(filename, width = 5, height = 6)

# Choose heatmap colors
hmap <- colorRampPalette(c("white", "#938fc2", "#40027e")) (n =100)

# Plot heatmap
readscores <- dba.plotHeatmap(dbObj, contrast = 1, correlations = FALSE,
                              scale = "row", colScheme = hmap)

dev.off()

print("Heatmap generated and saved in figures/fig3/")

