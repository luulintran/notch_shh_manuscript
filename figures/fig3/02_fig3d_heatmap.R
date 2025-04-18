# READ DBA OBJECT FOLLOWING DIFFBIND ANALYSIS FROM RDS FILE: -------------------
# Make sure it's the one after running diffbind analysis
dbObj <- readRDS(rds_dbObj)

# SAVE HEATMAP: ----------------------------------------------------------------
pdf(file.path(output_dir_figures, filename_heatmap), width = 5, height = 6)
readscores <- dba.plotHeatmap(dbObj, contrast = 1, correlations = FALSE,
                              scale = "row", colScheme = hmap)

dev.off()

print(paste0(filename_heatmap, " generated and saved in ", output_dir_figures))

