# READ DBA OBJECT FOLLOWING DIFFBIND ANALYSIS FROM RDS FILE: -------------------
# Make sure it's the one after running diffbind analysis
dbObj <- readRDS(rds_dbObj)

# MAKE PCA PLOT: ---------------------------------------------------------------
# Set file name and pdf settings for saving
pdf(file.path(output_dir_figures, filename_pca), width = 5, height = 5)

# Perform PCA on all accessible chromatin sites
dba.plotPCA(dbObj, DBA_CONDITION)

dev.off()

print(paste0(filename_pca, " generated and saved in ", output_dir_figures))

