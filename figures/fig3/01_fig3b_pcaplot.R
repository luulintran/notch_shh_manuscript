# Run after running 'analysis/02_diffbind_e16.R'
library(Diffbind)

# READ DBA OBJECT FOLLOWING DIFFBIND ANALYSIS FROM RDS FILE: -------------------
# Make sure it's the one after running diffbind analysis
dbObj.DB <- readRDS(
  "data/processed_data/atacseq_e16/r_objects/diffbind_dbObj.rds")

# MAKE PCA PLOT: ---------------------------------------------------------------
# Set file name and pdf settings for saving
filename <- "figures/fig3/fig3b_pcaplot.pdf"
pdf(filename, width = 5, height = 5)

# Perform PCA on all accessible chromatin sites
dba.plotPCA(dbObj, DBA_CONDITION)

dev.off()

print("PCA plot generated and saved in figures/fig3/")

