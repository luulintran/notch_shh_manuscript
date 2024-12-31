# SET UP
library(DiffBind)
library(tidyverse)
library(dplyr)
library(readr)

# Set working directory to project
# setwd("..")

# READ IN DATA:------------------------------------------------------------------------------------------
# Sample sheet containing paths to filtered bam files and narrowpeak files
samples <- read.csv("data/meta_data/atacseq_e16/SampleSheet.csv")

# CREATE DBA OBJECT: -------------------------------------------------------------------------------------
dbObj <- dba(sampleSheet=samples)

# COUNT READS: -------------------------------------------------------------------------------------------
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE, summits = 75)


# NORMALIZE DATA: -------------------------------------------------------------------------------------------
# Normalization is based on sequencing depth by default
dbObj <- dba.normalize(dbObj)

# Save the DBA object after normalizing read counts for later use (making normalized counts plots)
saveRDS(dbObj, file = "data/processed_data/atacseq_e16/r_objects/norm_read_counts_dbObj.rds")

# ESTABLISH MODEL DESIGN AND CONTRAST: ---------------------------------------------------------------------
# Make the "Control" condition the baseline or denominator of the default contrast.
dbObj <- dba.contrast(dbObj, 
                      reorderMeta = list(Condition = "Control"))

# PERFORM DIFFERENTIAL ANALYSIS: ---------------------------------------------------------------------
# By default, using DESeq2 but can change to edgeR
dbObj <- dba.analyze(dbObj)

# SAVE DBA OBJECT FOR LATER USE: -------------------------------------------------------------------------
# Save the DBA object after differential analysis
saveRDS(dbObj, file = "data/processed_data/atacseq_e16/r_objects/diffbind_dbObj.rds")

print("Diffbind analysis done. RDS files saved in r_objects/")


