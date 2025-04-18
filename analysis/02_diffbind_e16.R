# DEFINE FILES AND PATHS: ------------------------------------------------------
# Path to the gene raw counts matrix
input_file <- "data/meta_data/atacseq_e16/SampleSheet.csv"

# Path to output directory for r_objects
output_dir_robj <- "data/processed_data/atacseq_e16/r_objects"

# filename for RDS file containing normalized counts
filename_normcounts <- "norm_read_counts_dbObj.rds"

# filename for RDS file containing DBA object after differential analysis
filename_dbObj <- "diffbind_dbObj.rds"

# Make sure output directory exists, and if it doesn't, create one
if (!dir.exists(output_dir_robj)) {
  dir.create(output_dir_robj, recursive = TRUE)
}

# READ IN DATA:-----------------------------------------------------------------
# Sample sheet containing paths to filtered bam files and narrowpeak files
samples <- read.csv(input_file)

# CREATE DBA OBJECT: -----------------------------------------------------------
dbObj <- dba(sampleSheet=samples)

# COUNT READS: -----------------------------------------------------------------
dbObj <- dba.count(dbObj, 
                   bUseSummarizeOverlaps=TRUE, 
                   summits = 75)


# NORMALIZE DATA: --------------------------------------------------------------
# Normalization is based on sequencing depth by default
dbObj <- dba.normalize(dbObj)

# Save the DBA object after normalizing read counts for later use (making normalized counts plots)
saveRDS(dbObj, 
        file = file.path(output_dir_robj, 
                         filename_normcounts)
        )

# ESTABLISH MODEL DESIGN AND CONTRAST: -----------------------------------------
# Make the "Control" condition the baseline or denominator of the default contrast.
dbObj <- dba.contrast(dbObj, 
                      reorderMeta = list(Condition = "Control")
                      )

# PERFORM DIFFERENTIAL ANALYSIS: -----------------------------------------------
# By default, using DESeq2 but can change to edgeR
dbObj <- dba.analyze(dbObj)

# SAVE DBA OBJECT FOR LATER USE: -----------------------------------------------
# Save the DBA object after differential analysis
saveRDS(dbObj, file = file.path(output_dir_robj, 
                                filename_dbObj)
        )

print(paste0("Diffbind analysis done. RDS files saved in ", output_dir_robj))


