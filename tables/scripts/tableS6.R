# Run after running 'analysis/02_diffbind_e16.R'
# Requires .txt files resulting from HOMER findmotifsGenome.pl

# DEFINE FILES AND PATHS: ------------------------------------------------------
# input DBA object/ diffbind resultds rds file
rds_dbObj <- "data/processed_data/atacseq_e16/r_objects/diffbind_dbObj.rds"

# READ IN HOMER MOTIF INSTANCE FILES FROM HOMER findMotifsGenome.pl: -----------

RBPJ_CTRL <- read.csv(
  "data/processed_data/atacseq_e16/homer_output/RBPJ_kr_CTRL_findMotifs_inst.txt",  
  sep ="\t")
RBPJ_NICD <- read.csv(
  "data/processed_data/atacseq_e16/homer_output/RBPJ_kr_NICD_findMotifs_inst.txt",  
  sep ="\t")


# LOAD RDS FILE OF DIFFBIND DBA OBJECT: ----------------------------------------
dbObj <- readRDS(
  file = rds_dbObj)

# Create dba report with results including stats
dbObj.DB <- dba.report(dbObj)

# CONVERT SEQNAMES TO UCSC STYLE: ----------------------------------------------
# We want UCSC style and keep mm10 genome build consistent
# Convert to UCSC style from Ensembl

# Takes the seqnames of GRanges object and makes a vector of new seqnames that
# matches the UCSC style.
UCSC_newstyle <- mapSeqlevels(seqlevels(dbObj.DB), "UCSC") 

# Changes the seqnames in dbObj to UCSC style, so now seqnames column has "chr"
dbObj.DB_UCSC <- renameSeqlevels(dbObj.DB, UCSC_newstyle) 

# Store as a dataframe
out <- as.data.frame(dbObj.DB_UCSC)

# Add peak id as "peak_(rownumber)" for later use, 
# like for HOMER motif enrichment analysis
out <- out %>%
  mutate(peak_id = paste("peak_", seqnames, "_",row_number(), sep = ""))

# ASSIGN PEAKS TO NEAREST TSS: -------------------------------------------------
# Using ChIPseeker

# Get annotation data
# This is the UCSC annotation data for mm10
# The ATAC-seq data was mapped to mm10 so we want to keep it consistent:
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Annotate differential peaks to nearest TSS ***********************************
peak_gr <- GRanges(seqnames = out$seqnames,
                   ranges = IRanges(start = out$start, end = out$end))

# Annotate +/- 3 kb around the TSS
peak_anno <- annotatePeak(peak_gr, 
                          tssRegion = c(-3000,3000), 
                          TxDb = txdb, 
                          annoDb = "org.Mm.eg.db")
# Store annotated peaks as a dataframe
peak_anno_df <- as.data.frame(peak_anno)

# merge the out df (containing stats) and the 
# peak_anno_df (containing annotations) and match by seqnames, start, and end
merged_df <- merge(out, peak_anno_df, by = c("seqnames", "start", "end"))
merged_df <- subset(merged_df, select = -c(width.y, strand.y))

# Order df by FDR and rename df as peak_list
peak_list <- merged_df[order(merged_df$FDR), ]

# SEPARATE PEAK LIST BY CTRL (Fold < 0) and NICD (Fold > 0) ENRICHED: ----------

CTRL_enriched <- peak_list %>% 
  dplyr::filter(FDR < 0.05 & Fold < 0)

NICD_enriched <- peak_list %>% 
  dplyr::filter(FDR < 0.05 & Fold > 0)

# CLEAN UP DATAFRAMES: --------------------------------------------------------

#### Define functions for cleaning up and merging dataframes 

##### Function to remove unnecessary columns in diffbind dataframes ************
process_columns <- function(df) {
  # Remove unnecessary columns
  df <- df[, -c(4:8, 14:18, 20)]
  
  # Change order of columns
  #df <- df[, c(7:12, 1:6)]
  
  # Return modified dataframe
  return(df)
}


##### Function to clean up homer dataframes ************************************
process_dataframe <- function(df, tf_name) {
  # Rename the first column to "peak_id"
  colnames(df)[1] <- "peak_id"
  
  # Replace all values in the "Motif.Name" column with the TF name
  df$Motif.Name <- tf_name
  
  # Return cleaned up dataframe
  return(df)
}


### Clean up dataframes ********************************************************
# Apply the function to remove unnecessary columns 
# and reorder the remaining columns
CTRL_enriched <- process_columns(CTRL_enriched)
NICD_enriched <- process_columns(NICD_enriched)


# Apply function to replace column name of PositionID and values in Motif.Name

RBPJ_CTRL <- process_dataframe(RBPJ_CTRL, "RBPJ")

RBPJ_NICD <- process_dataframe(RBPJ_NICD, "RBPJ")


# Merge dataframes *************************************************************

##### Function to merge dataframes 
merge_dataframes <- function(df1, df2) {
  # Merge the data frames based on the "peak_id" column
  merged_df <- merge(df1, df2, by = "peak_id", all = TRUE)
  
  # Remove rows where the "Motif.Name" column is NA or empty
  merged_df <- merged_df[
    !is.na(merged_df$Motif.Name) & merged_df$Motif.Name != "", ]
  
  # Group by peak_id, keeping duplicates together
  grouped_df <- merged_df %>%
    group_by(peak_id) 
  
  # Return the grouped dataframe
  return(grouped_df)
}



### Merge dataframes into one. *************************************************
# Merge, remove NA's (peaks that didn't have the motif), 
# and group duplicate peaks together (peaks that contain the motif more than 
# once, depending on position)

rbpj_CTRL_motif_inst <- merge_dataframes(CTRL_enriched, RBPJ_CTRL)

rbpj_NICD_motif_inst <- merge_dataframes(NICD_enriched, RBPJ_NICD)


# SAVE TO FILES: ---------------------------------------------------------------
write.csv(rbpj_NICD_motif_inst, "tables/table_S6_rbpj_nicd.csv")
write.csv(rbpj_CTRL_motif_inst, "tables/table_S6_rbpj_ctrl.csv")


# Note: The individual csv files for Table S6 were combined into one excel file,
# 'Table_S6_rbpj_motif.xlxs', where each csv file is a tab. 
