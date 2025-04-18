# Run after running 'analysis/02_diffbind_e16.R'
# Requires .txt files resulting from HOMER findmotifsGenome.pl

# DEFINE FILES AND PATHS: ------------------------------------------------------
# input DBA object/ diffbind resultds rds file
rds_dbObj <- "data/processed_data/atacseq_e16/r_objects/diffbind_dbObj.rds"


# READ IN HOMER MOTIF INSTANCE FILES FROM HOMER findMotifsGenome.pl: -----------
TF_CTRL <- read.csv(
  "data/processed_data/atacseq_e16/homer_output/1NFIA_2NEUROD2_3MEF2C_14EOMES_CTRL_findMotifs_inst.txt",  
  sep ="\t")
TF_NICD <- read.csv(
  "data/processed_data/atacseq_e16/homer_output/1LHX1_2NF1_8SOX10_NICD_findMotifs_inst.txt",  
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

# SEPARATE MOTIF INSTANCE DFs BASED ON TF: -------------------------------------

# CTRL ENRICHED TFs: ***********************************************************
# Run this below to see what the motifs are within the Motif.Name column
#TF_CTRL <- split(TF_CTRL, TF_CTRL$Motif.Name)

list2env(split(TF_CTRL, TF_CTRL$Motif.Name), envir = .GlobalEnv)

NFIA_CTRL <- `1-TGCCAAGH,BestGuess:NFIA/MA0670.1/Jaspar(0.975)`
EOMES_CTRL <- `15-TTTCACAGCTCT,BestGuess:PB0013.1_Eomes_1/Jaspar(0.756)`
NEUROD2_CTRL <- `2-GCCATATG,BestGuess:NEUROD2/MA0668.1/Jaspar(0.935)`
MEF2C_CTRL <- `5-CTATTTTTAG,BestGuess:Mef2c(MADS)/GM12878-Mef2c-ChIP-Seq(GSE32465)/Homer(0.944)`

remove(`1-TGCCAAGH,BestGuess:NFIA/MA0670.1/Jaspar(0.975)`)
remove(`15-TTTCACAGCTCT,BestGuess:PB0013.1_Eomes_1/Jaspar(0.756)`)
remove(`2-GCCATATG,BestGuess:NEUROD2/MA0668.1/Jaspar(0.935)`)
remove(`5-CTATTTTTAG,BestGuess:Mef2c(MADS)/GM12878-Mef2c-ChIP-Seq(GSE32465)/Homer(0.944)`)

# NICD ENRICHED TFs: ***********************************************************
# Run this below to see what the motifs are within the Motif.Name column
#TF_NICD<- split(TF_NICD, TF_NICD$Motif.Name)

list2env(split(TF_NICD, TF_NICD$Motif.Name), envir = .GlobalEnv)

LHX1_NICD <- `1-BTAATTARNN,BestGuess:Lhx1(Homeobox)/EmbryoCarcinoma-Lhx1-ChIP-Seq(GSE70957)/Homer(0.989)`
NFI_NICD <- `2-BSCTGCCAAS,BestGuess:NF1-halfsite(CTF)/LNCaP-NF1-ChIP-Seq(Unpublished)/Homer(0.939)`
SOX10_NICD <- `9-GCAYTGTK,BestGuess:SOX10/MA0442.1/Jaspar(0.795)`

remove(`1-BTAATTARNN,BestGuess:Lhx1(Homeobox)/EmbryoCarcinoma-Lhx1-ChIP-Seq(GSE70957)/Homer(0.989)`)
remove(`2-BSCTGCCAAS,BestGuess:NF1-halfsite(CTF)/LNCaP-NF1-ChIP-Seq(Unpublished)/Homer(0.939)`)
remove(`9-GCAYTGTK,BestGuess:SOX10/MA0442.1/Jaspar(0.795)`)


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
NFIA_CTRL <- process_dataframe(NFIA_CTRL, "NFIA")
EOMES_CTRL <- process_dataframe(EOMES_CTRL, "EOMES")
NEUROD2_CTRL <- process_dataframe(NEUROD2_CTRL, "NEUROD2")
MEF2C_CTRL <- process_dataframe(MEF2C_CTRL, "MEF2C")

LHX1_NICD <- process_dataframe(LHX1_NICD, "LHX1")
NFI_NICD <- process_dataframe(NFI_NICD, "NFI-halfsite")
SOX10_NICD <- process_dataframe(SOX10_NICD, "SOX10")




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
nfia_CTRL_motif_inst <- merge_dataframes(CTRL_enriched, NFIA_CTRL)
eomes_CTRL_motif_inst <- merge_dataframes(CTRL_enriched, EOMES_CTRL)
neurod2_CTRL_motif_inst <- merge_dataframes(CTRL_enriched, NEUROD2_CTRL)
mef2c_CTRL_motif_inst <- merge_dataframes(CTRL_enriched, MEF2C_CTRL)

lhx1_NICD_motif_inst <- merge_dataframes(NICD_enriched, LHX1_NICD)
nfi_NICD_motif_inst <- merge_dataframes(NICD_enriched, NFI_NICD)
sox10_NICD_motif_inst <- merge_dataframes(NICD_enriched, SOX10_NICD)




# MAKE GENE LISTS: -------------------------------------------------------------
progenitor_genes <- c('Fabp7', 'Nes', 'Pax6', 'Slc1a3', 'Sox2', 
                      'Vim', 'Nr2e1', 'Hes1', 'Hes5', 'Ednrb', 
                      'Eomes', 'Ccne2', 'Clspn', 'Gins2', 'Pcna', 
                      'Atad2', 'Mcm7', 'Mcm3', 'Slbp', 'Gmnn', 
                      'Kiaa0101', 'Mcm10', 'Rad51', 'Cdc45', 'Exo1', 
                      'Hist1h4c', 'Cdk1', 'Hist1h1b', 'Hist1h1c', 
                      'Hist1h1e', 'Ube2c', 'Rrm2', 'Zwint', 'Hmgb2', 
                      'Ccna', 'Cdca5', 'Esco2', 'Aurkb', 
                      'Kif18b', 'Ckap2l', 'Hjurp', 'Cdca8', 'Ccnb1', 
                      'Cenpf', 'Cks2', 'Pttg1', 'Cdc20', 'Top2a', 
                      'Nusap1', 'Cenpa', 'Psrc1', 'Gas2l3', 'Plk1', 
                      'Kif20a','Dleu7', 'Ncald', 'Rfx4', 'Bcan')

RGC_genes <- c('Fabp7', 'Nes', 'Pax6', 'Slc1a3', 'Sox2', 'Vim', 
               'Nr2e1', 'Hes1', 'Hes5', 'Ednrb') 

IPC_genes <- c('Eomes', 'Sema3c', 'Neurod1', 'Neurog2', 'Sstr2', 
               'Gadd45g')

proliferative_genes <- c('Fabp7', 'Nes', 'Pax6', 'Slc1a3', 'Sox2', 
                         'Vim', 'Nr2e1')

neurogenic_genes <- c('Eomes', 'Neurog2', 'Tuba1a') 

neuronal_genes <- c('Map2', 'Mapt', 'Rbfox3', 'Tbr1', 'Tubb3', 'Neurod6', 
                    'Neurod2','Satb2', 'Gria2', 'Nrp1', 'Dab1', 
                    'Nrxn3', 'Neurod4')

newborn_neurons <- c('Foxg1', 'Neurod1', 'Unc5d', 'Rnd2', 'Rnd3', 'Dcx', 
                     'Pafah1b1', 'Cdk5') 

preOPC_genes <- c('Ascl1', 'Egfr', 'Egr1', 'Qk', 'Gas1', 
                  'Sall3', 'Gng12', 'Gsx2', 'Fam181b', 'Ccnd1')

OPC_genes <- c('Sox10', 'Pdgfra', 'Olig1', 'Olig2', 'Ascl1', 'Gng12', 
               'Cnp', 'Cspg4', 'Matn4', 'Brinp3', 
               'Lhfpl3', 'Cntn1')

OL_genes <- c('Mbp', 'Plp1', 'Mag', 'Cnp', 'Mog', 'Cldn11')

glial_lineage <- c('Sox8', 'Sox9', 'Nfia')

astrocyte_genes <- c('Aldh1l1', 'Fabp7', 'Aldoc', 'Hes5', 'Aqp4')


#  COMBINE GENE LISTS TO MAKE MORE BROAD CATAGORIES: ---------------------------
neurogenic_all <- c(IPC_genes, neurogenic_genes, neuronal_genes, 
                    newborn_neurons)
neurogenic_all <- unique(neurogenic_all)

gliogenic_all <- c(OPC_genes, preOPC_genes, astrocyte_genes, glial_lineage, OL_genes)
gliogenic_all <- unique(gliogenic_all)

progenitor_all <- c(progenitor_genes, RGC_genes, proliferative_genes)
progenitor_all <- unique(progenitor_all)

ol_all <- c(preOPC_genes, OPC_genes, OL_genes)
ol_all <- unique(ol_all)

# MAKE FUNCTION TO FILTER PEAKS BASED ON GENE LISTS: ---------------------------

# Filter peaks
#Function
filter_peaks_function <- function(peak_df, gene_list) {
  filtered_peaks <- peak_df[peak_df$SYMBOL %in% gene_list, ]
  return(filtered_peaks)
}

# ******************************************************************************
# Apply merge function. Now you should have a dataframe containing the peaks 
# that contain the motif
nfia_CTRL_motif_inst <- merge_dataframes(CTRL_enriched, NFIA_CTRL)
eomes_CTRL_motif_inst <- merge_dataframes(CTRL_enriched, EOMES_CTRL)
neurod2_CTRL_motif_inst <- merge_dataframes(CTRL_enriched, NEUROD2_CTRL)
mef2c_CTRL_motif_inst <- merge_dataframes(CTRL_enriched, MEF2C_CTRL)

lhx1_NICD_motif_inst <- merge_dataframes(NICD_enriched, LHX1_NICD)
nfi_NICD_motif_inst <- merge_dataframes(NICD_enriched, NFI_NICD)
sox10_NICD_motif_inst <- merge_dataframes(NICD_enriched, SOX10_NICD)

# Filter annotated peaks *******************************************************
neuron_nfia_CTRL_motif_inst <- filter_peaks_function(nfia_CTRL_motif_inst, 
                                                     neurogenic_all)
neuron_eomes_CTRL_motif_inst <- filter_peaks_function(eomes_CTRL_motif_inst, 
                                                      neurogenic_all)
neuron_neurod2_CTRL_motif_inst <- filter_peaks_function(neurod2_CTRL_motif_inst, 
                                                        neurogenic_all)
neuron_mef2c_CTRL_motif_inst <- filter_peaks_function(mef2c_CTRL_motif_inst, 
                                                      neurogenic_all)


glia_lhx1_NICD_motif_inst <- filter_peaks_function(lhx1_NICD_motif_inst, 
                                                   gliogenic_all)
glia_nfi_NICD_motif_inst <- filter_peaks_function(nfi_NICD_motif_inst, 
                                                  gliogenic_all)
glia_sox10_NICD_motif_inst <- filter_peaks_function(sox10_NICD_motif_inst, 
                                                    gliogenic_all)

# SAVE TO FILES: ---------------------------------------------------------------
write.csv(neuron_nfia_CTRL_motif_inst, "tables/table_S4_neuron_nfia_ctrl.csv")
write.csv(neuron_eomes_CTRL_motif_inst, "tables/table_S4_neuron_eomes_ctrl.csv")
write.csv(neuron_neurod2_CTRL_motif_inst, "tables/table_S4_neuron_neurod2_ctrl.csv")
write.csv(neuron_mef2c_CTRL_motif_inst, "tables/table_S4_neuron_mef2c_ctrl.csv")

write.csv(glia_lhx1_NICD_motif_inst, "tables/table_S5_glia_lhx1_nicd.csv")
write.csv(glia_nfi_NICD_motif_inst, "tables/table_S5_glia_nfi_nicd.csv")
write.csv(glia_sox10_NICD_motif_inst, "tables/table_S5_glia_sox10_nicd.csv")

# Note: The individual csv files for Table S4 were combined into one excel file,
# 'Table_S4_neurogenic_motifs.xlxs', where each csv file is a tab. 

# Note: The individual csv files for Table S5 were combined into one excel file,
# 'Table_S5_gliogenic_motifs.xlxs', where each csv file is a tab. 
