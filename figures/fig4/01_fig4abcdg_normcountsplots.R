# Run after running 'analysis/02_diffbind.R' and 'tables/scripts/tableS3.R'

# SET UP
library(DiffBind)
library(ChIPseeker)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)
library(dplyr)
library(readr)
library(ggplot2)
library(car)

# MAKE GENE LISTS: -------------------------------------------------------------
# Gene lists are based on cell type markers to determine 
# what cell identity progenitors are based on open chromatin
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
                      'Nusap1', 'Cenpa', 'Psrc1', 'Gas2l3', 'Plk1', 'Kif20a')

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

OPC_genes <- c('Sox10', 'Pdgfra', 'Olig1', 'Olig2', 'Ascl1', 'Gng12', 
               'Cnp', 'Cspg4', 'Matn4', 'Brinp3', 
               'Lhfpl3', 'Cntn1')

preOPC_genes <- c('Ascl1', 'Ccnd1', 'Dleu7', 'Egfr', 'Egr1', 'Qk', 'Gas1', 
                  'Sall3', 'Gng12', 'Ncald', 'Gsx2', 
                  'Fam181b', 'Rfx4', 'Bcan')

astrocyte_genes <- c('Aldh1l1', 'Fabp7', 'Aldoc', 'Hes5', 'Aqp4')

Shh_pathway_genes <- c('Gli1', 'Smo', 'Ptch1', 'Boc', 'Cdon', 'Gas1', 
                       'Gli2','Ptch2', 'Hhip', 'Gli3')



#  COMBINE GENE LISTS TO MAKE MORE BROAD CATAGORIES: ---------------------------
neurogenic_all <- c(IPC_genes, neurogenic_genes, neuronal_genes, 
                    newborn_neurons)

gliogenic_all <- c(OPC_genes, preOPC_genes, astrocyte_genes)

progenitor_all <- c(progenitor_genes, RGC_genes, proliferative_genes)

opc_all <- c(preOPC_genes, OPC_genes)

# RETRIEVE NORMALIZED READ COUNTS: ---------------------------------------------
# Load RDS file that contains dbObj following normalization
dbObj <- readRDS("data/processed_data/atacseq_e16/r_objects/norm_read_counts_dbObj.rds")

# Extract normalized read counts
norm_counts <- dba.peakset(dbObj, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)

# ANNOTATE PEAKS TO NEAREST GENE: ----------------------------------------------
# Load annotation database
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Retrieve peakset data from DBA object as dataframe ("DBA_DATA_FRAME")
peaks_gr <- dba.peakset(dbObj, DataType = DBA_DATA_FRAME, bRetrieve = TRUE)

# Convert to GRanges object, containing seqnames, start, end and the norm counts 
# for each replicate
peaks_gr <- GRanges(
  seqnames = peaks_gr$CHR,
  ranges = IRanges(start = peaks_gr$START, end = peaks_gr$END),
  CTRL1 = peaks_gr$CTRL1,
  CTRL2 = peaks_gr$CTRL2,
  CTRL3 = peaks_gr$CTRL3,
  NICD1 = peaks_gr$NICD1,
  NICD2 = peaks_gr$NICD2,
  NICD3 = peaks_gr$NICD3
)

# Map seq levels to UCSC style
UCSC_newstyle <- mapSeqlevels(seqlevels(peaks_gr), "UCSC") 
peaks_gr_UCSC <- renameSeqlevels(peaks_gr, UCSC_newstyle)


# Annotate peaks to the nearest gene (-3kb, +3kb from TSS) using Chipseeker
annotated_peaks <- annotatePeak(peaks_gr_UCSC, 
                                tssRegion = c(-3000,3000), 
                                TxDb = txdb, 
                                annoDb = "org.Mm.eg.db")

# Convert to a dataframe
# This dataframe has annotations, gene names, and the norm counts 
# for each replicate.
annotated_peaks <- as.data.frame(annotated_peaks)


annotated_peaks_subset <- subset(annotated_peaks, 
                                 select = -c(width, strand, annotation, 
                                           geneChr, geneStart, geneEnd, 
                                           geneLength, geneStrand, 
                                           transcriptId, distanceToTSS,
                                           geneId, ENSEMBL, SYMBOL, GENENAME))


# READ DIFFBIND RESULTS FILE CONTAINING SIGNIFICANT DIFFERENTIAL PEAKS: --------
diffbind_res_df <- read.csv("tables/table_S3_atacseq_e16_diffbind_results.csv")

# merge above dataframe with the new dataframe (merged_df) 
# with the normalized read counts. 

# Now I have a dataframe with only the significant peaks 
# (fold change > or < 0 and FDR < 0.05)
merged_df_enriched <- merge(diffbind_res_df, 
                            annotated_peaks_subset, 
                            by = c("seqnames", "start", "end"),
                            all.x = TRUE)


if (nrow(merged_df_enriched) == nrow(diffbind_res_df)) {
  print("Normalized counts and DiffBind results merged successfully")
} else {
  print("Warning: The merged dataframe does not have the same number of rows 
        as DiffBind results. Check for duplicates or missing values.")
}

# save over merged_df to use this in the violin plots below. 
merged_df <- merged_df_enriched


# Remove the unncessary columns and also Conc_NICD and Conc_Control 
# or these will be grepped with NICD and CTRL for making plots later
merged_df <- subset(merged_df, select = -c(width.x, strand.x, Conc, Conc_NICD,
                                           Conc_Control, p.value,
                                           geneChr, geneStart, geneEnd, 
                                           geneLength, geneStrand, 
                                           transcriptId)
                    )

merged_df <- merged_df[, c("seqnames", "start", "end", "peak_id", "Fold", "FDR",
                           "annotation", "geneId", "distanceToTSS", "ENSEMBL", 
                           "SYMBOL", "GENENAME", "CTRL1", "CTRL2", "CTRL3", 
                           "NICD1", "NICD2", "NICD3")]

# FILTER PEAKS FUNCTION: -------------------------------------------------------
# Make a function used for filtering peaks

filter_peaks_function <- function(peak_df, gene_list) {
  filtered_peaks <- peak_df[peak_df$SYMBOL %in% gene_list, ]
  return(filtered_peaks)
}


# NORMALIZED ACCESSIBILITY OF GENE CATAGORY LISTS: -----------------------------
# DEFINE FUNCTIONS FOR PLOTTING NORM ACCESSIBILITY *****************************

# Normalized counts function with log transformation ***************************
normalized_counts_function <- function(merged_df, gene_list) {
  normalized_counts <- merged_df %>% 
    dplyr::filter(SYMBOL %in% gene_list)
  
  # Calculate average read count for each gene for CTRL and NICD separately
  
  # indices for CTRL columns
  ctrl_cols <- grep("CTRL", colnames(normalized_counts)) 
  # indices for NICD columns
  nicd_cols <- grep("NICD", colnames(normalized_counts)) 
  
  normalized_counts$CTRL_avg <- rowMeans(normalized_counts[, ctrl_cols])
  normalized_counts$NICD_avg <- rowMeans(normalized_counts[, nicd_cols])
  
  # Apply log2 transformation with a pseudocount
  normalized_counts$CTRL_avg_log <- log2(normalized_counts$CTRL_avg + 1)
  normalized_counts$NICD_avg_log <- log2(normalized_counts$NICD_avg + 1)
  
  return(normalized_counts)
}

# Average read counts function for log-transformed values **********************
avg_read_counts_function <- function(normalized_counts) {
  avg_read_counts <- data.frame(
    Condition = rep(c("CTRL", "NICD"), 
                    each = nrow(normalized_counts)),
    Average_Read_Counts = c(normalized_counts$CTRL_avg_log, 
                            normalized_counts$NICD_avg_log)
  )
  return(avg_read_counts)
}

# Statistical testing function on log-transformed data *************************
stats_test_function <- function(normalized_counts) {
  # Test for normal distribution on log-transformed values
  shapiro_ctrl <- shapiro.test(normalized_counts$CTRL_avg_log)
  shapiro_nicd <- shapiro.test(normalized_counts$NICD_avg_log)
  
  # Test for equality of variances
  levene_test <- leveneTest(
    c(normalized_counts$CTRL_avg_log, 
      normalized_counts$NICD_avg_log),
    group = rep(c("CTRL", "NICD"), 
                c(length(normalized_counts$CTRL_avg_log), 
                  length(normalized_counts$NICD_avg_log)))
  )
  
  # Choose appropriate test based on normality and variance tests
  if (shapiro_ctrl$p.value > 0.05 && shapiro_nicd$p.value > 0.05) { # Normal
    if (levene_test$p.value > 0.05) { # Equal variance
      test_result <- t.test(normalized_counts$CTRL_avg_log, 
                            normalized_counts$NICD_avg_log, 
                            var.equal = TRUE) # Student's t-test
      test_used <- "Student's t-test"
    } else { # Unequal variance
      test_result <- t.test(normalized_counts$CTRL_avg_log, 
                            normalized_counts$NICD_avg_log, 
                            var.equal = FALSE) # Welch's t-test
      test_used <- "Welch's t-test"
    }
  } else { # Non-normal
    test_result <- wilcox.test(normalized_counts$CTRL_avg_log, 
                               normalized_counts$NICD_avg_log) # Wilcoxon test
    test_used <- "Wilcoxon rank-sum test"
  }
  
  mean_ctrl_log <- mean(normalized_counts$CTRL_avg_log)
  mean_nicd_log <- mean(normalized_counts$NICD_avg_log)
  
  # save stats result in this object
  stats_test_result <- list(
    test_used = test_used,
    p_value = test_result$p.value,
    mean_ctrl_log = mean_ctrl_log,
    mean_nicd_log = mean_nicd_log
  )
  
  return(stats_test_result)
}

# Plot function adjusted for log-transformed data ******************************
norm_accessibility_plot_function <- function(avg_counts, stats_test_result, plot_title) {
  # Create a label for the statistical test used stored in first item
  test_label <- paste("Test used:", stats_test_result[1])
  
  plot <- ggplot(avg_counts, aes(x = Condition, 
                                 y = Average_Read_Counts, 
                                 fill = Condition)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, 
                 position = position_dodge(width = 0.75), 
                 fill = "white") +
    scale_fill_manual(values = c("#bdb6d8", "#a3678c")) +
    labs(title = plot_title,
         x = "Condition",
         y = "log2(Normalized Read Counts + 1)") +
    theme(
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),  
      panel.background = element_blank(),
      axis.line.x = element_line(color = "black"), 
      axis.line.y = element_line(color = "black"),
      legend.position = "none"
    ) +
    annotate("text", 
             x = 1.5, 
             y = max(avg_counts$Average_Read_Counts) * 1.1,  
             label = paste("p-value:", 
                           signif(stats_test_result$p_value, digits = 3), 
                           "\n", test_label))
  
  return(plot)
}


# PROGENITOR NORM ACCESSIBILITY: -----------------------------------------------
norm_counts <- normalized_counts_function(merged_df, progenitor_all)
avg_counts <- avg_read_counts_function(norm_counts)
stats_test_result <- stats_test_function(norm_counts)

progenitor_plot <- norm_accessibility_plot_function(
  avg_counts, 
  stats_test_result, 
  "Normalized Accessibility of Progenitor Genes")

filename = 
  "figures/fig4/fig4a_normcountsplot_progenitor.pdf" 
pdf(filename, width = 5, height = 5)
print(progenitor_plot) 

dev.off()


# NEUROGENIC NORM ACCESSIBILITY: -----------------------------------------------
norm_counts <- normalized_counts_function(merged_df, neurogenic_all) 
avg_counts <- avg_read_counts_function(norm_counts)
stats_test_result <- stats_test_function(norm_counts)

neurogenic_plot <- norm_accessibility_plot_function(
  avg_counts, 
  stats_test_result, 
  "Normalized Accessibility of Neurogenic Genes") 

filename = 
  "figures/fig4/fig4b_normcountsplot_neurogenic.pdf"
pdf(filename, width = 5, height = 5)
print(neurogenic_plot) 

dev.off()


# GLIOGENIC NORM ACCESSIBILITY: ------------------------------------------------
# When calling function, change gene list name
norm_counts <- normalized_counts_function(merged_df, gliogenic_all) 
avg_counts <- avg_read_counts_function(norm_counts)
stats_test_result <- stats_test_function(norm_counts)


# change plot name
# change title
gliogenic_plot <- norm_accessibility_plot_function(
  avg_counts, 
  stats_test_result, 
  "Normalized Accessibility of Gliogenic Genes") 

# change filename
filename = 
  "figures/fig4/fig4c_normcountsplot_gliogenic.pdf" 
pdf(filename, width = 5, height = 5)
print(gliogenic_plot)

dev.off()


# OL LINEAGE NORM ACCESSIBILITY: ---------------------------------------------------------------------------------------
norm_counts <- normalized_counts_function(merged_df, opc_all)
avg_counts <- avg_read_counts_function(norm_counts)
stats_test_result <- stats_test_function(norm_counts)


# change plot name
opc_plot <- norm_accessibility_plot_function(
  avg_counts, 
  stats_test_result, 
  "Normalized Accessibility of OL Lineage Genes")

filename = 
  "figures/fig4/fig4d_normcountsplot_ol.pdf"
pdf(filename, width = 5, height = 5)
print(opc_plot)

dev.off()


# SHH PATHWAY NORM ACCESSIBILITY: ----------------------------------------------
norm_counts <- normalized_counts_function(merged_df, Shh_pathway_genes)
avg_counts <- avg_read_counts_function(norm_counts)
stats_test_result <- stats_test_function(norm_counts)


shh_plot <- norm_accessibility_plot_function(
  avg_counts, 
  stats_test_result, 
  "Normalized Accessibility of SHH Pathway Genes")

filename = 
  "figures/fig4/fig4g_normcountsplot_shh.pdf"
pdf(filename, width = 5, height = 5)
print(shh_plot)

dev.off()


# ------------------------------------------------------------------------------
print(
  "Normalized accessibility plots were generated and saved in figures/fig4"
  )