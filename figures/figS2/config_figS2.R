# Run 'analysis/03_deseq2_e17.R' and 'tables/scripts/tableS6.R' if you haven't already
#source("analysis/03_deseq2_e17.R")
#source("tables/scripts/tableS6.R")

# Run 'analysis/01_deseq2_e16.R' if you haven't already
#source("analysis/01_deseq2_e16.R")

# DEFINE FILES AND PATHS: ------------------------------------------------------
# input deseq2 results rds file
rds_deseq2_results_e17 <- "data/processed_data/rnaseq_e17/r_objects/deseq2_dds_e17.rds"
rds_deseq2_results_e16 <- "data/processed_data/rnaseq_e16/r_objects/deseq2_dds_e16.rds"


# input csv file containing deseq2 results
csv_deseq2_results <- "tables/table_S7_rnaseq_e17_deseq2_results.csv"

# output directory for figure S2
output_dir_figures <- "figures/figS2"

# output filename
filename_pca <-"figS2a_pcaplot.pdf"
filename_volcanoplot <- "figS2b_volcanoplot.pdf"
filename_upregulatedgo <- "figS2c_upregulatedgo_barplot.pdf"
filename_heatmap_e16 <- "figS2d_heatmap_tf_e16.pdf"
filename_heatmap_e17 <- "figS2d_heatmap_tf_e17.pdf"

# PCA plot colors
control_color <-  "#008080"
mutant_color <- "#dc9b06"

# Heatmap and volcano plot colors
downreg_color <- "#CDA2CC"
upreg_color <- "#FF6B35"

# GO term plot colors
upreg_go_color <- "#f48c67"

# VOLCANO PLOT; SELECTED GENE LIST AND THRESHOLD: ------------------------------
# Make a list of selected genes. Here, I want to show Notch genes, 
# a few neurogenic genes, and a few progenitor and gliogenic genes.
specific_genes <- c('Dll3', 'Dll1', 'Neurod4', 'Hey1', 
                    'Neurog2', 'Hey1', 'Hes5', 'Notch1', 
                    'Sox2', 'Sox9', 'Notch2')

# define significance threshold (padj 0.05)
alpha <- 0.05

# UPREGULATED GO PLOT; GO TERMS LIST: ------------------------------------------
GO_up_terms_list <- c('smoothened signaling pathway', 
                      'cilium assembly', 'cilium organization')

# HEATMAP; TRANSCRIPTION FACTORS LIST ------------------------------------------
# These TFs have motifs enriched in sites with decreased or 
# increased accessibility

# Heatmap with specific genes of TFs from atac-seq motif enrichment 
TFs <- c('Lhx1','Lhx2', 'Nfix', 'Nfia', 'Nfic', 'Sox2', 'Sox8', 'Sox9','Sox10', 
         'Neurog1', 'Neurog2', 'Neurod2', 'Mef2c', 'Mef2d', 'Eomes', 'Tbr1')

# ******************************************************************************
desired_order <- c('Lhx1','Lhx2', 'Nfix', 'Nfia', 'Nfic', 'Sox2', 'Sox8', 
                   'Sox9','Sox10', 'Neurog1', 'Neurog2', 'Neurod2', 'Mef2c', 
                   'Mef2d', 'Eomes', 'Tbr1')

print("Config file for figS2 ran successfully.")
