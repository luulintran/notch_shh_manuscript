# Run 'analysis/04_deseq2_e13-vs-e17.R' if you haven't already
#source("analysis/04_deseq2_e13-vs-e17.R")

# Run 'tables/scripts/tableS8.R' if you haven't already
# source("tables/scripts/tableS8.R")

# DEFINE FILES AND OUTPUT DIRS: ------------------------------------------------
# input deseq2 results rds file
rds_deseq2_results <- "data/processed_data/rnaseq_ctrl_e13/r_objects/deseq2_dds_e13-vs-e17.rds"

# input csv file containing deseq2 results
csv_deseq2_results <- "tables/table_S8_rnaseq_e13-vs-e17_deseq2_results.csv"

# output directory for figure S3
output_dir_figures <- "figures/figS3"

# output figure names
filename_pca <- "figS3a_pcaplot.pdf"
filename_volcanoplot <- "figS3b_volcanoplot.pdf"
filename_heatmap <- "figS3c_heatmap.pdf"

# PCA plot colors
control_color <- "#F4C430"
mutant_color <- "#C6A0F6"

# Volcano plot colors
downreg_color <- "#97cb5e"
upreg_color <- "#d14a96"

# heatmap colors
piyg_colors <- colorRampPalette(rev(brewer.pal(11, "PiYG")))(100)

# VOLCANO PLOT; DEFINE GENE LIST AND THRESHOLD: --------------------------------
# Make a list of selected genes. Here, I want to show Notch genes, 
# a few neurogenic genes, and a few progenitor and gliogenic genes.
specific_genes <- c('Dll3', 'Dll1', 'Neurod4', 'Hey1', 
                    'Neurog2', 'Hey1', 'Hes5', 'Notch1', 
                    'Sox2', 'Sox9', 'Notch2')

# define significance threshold (padj 0.05)
alpha <- 0.05


# HEATMAP; DEFINE GENE LISTS AND ORDER: ----------------------------------------
heatmap_genes_list <- c("Notch1", "Notch2", "Hes1", "Hes5", "Hey1", "Hes6",
                "Dll1", "Dll3", "Jag1", "Jag2")

print("Config file for figS3 ran successfully.")

