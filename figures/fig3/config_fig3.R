# Run 'analysis/02_diffbind_e16.R' and 'tables/scripts/tableS3.R' if you haven't already
# source("analysis/02_diffbind_e16.R")
# source("tables/scripts/tableS3.R")

# For tss enrichmentplot:
# Run scripts in 'preprocess/atacseq_e16_compute_matrix/'
# and put the matrix files in 'processed_data/atacseq_e16/deeptools_output'

# DEFINE FILES AND OUTPUT DIRS: ------------------------------------------------
# input DBA object rds file after running diffbind analysis
rds_dbObj <- "data/processed_data/atacseq_e16/r_objects/diffbind_dbObj.rds"

# Diffbind results csv file
csv_diffbind_results <- "tables/table_S3_atacseq_e16_diffbind_results.csv"

# Read in matrix files from deeptools

print("Reading in matrix files...")
# Unzip gz files and skip the JSON metadata line
ctrl_matrix <- fread(
  "zcat < data/processed_data/atacseq_e16/deeptools_output/computematrix/CTRL.matrix.gz", 
  skip = 1)

nicd_matrix <- fread(
  "zcat < data/processed_data/atacseq_e16/deeptools_output/computematrix/NICD.matrix.gz",
  skip = 1) 

print("Matrix files unzipped")

# output directory for figure 3
output_dir_figures <- "figures/fig3"

# output figure names
filename_pca <- "fig3b_pcaplot.pdf"
filename_heatmap <- "fig3d_heatmap.pdf"
filename_logfoldchange <- "fig3i_log2foldchangeplot.pdf"
filename_tssenrichment <- "fig3c_tssenrichmentplot.pdf"

# Plot colors
control_color <- "#AFA2C2"
mutant_color <- "#833D66"

# HEATMAP; DIFF ACCESSIBLE SITES: ----------------------------------------------
# Choose heatmap colors
hmap <- colorRampPalette(c("white", "#938fc2", "#40027e")) (n =100)


# LOG2 FOLD CHANGE PLOT: -------------------------------------------------------
# Choose barplot colors
bar_colors_midnights <- c("#121D27", "#5A658B", "#6F86A2", "#85A7BA", "#AA9EB6")

# make lists of selected gene names
selected_neurogenic <- c('Neurod2','Eomes')
selected_progenitor <- c('Sox2', 'Pax6', 'Hes1', 'Hes5')
selected_gliogenic <- c("Olig2", "Egfr", "Ascl1")

desired_order <- c("Hes1", "Hes5", "Sox2", "Pax6", "Olig2", "Egfr", "Ascl1", "Neurod2", "Eomes")

print("Config file for fig3 ran successfully.")
