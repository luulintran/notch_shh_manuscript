# Run 'analysis/03_deseq2_e17.R' and 'tables/scripts/tableS6.R' if you haven't already
#source("analysis/03_deseq2_e17.R")
#source("tables/scripts/tableS6.R")

# DEFINE FILES AND OUTPUT DIRS: ------------------------------------------------
# input rds file containing deseq2 results
rds_deseq2_results <- "data/processed_data/rnaseq_e17/r_objects/deseq2_dds_e17.rds"

# input csv file containing deseq2 results
csv_deseq2_results <- "tables/table_S7_rnaseq_e17_deseq2_results.csv"

# output directory for figure 5
output_dir_figures <-  "figures/fig5"

# output figure names
filename_upregulatedgo <- "fig5f_upregulatedgo_barplot_smo.pdf"
filename_downregulatedgo <- "fig5g_downregulatedgo_barplot.pdf"
filename_heatmap <- "fig5h_heatmap.pdf"

# Plot colors
# (downregulated)
control_color = "#ca8bc9"

# (upregulated)
mutant_color = "#f48c67"

# Heatmap colors
downreg_color <- "#ca8bc9"
upreg_color <- "#FF6B35"

# GO TERMS PLOTS: --------------------------------------------------------------
# Define the GO terms of interest
GO_up_specific_terms <- c('smoothened signaling pathway',
                    'gliogenesis', 
                    'glial cell differentiation', 
                    'glial cell development', 
                    'oligodendrocyte development')

# Define the GO terms of interest
GO_down_specific_terms <- c('regulation of neurogenesis', 
                     'dendrite development', 
                     'positive regulation of neurogenesis', 
                     'regulation of synapse structure or activity', 
                     'regulation of neuron differentiation')

# HEATMAP; SELECTED SHH GENES: -------------------------------------------------
# MAKE LIST OF SELECTED SHH PATHWAY GENES: -------------------------------------
Shh_gene_list <- c("Smo", "Boc", "Cdon", 
                   "Gas1", "Gli1", "Gli2", "Gli3", 
                   "Sufu", "Disp1", "Iqce", "Efcab7", 
                   "Ptch1", "Ptch2", "Hhip", "Sufu", 
                   "Gsk3b", "Ck1", "Pcaf", "Cul3", "Kif7", "Evc", "Evc2")

desired_order <- c("Evc2", "Kif7", "Gas1", "Boc", "Smo", 
                   "Iqce", "Evc", "Gli1", "Ptch1", "Ptch2", "Hhip")

print("Config file for fig5 ran successfully.")