# Run 'analysis/01_deseq2_e16.R' if you haven't already
#source("analysis/01_deseq2_e16.R")

# Run 'tables/scripts/tableS1.R' if you haven't already
# source("tables/scripts/tableS1.R")

# DEFINE FILES AND OUTPUT DIRS: ------------------------------------------------
# input rds file containing deseq2 results
rds_deseq2_results <- "data/processed_data/rnaseq_e16/r_objects/deseq2_dds_e16.rds"

# input csv file containing deseq2 results
csv_deseq2_results <- "tables/table_S1_rnaseq_e16_deseq2_results.csv"

# output directory for figure 2
output_dir_figures <- "figures/fig2"

# output figure names
filename_pca <- "fig2b_pcaplot.pdf"
filename_volcanoplot <- "fig2c_volcanoplot.pdf"
filename_heatmap <- "fig2d_heatmap.pdf"

# Plot colors
control_color = "#76BAE0"
mutant_color = "#B8396B"

# VOLCANO PLOT; DEFINE GENE LIST AND THRESHOLD: --------------------------------
# Make a list of selected genes. Here, I want to show Notch genes, 
# a few neurogenic genes, and a few progenitor and gliogenic genes.
specific_genes <- c('Dll3', 'Dll1', 'Neurod4', 'Eomes', 'Hey1', 
                    'Neurog2', 'Hes1', 'Notch1', 
                    'Nes', 'Sox2', 'Sox9')

# define significance threshold (padj 0.05)
alpha <- 0.05


# HEATMAP; DEFINE GENE LISTS AND ORDER: ----------------------------------------
progenitor <- c('Nes', 'Sox2', 'Vim', 'Nr2e1', 'Hes1', 'Ednrb', 'Dleu7', 
                'Ncald','Rfx4', 'Bcan')

glia <- c('Ccnd1', 'Aldh1l1', 'Cntn1', 'Olig1', 'Olig2', 'Pdgfra', 'Sox10')

neuron <- c('Eomes', 'Neurod1', 'Neurog2', 'Neurog1', 'Btg2', 'Tbr1', 'Satb2',
            'Neurod4')


# Combine progenitor, neurogenic, and gliogenic gene lists for heatmap
combined_gene_list <- c(progenitor, glia, neuron)

# define order of genes on heatmap
desired_order <- c(progenitor, glia, neuron)

print("Config file for fig2 ran successfully.")

