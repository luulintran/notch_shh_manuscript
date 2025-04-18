# Run 'analysis/02_diffbind.R' and 'tables/scripts/tableS3.R' if you haven't already
#source("analysis/02_diffbind_e16.R")
#source("tables/scripts/tableS3.R")

# DEFINE FILES AND OUTPUT DIRS: ------------------------------------------------
# input rds file containing diffbind norm read counts results
rds_normcounts <- "data/processed_data/atacseq_e16/r_objects/norm_read_counts_dbObj.rds"

# input diffbind results as table
csv_diffbind_results <- "tables/table_S3_atacseq_e16_diffbind_results.csv"

# output directory for figure 4
output_dir_figures <- "figures/fig4"

# output filename for log2 fold change plot
filename_logfoldchange <- "fig4h_log2foldchangeplot_shh.pdf"

# output figure name for each norm read counts plot
filename_progenitor <- "fig4a_normcountsplot_progenitor_041825.pdf"
filename_neuro <- "fig4b_normcountsplot_neurogenic_041825.pdf"
filename_glia <- "fig4c_normcountsplot_gliogenic_041825.pdf" 
filename_ol <- "fig4d_normcountsplot_ol_041825.pdf"
filename_shh <- "fig4g_normcountsplot_shh_041825.pdf"

# Plot colors
control_color = "#bdb6d8"
mutant_color = "#a3678c"


# NORM READ COUNTS PLOTS; MAKE GENE LISTS: -------------------------------------
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

Shh_pathway_genes <- c('Gli1', 'Smo', 'Ptch1', 'Boc', 'Cdon', 'Gas1', 
                       'Gli2','Ptch2', 'Hhip', 'Gli3')

#  COMBINE GENE LISTS TO MAKE MORE BROAD CATAGORIES: ---------------------------
# combine gene lists and then remove any repeats/ duplicated genes with unique()
neurogenic_all <- c(IPC_genes, neurogenic_genes, neuronal_genes, 
                    newborn_neurons)
neurogenic_all <- unique(neurogenic_all)

gliogenic_all <- c(OPC_genes, preOPC_genes, astrocyte_genes, glial_lineage, OL_genes)
gliogenic_all <- unique(gliogenic_all)

progenitor_all <- c(progenitor_genes, RGC_genes, proliferative_genes)
progenitor_all <- unique(progenitor_all)

ol_all <- c(preOPC_genes, OPC_genes, OL_genes)
ol_all <- unique(ol_all)

# LOG2 FOLD CHANGE PLOT; SELECTED SHH PATHWAY GENES: ---------------------------
# Selected SHH genes
selected_shh_pos <- c('Smo', 'Gli1', 'Gli2', 'Gas1', 'Cdon', 'Boc')
selected_shh_neg <- c('Ptch1', 'Ptch2', 'Gli3', 'Hhip')

# Bar plot colors
bar_colors_midnights <- c("#121D27", "#5A658B", "#6F86A2", "#85A7BA", "#AA9EB6")

print("Config file for fig4 ran successfully.")
