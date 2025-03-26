# Run after running 'analysis/01_deseq2_e16.R' 

library(org.Mm.eg.db)
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(readr)
library(dplyr)
library(pheatmap)

# READ DDS OBJECT FOLLOWING DESEQ2 ANALYSIS FROM RDS FILE: --------------------
dds <- readRDS("data/processed_data/rnaseq_e16/r_objects/deseq2_dds_e16.rds")
# Store results in res
res <- results(dds)

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
neurogenic_all <- c(IPC_genes, neurogenic_genes, neuronal_genes, 
                    newborn_neurons)

gliogenic_all <- c(OPC_genes, preOPC_genes, astrocyte_genes, glial_lineage, OL_genes)

progenitor_all <- c(progenitor_genes, RGC_genes, proliferative_genes)

ol_all <- c(preOPC_genes, OPC_genes, OL_genes)

# ANNOTATE DESEQ2 RESULTS WITH GENE SYMBOLS AND ENTREZ IDS: --------------------
ensembl_ids <- rownames(res)

# annotate with gene symobols using org.Mm.eg.db package
res$symbol <- mapIds(org.Mm.eg.db,
                     keys=ensembl_ids,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
# annotate with entrez ids
res$entrez <- mapIds(org.Mm.eg.db,
                     keys=ensembl_ids,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

# Re-order res based on padj
resOrdered <- res[order(res$padj),]

# Move symbol and entrez columns to the front
resOrdered <- resOrdered[, c("symbol", 
                             "entrez", 
                             setdiff(colnames(resOrdered), 
                                     c("symbol", "entrez")))
                         ]


# PREPARE DATA FOR HEATMAP SHOWING SPECIFIC GENE SET: --------------------------

# Combine progenitor, neurogenic, and gliogenic gene lists for heatmap
combined_gene_list <- c(progenitor_all, gliogenic_all, neurogenic_all)

# Remove NAs from padj before filtering for significant genes with padj < 0.05
sig_res <- resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.05, ]

# Extract rownames/ensembl id's of significant DEGs
sig_genes <- rownames(sig_res)

# Map ensembl id's to symbols
sig_symbols <- mapIds(org.Mm.eg.db,
                      keys = sig_genes,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")

# Filter sig_genes list for genes in combined_gene list by symbol
selected_genes <- sig_genes[sig_symbols %in% combined_gene_list]

# EXTRACT TRANSFORMED VALUES: **************************************************
vsd <- vst(dds, blind=FALSE)

# Subset assay matrix for the selected genes
mat <- assay(vsd)[selected_genes, ]

# Update row names with gene symbols for the selected genes
rownames(mat) <- sig_symbols[sig_symbols %in% combined_gene_list]

# Center the data
mat <- mat - rowMeans(mat)

# Make Heatmap in specific order ***********************************************
progenitor <- c('Nes', 'Sox2', 'Vim', 'Nr2e1', 'Hes1', 'Ednrb', 'Dleu7', 
                'Ncald','Rfx4', 'Bcan')
progenitor <- c('Hes1', 'Dleu7', 'Rfx4', 'Sox2', 'Bcan', 'Ncald', 'Nes', 'Ednrb', 'Vim', 'Nr2e1')
glia <- c('Ccnd1', 'Aldh1l1', 'Cntn1', 'Olig1', 'Olig2', 'Pdgfra', 'Sox10', 'Sox8', 'Sox9', 'Aqp4','Mbp')
glia <- c('Sox8', 'Sox9', 'Aqp4', 'Olig2', 'Ccnd1', 'Aldh1l1', 'Cntn1', 'Sox10', 'Pdgfra', 'Olig2', 'Mbp')
neuron <- c('Eomes', 'Neurod1', 'Neurog2', 'Neurog1', 'Btg2', 'Tbr1', 'Satb2',
            'Neurod4')

desired_order <- c(progenitor, glia, neuron)

# Only include genes present in the matrix
ordered_genes <- desired_order[desired_order %in% rownames(mat)]

# Subset the matrix to include only the ordered genes
mat_ordered <- mat[match(ordered_genes, rownames(mat)), ]

# MAKE HEATMAP WITH ORDERED GENES: ---------------------------------------------
ordered_combined_gene_list_heatmap <- pheatmap(
  mat_ordered,
  # Set to FALSE to use specified order
  cluster_rows = FALSE,  
  cluster_cols = TRUE, 
  # Show gene symbols
  show_rownames = TRUE,  
  annotation_col = as.data.frame(
    colData(vsd)[, "condition", drop=FALSE]),
  color = colorRampPalette(c("#76BAE0", "white", "#B8396B"))(50))

filename = "figures/fig2/fig2d_heatmap_ordered_032625.pdf"
pdf(filename, width = 6, height = 8)
print(ordered_combined_gene_list_heatmap)

dev.off()

print("Heatmap with specific genes generated and saved in figures/fig2/")