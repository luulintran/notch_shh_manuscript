# Run after running 'analysis/02_diffbind_e16.R' and 'tables/scripts/tableS3.R'

# SET UP
library(tidyverse)
library(dplyr)
library(readr)
library(ggplot2)

# READ DIFFBIND RESULTS FILE WITH GENE ANNOTATIONS AND STATS RESULTS: ----------
diffbind_res_df <- read.csv("tables/table_S3_atacseq_e16_diffbind_results.csv")

# SEPARATE BETWEEN CTRL-ENRICHED (Fold < 0) and NICD-ENRICHED (Fold > 0): ------
# This dataframe has everything including the stats info
CTRL_enrich <- diffbind_res_df %>% 
  dplyr::filter(FDR < 0.05 & Fold < 0)

# This dataframe has everything including the stats info
NICD_enrich <- diffbind_res_df %>% 
  dplyr::filter(FDR < 0.05 & Fold > 0)

# PREPARE DATA FOR PLOTTING FOLD CHANGE: ---------------------------------------

# Extract fold change values and group by gene symbol
#CTRL
CTRL_enriched_fold <- CTRL_enrich[, c("SYMBOL", "Fold")]

CTRL_enriched_gene_logfold <- CTRL_enriched_fold %>%
  group_by(SYMBOL) %>%
  #get the mean of log fold changes
  summarise(logFoldChange = mean(Fold, na.rm = TRUE)) 


#NICD
NICD_enriched_fold <- NICD_enrich[, c("SYMBOL", "Fold")]

NICD_enriched_gene_logfold <- NICD_enriched_fold %>%
  group_by(SYMBOL) %>%
  summarise(logFoldChange = mean(Fold, na.rm = TRUE))

# MAKE GENE LISTS: -------------------------------------------------------------
# make lists of selected gene names
selected_neurogenic <- c('Neurod2','Eomes')

selected_progenitor <- c('Sox2', 'Pax6', 'Hes1', 'Hes5')

selected_gliogenic <- c("Olig2", "Egfr", "Ascl1")

# Filter lists for the selected neurogenic genes in CTRL: ----------------------
selected_CTRL_fold_neurogenic <- CTRL_enriched_fold %>%
  dplyr:: filter(SYMBOL %in% selected_neurogenic)

selected_CTRL_fold_neurogenic_datasumm <- selected_CTRL_fold_neurogenic %>%
  group_by(SYMBOL) %>%
  summarise(mean_fold = mean(Fold),
            sd_fold = sd(Fold),
            n_fold = n(),
            SE_fold = sd(Fold)/ sqrt(n()))


# Filter lists for the selected progenitor genes in NICD: ----------------------
selected_NICD_fold_progenitor <- NICD_enriched_fold %>%
  dplyr:: filter(SYMBOL %in% selected_progenitor)

selected_NICD_fold_progenitor_datasumm <- selected_NICD_fold_progenitor %>%
  group_by(SYMBOL) %>%
  summarise(mean_fold = mean(Fold),
            sd_fold = sd(Fold),
            n_fold = n(),
            SE_fold = sd(Fold)/ sqrt(n()))


# Filter lists for the selected gliogenic genes in NICD: -----------------------
selected_NICD_fold_gliogenic <- NICD_enriched_fold %>%
  dplyr:: filter(SYMBOL %in% selected_gliogenic)

selected_NICD_fold_gliogenic_datasumm <- selected_NICD_fold_gliogenic %>%
  group_by(SYMBOL) %>%
  summarise(mean_fold = mean(Fold),
            sd_fold = sd(Fold),
            n_fold = n(),
            SE_fold = sd(Fold)/ sqrt(n()))
head(selected_NICD_fold_gliogenic_datasumm)

# Combine filtered dataframes: -------------------------------------------------
selected_df <- rbind(selected_CTRL_fold_neurogenic, 
                     selected_NICD_fold_progenitor, 
                     selected_NICD_fold_gliogenic)
selected_df_datasumm <- rbind(selected_CTRL_fold_neurogenic_datasumm, 
                              selected_NICD_fold_progenitor_datasumm, 
                              selected_NICD_fold_gliogenic_datasumm)
selected_df$SYMBOL <- factor(selected_df$SYMBOL, 
                             levels = c("Hes1", "Hes5", "Sox2", "Pax6", 
                                        "Olig2", "Egfr", "Ascl1", "Neurod2", 
                                        "Eomes"))

# Choose barplot colors: -------------------------------------------------------
bar_colors_midnights <- c("#121D27", "#5A658B", "#6F86A2", "#85A7BA", "#AA9EB6")

#Repeat the colors to match the number of unique genes
num_genes <- length(unique(selected_df$SYMBOL))
bar_colors_cycle <- rep(bar_colors_midnights, length.out = num_genes)

# PLOT LOG2 FOLD CHANGE: -------------------------------------------------------
cellfategenes_log2foldchange_plot <- ggplot(
  selected_df_datasumm, aes(x= SYMBOL, y= mean_fold)) + 
  geom_col(width= 0.8, 
           aes(fill= SYMBOL),
           color= 'black',
           size= 0.5) + 
  geom_errorbar(aes(ymin = mean_fold - SE_fold, 
                    ymax = mean_fold + SE_fold), 
                width=0.3,
                size= 0.5) +
  geom_point(data= selected_df, 
             aes(x=SYMBOL, 
                 y=Fold,
                 fill= SYMBOL),
             binaxis='y', 
             stackdir='center', 
             dotsize=0.25, 
             show.legend = FALSE) +
  scale_fill_manual(values = bar_colors_cycle) +
  ylab("log2 Fold Change") +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color = "black"), 
    axis.line.y = element_line(color = "black"),  
    axis.text.x = element_text(angle = 90, vjust = 1, hjust =1) #add line at 0
  ) +
  geom_hline(yintercept = 0, color = "black") + 
  theme(legend.position = "none")

filename = "figures/fig3/fig3i_log2foldchangeplot.pdf"
pdf(filename, width = 3, height = 3)
print(cellfategenes_log2foldchange_plot)

dev.off()

print("Log2 fold change plot generated and saved in figures/fig3")