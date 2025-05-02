# READ DIFFBIND RESULTS FILE WITH GENE ANNOTATIONS AND STATS RESULTS: ----------
diffbind_res_df <- read.csv(csv_diffbind_results)

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

# SUMMARIZE DATA FOR SHH GENES: ------------------------------------------------
#shh positive CTRL
selected_CTRL_fold_shh_pos <- CTRL_enriched_fold %>%
  dplyr:: filter(SYMBOL %in% selected_shh_pos)

selected_CTRL_fold_shh_pos_datasumm <- selected_CTRL_fold_shh_pos %>%
  group_by(SYMBOL) %>%
  summarise(mean_fold = mean(Fold),
            sd_fold = sd(Fold),
            n_fold = n(),
            SE_fold = sd(Fold)/ sqrt(n()))


#shh negative CTRL
selected_CTRL_fold_shh_neg <- CTRL_enriched_fold %>%
  dplyr:: filter(SYMBOL %in% selected_shh_neg)

selected_CTRL_fold_shh_neg_datasumm <- selected_CTRL_fold_shh_neg %>%
  group_by(SYMBOL) %>%
  summarise(mean_fold = mean(Fold),
            sd_fold = sd(Fold),
            n_fold = n(),
            SE_fold = sd(Fold)/ sqrt(n()))


#shh positive NICD
selected_NICD_fold_shh_pos <- NICD_enriched_fold %>%
  dplyr:: filter(SYMBOL %in% selected_shh_pos)

selected_NICD_fold_shh_pos_datasumm <- selected_NICD_fold_shh_pos %>%
  group_by(SYMBOL) %>%
  summarise(mean_fold = mean(Fold),
            sd_fold = sd(Fold),
            n_fold = n(),
            SE_fold = sd(Fold)/ sqrt(n()))

#shh negative NICD
selected_NICD_fold_shh_neg <- NICD_enriched_fold %>%
  dplyr:: filter(SYMBOL %in% selected_shh_neg)

selected_NICD_fold_shh_neg_datasumm <- selected_NICD_fold_shh_neg %>%
  group_by(SYMBOL) %>%
  summarise(mean_fold = mean(Fold),
            sd_fold = sd(Fold),
            n_fold = n(),
            SE_fold = sd(Fold)/ sqrt(n()))

# Combine dfs: -----------------------------------------------------------------
selected_df <- rbind(selected_CTRL_fold_shh_pos, 
                     selected_CTRL_fold_shh_neg, 
                     selected_NICD_fold_shh_pos, 
                     selected_NICD_fold_shh_neg)
selected_df_datasumm <- rbind(selected_CTRL_fold_shh_pos_datasumm, 
                              selected_CTRL_fold_shh_neg_datasumm, 
                              selected_NICD_fold_shh_pos_datasumm, 
                              selected_NICD_fold_shh_neg_datasumm)

# Change order of genes
selected_df$SYMBOL <- factor(selected_df$SYMBOL, levels = c('Ptch1', 'Boc',
                                                            'Gli3', 'Gli2',
                                                            'Gas1', 'Cdon',
                                                            'Gli1', 'Smo',
                                                            'Hhip'))

# Bar plot colors: -------------------------------------------------------------
#Repeat the colors to match the number of unique genes
num_genes <- length(unique(selected_df$SYMBOL))
bar_colors_cycle <- rep(bar_colors_midnights, length.out = num_genes)

# PLOT SHH LOG2FOLD CHANGE: ----------------------------------------------------
shh_log2foldchange_plot <- ggplot(selected_df_datasumm, aes(x= SYMBOL, 
                                                            y= mean_fold)) + 
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
             dotsize=0.8, 
             show.legend = FALSE) +
  scale_fill_manual(values = bar_colors_cycle) +
  ylab("log2 Fold Change") +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 1, hjust =1)
  ) +
  geom_hline(yintercept = 0, color = "black") + 
  theme(legend.position = "none")

# SAVE PLOT: -------------------------------------------------------------------
pdf(file.path(output_dir_figures, filename_logfoldchange), width = 5, height = 5)
print(shh_log2foldchange_plot)

dev.off()

print(
  paste0(
    filename_logfoldchange, 
    " generated and saved in ",output_dir_figures
    )
  )
