# Run after running 'analysis/02_diffbind_e16.R' and 'tables/scripts/tableS3.R'

library(readr)
library(tidyverse)
library(dplyr)

# READ IN DIFFBIND RESULTS: ----------------------------------------------------
diffbind_res_df <- read.csv("tables/table_S3_atacseq_e16_diffbind_results.csv")

# Create dataframes for CTRL (Fold < 0) keeping only significant peaks 
# (FDR < 0.05)
CTRL_enrich <- diffbind_res_df %>% 
  dplyr::filter(FDR < 0.05 & Fold < 0) %>% 
  dplyr::select(seqnames, start, end, peak_id)

# Create dataframes for NICD keeping only significant peaks (FDR < 0.05)
NICD_enrich <- diffbind_res_df %>% 
  dplyr::filter(FDR < 0.05 & Fold > 0) %>% 
  dplyr::select(seqnames, start, end, peak_id)

# Write to bed files
write.table(
  CTRL_enrich, 
  file="data/processed_data/atacseq_e16/diffbind_results_bed/CTRL_e16_enriched.bed", 
  sep="\t", quote=F, row.names=F, col.names=F)

write.table(
  NICD_enrich, 
  file="data/processed_data/atacseq_e16/diffbind_results_bed/NICD_e16_enriched.bed", 
  sep="\t", quote=F, row.names=F, col.names=F)

print("Bed files saved in data/processed_data/atacseq_e16/diffbind_results_bed")
