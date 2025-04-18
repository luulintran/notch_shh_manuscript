# DEFINE FILES AND PATHS: ------------------------------------------------------
# Path to the gene raw counts matrix
input_cts <- "data/processed_data/rnaseq_e17/raw_counts/rnaseq_e17_raw_counts.tsv"

# Path to the sample sheet file
input_coldata <- "data/meta_data/rnaseq_e17/coldata.csv"

# Path to output directory for r_objects
output_dir_robj <- "data/processed_data/rnaseq_e17/r_objects"

# filename for RDS file
filename <- "deseq2_dds_e17.rds"

# Make sure output directory exists, and if it doesn't, create one
if (!dir.exists(output_dir_robj)) {
  dir.create(output_dir_robj, recursive = TRUE)
}

# READ IN DATA:-----------------------------------------------------------------
# raw gene counts data from star_salmon; remove gene_name column
cts <- read.csv(input_cts, 
                sep= "\t", 
                row.names = 1) # makes the rownames the first col
cts <- as.data.frame(cts[, -c(1)]) # remove the gene_name col

# sample metadata
coldata <- read.csv(input_coldata, 
                    row.names = 1)

# CREATE DESEQ2 DATASET:--------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = round(cts), 
                              colData = coldata, 
                              design = ~ condition)

dds$condition <- relevel(dds$condition, 
                         ref = "control")

# RUN DESEQ2 ANALYSIS: ---------------------------------------------------------
dds <- DESeq(dds)

# SAVE DESEQ2 OBJECT TO RDS FOR LATER: -----------------------------------------
saveRDS(dds, 
        file = file.path(output_dir_robj, 
                         filename)
)

print(paste0("DESeq2 analysis done. RDS file saved in ", output_dir_robj))
