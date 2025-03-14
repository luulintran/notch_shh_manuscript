# setup_environment.R

# Load renv for package management
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
library(renv)

# Restore the environment from the renv.lock file
if (file.exists("renv.lock")) {
  cat("\nRestoring environment from renv.lock...\n")
  renv::restore()
} else {
  stop("renv.lock file not found. Make sure it exists in the project directory.")
}

# Set up directory structure
required_directories <- c(
  "data/raw_data",
  "data/processed_data",
  "data/meta_data",
  "data/motif_files",
  "data/processed_data/atacseq_e16",
  "data/processed_data/rnaseq_e16",
  "data/processed_data/rnaseq_e17",
  "data/processed_data/atacseq_e16/deeptools_output",
  "data/processed_data/atacseq_e16/diffbind_results_bed",
  "data/processed_data/atacseq_e16/filtered_bams",
  "data/processed_data/atacseq_e16/filtered_macs2",
  "data/processed_data/atacseq_e16/r_objects",
  "data/processed_data/homer_output",
  "data/processed_data/rnaseq_e16/normalized_counts",
  "data/processed_data/rnaseq_e16/r_objects",
  "data/processed_data/rnaseq_e17/normalized_counts",
  "data/processed_data/rnaseq_e17/r_objects/",
  "analysis",
  "figures",
  "preprocess",
  "setup",
  "tables"
)

for (dir in required_directories) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat(paste0("Created directory: ", dir, "\n"))
  }
}

cat("\nEnvironment setup complete.\n")