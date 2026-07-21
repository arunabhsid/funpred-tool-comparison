##############################################################
# 05_mimosa_prep.R  -  prepare inputs for MIMOSA2
##############################################################
#
# We run MIMOSA2 through the free web app:
# (https://borenstein-lab.github.io/MIMOSA2shiny/)
# This script only formats and saves the input files, then prints the manual
# upload steps. 
##############################################################

# config.R and 02_load_and_harmonize.R must have been run first
if (!exists("all_data")) {
  stop("Run config.R and 02_load_and_harmonize.R first (all_data not found).")
}

library(data.table)

for (dataset_name in names(all_data)) {

  cat("\n=== Preparing MIMOSA2 files for:", dataset_name, "===\n")
  dataset <- all_data[[dataset_name]]

  # MIMOSA2 expects features in rows and samples in columns, with the first
  # column named "KO" (KO tables) or "KEGG" (metabolites). Our matrices are
  # samples x features, so we transpose them and move the IDs into a column.
  mimosa_dir <- file.path(dataset$res_dir, "mimosa2_input")
  dir.create(mimosa_dir, recursive = TRUE, showWarnings = FALSE)

  # PICRUSt2 KOs
  pic_t <- t(dataset$pic_common)
  ko_picrust <- data.frame(KO = rownames(pic_t), pic_t, check.names = FALSE, row.names = NULL)

  # Tax4Fun2 KOs
  t4f_t <- t(dataset$t4f_common)
  ko_tax4fun <- data.frame(KO = rownames(t4f_t), t4f_t, check.names = FALSE, row.names = NULL)

  # Metabolites
  mtb_t <- t(dataset$mtb)
  metabolites <- data.frame(KEGG = rownames(mtb_t), mtb_t, check.names = FALSE, row.names = NULL)

  # match all three files on the same samples
  stopifnot(
    identical(colnames(ko_picrust)[-1], colnames(metabolites)[-1]),
    identical(colnames(ko_tax4fun)[-1], colnames(metabolites)[-1])
  )

  # save
  pic_fp <- file.path(mimosa_dir, "ko_picrust2_for_mimosa.tsv")
  t4f_fp <- file.path(mimosa_dir, "ko_tax4fun2_for_mimosa.tsv")
  met_fp <- file.path(mimosa_dir, "metabolites_kegg_for_mimosa.tsv")
  fwrite(ko_picrust,  pic_fp, sep = "\t")
  fwrite(ko_tax4fun,  t4f_fp, sep = "\t")
  fwrite(metabolites, met_fp, sep = "\t")
  cat("  files saved to:", mimosa_dir, "\n")

  # MANUAL STEP: run MIMOSA2 on the web app for both imputs (PICRUSt2 and Tax4Fun2) 
  cat("\n  MANUAL STEP: run MIMOSA2 at https://borenstein-lab.github.io/MIMOSA2shiny/\n")
  cat("  Parameters for both runs: Data type 'Metagenome: Total KO abundances';",
      "Metabolite ID 'KEGG Compound IDs';",
      "Metabolic model 'KEGG metabolic model'; Rank-based ON.\n\n")

  cat("  Run A - PICRUSt2:\n")
  cat("    microbiome file:", pic_fp, "\n")
  cat("    metabolite file:", met_fp, "\n")
  cat("    download the 'Model Summaries' table and save as:\n")
  cat("     ", file.path(mimosa_dir, "mimosa_pic_model_summaries.tsv"), "\n\n")

  cat("  Run B - Tax4Fun2:\n")
  cat("    microbiome file:", t4f_fp, "\n")
  cat("    metabolite file:", met_fp, "\n")
  cat("    download the 'Model Summaries' table and save as:\n")
  cat("     ", file.path(mimosa_dir, "mimosa_t4f_model_summaries.tsv"), "\n")
}

cat("\n*** MIMOSA2 input files prepared. Complete the manual uploads for every dataset, then run script 06. ***\n")
