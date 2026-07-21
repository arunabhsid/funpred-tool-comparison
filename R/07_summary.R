##############################################################
# 07_summary.R  -  cross-dataset summary table 
##############################################################
#
# Takes the per-dataset results from the earlier analyses into one overview
# table: dataset sizes, how often each tool won in each analysis, and the
# MIMOSA MGM counts and median model fit.
#
##############################################################

# config.R and 02_load_and_harmonize.R must have been run first
if (!exists("all_data")) {
  stop("Run config.R and 02_load_and_harmonize.R first (all_data not found).")
}

library(data.table)


summary_rows <- list()

for (dataset_name in names(all_data)) {

  dataset     <- all_data[[dataset_name]]
  results_dir <- dataset$res_dir

  # dataset dimensions (from all_data)
  summary_row <- list(
    Dataset        = dataset_name,
    n_samples      = dataset$n_samples,
    n_metabolites  = ncol(dataset$mtb),
    n_KO_picrust   = ncol(dataset$pic_full),
    n_KO_tax4fun   = ncol(dataset$t4f_full),
    maxcorr_median_pic  = NA_real_, maxcorr_median_t4f = NA_real_,
    maxcorr_win_pic     = NA_integer_, maxcorr_win_t4f = NA_integer_,
    targeted_n_tested   = NA_integer_,
    targeted_win_pic    = NA_integer_, targeted_win_t4f = NA_integer_,
    mimosa_mgm_both     = NA_integer_, mimosa_mgm_pic_only = NA_integer_,
    mimosa_mgm_t4f_only = NA_integer_,
    mimosa_median_r2_pic = NA_real_, mimosa_median_r2_t4f = NA_real_
  )

  # max KO-metabolite correlation (script 03) 
  maxcorr_fp <- file.path(results_dir, "max_ko_metabolite_correlation.tsv")
  if (file.exists(maxcorr_fp)) {
    maxcorr <- fread(maxcorr_fp)
    summary_row$maxcorr_median_pic <- round(median(maxcorr$max_rho_picrust2, na.rm = TRUE), 3)
    summary_row$maxcorr_median_t4f <- round(median(maxcorr$max_rho_tax4fun2, na.rm = TRUE), 3)
    summary_row$maxcorr_win_pic    <- sum(maxcorr$winner == "PICRUSt2", na.rm = TRUE)
    summary_row$maxcorr_win_t4f    <- sum(maxcorr$winner == "Tax4Fun2", na.rm = TRUE)
  }

  # targeted correlation (script 04) 
  targeted_fp <- file.path(results_dir, "targeted_correlation.tsv")
  if (file.exists(targeted_fp)) {
    targeted <- fread(targeted_fp)
    summary_row$targeted_n_tested <- nrow(targeted)
    summary_row$targeted_win_pic  <- sum(targeted$winner == "PICRUSt2", na.rm = TRUE)
    summary_row$targeted_win_t4f  <- sum(targeted$winner == "Tax4Fun2", na.rm = TRUE)
  }

  # MIMOSA (script 06): merge the two tools and reclassify with the same rule 
  pic_results_fp <- file.path(results_dir, "mimosa_pic_results.csv")
  t4f_results_fp <- file.path(results_dir, "mimosa_t4f_results.csv")
  if (file.exists(pic_results_fp) && file.exists(t4f_results_fp)) {
    mimosa_pic <- fread(pic_results_fp)
    names(mimosa_pic) <- c("metabolite_cid", "r2_pic", "p_pic", "q_pic", "slope_pic")
    mimosa_t4f <- fread(t4f_results_fp)
    names(mimosa_t4f) <- c("metabolite_cid", "r2_t4f", "p_t4f", "q_t4f", "slope_t4f")
    both <- merge(mimosa_pic, mimosa_t4f, by = "metabolite_cid")

    mgm_pic <- both$q_pic < 0.05
    mgm_t4f <- both$q_t4f < 0.05
    summary_row$mimosa_mgm_both      <- sum(mgm_pic & mgm_t4f,  na.rm = TRUE)
    summary_row$mimosa_mgm_pic_only  <- sum(mgm_pic & !mgm_t4f, na.rm = TRUE)
    summary_row$mimosa_mgm_t4f_only  <- sum(!mgm_pic & mgm_t4f, na.rm = TRUE)
    summary_row$mimosa_median_r2_pic <- round(median(both$r2_pic, na.rm = TRUE), 4)
    summary_row$mimosa_median_r2_t4f <- round(median(both$r2_t4f, na.rm = TRUE), 4)
  }

  summary_rows[[dataset_name]] <- as.data.frame(summary_row, stringsAsFactors = FALSE)
}

summary_table <- do.call(rbind, summary_rows)
fwrite(summary_table, "results/Summary_table_all_datasets.tsv", sep = "\t")
cat("saved: results/Summary_table_all_datasets.tsv\n\n")
print(summary_table)

cat("\n*** Summary complete. Manuscript figures are in plots/final/. ***\n")
