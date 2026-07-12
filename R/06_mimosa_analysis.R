##############################################################
# 06_mimosa_analysis.R  -  analyse MIMOSA2 results
##############################################################
#
# After running MIMOSA2 on the web app (script 05) and downloading its
# "Model Summaries" tables, this script reads them, classifies each metabolite
# as microbiome-governed (MGM) or not, and makes the figures
##############################################################

# config.R and 02_load_and_harmonize.R must have been run first
if (!exists("all_data")) 
stop("Run config.R and 02_load_and_harmonize.R first (all_data not found).")

# libraries
library(data.table)
library(ggplot2)

mimosa_results <- list()

for (dataset_name in names(all_data)) {

  cat("\n=== MIMOSA2 results for", dataset_name, "===\n")
  dataset    <- all_data[[dataset_name]]
  mimosa_dir <- file.path(dataset$res_dir, "mimosa2_input")

  # 1. check if the two downloaded MIMOSA2 summaries are present
  pic_summary_fp <- file.path(mimosa_dir, "mimosa_pic_model_summaries.tsv")
  t4f_summary_fp <- file.path(mimosa_dir, "mimosa_t4f_model_summaries.tsv")
  if (!file.exists(pic_summary_fp) || !file.exists(t4f_summary_fp)) {
    cat("  MIMOSA2 output not found for", dataset_name, "- do the manual step in script 05 first.\n")
    cat("    expected:", pic_summary_fp, "\n")
    cat("             ", t4f_summary_fp, "\n\n")
    next
  }

  # 2. read and tidy. MIMOSA2 columns: compound, Rsq, PVal, Slope. Add a BH-corrected q-value
  pic_summary <- fread(pic_summary_fp)
  t4f_summary <- fread(t4f_summary_fp)

  mimosa_pic <- data.frame(
    metabolite_cid = pic_summary$compound,
    r2             = pic_summary$Rsq,
    p_value        = pic_summary$PVal,
    q_value        = p.adjust(pic_summary$PVal, method = "BH"),
    slope          = pic_summary$Slope
  )
  mimosa_t4f <- data.frame(
    metabolite_cid = t4f_summary$compound,
    r2             = t4f_summary$Rsq,
    p_value        = t4f_summary$PVal,
    q_value        = p.adjust(t4f_summary$PVal, method = "BH"),
    slope          = t4f_summary$Slope
  )

  fwrite(mimosa_pic, file.path(dataset$res_dir, "mimosa_pic_results.csv"))
  fwrite(mimosa_t4f, file.path(dataset$res_dir, "mimosa_t4f_results.csv"))

  # 3. merge the two tools by metabolite and then classify
  pic_named <- mimosa_pic
  names(pic_named) <- c("metabolite_cid", "r2_pic", "p_pic", "q_pic", "slope_pic")
  t4f_named <- mimosa_t4f
  names(t4f_named) <- c("metabolite_cid", "r2_t4f", "p_t4f", "q_t4f", "slope_t4f")
  both <- merge(pic_named, t4f_named, by = "metabolite_cid")
  both$Dataset <- dataset_name

  # an MGM is a metabolite with FDR < 0.05 for that tool
  both$mgm_pic <- both$q_pic < 0.05
  both$mgm_t4f <- both$q_t4f < 0.05

  # which tool(s) called it an MGM. which() drops NA q-values, so a metabolite
  # with a missing q-value falls through to "Neither" 
  both$mgm_cat <- "Neither"
  both$mgm_cat[which(both$mgm_pic & both$mgm_t4f)]  <- "Both MGM"
  both$mgm_cat[which(both$mgm_pic & !both$mgm_t4f)] <- "PICRUSt2-only MGM"
  both$mgm_cat[which(!both$mgm_pic & both$mgm_t4f)] <- "Tax4Fun2-only MGM"

  # short summary
  cat("  metabolites evaluated:", nrow(both), "\n")
  print(table(both$mgm_cat))
  cat("  median R2 - PICRUSt2:", round(median(both$r2_pic, na.rm = TRUE), 4),
      "| Tax4Fun2:", round(median(both$r2_t4f, na.rm = TRUE), 4), "\n")

  # per-dataset MGM table (the non-Neither metabolites)
  mgm_rows <- both[both$mgm_cat != "Neither", ]
  mgm_rows$mgm_cat <- factor(mgm_rows$mgm_cat,
                        levels = c("Both MGM", "PICRUSt2-only MGM", "Tax4Fun2-only MGM"))
  mgm_rows <- mgm_rows[order(mgm_rows$mgm_cat), ]
  mgm_table <- data.frame(
    Dataset       = mgm_rows$Dataset,
    KEGG_compound = mgm_rows$metabolite_cid,
    MGM_category  = as.character(mgm_rows$mgm_cat),
    R2_PICRUSt2   = round(mgm_rows$r2_pic, 3),
    q_PICRUSt2    = signif(mgm_rows$q_pic, 3),
    R2_Tax4Fun2   = round(mgm_rows$r2_t4f, 3),
    q_Tax4Fun2    = signif(mgm_rows$q_t4f, 3)
  )
  fwrite(mgm_table, file.path(dataset$res_dir, "MIMOSA_MGM_table.csv"))

  mimosa_results[[dataset_name]] <- both
  cat("  done:", dataset_name, "\n")
}

##############################################
# Figures + Supplementary Table S1
##############################################

if (length(mimosa_results) > 0) {

  cat("\n=== Figures ===\n")

  dir.create("../plots", showWarnings = FALSE)

  mimosa_all <- do.call(rbind, mimosa_results)
  # fix the category order for the legend 
  mimosa_all$mgm_cat <- factor(mimosa_all$mgm_cat, levels = names(mgm_colors))

  # Figure 4: proportion of metabolites in each MGM category, per dataset
  mgm_counts <- as.data.frame(table(Dataset = mimosa_all$Dataset, mgm_cat = mimosa_all$mgm_cat))
  names(mgm_counts)[names(mgm_counts) == "Freq"] <- "n"
  totals <- tapply(mgm_counts$n, mgm_counts$Dataset, sum)
  mgm_counts$pct <- mgm_counts$n / totals[as.character(mgm_counts$Dataset)]

  fig4 <- ggplot(mgm_counts, aes(x = Dataset, y = pct, fill = mgm_cat)) +
    geom_col(width = 0.6) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = mgm_colors, name = NULL) +
    labs(x = NULL, y = "Proportion of metabolites") +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")

  # MANUSCRIPT Figure 4
  jpeg("../plots/mimosa_mgm_proportions.jpg", height = 5, width = 8, units = "in", res = 600)
  print(fig4)
  dev.off()

  cat("  saved: ../plots/mimosa_mgm_proportions.jpg\n")

  # Figure 5: MIMOSA2 model fit (R2) for PICRUSt2 vs Tax4Fun2 per dataset
  fig5 <- ggplot(mimosa_all, aes(x = r2_pic, y = r2_t4f, colour = mgm_cat)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
    geom_point(alpha = 0.6, size = 2) +
    scale_colour_manual(values = mgm_colors, name = NULL) +
    coord_fixed() +
    facet_wrap(~ Dataset) +
    labs(x = expression(R^2 * " (PICRUSt2)"), y = expression(R^2 * " (Tax4Fun2)")) +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(face = "bold"))

  # MANUSCRIPT Figure 5
  jpeg("../plots/mimosa_r2_scatter.jpg", height = 6, width = 9, units = "in", res = 600)
  print(fig5)
  dev.off()

  cat("  saved: ../plots/mimosa_r2_scatter.jpg\n")

  # Supplementary Table S1: the full MGM list across datasets
  mgm_all <- mimosa_all[mimosa_all$mgm_cat != "Neither", ]
  mgm_all$mgm_cat <- factor(mgm_all$mgm_cat,
                            levels = c("Both MGM", "PICRUSt2-only MGM", "Tax4Fun2-only MGM"))
  mgm_all <- mgm_all[order(mgm_all$Dataset, mgm_all$mgm_cat), ]
  supp_table_s1 <- data.frame(
    Dataset       = mgm_all$Dataset,
    KEGG_compound = mgm_all$metabolite_cid,
    MGM_category  = as.character(mgm_all$mgm_cat),
    R2_PICRUSt2   = round(mgm_all$r2_pic, 3),
    q_PICRUSt2    = signif(mgm_all$q_pic, 3),
    R2_Tax4Fun2   = round(mgm_all$r2_t4f, 3),
    q_Tax4Fun2    = signif(mgm_all$q_t4f, 3)
  )
  # MANUSCRIPT Supplementary Table S1
  fwrite(supp_table_s1, "supplementary_table_s1_mgm_list.csv")
  cat("  saved: supplementary_table_s1_mgm_list.csv\n")
}
cat("\n*** Done. Proceed to script 07 (summary/combine). ***\n")
