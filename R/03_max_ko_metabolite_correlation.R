##############################################################
# 03_max_ko_metabolite_correlation.R
##############################################################
#
# For each metabolite, correlate it with every predicted KO common to
# both tools and keep the largest value, max |rho|. This is a hypothesis-free
# measure of the strongest signal each tool recovers for that metabolite; the
# tool with the higher max |rho| is the winner.
##############################################################

# libraries
library(data.table)
library(ggplot2)

# config.R and 02_load_and_harmonize.R must have been run first
if (!exists("all_data")) {
  stop("Run config.R and 02_load_and_harmonize.R first (all_data not found).")
}

###############################
# run analysis for each dataset
###############################

maxcorr_results <- list()

for (ds_name in names(all_data)) {

  cat("\n=== Max KO-metabolite correlation:", ds_name, "===\n")
  d <- all_data[[ds_name]]

  # 1. transform- CLR for the KO predictions (they are compositional). We use the
  #    shared KOs so both tools are judged on the same set. 
  pic_clr <- clr_transform(d$pic_common)
  t4f_clr <- clr_transform(d$t4f_common)

  # 2. for each metabolite, get its best |rho| against any KO, for each tool
  metabs <- colnames(d$mtb)
  cat("  testing", length(metabs), "metabolites...\n")

  max_rho_pic <- numeric(length(metabs))
  max_rho_t4f <- numeric(length(metabs))

  for (i in seq_along(metabs)) {
    y <- d$mtb[, metabs[i]]
    max_rho_pic[i] <- max_abs_rho(y, pic_clr)
    max_rho_t4f[i] <- max_abs_rho(y, t4f_clr)
  }

  # 3. collect the results and note which tool scored higher for each metabolite
  res <- data.table(
    Dataset          = ds_name,
    KEGG             = metabs,
    max_rho_picrust2 = max_rho_pic,
    max_rho_tax4fun2 = max_rho_t4f
  )
  res <- res[!is.na(max_rho_picrust2) | !is.na(max_rho_tax4fun2)]
  res$winner <- ifelse(res$max_rho_picrust2 > res$max_rho_tax4fun2, "PICRUSt2",
              ifelse(res$max_rho_tax4fun2 > res$max_rho_picrust2, "Tax4Fun2", "Tie"))

  # write to file (needed to produce the figures)
  fwrite(res, file.path(d$res_dir, "max_ko_metabolite_correlation.tsv"), sep = "\t")

  # short summary output
  med_pic <- median(res$max_rho_picrust2, na.rm = TRUE)
  med_t4f <- median(res$max_rho_tax4fun2, na.rm = TRUE)
  cat("  median max |rho| - PICRUSt2:", round(med_pic, 3),
      "| Tax4Fun2:", round(med_t4f, 3), "\n")
  cat("  winner counts:\n")
  print(table(res$winner))

  # Store the table
  maxcorr_results[[ds_name]] <- res
  cat("  done:", ds_name, "\n")
}

###################
# Figures (all datasets together)
###################

cat("\n=== Figures ===\n")
dir.create("plots/final", recursive = TRUE, showWarnings = FALSE)

# Create one combined table (per-dataset tables into one)
maxcorr_all <- rbindlist(maxcorr_results)

# build the long table: one block per tool, then stack them
pic <- data.frame(Dataset = maxcorr_all$Dataset, KEGG = maxcorr_all$KEGG,
                  Method = "PICRUSt2", rho = maxcorr_all$max_rho_picrust2)
t4f <- data.frame(Dataset = maxcorr_all$Dataset, KEGG = maxcorr_all$KEGG,
                  Method = "Tax4Fun2", rho = maxcorr_all$max_rho_tax4fun2)
maxcorr_long <- rbind(pic, t4f)
maxcorr_long <- maxcorr_long[is.finite(maxcorr_long$rho), ]

# Figure 1A: max |rho| by tool
fig1a <- ggplot(maxcorr_long, aes(x = Method, y = rho, fill = Method)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.25) +
  facet_wrap(~ Dataset, scales = "free_y") +
  scale_fill_manual(values = tool_colors) +
  labs(y = "Max |Spearman rho| per metabolite", x = "Prediction method") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

# Figure 1B: how often each tool wins (PICRUSt2 left, Tax4Fun2 right)
wins <- maxcorr_all[maxcorr_all$winner != "Tie" & !is.na(maxcorr_all$winner), ]
win_counts <- as.data.frame(table(Dataset = wins$Dataset, winner = wins$winner))
names(win_counts)[names(win_counts) == "Freq"] <- "N"
win_counts <- win_counts[win_counts$N > 0, ]
win_counts$signed_N   <- ifelse(win_counts$winner == "PICRUSt2", -win_counts$N, win_counts$N)
win_counts$label_side <- ifelse(win_counts$signed_N < 0, 1.2, -0.2)   # push labels outside the bars

# rev factor matching the boxplot facet order
win_counts$Dataset <- factor(win_counts$Dataset, levels = rev(levels(factor(win_counts$Dataset))))

fig1b <- ggplot(win_counts, aes(x = signed_N, y = Dataset, fill = winner)) +
  geom_col(width = 0.55) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_text(aes(label = N, hjust = label_side), size = 2.6) +
  scale_fill_manual(values = tool_colors) +
  scale_x_continuous(labels = abs, expand = expansion(mult = 0.18)) +
  labs(x = "Number of metabolites", y = NULL, fill = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", axis.text.y = element_text(size = 8),
        panel.grid.major.y = element_blank())
# combine
fig1 <- cowplot::plot_grid(fig1a, fig1b, labels = c("A", "B"), ncol = 1, rel_heights = c(4.5, 1))

# MANUSCRIPT Figure 1
jpeg("plots/final/fig1_max_correlation.jpg", height = 7, width = 8, units = "in", res = 600)
print(fig1)
dev.off()
cat("  saved: plots/final/fig1_max_correlation.jpg\n")

cat("\n=== Done. Proceed to script 04 ===\n")
