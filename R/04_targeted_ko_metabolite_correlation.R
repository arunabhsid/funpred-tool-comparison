##############################################################
# 04_targeted_ko_metabolite_correlation.R
##############################################################
#
# For each metabolite, correlate the mean predicted abundance of the KOs that
# KEGG links to it with its measured level (separately for PICRUSt2 and 
# Tax4Fun2) using each tool's full KO set and raw abundances. The tool with 
# the higher absolute Spearman correlation is the winner.
#
# Part A builds the compound-to-KO mapping from the KEGG API (run only once).
# Part B runs the correlations for every dataset and makes the figures.
##############################################################

# config.R and 02_load_and_harmonize.R must have been run first
if (!exists("all_data")) {
  stop("Run config.R and 02_load_and_harmonize.R first (all_data not found).")
  }

library(data.table)
library(KEGGREST)
library(ggplot2)

##############################################
# PART A: build the compound -> KO mapping
##############################################
# This queries the KEGG API (it is slow: ~30-60 min for ~200
# compounds) and saves the result. If the file already exists we skip
# to Part B.
# The compound -> KO lookup below was inspired by a KEGGREST module -> KO
# mapping approach (https://www.biostars.org/p/496921/)

# look up the KOs (enzymes) linked to one KEGG compound, through its reactions
get_kos_for_compound <- function(compound_id) {
  reactions <- tryCatch(keggLink("reaction", paste0("cpd:", compound_id)), error = function(e) NULL)
  reactions <- unique(unlist(reactions))
  if (length(reactions) == 0) return(character(0))

  KOs <- character(0)
  for (reaction_id in reactions) {
    Sys.sleep(1)                                   # pause between requests 
    KO_result <- tryCatch(keggLink("ko", reaction_id), error = function(e) NULL)
    KOs <- c(KOs, KO_result)
  }
  unique(sub("ko:", "", KOs))                      # "ko:K00001" -> "K00001"
}
# check
if (file.exists(kegg_mapping_fp)) {
  cat("KEGG mapping available:", kegg_mapping_fp, "\n")
  cat("(delete that file and re-run to rebuild it)\n\n")
} else {

  cat("Building the KEGG compound -> KO mapping (this queries KEGG)...\n")

  # every unique KEGG compound measured across the four datasets
  all_compounds <- character(0)
  for (d in all_data) {
  all_compounds <- c(all_compounds, colnames(d$mtb))
  }
  all_compounds <- unique(all_compounds)
  cat("  unique compounds to look up:", length(all_compounds), "\n")

  # look up the KOs for each compound
  compound_to_ko <- list()
  for (i in seq_along(all_compounds)) {
    compound_id <- all_compounds[i]
    if (i %% 25 == 0) {
      cat("  ", i, "/", length(all_compounds), "\n") # check check
    }
    Sys.sleep(1)
    compound_to_ko[[compound_id]] <- get_kos_for_compound(compound_id)
  }

  # save one row per compound, its KOs joined by commas
  rows <- list()
  for (compound_id in names(compound_to_ko)) {
    KOs <- compound_to_ko[[compound_id]]
    if (length(KOs) == 0) next
    rows[[compound_id]] <- data.frame(Compound = compound_id, KOs = paste(KOs, collapse = ","))
  }

  mapping <- do.call(rbind, rows)
  fwrite(mapping, kegg_mapping_fp, sep = "\t")
  cat("  saved mapping for", nrow(mapping), "compounds ->", kegg_mapping_fp, "\n\n")
}


##############################################
# PART B: targeted correlations per dataset
##############################################

# load the mapping and expand it to one row per Compound-KO pair
map_wide <- fread(kegg_mapping_fp) # i saved it in main dir

pair_rows <- list()
for (i in seq_len(nrow(map_wide))) {
  KOs <- strsplit(map_wide$KOs[i], ",")[[1]]       # split the comma-separated KOs (strsplit always returns a list- get the vector)
  KOs <- trimws(KOs)
  KOs <- KOs[KOs != "" & KOs != "NA"]
  if (length(KOs) == 0) next
  pair_rows[[i]] <- data.frame(Compound = map_wide$Compound[i], KO = KOs)
}
pairs <- do.call(rbind, pair_rows) # combine into one long table

# clean the IDs and keep only proper KEGG codes
pairs$KO       <- clean_ko_ids(pairs$KO)
pairs$Compound <- clean_compound_ids(pairs$Compound)
pairs <- pairs[grepl("^K[0-9]{5}$", pairs$KO) & grepl("^C[0-9]{5}$", pairs$Compound), ]
pairs <- unique(pairs)
cat("mapping:", nrow(pairs), "compound-KO pairs |",
    length(unique(pairs$Compound)), "compounds |",
    length(unique(pairs$KO)), "KOs\n\n")

targeted_results <- list()

for (dataset_name in names(all_data)) {

  cat("--- Targeted correlation:", dataset_name, "---\n")
  dataset   <- all_data[[dataset_name]]
  mtb <- dataset$mtb
  pic <- dataset$pic_full
  t4f <- dataset$t4f_full

  # compounds we can test: measured here AND present in the mapping
  testable <- intersect(colnames(mtb), unique(pairs$Compound))
  cat("  testable compounds:", length(testable), "\n")
  if (length(testable) == 0) { cat("  none testable, skipping\n\n"); next }

  # one slot per testable compound
  n <- length(testable)
  nKO_pic <- integer(n); rho_pic <- rep(NA_real_, n); p_pic <- rep(NA_real_, n)
  nKO_t4f <- integer(n); rho_t4f <- rep(NA_real_, n); p_t4f <- rep(NA_real_, n)

  for (i in seq_along(testable)) {
    compound_id <- testable[i]
    metabolite   <- mtb[, compound_id]                              # measured metabolite
    mapped_KOs <- pairs$KO[pairs$Compound == compound_id]
    pic_KOs <- intersect(mapped_KOs, colnames(pic))     # mapped KOs present in PICRUSt2
    t4f_KOs <- intersect(mapped_KOs, colnames(t4f))     # mapped KOs present in Tax4Fun2
    nKO_pic[i] <- length(pic_KOs)
    nKO_t4f[i] <- length(t4f_KOs)

    # PICRUSt2: average the mapped KOs, then correlate with the metabolite
    if (length(pic_KOs) > 0) {
      ko_mean <- rowMeans(pic[, pic_KOs, drop = FALSE])
      spearman_result <- spearman_with_p(ko_mean, metabolite)
      rho_pic[i] <- spearman_result["rho"]; p_pic[i] <- spearman_result["p"]
    }
    # Tax4Fun2: same like PICRUST2
    if (length(t4f_KOs) > 0) {
      ko_mean <- rowMeans(t4f[, t4f_KOs, drop = FALSE])
      spearman_result <- spearman_with_p(ko_mean, metabolite)
      rho_t4f[i] <- spearman_result["rho"]; p_t4f[i] <- spearman_result["p"]
    }
  }

  res <- data.frame(
    Dataset  = dataset_name, Compound = testable,
    nKO_pic  = nKO_pic, nKO_t4f = nKO_t4f,
    rho_pic  = rho_pic, p_pic = p_pic,
    rho_t4f  = rho_t4f, p_t4f = p_t4f
  )

  # multiple-testing correction (Benjamini-Hochberg) per tool
  res$fdr_pic <- p.adjust(res$p_pic, method = "BH")
  res$fdr_t4f <- p.adjust(res$p_t4f, method = "BH")

  # winner: the tool has the larger |rho|. if only one tool could be
  # tested for a compound it wins by default; if neither, it is a tie.
  abs_rho_pic <- abs(res$rho_pic)
  abs_rho_t4f <- abs(res$rho_t4f)
  res$winner <- "Tie"
  res$winner[is.finite(abs_rho_pic) & !is.finite(abs_rho_t4f)] <- "PICRUSt2"
  res$winner[!is.finite(abs_rho_pic) & is.finite(abs_rho_t4f)] <- "Tax4Fun2"
  res$winner[is.finite(abs_rho_pic) & is.finite(abs_rho_t4f) & abs_rho_pic > abs_rho_t4f]  <- "PICRUSt2"
  res$winner[is.finite(abs_rho_pic) & is.finite(abs_rho_t4f) & abs_rho_t4f > abs_rho_pic]  <- "Tax4Fun2"
  
  # keep compounds with at least one testable KO
  res <- res[res$nKO_pic >= 1 | res$nKO_t4f >= 1, ]

  fwrite(res, file.path(dataset$res_dir, "targeted_correlation.tsv"), sep = "\t")

  cat("  PICRUSt2 wins:", sum(res$winner == "PICRUSt2"),
      "| Tax4Fun2 wins:", sum(res$winner == "Tax4Fun2"),
      "| ties:", sum(res$winner == "Tie"), "\n\n")

  targeted_results[[dataset_name]] <- res
}

##############################################
# Figures
##############################################

cat("=== Figures ===\n")
dir.create("plots/final", recursive = TRUE, showWarnings = FALSE)

targeted_all <- rbindlist(targeted_results)

# Figure 2A: signed rho, PICRUSt2 vs Tax4Fun2
scatter_dat <- targeted_all[is.finite(targeted_all$rho_pic) & is.finite(targeted_all$rho_t4f), ]
scatter_dat$higher <- ifelse(abs(scatter_dat$rho_t4f) > abs(scatter_dat$rho_pic),
                             "Tax4Fun2 higher", "PICRUSt2 higher")
fig2a <- ggplot(scatter_dat, aes(x = rho_pic, y = rho_t4f, colour = higher)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.6, size = 2) +
  facet_wrap(~ Dataset) +
  scale_colour_manual(values = c("PICRUSt2 higher" = tool_colors[["PICRUSt2"]],
                                 "Tax4Fun2 higher" = tool_colors[["Tax4Fun2"]])) +
  labs(x = "Spearman rho (PICRUSt2)", y = "Spearman rho (Tax4Fun2)", colour = "Higher correlation") +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

# Figure 2B: how often each tool wins
wins <- targeted_all[targeted_all$winner != "Tie" & !is.na(targeted_all$winner), ]
win_counts <- as.data.frame(table(Dataset = wins$Dataset, winner = wins$winner))
names(win_counts)[names(win_counts) == "Freq"] <- "N"
win_counts <- win_counts[win_counts$N > 0, ]
win_counts$signed_N   <- ifelse(win_counts$winner == "PICRUSt2", -win_counts$N, win_counts$N)
win_counts$label_side <- ifelse(win_counts$signed_N < 0, 1.2, -0.2)   # push labels outside the bars

# rev factor matching the scatter facet order
win_counts$Dataset <- factor(win_counts$Dataset, levels = rev(levels(factor(win_counts$Dataset))))

fig2b <- ggplot(win_counts, aes(x = signed_N, y = Dataset, fill = winner)) +
  geom_col(width = 0.55) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_text(aes(label = N, hjust = label_side), size = 2.6) +
  scale_fill_manual(values = tool_colors) +
  scale_x_continuous(labels = abs, expand = expansion(mult = 0.18)) +
  labs(x = "Number of metabolites", y = NULL, fill = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", axis.text.y = element_text(size = 8),
        panel.grid.major.y = element_blank())

# combine both into one two-panel Figure 2 (A = scatter, B = win-count strip)
fig2 <- cowplot::plot_grid(fig2a, fig2b, labels = c("A", "B"), ncol = 1, rel_heights = c(4.5, 1))

# MANUSCRIPT Figure 2
jpeg("plots/final/fig2_targeted_correlation.jpg", height = 8, width = 9, units = "in", res = 600)
print(fig2)
dev.off()
cat("  saved: plots/final/fig2_targeted_correlation.jpg\n")

# Summary for manuscript
for (nm in names(targeted_results)) {
  r <- targeted_results[[nm]]
  cat(nm, "| PIC:", round(median(abs(r$rho_pic), na.rm = TRUE), 3),
      "| T4F:", round(median(abs(r$rho_t4f), na.rm = TRUE), 3), "\n")
}

cat("\n=== Done. Proceed to script 05 (MIMOSA prep). ===\n")
