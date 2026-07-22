##############################################################
# 08_shotgun_analysis.R  -  shotgun metagenomics reference
##############################################################
#
# Dataset: Franzosa et al. 2019 (220 IBD samples)- shotgun KOs from HUMAnN3
# with paired faecal metabolomics.
#
# Samples are matched to the metabolomics by SRA accession (1). The 43 UMCG
# validation samples use different IDs in the two sources and are left out,
# leaving 177 matched samples.
#
# We run the same two analyses as before on the observed KOs - the maximum and
# the targeted KO-metabolite correlation - and compare them with the 16S
# results
##############################################################

source("R/config.R")

# libraries
library(data.table)
library(ggplot2)

# paths
shotgun_data_dir <- "data/FRANZOSA_IBD_2019"
shotgun_res_dir  <- "results/FRANZOSA_IBD_2019_shotgun"
shotgun_name     <- "FRANZOSA_IBD_2019"
dir.create(shotgun_res_dir, recursive = TRUE, showWarnings = FALSE) 

###################
# 1: match the HUMAnN sample IDs to the study IDs
###################
# HUMAnN's columns are SRR sequencing accessions, but the
# metabolomics uses study IDs like "PRISM.7122". We match them with the SRA run
# table (franzosa_runinfo.tsv, downloaded from ENA for BioProject PRJNA400072),
# which gives the study ID (sample_alias) behind each SRR accession.
#
# The alias is written a little differently from the curated study IDs, so we
# rename it to match ("PRISM_7941" -> "PRISM.7941", "LLDeep_0012" ->
# "Validation.LLDeep_0012")

cat("=== 1: matching sample IDs ===\n")

runinfo_fp <- file.path(shotgun_data_dir, "franzosa_runinfo.tsv")
humann_fp  <- file.path(shotgun_data_dir, "humann_kos.tsv")
if (!file.exists(runinfo_fp)) stop("franzosa_runinfo.tsv not found in ", shotgun_data_dir)
if (!file.exists(humann_fp))  stop("humann_kos.tsv not found in ", shotgun_data_dir)

runinfo <- fread(runinfo_fp)

# rewrite the SRA alias into the curated study-ID style
study_id <- runinfo$sample_alias
study_id <- sub("^PRISM_",  "PRISM.",             study_id)
study_id <- sub("^LLDeep_", "Validation.LLDeep_", study_id)

# lookup table: SRR accession -> study ID
id_map <- setNames(study_id, runinfo$run_accession)
cat("  run table:", length(id_map), "SRR accessions\n\n")

###################
# 2: load and remap the HUMAnN KO table
###################
# HUMAnN puts KO IDs in the first column and samples across. TODO: drop the
# stratified per-species rows, if any (they contain "|") and the UNMAPPED/UNGROUPED
# rows and rename the sample columns to study IDs, then transpose to samples x KOs

cat("=== 2: HUMAnN KO table ===\n")

ko_humann <- fread(humann_fp)
names(ko_humann)[1] <- "KO"

stratified <- grepl("\\|", ko_humann$KO)
special    <- ko_humann$KO %in% c("UNMAPPED", "UNGROUPED")
ko_humann <- ko_humann[!stratified & !special, ]
cat("  KOs after removing stratified/UNMAPPED/UNGROUPED:", nrow(ko_humann), "\n")

# rename each sample column to its study ID, matched by SRR accession. 
#The column names look like "SRR6468499_multi_omic...", so
# we pull out the SRR first
srr <- sub("_.*", "", names(ko_humann)[-1]) # SRR from each column name
names(ko_humann)[-1] <- id_map[srr]

ko_humann$KO <- clean_ko_ids(ko_humann$KO)

# transpose to a samples x KOs matrix
ko_matrix <- as.matrix(ko_humann[, -1, with = FALSE])
rownames(ko_matrix) <- ko_humann$KO
shotgun_ko_all <- t(ko_matrix)     # rows = samples, cols = KOs
cat("  shotgun KO matrix:", nrow(shotgun_ko_all), "samples x", ncol(shotgun_ko_all), "KOs\n\n")


###################
# 3: metabolomics
###################
# Same reshaping as the 16S pipeline: relabel each metabolite by its KEGG
# compound ID, drop those without one, average duplicates -> samples x compounds.

cat("=== 3: metabolomics ===\n")

mtb <- fread(file.path(shotgun_data_dir, "mtb.tsv"))
mtb_map <- fread(file.path(shotgun_data_dir, "mtb.map.tsv"))

long <- melt(mtb, id.vars = "Sample", variable.name = "Compound", value.name = "abundance")
long <- merge(long, unique(mtb_map[, c("Compound", "KEGG")]), by = "Compound", all.x = TRUE)
long <- long[!is.na(KEGG) & KEGG != "" & KEGG != "NA"]
wide <- dcast(long, Sample ~ KEGG, value.var = "abundance", fun.aggregate = mean)

mtb_matrix <- as.matrix(wide[, setdiff(names(wide), "Sample"), with = FALSE])
rownames(mtb_matrix) <- wide$Sample
colnames(mtb_matrix) <- clean_compound_ids(colnames(mtb_matrix))
cat("  metabolomics:", nrow(mtb_matrix), "samples x", ncol(mtb_matrix), "KEGG compounds\n\n")


###################
# 4: keep shared samples and clean
###################

common_samples <- intersect(rownames(shotgun_ko_all), rownames(mtb_matrix))
cat("=== 4: shared samples:", length(common_samples), "===\n")
if (length(common_samples) < 10) stop("Too few overlapping samples (", length(common_samples), ").")

shotgun_ko  <- to_clean_matrix(shotgun_ko_all[common_samples, , drop = FALSE])
shotgun_mtb <- to_clean_matrix(mtb_matrix[common_samples, , drop = FALSE])
cat("  after cleaning - KOs:", ncol(shotgun_ko), "| metabolites:", ncol(shotgun_mtb), "\n\n")


###################
# 5: maximum KO-metabolite correlation (shotgun)
###################
# For each metabolite, the strongest correlation with any observed KO. KOs are
# CLR-transformed first (they are compositional); the metabolites are left raw
# (Spearman ranks them internally).

cat("=== 5: maximum KO-metabolite correlation ===\n")

shotgun_clr <- clr_transform(shotgun_ko)
metabolites <- colnames(shotgun_mtb)
cat("  testing", length(metabolites), "metabolites...\n")

max_rho_shotgun <- numeric(length(metabolites))
for (i in seq_along(metabolites)) {
  max_rho_shotgun[i] <- max_abs_rho(shotgun_mtb[, metabolites[i]], shotgun_clr)
}

max_corr_shotgun <- data.frame(
  Dataset = shotgun_name,
  KEGG    = metabolites,
  max_rho_shotgun = max_rho_shotgun
)
max_corr_shotgun <- max_corr_shotgun[!is.na(max_corr_shotgun$max_rho_shotgun), ]

median_max_shotgun <- median(max_corr_shotgun$max_rho_shotgun, na.rm = TRUE)
cat("  median max |rho| (shotgun):", round(median_max_shotgun, 3), "\n")
fwrite(max_corr_shotgun, file.path(shotgun_res_dir, "shotgun_max_correlation.tsv"), sep = "\t")
cat("  done\n\n")


###################
# 6: targeted KO-metabolite correlation (shotgun)
###################
# For each metabolite, average the observed abundance of the KOs KEGG links to
# it and correlate that with the metabolite. Uses the mapping built in script 04.

cat("=== 6: targeted KO-metabolite correlation ===\n")

if (!file.exists(kegg_mapping_fp)) {
  stop("KEGG mapping not found: ", kegg_mapping_fp, " (build it in script 04).")
  }

# expand the mapping to one Compound-KO pair per row (same as script 04)
map_wide <- fread(kegg_mapping_fp)
pair_rows <- list()
for (i in seq_len(nrow(map_wide))) {
  KOs <- strsplit(map_wide$KOs[i], ",")[[1]]
  KOs <- trimws(KOs)
  KOs <- KOs[KOs != "" & KOs != "NA"]
  if (length(KOs) == 0) next
  pair_rows[[i]] <- data.frame(Compound = map_wide$Compound[i], KO = KOs)
}
pairs <- do.call(rbind, pair_rows)
pairs$KO       <- clean_ko_ids(pairs$KO)
pairs$Compound <- clean_compound_ids(pairs$Compound)
pairs <- pairs[grepl("^K[0-9]{5}$", pairs$KO) & grepl("^C[0-9]{5}$", pairs$Compound), ]
pairs <- unique(pairs)

testable <- intersect(colnames(shotgun_mtb), unique(pairs$Compound))
cat("  testable compounds:", length(testable), "\n")

targeted_shotgun <- NULL
if (length(testable) > 0) {

  n <- length(testable)
  nKO_shotgun <- integer(n)
  rho_shotgun <- rep(NA_real_, n)
  p_shotgun   <- rep(NA_real_, n)

  for (i in seq_along(testable)) {
    compound_id <- testable[i]
    metabolite  <- shotgun_mtb[, compound_id]
    mapped_KOs  <- pairs$KO[pairs$Compound == compound_id]
    shotgun_KOs <- intersect(mapped_KOs, colnames(shotgun_ko))
    nKO_shotgun[i] <- length(shotgun_KOs)

    if (length(shotgun_KOs) > 0) {
      ko_mean <- rowMeans(shotgun_ko[, shotgun_KOs, drop = FALSE])
      spearman_result <- spearman_with_p(ko_mean, metabolite)
      rho_shotgun[i] <- spearman_result["rho"]
      p_shotgun[i]   <- spearman_result["p"]
    }
  }

  targeted_shotgun <- data.frame(
    Dataset  = shotgun_name, Compound = testable,
    nKO_shotgun = nKO_shotgun, rho_shotgun = rho_shotgun, p_shotgun = p_shotgun
  )
  targeted_shotgun$fdr_shotgun <- p.adjust(targeted_shotgun$p_shotgun, method = "BH")
  targeted_shotgun <- targeted_shotgun[targeted_shotgun$nKO_shotgun >= 1, ]

  median_targeted_shotgun <- median(abs(targeted_shotgun$rho_shotgun), na.rm = TRUE)
  fwrite(targeted_shotgun, file.path(shotgun_res_dir, "shotgun_targeted_correlation.tsv"), sep = "\t")
  cat("  compounds with >=1 mapped KO:", nrow(targeted_shotgun),
      "| median |rho|:", round(median_targeted_shotgun, 3), "\n\n")
} else {
  cat("  no testable compounds - skipping\n\n")
}


###################
# 7: compare shotgun with the 16S predictions (figures)
###################

cat("=== 7: shotgun vs 16S ===\n")
dir.create("plots/final", showWarnings = FALSE)

### maximum KO-metabolite correlation -> Supplementary Figure S1 
max_files <- list.files("results", pattern = "max_ko_metabolite_correlation.tsv",
                        recursive = TRUE, full.names = TRUE)
max_files <- max_files[!grepl("shotgun", max_files)]

if (length(max_files) > 0) {
  # read every 16S max-correlation file and stack them
  max_16s_list <- list()
  for (file_path in max_files) {
    max_16s_list[[file_path]] <- fread(file_path)
  }
  max_16s <- do.call(rbind, max_16s_list)

  # one block per tool (16S) and one for shotgun, then stack
  pic_block     <- data.frame(Dataset = paste0(max_16s$Dataset, " (16S)"), Method = "PICRUSt2", rho = max_16s$max_rho_picrust2)
  t4f_block     <- data.frame(Dataset = paste0(max_16s$Dataset, " (16S)"), Method = "Tax4Fun2", rho = max_16s$max_rho_tax4fun2)
  shotgun_block <- data.frame(Dataset = paste0(shotgun_name, " (Shotgun)"), Method = "Shotgun",  rho = max_corr_shotgun$max_rho_shotgun)
  max_combined  <- rbind(pic_block, t4f_block, shotgun_block)
  max_combined  <- max_combined[is.finite(max_combined$rho), ]

  fig_max <- ggplot(max_combined, aes(x = Method, y = rho, fill = Method)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.25) +
    facet_wrap(~ Dataset, scales = "free_y") +
    scale_fill_manual(values = tool_colors) +
    labs(y = "Max |Spearman rho| per metabolite", x = "Method") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none", strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1))

  # MANUSCRIPT Supplementary Figure S1
  jpeg("plots/final/figS1_shotgun_vs_16s_max_correlation.jpg", height = 6, width = 11, units = "in", res = 600)
  print(fig_max)
  dev.off()
  cat("  saved: plots/final/figS1_shotgun_vs_16s_max_correlation.jpg\n")
} else {
  cat("  no 16S maximum-correlation results found - run the 16S pipeline first\n")
}

### targeted KO-metabolite correlation -> Supplementary Figure S2 
targeted_files <- list.files("results", pattern = "targeted_correlation.tsv",
                             recursive = TRUE, full.names = TRUE)
targeted_files <- targeted_files[!grepl("shotgun", targeted_files)]

if (length(targeted_files) > 0 && !is.null(targeted_shotgun)) {
  # read every 16S targeted file and stack them
  targeted_16s_list <- list()
  for (file_path in targeted_files) {
    targeted_16s_list[[file_path]] <- fread(file_path)
  }
  targeted_16s <- do.call(rbind, targeted_16s_list)

  pic_block     <- data.frame(Dataset = paste0(targeted_16s$Dataset, " (16S)"), Method = "PICRUSt2", rho = abs(targeted_16s$rho_pic))
  t4f_block     <- data.frame(Dataset = paste0(targeted_16s$Dataset, " (16S)"), Method = "Tax4Fun2", rho = abs(targeted_16s$rho_t4f))
  shotgun_block <- data.frame(Dataset = paste0(shotgun_name, " (Shotgun)"), Method = "Shotgun",  rho = abs(targeted_shotgun$rho_shotgun))
  targeted_combined <- rbind(pic_block, t4f_block, shotgun_block)
  targeted_combined <- targeted_combined[is.finite(targeted_combined$rho), ]

  fig_targeted <- ggplot(targeted_combined, aes(x = Method, y = rho, fill = Method)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.25) +
    facet_wrap(~ Dataset, scales = "free_y") +
    scale_fill_manual(values = tool_colors) +
    labs(y = "|Spearman rho|", x = "Method") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none", strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1))

  # MANUSCRIPT Supplementary Figure S2
  jpeg("plots/final/figS2_shotgun_vs_16s_targeted.jpg", height = 6, width = 11, units = "in", res = 600)
  print(fig_targeted)
  dev.off()
  cat("  saved: plots/final/figS2_shotgun_vs_16s_targeted.jpg\n")
} else {
  cat("  no 16S targeted results (or no shotgun targeted) - skipping Supplementary Figure S2\n")
}


###################
# 8: summary
###################

cat("\n=== Shotgun analysis complete ===\n")
cat("  dataset:", shotgun_name, "| samples:", length(common_samples),
    "| KOs:", ncol(shotgun_ko), "| metabolites:", ncol(shotgun_mtb), "\n")
cat("  median max |rho| (shotgun):", round(median_max_shotgun, 3), "\n")
if (!is.null(targeted_shotgun)) {
  cat("  median targeted |rho| (shotgun):",
      round(median(abs(targeted_shotgun$rho_shotgun), na.rm = TRUE), 3), "\n")
}
cat("  figures: ../plots/shotgun_vs_16s_max_correlation.jpg (Supplementary Figure S1),",
    "shotgun_vs_16s_targeted.jpg (Supplementary Figure S2)\n")
