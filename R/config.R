# =============================================================
# config.R - run this first
# =============================================================
#
# Manually edit this file once. Point the wd and dataset names

############### 
# 1. Datasets
############### 

# Details:
#   name        short label, used in Muller et al.
#   data_dir    local dir, containing datasets
#   biom_fp     ASV feature table, BIOM format (input to PICRUSt2 / Tax4Fun2)
#   fasta_fp    ASV representative sequences, FASTA
#   mtb_fp      metabolite abundance table, TSV, one row per sample
#   mtb_map_fp  metabolite-to-KEGG table, TSV, with Compound and KEGG columns
#   res_dir     results are written here
#
# mtb.tsv needs a "Sample" column followed by one column per metabolite.
# mtb.map.tsv needs at least "Compound" and "KEGG" columns.

datasets <- list(
  
  list(
    name       = "KIM_ADENOMAS_2020",
    data_dir   = "data/KIM_ADENOMAS_2020",
    biom_fp    = "data/KIM_ADENOMAS_2020/feature-table.biom",
    fasta_fp   = "data/KIM_ADENOMAS_2020/dna-sequences.fasta",
    mtb_fp     = "data/KIM_ADENOMAS_2020/mtb.tsv",
    mtb_map_fp = "data/KIM_ADENOMAS_2020/mtb.map.tsv",
    res_dir    = "results/KIM_ADENOMAS_2020"
  ),
  
  list(
    name        = "KANG_AUTISM_2018",
    data_dir    = "data/KANG_AUTISM_2018",
    biom_fp     = "data/KANG_AUTISM_2018/feature-table.biom",
    fasta_fp    = "data/KANG_AUTISM_2018/dna-sequences.fasta",
    mtb_fp      = "data/KANG_AUTISM_2018/mtb.tsv",
    mtb_map_fp  = "data/KANG_AUTISM_2018/mtb.map.tsv",
    metadata_fp = "data/KANG_AUTISM_2018/metadata.tsv",
    res_dir     = "results/KANG_AUTISM_2018"
  ),
  
  list(
    name       = "JACOBS_IBD_FAMILIES_2016",
    data_dir   = "data/JACOBS_IBD_FAMILIES_2016",
    biom_fp    = "data/JACOBS_IBD_FAMILIES_2016/feature-table.biom",
    fasta_fp   = "data/JACOBS_IBD_FAMILIES_2016/dna-sequences.fasta",
    mtb_fp     = "data/JACOBS_IBD_FAMILIES_2016/mtb.tsv",
    mtb_map_fp = "data/JACOBS_IBD_FAMILIES_2016/mtb.map.tsv",
    res_dir    = "results/JACOBS_IBD_FAMILIES_2016"
  ),
  
  
  list(
    name       = "HE_INFANTS_MFGM_2019",
    data_dir   = "data/HE_INFANTS_MFGM_2019",
    biom_fp    = "data/HE_INFANTS_MFGM_2019/feature-table.biom",
    fasta_fp   = "data/HE_INFANTS_MFGM_2019/dna-sequences.fasta",
    mtb_fp     = "data/HE_INFANTS_MFGM_2019/mtb.tsv",
    mtb_map_fp = "data/HE_INFANTS_MFGM_2019/mtb.map.tsv",
    res_dir    = "results/HE_INFANTS_MFGM_2019"
  )
)

##############################
# 2. Tax4Fun2 reference data
##############################

# Download this once; it's the same for every dataset

if (!dir.exists("downloads/Tax4Fun2_ReferenceData_v2")) {
  system("mkdir -p downloads")
  system("wget -P downloads https://zenodo.org/records/10035668/files/Tax4Fun2_ReferenceData_v2.tar.gz")
  system("tar -xzvf downloads/Tax4Fun2_ReferenceData_v2.tar.gz -C downloads")
}

# This gives:
# ~/Documents/funpred-tool-comparison/downloads/Tax4Fun2_ReferenceData_v2/
#     ├── Ref99NR/
#     ├── Ref100NR/
#     ├── (other mapping / genome files)

# Install BLAST+ on your system (Ubuntu):
# BLAST+ (used by Tax4Fun2). If it's not on your PATH, print how to install it.
if (Sys.which("blastn") == "" || Sys.which("makeblastdb") == "") {
  cat("BLAST+ not found. Install it once on your system (Ubuntu):\n",
      "  sudo apt-get update\n",
      "  sudo apt-get install -y ncbi-blast+\n",
      "Check with:\n",
      "  which blastn; which makeblastdb\n",
      "  blastn -version; makeblastdb -version\n")
}

tax4fun2_ref_dir <- "downloads/Tax4Fun2_ReferenceData_v2"


#################################
# 3. KEGG compound-to-KO mapping
#################################
# Links each metabolite to the enzymes (KOs) that act on it, via KEGG
# reactions. It is built once in script 05 and reused everywhere.

kegg_mapping_fp <- "compound_to_ko_mapping.tsv"


###################
# Packages used
###################

# CRAN
cran_pkgs <- c("data.table", "ggplot2", "ggrepel", "dplyr",
               "tidyr", "readr", "cowplot", "scales")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

# Bioconductor (installed through BiocManager)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("biomformat", "KEGGREST")) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p)
}

# Tax4Fun2 (v1.1.5, the version from Wemheuer et al. 2020) can be 
# installed fromsource. I could not find the original repo; See:
# github.com/fjossandon/Tax4Fun2 for instructions to install
# wget https://github.com/bwemheu/Tax4Fun2/releases/download/1.1.5/Tax4Fun2_1.1.5.tar.gz
if (!requireNamespace("Tax4Fun2", quietly = TRUE)) {
  install.packages("path/to/Tax4Fun2_1.1.5.tar.gz", repos = NULL, type = "source")
}

# PICRUSt2 installation instructions
cat("Run PICRUSt2 from the terminal, not R.\n",
    "If not installed already, install it once with:\n",
    "  conda create -n picrust2 -c bioconda -c conda-forge picrust2=2.4.1\n",
    "  conda activate picrust2\n",
    "Then run the prediction from the terminal (script 02 prints the command).\n")


###################
# Helper functions 
###################

# Clean KO IDs, e.g. "KO:K00001" -> "K00001".
clean_ko_ids <- function(v) {
  v <- toupper(trimws(v))
  sub("^KO:", "K", v)
}

# Clean KEGG compound IDs, e.g. "CPD:C00022" -> "C00022".
clean_compound_ids <- function(v) {
  v <- toupper(trimws(v))
  sub("^CPD:", "C", v)
}

# clr transform
clr_transform <- function(M, pseudocount = 1e-6) {
  L <- log(M + pseudocount)
  L - rowMeans(L, na.rm = TRUE)
}

# Spearman correlation
spearman_with_p <- function(x, y, min_obs = 10) {
  ok <- is.finite(x) & is.finite(y)
  # if too few usable samples, skip test and return NA
  if (sum(ok) < min_obs) return(c(rho = NA_real_, p = NA_real_))
  ct <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  c(rho = as.numeric(ct$estimate), p = ct$p.value)
}

# Remove zero variance columns
to_clean_matrix <- function(M) {
  M2 <- matrix(as.numeric(M), nrow = nrow(M), ncol = ncol(M),
               dimnames = dimnames(M))
  sds <- apply(M2, 2, sd, na.rm = TRUE) # standard deviation of each column
  M2[, sds > 0, drop = FALSE] # keep only columns whose sd > 0
}

############
# Summary 
############

cat("\nConfig file loaded.\n")
cat("Datasets:", length(datasets), "\n")
for (ds in datasets) cat("  -", ds$name, "\n")
cat("\n")