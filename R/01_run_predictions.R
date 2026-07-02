##############################################################
# 01_run_predictions.R  -  predict KO abundances from 16S data
##############################################################
#
# For each dataset, generate predicted KO tables:
#   - PICRUSt2   run in the terminal
#   - Tax4Fun2   run here in R 
##############################################################

source("R/config.R")
library(data.table)
library(biomformat)
library(Matrix)
library(Tax4Fun2)

###############################################################
# Prepare HE_INFANTS: keep only the 12-month timepoint
###############################################################
# HE_INFANTS is longitudinal - the same infants are sampled at 
# 2, 4, 6 and 12 months. We keep only one timepoint.

# pull the HE entry out of the config list
for (d in datasets) {
  if (d$name == "HE_INFANTS_MFGM_2019") {
    he <- d
  }
}

raw_dir <- "data/HE_raw"
out_dir <- he$data_dir

# skip if we've already prepared it from previous run
if (file.exists(he$mtb_fp) && file.exists(he$biom_fp) && file.exists(he$fasta_fp)) {
  cat("HE_INFANTS already prepared. Skipping.\n\n")
} else {
  cat("Preparing HE_INFANTS (filtering to 12 months)...\n")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # 1. which samples are at 12 months?
  meta <- fread(file.path(raw_dir, "metadata.tsv"))
  meta_12 <- meta[Age == 12]
  keep_samples <- meta_12$Sample
  cat("samples at 12 months:", length(keep_samples), "\n")
  fwrite(meta_12, file.path(out_dir, "metadata.tsv"), sep = "\t")

  # 2. metabolite table: keep those samples
  mtb <- fread(file.path(raw_dir, "mtb.tsv"))
  mtb_12 <- mtb[Sample %in% keep_samples]
  fwrite(mtb_12, file.path(out_dir, "mtb.tsv"), sep = "\t")

  # 3. compound map: just copy it across (nothing to filter)
  #    the raw file is named mtb_map.tsv. we save it as mtb.map.tsv
  file.copy(file.path(raw_dir, "mtb_map.tsv"),
            file.path(out_dir, "mtb.map.tsv"), overwrite = TRUE)

  # 4. BIOM (ASV counts): keep those samples, then drop ASVs that are now empty
  bt     <- read_biom(file.path(raw_dir, "feature-table.biom"))
  counts <- as(biom_data(bt), "matrix")        # ASVs in rows, samples in columns

  # the BIOM sample names carry a study prefix like "12021."; strip it so they
  # match the metadata names ("12021.BF.13.1" -> "BF.13.1")
  colnames(counts) <- sub("^[0-9]+\\.", "", colnames(counts))

  counts_12 <- counts[, colnames(counts) %in% keep_samples, drop = FALSE]
  counts_12 <- counts_12[rowSums(counts_12) > 0, , drop = FALSE]
  cat("ASVs kept (non-zero at 12 months):", nrow(counts_12), "of", nrow(counts), "\n")

  write_biom(make_biom(counts_12), file.path(out_dir, "feature-table.biom"))

  # 5. FASTA: keep the sequences for the ASVs we retained. Each line belongs to
  #    the most recent header above it, so cumsum() over the header positions
  #    numbers every line by its record; we then keep the lines whose record is
  #    one we want. This works whether a sequence sits on one line or wraps
  keep_asvs   <- rownames(counts_12)
  fa          <- readLines(file.path(raw_dir, "dna-sequences.fasta"))
  is_header   <- grepl("^>", fa)
  ids         <- sub("^>", "", fa[is_header])         # one per record, in order
  record      <- cumsum(is_header)                    # record number of each line
  keep_record <- record %in% which(ids %in% keep_asvs)
  writeLines(fa[keep_record], file.path(out_dir, "dna-sequences.fasta"))

  cat("  done:", length(keep_samples), "samples,", nrow(counts_12), "ASVs,",
      ncol(mtb_12) - 1, "metabolites\n\n")
}

###############################################################
# Run predictions for each dataset
###############################################################

for (ds in datasets) {

  cat("\n=== Loading", ds$name, "===\n")
  dir.create(ds$res_dir, recursive = TRUE)

  ##########################################
  # PICRUSt2 (run in the terminal, not in R)
  ##########################################
  # We don't run PICRUSt2 from here. we just check whether its output exists.

  picrust2_ko_fp <- file.path(ds$res_dir,
    "picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz")

  if (!file.exists(picrust2_ko_fp)) {
    cat("  PICRUSt2 output not found. Run this in your terminal, then re-run current script:\n\n")
    cat("    picrust2_pipeline.py \\\n")
    cat("      -s", ds$fasta_fp, "\\\n")
    cat("      -i", ds$biom_fp, "\\\n")
    cat("      -o", file.path(ds$res_dir, "picrust2_out"), "\\\n")
    cat("      -p 4\n\n")
    next
  } else {
  cat("  PICRUSt2 output found.\n")
  }

##########################################
  # Tax4Fun2 (runs in R, calls BLAST)
##########################################
  
  # skip if already done

  t4f_tmp  <- file.path(ds$res_dir, "Tax4Fun2_blast")
  existing <- list.files(t4f_tmp, recursive = TRUE, full.names = TRUE,
                         pattern = "(?i)(Functional|Profile|UProC|Predicted|KO)")
  existing <- existing[!grepl("backup", existing, ignore.case = TRUE)]
  if (length(existing) > 0) {
    cat("Tax4Fun2 output found.\n")
    cat("  === done:", ds$name, "===\n")
    next
  }

  cat("Running Tax4Fun2...\n")
  dir.create(t4f_tmp, recursive = TRUE)

  # Tax4Fun2 expects blastn and makeblastdb inside its reference folder, not
  # wherever they're actually installed. So we just copy them there
  blastn_path      <- system("command -v blastn", intern = TRUE)
  makeblastdb_path <- system("command -v makeblastdb", intern = TRUE)
  if (length(blastn_path) == 0 || length(makeblastdb_path) == 0) {

    stop("BLAST+ not found. Install it: sudo apt-get install -y ncbi-blast+")
  }

  blast_dir <- file.path(path.expand(tax4fun2_ref_dir), "blast_bin", "bin")
  dir.create(blast_dir, recursive = TRUE)
  file.copy(blastn_path,      file.path(blast_dir, "blastn"),overwrite = TRUE) # also check symlink
  file.copy(makeblastdb_path, file.path(blast_dir, "makeblastdb"), overwrite = TRUE)

  # Tax4Fun2 wants the ASV table as a tab-delimited OTU file
  bt <- read_biom(ds$biom_fp)
  mx <- as.matrix(biom_data(bt)) # ASVs in rows, samples in columns

  otu <- data.frame(OTU = rownames(mx), mx, check.names = FALSE)
  otu_fp <- file.path(ds$res_dir, "otu_table.txt")
  write.table(otu, otu_fp, sep = "\t", quote = FALSE, row.names = FALSE)

  # BLAST the sequences against the reference, then make the predictions
  runRefBlast(
    path_to_otus           = ds$fasta_fp,
    path_to_reference_data = tax4fun2_ref_dir,
    path_to_temp_folder    = t4f_tmp,
    database_mode          = "Ref99NR",
    use_force              = TRUE,
    num_threads            = 4
  )
  makeFunctionalPrediction(
    path_to_otu_table        = otu_fp,
    path_to_reference_data   = tax4fun2_ref_dir,
    path_to_temp_folder      = t4f_tmp,
    database_mode            = "Ref99NR",
    normalize_by_copy_number = TRUE
  )

  # find the KO prediction file Tax4Fun2 wrote (its name varies a bit run to run)
  all_files <- list.files(t4f_tmp, recursive = TRUE, full.names = TRUE)
  hits <- grep("Functional|Profile|UProC|Predicted|KO", basename(all_files), ignore.case = TRUE)
  ko_files <- all_files[hits]   # every file whose name matched
  ko_file  <- ko_files[1]       # just take the first one

  cat("  Tax4Fun2 output:", ko_file, "\n")
  cat("  === Done:", ds$name, "===\n")
}

cat("\nAll datasets processed. Proceed to script 02.\n")
