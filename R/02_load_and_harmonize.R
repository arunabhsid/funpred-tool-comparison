##############################################################
# 02_load_and_harmonize.R  
##############################################################
#
# For each dataset we read three tables (PICRUSt2 KOs, Tax4Fun2 KOs &
# metabolites). Then, turn each into a matrix with samples in the rows, 
# and keep only the common samples. Everything is stored in one list,
# all_data, that the later scripts read from.
#
# For KANG the KO tables use raw sequencing IDs while the metabolomics 
# uses study IDs, so we rename the KO columns using its metadata.tsv 
# file.
# 
# The other three datasets already match.

##############################################################

source("R/config.R")
library(data.table)

###################
# Helper functions
###################

# Line up a KO table's sample names with the metabolomics. Two fixes:
#   - strip a leading "X" that R adds to names starting with a digit
#     ("X12b" -> "12b"); Tax4Fun2 sometimes writes these.
#   - some studies name the same samples differently in the 16S data and in the
#     metabolomics. If the dataset comes with a lookup file (metadata.tsv, which
#     has a RawSampleID column and a Sample column side by side, like a use it 
#     to rename the KO columns from the raw sequencing IDs over to the matching
#     study sample IDs. Only KANG needs this.

fix_sample_names <- function(ko, ds) {

  names(ko) <- sub("^X([0-9])", "\\1", names(ko)) 

  if (!is.null(ds$metadata_fp) && file.exists(ds$metadata_fp)) {
    meta   <- fread(ds$metadata_fp, select = c("Sample", "RawSampleID"))
    id_map <- setNames(meta$Sample, meta$RawSampleID)   # raw ID -> study ID
    hit    <- names(ko) %in% names(id_map)
    names(ko)[hit] <- id_map[names(ko)[hit]]
    cat("KANG: renamed", sum(hit), "sample columns using the lookup table\n")
  }
  ko
}

# Find the KO output file Tax4Fun2 wrote (its name varies a little run to run).
find_tax4fun2_output <- function(ds) {
  files <- list.files(file.path(ds$res_dir, "Tax4Fun2_blast"),
                      recursive = TRUE, full.names = TRUE)
  files <- files[!grepl("backup", files, ignore.case = TRUE)]
  hits  <- grep("Functional|Profile|UProC|Predicted|KO", basename(files), ignore.case = TRUE)
  files[hits][1]
}

# Turn a KO table (a "KO" column plus one column per sample) into a matrix
# with samples in the rows and KOs in the columns.
ko_table_to_matrix <- function(ko, ds) {
  names(ko)[1] <- "KO"
  ko <- fix_sample_names(ko, ds)

  # Tax4Fun2 sometimes reads a sample column in as text; make those numbers
  for (col in setdiff(names(ko), "KO")) {
    if (is.character(ko[[col]])) ko[[col]] <- as.numeric(ko[[col]])
  }

  m <- as.matrix(ko[, -1, with = FALSE])   # rows = KOs, cols = samples
  rownames(m) <- ko$KO
  m <- t(m)                                # flip: rows = samples, cols = KOs
  colnames(m) <- clean_ko_ids(colnames(m))
  m
}

# Convert wide metabolite table into a matrix of samples x KEGG compounds.
# Relabel each metabolite by its KEGG compound ID, so the metabolites can later
# be matched to the enzymes that act on them. Metabolites without a KEGG ID are
# dropped, and where two or more metabolites map to the same KEGG compound their
# abundances are averaged into one value.

metabolites_to_matrix <- function(mtb, mtb_map) {
  long <- melt(mtb, id.vars = "Sample", variable.name = "Compound", value.name = "abundance")
  long <- merge(long, unique(mtb_map[, .(Compound, KEGG)]), by = "Compound", all.x = TRUE)
  long <- long[!is.na(KEGG) & KEGG != ""]
  wide <- dcast(long, Sample ~ KEGG, value.var = "abundance", fun.aggregate = mean)

  m <- as.matrix(wide[, -"Sample"])
  rownames(m) <- wide$Sample
  colnames(m) <- clean_compound_ids(colnames(m))
  m
}

###################
# Load each dataset
###################

all_data <- list()

for (ds in datasets) {

  cat("\n=== Loading", ds$name, "===\n")

  # 1. read the three tables 
  mtb     <- fread(ds$mtb_fp)
  mtb_map <- fread(ds$mtb_map_fp)
  pic_fp  <- file.path(ds$res_dir, "picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz")
  ko_pic  <- fread(pic_fp)
  ko_t4f  <- fread(find_tax4fun2_output(ds))

  # 2. convert into a matrix with samples in the rows
  metab <- metabolites_to_matrix(mtb, mtb_map)
  pic   <- ko_table_to_matrix(ko_pic, ds)
  t4f   <- ko_table_to_matrix(ko_t4f, ds)
  cat("  PICRUSt2:", ncol(pic), "KOs | Tax4Fun2:", ncol(t4f),
      "KOs | metabolites:", ncol(metab), "compounds\n")

  # 3. keep only the samples that appear in all three
  common <- Reduce(intersect, list(rownames(pic), rownames(t4f), rownames(metab)))
  if (length(common) == 0) {
    cat("  no samples in common - check the IDs. Skipping.\n")
    next
  }
  cat("  samples in common:", length(common), "\n")

  # 4. subset to those samples, make numeric, drop columns that never vary
  pic   <- to_clean_matrix(pic[common, , drop = FALSE])
  t4f   <- to_clean_matrix(t4f[common, , drop = FALSE])
  metab <- to_clean_matrix(metab[common, , drop = FALSE])

  # 5. KOs that are common to both PICRUSt & Tax4Fun 
  common_kos <- intersect(colnames(pic), colnames(t4f))
  cat("  KOs shared by both tools:", length(common_kos), "\n")

  # 6. save
  all_data[[ds$name]] <- list(
    name       = ds$name,
    res_dir    = ds$res_dir,
    pic_full   = pic,                              # PICRUSt2, all KOs
    t4f_full   = t4f,                              # Tax4Fun2, all KOs
    pic_common = pic[, common_kos, drop = FALSE],  # PICRUSt2, shared KOs only
    t4f_common = t4f[, common_kos, drop = FALSE],  # Tax4Fun2, shared KOs only
    mtb        = metab,                            # metabolites
    n_samples  = length(common)
  )
  cat("  done:", ds$name, "\n")
}

cat("\nAll datasets loaded. Proceed to script 03.\n")
cat("Ready:", paste(names(all_data), collapse = ", "), "\n\n")
