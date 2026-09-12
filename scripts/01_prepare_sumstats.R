#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ------------------------------------------------------------
# Argument parser
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  
  if (length(hit) == 0) {
    return(default)
  }
  
  if (hit[1] == length(args)) {
    stop("Missing value for ", flag)
  }
  
  args[hit[1] + 1]
}

disease <- tolower(get_arg("--disease"))
input   <- get_arg("--input")
object  <- get_arg("--object")
out_dir <- get_arg("--out-dir", "outputs/sumstats")

if (is.null(disease)) {
  stop("Required argument: --disease")
}

if (is.null(input)) {
  stop("Required argument: --input")
}

if (!disease %in% c("t2d")) {
  stop(
    "This workflow currently supports T2D only."
  )
}

# ------------------------------------------------------------
# Load summary statistics
# ------------------------------------------------------------

load_sumstats <- function(path, object = NULL) {
  
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    return(as.data.table(readRDS(path)))
  }
  
  if (grepl("\\.(rdata|rds)$", path, ignore.case = TRUE)) {
    
    env <- new.env(parent = emptyenv())
    loaded <- load(path, envir = env)
    
    if (is.null(object)) {
      stop(
        "For .RData/.rda input, specify the object with --object. ",
        "Objects found: ", paste(loaded, collapse = ", ")
      )
    }
    
    if (!object %in% loaded) {
      stop(
        "Object '", object, "' not found. Objects found: ",
        paste(loaded, collapse = ", ")
      )
    }
    
    return(as.data.table(env[[object]]))
  }
  
  stop("Input must be .rds, .RData, or .rda")
}

ss <- load_sumstats(input, object)

cat("Disease: ", disease, "\n", sep = "")
cat("Input variants: ", format(nrow(ss), big.mark = ","), "\n", sep = "")

# ------------------------------------------------------------
# T2D configuration
#
# No new allele harmonization / ambiguous-SNP filtering added.
# ------------------------------------------------------------

if (disease == "t2d") {
  
  required <- c(
    "chromosome",
    "base_pair_location",
    "effect_allele",
    "other_allele",
    "odds_ratio",
    "p_value"
  )
  
  missing_cols <- setdiff(required, names(ss))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required T2D columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  thresholds <- c(
    p1e5  = 1e-5,
    p1e3  = 1e-3,
    p0.01 = 0.01,
    p0.05 = 0.05,
    p0.1  = 0.1,
    p0.5  = 0.5,
    p1    = 1
  )
  
  # Match original representation
  ss[, chromosome := as.character(chromosome)]
  
  ss[, variant_id := paste0(
    "c", chromosome, "_",
    base_pair_location, "_snv_",
    effect_allele, "_",
    other_allele
  )]
  
  # Original effect-size transformation
  ss[, beta := log(odds_ratio)]
  
  n_input <- nrow(ss)
  
  # ----------------------------------------------------------
  # Original filtering
  # ----------------------------------------------------------
  
  score_base <- ss[
    chromosome %in% as.character(1:22) &
      !is.na(base_pair_location) &
      !is.na(effect_allele) &
      effect_allele %in% c("A", "C", "G", "T") &
      other_allele %in% c("A", "C", "G", "T") &
      !is.na(p_value) &
      p_value > 0 &
      p_value <= 1,
    .(
      chromosome,
      variant_id,
      effect_allele,
      beta = as.numeric(beta),
      p_value = as.numeric(p_value)
    )
  ]
  
  n_after_filter <- nrow(score_base)
  
  # Original duplicate handling:
  # sort by ID and p-value, retain smallest p-value per ID
  setorder(score_base, variant_id, p_value)
  
  n_duplicate_rows <- sum(duplicated(score_base$variant_id))
  
  score_base <- score_base[
    !duplicated(variant_id)
  ]
  
  n_final <- nrow(score_base)
}

# ------------------------------------------------------------
# Output directories
# ------------------------------------------------------------

disease_dir <- file.path(out_dir, disease)
score_dir   <- file.path(disease_dir, "score_files")

dir.create(
  score_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Save cleaned score base
# ------------------------------------------------------------

base_file <- file.path(
  disease_dir,
  paste0(disease, "_score_base.tsv.gz")
)

fwrite(
  score_base,
  base_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat(
  "Wrote score base: ",
  base_file,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# Write chromosome × threshold PLINK score files
# ------------------------------------------------------------

manifest <- list()

i <- 1L

for (chr in 1:22) {
  
  for (label in names(thresholds)) {
    
    threshold <- thresholds[[label]]
    
    score <- score_base[
      chromosome == as.character(chr) &
        p_value <= threshold,
      .(
        variant_id,
        effect_allele,
        beta
      )
    ]
    
    score_file <- file.path(
      score_dir,
      paste0(
        disease,
        "_chr", chr,
        "_", label,
        ".score.tsv"
      )
    )
    
    fwrite(
      score,
      score_file,
      sep = "\t",
      quote = FALSE
    )
    
    manifest[[i]] <- data.table(
      disease = disease,
      chromosome = chr,
      threshold_label = label,
      p_threshold = threshold,
      n_variants = nrow(score),
      score_file = score_file
    )
    
    i <- i + 1L
    
    cat(
      disease,
      " chr", chr,
      " ", label,
      ": ",
      format(nrow(score), big.mark = ","),
      " variants\n",
      sep = ""
    )
  }
}

manifest <- rbindlist(manifest)

# ------------------------------------------------------------
# QC summaries
# ------------------------------------------------------------

qc <- data.table(
  disease = disease,
  n_input = n_input,
  n_after_original_filters = n_after_filter,
  n_duplicate_variant_ids_removed = n_duplicate_rows,
  n_score_base = n_final
)

fwrite(
  qc,
  file.path(disease_dir, "summary_stats_qc.tsv"),
  sep = "\t"
)

fwrite(
  manifest,
  file.path(disease_dir, "score_file_manifest.tsv"),
  sep = "\t"
)

cat("\nSummary\n")
print(qc)

cat(
  "\nGenerated ",
  nrow(manifest),
  " chromosome × threshold score files.\n",
  sep = ""
)
