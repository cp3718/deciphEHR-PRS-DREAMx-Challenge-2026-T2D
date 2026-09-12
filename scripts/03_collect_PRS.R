#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# Arguments
# ============================================================

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

prs_dir <- get_arg(
  "--prs-dir",
  "outputs/prs/t2d"
)

out_dir <- get_arg(
  "--out-dir",
  "outputs/prs/t2d/combined"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# P-value thresholds
#
# Must match 01_prepare_sumstats.R and 02_generate_prs.sbatch.
# ============================================================

thresholds <- c(
  "1e5",
  "1e3",
  "0.01",
  "0.05",
  "0.1",
  "0.5",
  "1"
)

# ============================================================
# Combine chromosome-level PLINK scores
#
# PLINK reports SCORE1_AVG for each chromosome.
#
# To recover the genome-wide average score:
#
#   chromosome numerator =
#       SCORE1_AVG * ALLELE_CT
#
#   genome-wide PRS =
#       sum(chromosome numerators) /
#       sum(chromosome ALLELE_CT)
#
# This is the aggregation used in the reconstructed challenge
# workflow and was validated against the historical .sscore
# files.
# ============================================================

all_thresholds <- vector(
  "list",
  length(thresholds)
)

names(all_thresholds) <- thresholds

for (threshold in thresholds) {

  cat(
    "\nCombining threshold:",
    threshold,
    "\n"
  )

  chr_scores <- vector(
    "list",
    22
  )

  reference_ids <- NULL

  for (chr in 1:22) {

    score_file <- file.path(
      prs_dir,
      paste0(
        "t2d_chr",
        chr,
        "_p",
        threshold,
        "_prs.sscore"
      )
    )

    if (!file.exists(score_file)) {
      stop(
        "Missing file: ",
        score_file
      )
    }

    x <- fread(
      score_file
    )

    # PLINK 2 normally writes the sample identifier as #IID.
    if ("#IID" %in% names(x)) {
      setnames(
        x,
        "#IID",
        "IID"
      )
    }

    required_cols <- c(
      "IID",
      "ALLELE_CT",
      "SCORE1_AVG"
    )

    missing_cols <- setdiff(
      required_cols,
      names(x)
    )

    if (length(missing_cols) > 0) {
      stop(
        "Missing required columns in ",
        score_file,
        ": ",
        paste(
          missing_cols,
          collapse = ", "
        )
      )
    }

    x <- x[
      ,
      .(
        IID = as.character(IID),
        ALLELE_CT = as.numeric(ALLELE_CT),
        SCORE1_AVG = as.numeric(SCORE1_AVG)
      )
    ]

    # --------------------------------------------------------
    # Integrity checks
    # --------------------------------------------------------

    if (anyDuplicated(x$IID)) {
      stop(
        "Duplicate IID values in ",
        score_file
      )
    }

    if (anyNA(x$ALLELE_CT) ||
        anyNA(x$SCORE1_AVG)) {

      stop(
        "Missing ALLELE_CT or SCORE1_AVG values in ",
        score_file
      )
    }

    if (any(!is.finite(x$ALLELE_CT)) ||
        any(!is.finite(x$SCORE1_AVG))) {

      stop(
        "Non-finite score values in ",
        score_file
      )
    }

    if (any(x$ALLELE_CT <= 0)) {
      stop(
        "Non-positive ALLELE_CT detected in ",
        score_file
      )
    }

    # Check that each chromosome contains the same individuals.
    current_ids <- sort(x$IID)

    if (is.null(reference_ids)) {

      reference_ids <- current_ids

    } else if (!identical(
      reference_ids,
      current_ids
    )) {

      stop(
        "Sample IDs differ across chromosomes at threshold ",
        threshold,
        ". Problem detected at chromosome ",
        chr,
        "."
      )
    }

    # --------------------------------------------------------
    # Recover chromosome-level numerator
    # --------------------------------------------------------

    x[
      ,
      weighted_score :=
        SCORE1_AVG * ALLELE_CT
    ]

    chr_scores[[chr]] <- x[
      ,
      .(
        IID,
        ALLELE_CT,
        weighted_score
      )
    ]

    cat(
      "  Chromosome ",
      chr,
      ": ",
      nrow(x),
      " individuals\n",
      sep = ""
    )
  }

  # ==========================================================
  # Genome-wide aggregation
  # ==========================================================

  genome <- rbindlist(
    chr_scores
  )

  genome <- genome[
    ,
    .(
      total_allele_ct =
        sum(ALLELE_CT),

      total_weighted_score =
        sum(weighted_score)
    ),
    by = IID
  ]

  if (any(genome$total_allele_ct <= 0)) {
    stop(
      "Non-positive genome-wide allele count at threshold ",
      threshold
    )
  }

  genome[
    ,
    prs :=
      total_weighted_score /
      total_allele_ct
  ]

  if (anyNA(genome$prs) ||
      any(!is.finite(genome$prs))) {

    stop(
      "Invalid genome-wide PRS at threshold ",
      threshold
    )
  }

  genome[
    ,
    threshold := threshold
  ]

  all_thresholds[[threshold]] <- genome[
    ,
    .(
      IID,
      threshold,
      prs,
      total_allele_ct
    )
  ]

  cat(
    "  Genome-wide PRS complete: ",
    nrow(genome),
    " individuals\n",
    sep = ""
  )
}

# ============================================================
# Long-format output
# ============================================================

prs_long <- rbindlist(
  all_thresholds
)

long_file <- file.path(
  out_dir,
  "t2d_prs_thresholds_long.tsv"
)

fwrite(
  prs_long,
  long_file,
  sep = "\t"
)

# ============================================================
# Wide-format output for downstream modeling
# ============================================================

prs_wide <- dcast(
  prs_long,
  IID ~ threshold,
  value.var = "prs"
)

setnames(
  prs_wide,
  old = thresholds,
  new = paste0(
    "prs_p",
    thresholds
  )
)

wide_file <- file.path(
  out_dir,
  "t2d_prs_thresholds_wide.tsv"
)

fwrite(
  prs_wide,
  wide_file,
  sep = "\t"
)

# ============================================================
# Final summary
# ============================================================

cat(
  "\n----------------------------------------\n"
)

cat(
  "T2D PRS collection complete\n"
)

cat(
  "----------------------------------------\n"
)

cat(
  "Individuals:",
  uniqueN(prs_long$IID),
  "\n"
)

cat(
  "Thresholds:",
  uniqueN(prs_long$threshold),
  "\n"
)

cat(
  "Long-format output:",
  long_file,
  "\n"
)

cat(
  "Wide-format output:",
  wide_file,
  "\n\n"
)
