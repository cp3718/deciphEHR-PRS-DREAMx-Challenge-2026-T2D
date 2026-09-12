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

  if (!length(hit)) {
    return(default)
  }

  if (hit[1] == length(args)) {
    stop("Missing value for ", flag)
  }

  args[hit[1] + 1]
}

prs_file <- get_arg("--prs")

out_dir <- get_arg(
  "--out-dir",
  "outputs/ensemble/t2d"
)

if (is.null(prs_file)) {
  stop("Required: --prs")
}

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Load threshold-specific PRSs
# ============================================================

prs <- fread(prs_file)

required_cols <- c(
  "IID",
  "prs_p1",
  "prs_p0.5",
  "prs_p0.1"
)

missing_cols <- setdiff(
  required_cols,
  names(prs)
)

if (length(missing_cols)) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

if (anyDuplicated(prs$IID)) {
  stop("Duplicate IID values detected in PRS table.")
}

score_cols <- c(
  "prs_p1",
  "prs_p0.5",
  "prs_p0.1"
)

if (
  anyNA(prs[, ..score_cols]) ||
  any(
    !is.finite(
      as.matrix(
        prs[, ..score_cols]
      )
    )
  )
) {
  stop(
    "Missing or non-finite threshold-specific PRS values detected."
  )
}

# ============================================================
# Historical ensemble exploration
#
# During the challenge, several weighted combinations of the
# best-performing broad-threshold PRSs were explored
# empirically.
#
# These weights were heuristic rather than learned by a
# formal optimization procedure.
# ============================================================

# 50% P<=1 + 50% P<=0.5
prs[
  ,
  ens_50_50 :=
    0.50 * prs_p1 +
    0.50 * prs_p0.5
]

# 70% P<=1 + 30% P<=0.5
prs[
  ,
  ens_70_30 :=
    0.70 * prs_p1 +
    0.30 * prs_p0.5
]

# 60% P<=1 + 30% P<=0.5 + 10% P<=0.1
prs[
  ,
  ens_60_30_10 :=
    0.60 * prs_p1 +
    0.30 * prs_p0.5 +
    0.10 * prs_p0.1
]

# Final historical ensemble:
#
# 60% P<=1
# 25% P<=0.5
# 15% P<=0.1
prs[
  ,
  ens_60_25_15 :=
    0.60 * prs_p1 +
    0.25 * prs_p0.5 +
    0.15 * prs_p0.1
]

# ============================================================
# QC
# ============================================================

ensemble_cols <- c(
  "ens_50_50",
  "ens_70_30",
  "ens_60_30_10",
  "ens_60_25_15"
)

if (
  anyNA(prs[, ..ensemble_cols]) ||
  any(
    !is.finite(
      as.matrix(
        prs[, ..ensemble_cols]
      )
    )
  )
) {
  stop(
    "Missing or non-finite ensemble values detected."
  )
}

# ============================================================
# Save
#
# Scores remain on their raw scale here.
#
# Standardization is intentionally deferred to model fitting,
# where mean and SD must be estimated from the appropriate
# training cohort and then reused for held-out samples.
# ============================================================

out_file <- file.path(
  out_dir,
  "t2d_ensemble_scores.tsv"
)

fwrite(
  prs[
    ,
    c(
      "IID",
      "prs_p1",
      "prs_p0.5",
      "prs_p0.1",
      ensemble_cols
    ),
    with = FALSE
  ],
  out_file,
  sep = "\t"
)

# ============================================================
# Summary
# ============================================================

cat(
  "\n----------------------------------------\n"
)

cat(
  "T2D ensemble construction complete\n"
)

cat(
  "----------------------------------------\n"
)

cat(
  "Individuals:",
  nrow(prs),
  "\n"
)

cat(
  "Final ensemble:",
  "0.60 * prs_p1 + 0.25 * prs_p0.5 + 0.15 * prs_p0.1\n"
)

cat(
  "Output:",
  out_file,
  "\n"
)
