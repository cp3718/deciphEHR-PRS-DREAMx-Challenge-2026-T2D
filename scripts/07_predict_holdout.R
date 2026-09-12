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

prs_dir <- get_arg("--prs-dir")
model_file <- get_arg("--model")
covar_file <- get_arg("--covariates")

out_dir <- get_arg(
  "--out-dir",
  "outputs/final_prediction/t2d"
)

if (
  is.null(prs_dir) ||
  is.null(model_file) ||
  is.null(covar_file)
) {
  stop(
    "Required: --prs-dir --model --covariates"
  )
}

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Helper: combine chromosome-level scores for one threshold
#
# Genome-wide PRS =
#
#   sum(SCORE1_AVG * ALLELE_CT) /
#   sum(ALLELE_CT)
#
# This is identical to the aggregation used for the historical
# challenge submission.
# ============================================================

combine_threshold <- function(threshold) {

  chr_list <- vector(
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
        "Missing score file: ",
        score_file
      )
    }

    x <- fread(
      score_file
    )

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

    if (length(missing_cols)) {
      stop(
        "Missing columns in ",
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

    if (
      anyNA(x$ALLELE_CT) ||
      anyNA(x$SCORE1_AVG)
    ) {
      stop(
        "Missing score values in ",
        score_file
      )
    }

    if (
      any(!is.finite(x$ALLELE_CT)) ||
      any(!is.finite(x$SCORE1_AVG))
    ) {
      stop(
        "Non-finite score values in ",
        score_file
      )
    }

    if (any(x$ALLELE_CT <= 0)) {
      stop(
        "Non-positive ALLELE_CT in ",
        score_file
      )
    }

    # Ensure every chromosome contains the same subjects.
    current_ids <- sort(
      x$IID
    )

    if (is.null(reference_ids)) {

      reference_ids <- current_ids

    } else if (!identical(
      reference_ids,
      current_ids
    )) {

      stop(
        "Sample IDs differ across chromosomes ",
        "for threshold ",
        threshold,
        ". Problem detected at chromosome ",
        chr,
        "."
      )
    }

    # --------------------------------------------------------
    # Recover chromosome numerator
    # --------------------------------------------------------

    x[
      ,
      weighted_score :=
        SCORE1_AVG * ALLELE_CT
    ]

    chr_list[[chr]] <- x[
      ,
      .(
        IID,
        ALLELE_CT,
        weighted_score
      )
    ]
  }

  # ==========================================================
  # Genome-wide aggregation
  # ==========================================================

  z <- rbindlist(
    chr_list
  )

  z <- z[
    ,
    .(
      total_allele_ct =
        sum(ALLELE_CT),

      total_weighted_score =
        sum(weighted_score)
    ),
    by = IID
  ]

  if (any(z$total_allele_ct <= 0)) {
    stop(
      "Invalid genome-wide allele count ",
      "for threshold ",
      threshold
    )
  }

  z[
    ,
    prs :=
      total_weighted_score /
      total_allele_ct
  ]

  if (
    anyNA(z$prs) ||
    any(!is.finite(z$prs))
  ) {
    stop(
      "Invalid genome-wide PRS for threshold ",
      threshold
    )
  }

  z[
    ,
    .(
      IID,
      prs
    )
  ]
}

# ============================================================
# 1. Combine the three threshold-specific holdout PRSs used
#    by the submitted ensemble
# ============================================================

cat(
  "Combining P <= 0.1...\n"
)

p01 <- combine_threshold(
  "0.1"
)

setnames(
  p01,
  "prs",
  "prs_p0.1"
)

cat(
  "Combining P <= 0.5...\n"
)

p05 <- combine_threshold(
  "0.5"
)

setnames(
  p05,
  "prs",
  "prs_p0.5"
)

cat(
  "Combining P <= 1...\n"
)

p1 <- combine_threshold(
  "1"
)

setnames(
  p1,
  "prs",
  "prs_p1"
)

prs <- Reduce(
  function(x, y) {
    merge(
      x,
      y,
      by = "IID",
      all = FALSE
    )
  },
  list(
    p01,
    p05,
    p1
  )
)

cat(
  "Holdout individuals with all three PRSs:",
  nrow(prs),
  "\n"
)

# ============================================================
# 2. Historical final ensemble
#
# Final challenge weighting:
#
#   60% P <= 1
#   25% P <= 0.5
#   15% P <= 0.1
# ============================================================

prs[
  ,
  ens_60_25_15 :=
    0.60 * prs_p1 +
    0.25 * prs_p0.5 +
    0.15 * prs_p0.1
]

# ============================================================
# 3. Load final fitted model and training preprocessing
# ============================================================

bundle <- readRDS(
  model_file
)

required_bundle <- c(
  "model",
  "prs_mean",
  "prs_sd",
  "sex_levels"
)

missing_bundle <- setdiff(
  required_bundle,
  names(bundle)
)

if (length(missing_bundle)) {
  stop(
    "Model bundle is missing: ",
    paste(
      missing_bundle,
      collapse = ", "
    )
  )
}

if (
  !is.finite(bundle$prs_mean) ||
  !is.finite(bundle$prs_sd) ||
  bundle$prs_sd <= 0
) {
  stop(
    "Invalid training-derived PRS standardization parameters."
  )
}

cat(
  "Training PRS mean:",
  format(
    bundle$prs_mean,
    digits = 12
  ),
  "\n"
)

cat(
  "Training PRS SD:",
  format(
    bundle$prs_sd,
    digits = 12
  ),
  "\n"
)

# ============================================================
# 4. Apply TRAINING-derived standardization to holdout
#
# The holdout cohort must never be standardized using its own
# mean or SD.
# ============================================================

prs[
  ,
  prs_std :=
    (
      ens_60_25_15 -
        bundle$prs_mean
    ) /
      bundle$prs_sd
]

# ============================================================
# 5. Load holdout covariates
# ============================================================

covar <- fread(
  covar_file
)

required_covars <- c(
  "sampleid",
  "sex",
  "age_sd",
  paste0(
    "PC",
    1:10,
    "_sd"
  )
)

missing_covars <- setdiff(
  required_covars,
  names(covar)
)

if (length(missing_covars)) {
  stop(
    "Missing holdout covariates: ",
    paste(
      missing_covars,
      collapse = ", "
    )
  )
}

if (anyDuplicated(covar$sampleid)) {
  stop(
    "Duplicate sampleid values in holdout covariate table."
  )
}

# Preserve original holdout order.
covar[
  ,
  original_order := .I
]

hold <- merge(
  covar,
  prs,
  by.x = "sampleid",
  by.y = "IID",
  all.x = TRUE,
  sort = FALSE
)

setorder(
  hold,
  original_order
)

cat(
  "Holdout covariate rows:",
  nrow(covar),
  "\n"
)

cat(
  "Matched holdout rows:",
  sum(
    !is.na(
      hold$ens_60_25_15
    )
  ),
  "\n"
)

if (
  anyNA(
    hold$ens_60_25_15
  )
) {
  stop(
    "Not all holdout individuals matched to PRS scores."
  )
}

# ============================================================
# 6. Match training factor levels
# ============================================================

hold[
  ,
  sex := factor(
    sex,
    levels = bundle$sex_levels
  )
]

if (anyNA(hold$sex)) {
  stop(
    "Holdout contains sex values absent from ",
    "model-training levels."
  )
}

prediction_vars <- c(
  "prs_std",
  "age_sd",
  "sex",
  paste0(
    "PC",
    1:10,
    "_sd"
  )
)

if (
  any(
    !complete.cases(
      hold[
        ,
        ..prediction_vars
      ]
    )
  )
) {
  stop(
    "Missing PRS or covariate values detected ",
    "in holdout prediction data."
  )
}

# ============================================================
# 7. Predict final T2D probability
# ============================================================

hold[
  ,
  T2D_probability :=
    as.numeric(
      predict(
        bundle$model,
        newdata = hold,
        type = "response"
      )
    )
]

if (
  anyNA(hold$T2D_probability) ||
  any(
    !is.finite(
      hold$T2D_probability
    )
  )
) {
  stop(
    "Invalid holdout prediction values generated."
  )
}

# ============================================================
# 8. QC summary
# ============================================================

cat(
  "\nPrediction summary:\n"
)

print(
  summary(
    hold$T2D_probability
  )
)

cat(
  "\nStandardized PRS summary:\n"
)

print(
  summary(
    hold$prs_std
  )
)

# ============================================================
# 9. Save detailed prediction table
# ============================================================

detailed_file <- file.path(
  out_dir,
  "t2d_holdout_predictions_detailed.tsv"
)

fwrite(
  hold[
    ,
    .(
      sampleid,
      prs_p0.1,
      prs_p0.5,
      prs_p1,
      ens_60_25_15,
      prs_std,
      T2D_probability
    )
  ],
  detailed_file,
  sep = "\t"
)

# ============================================================
# 10. Save clean challenge-style prediction file
# ============================================================

prediction_file <- file.path(
  out_dir,
  "t2d_holdout_predictions.tsv"
)

fwrite(
  hold[
    ,
    .(
      sampleid,
      T2D = T2D_probability
    )
  ],
  prediction_file,
  sep = "\t"
)

cat(
  "\nDetailed output:",
  detailed_file,
  "\n"
)

cat(
  "Prediction output:",
  prediction_file,
  "\n"
)
