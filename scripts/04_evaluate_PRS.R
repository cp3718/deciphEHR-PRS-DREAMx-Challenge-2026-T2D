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
pheno_file <- get_arg("--phenotype")

out_dir <- get_arg(
  "--out-dir",
  "outputs/evaluation/t2d"
)

if (is.null(prs_file) || is.null(pheno_file)) {
  stop("Required: --prs --phenotype")
}

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Load PRSs and phenotype/covariates
# ============================================================

prs <- fread(prs_file)
pheno <- fread(pheno_file)

if (!"IID" %in% names(prs)) {
  stop("PRS table must contain IID.")
}

if (!"sampleid" %in% names(pheno)) {
  stop("Phenotype table must contain sampleid.")
}

required_pheno <- c(
  "sampleid",
  "T2D",
  "age_sd",
  "sex",
  paste0("PC", 1:10, "_sd")
)

missing_pheno <- setdiff(
  required_pheno,
  names(pheno)
)

if (length(missing_pheno)) {
  stop(
    "Missing phenotype/covariate columns: ",
    paste(missing_pheno, collapse = ", ")
  )
}

dt <- merge(
  pheno,
  prs,
  by.x = "sampleid",
  by.y = "IID"
)

# ============================================================
# Eligible labeled analysis set
#
# This step is an in-sample threshold diagnostic, not an
# out-of-sample validation analysis.
# ============================================================

dt <- dt[
  T2D %in% c(0, 1) &
    sex %in% c("Male", "Female")
]

dt[
  ,
  sex := factor(sex)
]

prs_cols <- grep(
  "^prs_p",
  names(dt),
  value = TRUE
)

if (!length(prs_cols)) {
  stop("No PRS columns beginning with 'prs_p' were found.")
}

# Keep the seven historical thresholds in explicit order.
expected_prs_cols <- c(
  "prs_p1e5",
  "prs_p1e3",
  "prs_p0.01",
  "prs_p0.05",
  "prs_p0.1",
  "prs_p0.5",
  "prs_p1"
)

missing_prs <- setdiff(
  expected_prs_cols,
  prs_cols
)

if (length(missing_prs)) {
  stop(
    "Missing expected PRS columns: ",
    paste(missing_prs, collapse = ", ")
  )
}

prs_cols <- expected_prs_cols

# ============================================================
# Covariates
# ============================================================

covariates <- c(
  "age_sd",
  "sex",
  paste0("PC", 1:10, "_sd")
)

# ============================================================
# Threshold diagnostics
# ============================================================

results <- vector(
  "list",
  length(prs_cols)
)

names(results) <- prs_cols

for (col in prs_cols) {

  model_vars <- c(
    "T2D",
    col,
    covariates
  )

  x <- copy(
    dt[
      complete.cases(dt[, ..model_vars])
    ]
  )

  if (!nrow(x)) {
    stop(
      "No complete observations available for ",
      col
    )
  }

  # Standardize the PRS within the analysis cohort so the
  # coefficient represents the effect per 1 SD increase.
  x[
    ,
    prs_std :=
      as.numeric(
        scale(get(col))
      )
  ]

  model <- glm(
    T2D ~
      prs_std +
      age_sd +
      sex +
      PC1_sd +
      PC2_sd +
      PC3_sd +
      PC4_sd +
      PC5_sd +
      PC6_sd +
      PC7_sd +
      PC8_sd +
      PC9_sd +
      PC10_sd,
    data = x,
    family = binomial()
  )

  sm <- summary(model)$coefficients

  results[[col]] <- data.table(
    prs = col,

    n = nobs(model),

    n_cases =
      sum(x$T2D == 1),

    n_controls =
      sum(x$T2D == 0),

    case_mean_prs_z =
      mean(
        x[T2D == 1, prs_std]
      ),

    control_mean_prs_z =
      mean(
        x[T2D == 0, prs_std]
      ),

    prs_beta =
      sm["prs_std", "Estimate"],

    prs_se =
      sm["prs_std", "Std. Error"],

    prs_or_per_sd =
      exp(
        sm["prs_std", "Estimate"]
      ),

    prs_p =
      sm["prs_std", "Pr(>|z|)"],

    AIC =
      AIC(model),

    null_deviance =
      model$null.deviance,

    residual_deviance =
      model$deviance
  )
}

results <- rbindlist(
  results
)

# ============================================================
# Save
# ============================================================

output_file <- file.path(
  out_dir,
  "t2d_threshold_diagnostics.tsv"
)

fwrite(
  results,
  output_file,
  sep = "\t"
)

# ============================================================
# Summary
# ============================================================

cat(
  "\n----------------------------------------\n"
)

cat(
  "T2D threshold diagnostics\n"
)

cat(
  "----------------------------------------\n"
)

print(
  results[
    order(AIC)
  ]
)

cat(
  "\nOutput:",
  output_file,
  "\n"
)
