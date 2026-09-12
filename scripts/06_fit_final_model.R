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

ensemble_file <- get_arg("--ensemble")
pheno_file    <- get_arg("--phenotype")

out_dir <- get_arg(
  "--out-dir",
  "outputs/final_model/t2d"
)

if (is.null(ensemble_file) || is.null(pheno_file)) {
  stop("Required: --ensemble --phenotype")
}

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Load data
# ============================================================

prs <- fread(ensemble_file)
pheno <- fread(pheno_file)

required_prs <- c(
  "IID",
  "ens_60_25_15"
)

required_pheno <- c(
  "sampleid",
  "T2D",
  "sex",
  "age_sd",
  paste0("PC", 1:10, "_sd")
)

missing_prs <- setdiff(
  required_prs,
  names(prs)
)

missing_pheno <- setdiff(
  required_pheno,
  names(pheno)
)

if (length(missing_prs)) {
  stop(
    "Missing PRS columns: ",
    paste(missing_prs, collapse = ", ")
  )
}

if (length(missing_pheno)) {
  stop(
    "Missing phenotype/covariate columns: ",
    paste(missing_pheno, collapse = ", ")
  )
}

if (anyDuplicated(prs$IID)) {
  stop("Duplicate IID values detected in ensemble table.")
}

if (anyDuplicated(pheno$sampleid)) {
  stop("Duplicate sampleid values detected in phenotype table.")
}

# ============================================================
# Merge final labeled training cohort with ensemble PRS
# ============================================================

train <- merge(
  pheno,
  prs[
    ,
    .(
      IID,
      ens_60_25_15
    )
  ],
  by.x = "sampleid",
  by.y = "IID"
)

cat(
  "Phenotype rows:",
  nrow(pheno),
  "\n"
)

cat(
  "Matched PRS rows:",
  nrow(train),
  "\n"
)

# ============================================================
# Define final model analysis set
#
# Standardization parameters must be estimated from the same
# training subjects used to fit the model.
# ============================================================

model_vars <- c(
  "T2D",
  "ens_60_25_15",
  "age_sd",
  "sex",
  paste0("PC", 1:10, "_sd")
)

train <- train[
  T2D %in% c(0, 1) &
    sex %in% c("Female", "Male")
]

train <- train[
  complete.cases(
    train[
      ,
      ..model_vars
    ]
  )
]

if (!nrow(train)) {
  stop("No complete observations available for final model.")
}

if (
  any(
    !is.finite(
      train$ens_60_25_15
    )
  )
) {
  stop("Non-finite ensemble PRS values detected.")
}

# Explicit reference level reproduces sexMale coefficient.
train[
  ,
  sex := factor(
    sex,
    levels = c(
      "Female",
      "Male"
    )
  )
]

cat(
  "Final model N:",
  nrow(train),
  "\n"
)

cat(
  "Cases:",
  sum(train$T2D == 1),
  "\n"
)

cat(
  "Controls:",
  sum(train$T2D == 0),
  "\n"
)

# ============================================================
# Training-derived PRS standardization
#
# IMPORTANT:
# Mean and SD are estimated only from the final labeled
# training cohort. The same parameters are later applied to
# holdout subjects.
# ============================================================

prs_mean <- mean(
  train$ens_60_25_15
)

prs_sd <- sd(
  train$ens_60_25_15
)

if (
  !is.finite(prs_mean) ||
  !is.finite(prs_sd) ||
  prs_sd <= 0
) {
  stop("Invalid training-derived PRS mean/SD.")
}

train[
  ,
  prs_std :=
    (
      ens_60_25_15 -
        prs_mean
    ) /
      prs_sd
]

cat(
  "PRS training mean:",
  format(prs_mean, digits = 12),
  "\n"
)

cat(
  "PRS training SD:",
  format(prs_sd, digits = 12),
  "\n"
)

# ============================================================
# Final logistic model
#
# Historical final model:
#
# T2D ~ standardized ensemble PRS
#       + age
#       + sex
#       + PCs 1-10
# ============================================================

formula_final <- T2D ~
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
  PC10_sd

model <- glm(
  formula_final,
  data = train,
  family = binomial(
    link = "logit"
  )
)

cat(
  "\n----------------------------------------\n"
)

cat(
  "Final T2D logistic model\n"
)

cat(
  "----------------------------------------\n"
)

print(
  summary(model)
)

# ============================================================
# Save model and preprocessing parameters
#
# The bundle contains everything required to reproduce
# prediction on an independent holdout cohort.
# ============================================================

model_bundle <- list(

  model = model,

  prs_mean = prs_mean,

  prs_sd = prs_sd,

  sex_levels =
    levels(train$sex),

  ensemble =
    "0.60*prs_p1 + 0.25*prs_p0.5 + 0.15*prs_p0.1",

  formula =
    formula_final,

  n_training =
    nrow(train),

  n_cases =
    sum(train$T2D == 1),

  n_controls =
    sum(train$T2D == 0)
)

model_file <- file.path(
  out_dir,
  "t2d_final_model.rds"
)

saveRDS(
  model_bundle,
  model_file
)

# ============================================================
# Human-readable coefficient table
# ============================================================

sm <- summary(model)$coefficients

coef_dt <- data.table(
  term = rownames(sm),
  estimate = sm[, "Estimate"],
  std_error = sm[, "Std. Error"],
  z_value = sm[, "z value"],
  p_value = sm[, "Pr(>|z|)"]
)

coef_file <- file.path(
  out_dir,
  "t2d_final_model_coefficients.tsv"
)

fwrite(
  coef_dt,
  coef_file,
  sep = "\t"
)

# ============================================================
# Final summary
# ============================================================

cat(
  "\nSaved model bundle:",
  model_file,
  "\n"
)

cat(
  "Saved coefficient table:",
  coef_file,
  "\n"
)
