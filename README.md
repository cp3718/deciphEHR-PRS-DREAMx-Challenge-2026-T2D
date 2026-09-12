# deciphEHR PRS DREAMx Challenge 2026 – T2D Workflow

This repository contains a reconstruction of the Type 2 Diabetes (T2D) polygenic risk score workflow used for our challenge submission.

The repository currently includes the submitted analysis workflow only.

## Workflow

The T2D workflow consisted of:

1. **GWAS summary-statistic preparation**
   - Construct challenge-compatible variant IDs.
   - Convert GWAS odds ratios to log-odds effect sizes.
   - Remove invalid variants and duplicate variant IDs.
   - Generate score files at seven GWAS P-value thresholds:

   `1e-5, 1e-3, 0.01, 0.05, 0.1, 0.5, 1`

2. **PRS calculation**
   - Calculate chromosome-specific PRSs using PLINK 2.
   - No LD clumping or pruning was applied.
   - Combine chromosome-level scores into genome-wide PRSs for each threshold.

3. **Threshold diagnostics**
   - Compare threshold-specific PRSs in the labeled training cohort using logistic regression with age, sex, and the first 10 ancestry principal components as covariates.

4. **PRS ensemble**
   - Several weighted combinations of the broader-threshold PRSs were explored empirically.
   - The final ensemble selected was:

   \[
   PRS_{ensemble} =
   0.60\,PRS_{P\leq1}
   + 0.25\,PRS_{P\leq0.5}
   + 0.15\,PRS_{P\leq0.1}
   \]

5. **Final model**
   - Standardize the ensemble PRS using the labeled training cohort.
   - Fit a logistic regression model including:
     - ensemble PRS
     - age
     - sex
     - PCs 1–10

6. **Holdout prediction**
   - Apply the training-derived PRS mean and standard deviation to the holdout cohort.
   - Generate final T2D probabilities using the fitted logistic model.

## Scripts

```text
scripts/
├── 01_prepare_sumstats.R
├── 02_generate_prs.sbatch
├── 03_collect_prs.R
├── 04_threshold_diagnostics.R
├── 05_build_ensemble.R
├── 06_fit_final_model.R
└── 07_predict_holdout.R
