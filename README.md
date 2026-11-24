# ABCD suicidal thoughts and behaviors: analysis scripts

This repository organizes the analyses reported in the manuscript into four R scripts. The scripts use prepared ABCD research data, which are not distributed in this repository. Access to participant-level data is controlled by ABCD.

## Setup

Open R with this repository as the working directory. Set `ABCD_PREPARED_DIR` to the directory containing `abcd_SuicideData.RData` and `abcd_ImagingData.RData`, and `ABCD_DATA_DIR` to the directory containing `abcd_covariables_full.csv`. Use local paths appropriate for your computer:

```r
Sys.setenv(ABCD_PREPARED_DIR = '/path/to/prepared_data', ABCD_DATA_DIR = '/path/to/ABCD')
source('01_data_and_clustering.R')
```

The scripts require R packages including `factoextra`, `chisq.posthoc.test`, `reshape2`, `stringr`, `lme4`, `lmerTest`, and `emmeans`; the prepared object may require `ggpubr` to restore all attributes. Install missing packages separately. Clustering recomputation is optional and runs only when `RUN_CLUSTER=1`.

## Scripts

- `01_data_and_clustering.R`: loads the prepared workspace and displays the original five profile counts. Contains the original optional K-means call.
- `02_behavior_environment_and_prs.R`: original behavioral/environmental comparisons and six psychiatric PRS ANOVAs; run after script 01 in the same R session.
- `03_neuroimaging.R`: structural MRI, resting-state connectivity and diffusion models, extracted from the original scripts; run after script 01.
- `04_longitudinal.R`: reporter-specific mixed models and direct change contrasts for the revision; run after script 01. This is a sensitivity analysis, not a validated replacement for the originally reported analysis.

Review code and output before reporting. The scripts retain the original variable indexing and expect the original prepared workspace layout. Do not upload raw ABCD data, derived individual records, local results or genetic files.

## Reconciliation notes

- The present Supplementary Tables S10–S11 ADHD-PRS results reproduce from the saved `prs_prsice` object, not the saved `prs_prscs` object; the manuscript describes PRS-CS. Script 02 follows the numerical table source; method provenance must be resolved before submission. The legacy PRS ANOVA does not include genetic principal components or batch.
- Supplementary Table S7 has 68 cortical thickness, 68 surface area, 14 subcortical variables and whole-brain volume. The original structural block selected 12 subcortical variables; two nucleus accumbens measures were added to match the table. A single cortical-area model was checked against S7, but the complete structural, functional and diffusion analyses have not been independently reproduced.
- The original longitudinal score sums 49 items per reporter, including NSSI-labelled items. Verify item coding, exclusions, missingness and validity before using script 04's output in the manuscript. An all-missing wave is treated as missing in script 04.
- The multiple-testing family, environmental covariate models, family nesting, exact gap-statistic selection rule and tables S13–S14 require further audit. New calculations should not be represented as the original analysis until compared with the archived results.
- These scripts are organized research code, not a fully validated end-to-end reproduction pipeline. Some original intermediate variables, R packages and data formats must be checked in the source workspace.
