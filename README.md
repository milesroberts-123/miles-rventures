# rventures

An R package of convenience functions for use across RStudio projects.

## Contents

- [Installation](#installation)
- [Usage](#usage)
- [Functions](#functions)
  - [Plotting](#plotting)
  - [Statistics](#statistics)
  - [Population genetics](#population-genetics)
  - [Temporal-replicate toolkit](#temporal-replicate-toolkit)
  - [Dataset structures](#dataset-structures)
  - [Gene expression](#gene-expression)
  - [File management](#file-management)
- [Development](#development)

## Installation

From GitHub (requires the `devtools` package):

```r
devtools::install_github("milesroberts-123/miles-rventures")
```

Or from a local clone:

```r
devtools::install()
```

## Usage

```r
library(rventures)
```

## Functions

### Plotting

- `save_plot()` — save a ggplot as PDF/PNG with by-date/by-analysis symlinks
- `plot_groups_pdf()` — split a data frame by grouping columns and write one plot per group to a multi-page PDF
- `plot_manhattan()` — GWAS Manhattan plot with per-chromosome x-axis and optional SNP highlighting
- `plot_paf_dotplot()` — dotplot of aligned feature (PAF) coordinates, query vs. target
- `plot_var_cov_matrix()` — upper-triangle heatmap of a variance-covariance matrix
- `plot_predictions_truth_scatter()` — predictions vs. truth scatter with identity line, fit, and correlation stats

### Statistics

- `corr_ci_autocorr()` — Pearson correlation CI adjusted for AR(1) autocorrelation
- `lm_sim()`, `multicol_sim()` — linear-model simulations
- `train_abc()` — ABC cross-validation: predict each validation simulation from training simulations via `abc` (requires the suggested `abc` package)

### Population genetics

- `fitfreq()`, `WF_sel()` — selection/drift allele-frequency simulations
- `fc()` — standardized variance in allele-frequency change (Waples 1989)
- `waples_ne()` — temporal Ne estimate corrected for selfing
- `hill_weir_r2()` — expected LD (r^2) between loci a given distance apart (Hill & Weir 1988)

### Temporal-replicate toolkit

- `freq_increments()` — allele-frequency changes between adjacent time points
- `rm_na_after_na()` — NA-truncate a time series at its first NA
- `arcsin_sqrt()` — variance-stabilizing arcsine-square-root transform
- `sum_of_het()`, `geom_pairwise_mean()`, `rolling_matrix_sum()` — small matrix/vector helpers
- `sign_permute_increments()` — sign randomization (window/cell/genome) for null tests
- `covmat_from_pmat()` — covariance matrix of allele-frequency changes, with sample-size correction and heterozygosity standardization
- `standard_cov_by_het()` — standardize a covariance matrix by heterozygosity
- `g_prime()` — linked-selection G statistic corrected for selection-inflated variance
- `gt_from_covmat()` — G ratios accumulated over nested top-left squares of a covariance matrix
- `conv_cor_wn_env()` — convergence correlation among replicates within an environment
- `replicate_gt()` — variance/covariance partition by replicate and time labels
- `covmat_pop_pair()` — mean standardized covariance between a pair of populations
- `feder_t_test()` — Feder-style one-sample t-test of standardized allele-frequency changes, pooled across time points and replicates per population
- `fit_af_glm()` — quasibinomial GLM of allele frequency over time, with optional shared-baseline observation and fixed-intercept offset

### Dataset structures

- `freq_matrix()` — validated L x S allele-frequency matrix
- `snp_coords()` — validated L x 2 chromosome/position table
- `p0_vec()` — validated initial allele-frequency vector
- `sample_info()` — validated S x 4 population/time/replicate/sample-size metadata
- `validate_af_dataset()` — cross-check the four objects for consistent dimensions
- `extract_samples()` — subset a freq_matrix by population/time/replicate with rebuilt column names
- `switch_tracked_allele()` — randomly switch which allele is tracked in a freq_matrix/p0_vec pair, using a shared variant subset
- `filter_fixations()` — drop variants starting near the absorbing boundaries (within tol of 0/1, or NA) and mark later fixation/loss entries NA, ready for `rm_na_after_na()`

### Gene expression

- `expr_matrix()` — validated G x S gene-expression matrix (G genes x S samples)
- `expr_sample_meta()` — validated sample metadata table (S samples x any number of annotation columns)
- `tau()` — cell-type specificity index of an expression vector (0 = ubiquitous, 1 = single sample), with NA-handling modes
- `mean_no_zeros()` — mean expression conditioned on expression above a threshold

### File management

- `read_paf()` — read a PAF alignment file
- `bind_files()` — read a list of files with `fread` and bind their tables row-wise, with a source-file column
- `append_table()` — append a data frame to a file
- `check_file_exists()` — delete a file left over from a previous run, if present
- `ensure_parent_dir()` — create the parent folder structure of a target file path, if missing
- `csv_to_parquet()` — one-time conversion of a large CSV to a chunked parquet dataset (skipped if the directory exists), returned as a lazy `arrow::Dataset`

## Development

After editing `R/` files, regenerate docs and NAMESPACE:

```r
devtools::document()
```