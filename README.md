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

### Dataset structures

- `freq_matrix()` — validated L x S allele-frequency matrix
- `snp_coords()` — validated L x 2 chromosome/position table
- `p0_vec()` — validated initial allele-frequency vector
- `sample_info()` — validated S x 3 population/time/replicate metadata
- `validate_af_dataset()` — cross-check the four objects for consistent dimensions
- `extract_samples()` — subset a freq_matrix by population/time/replicate with rebuilt column names

### File management

- `read_paf()` — read a PAF alignment file
- `append_table()` — append a data frame to a file
- `check_file_exists()` — stop with an informative error if a file is missing

## Development

After editing `R/` files, regenerate docs and NAMESPACE:

```r
devtools::document()
```