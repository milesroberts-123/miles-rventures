# rventures

An R package of convenience functions for use across RStudio projects.

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

Functions include:

- `save_plot()` — save a ggplot as PDF/PNG with by-date/by-analysis symlinks
- `plot_groups_pdf()` — split a data frame by grouping columns and write one plot per group to a multi-page PDF
- `corr_ci_autocorr()` — Pearson correlation CI adjusted for AR(1) autocorrelation
- `fitfreq()`, `WF.sel()` — population-genetics selection/drift simulations
- `lm_sim()`, `multicol_sim()` — linear-model simulations
- `append_table()`, `check_file_exists()` — small file utilities

## Development

After editing `R/` files, regenerate docs and NAMESPACE:

```r
devtools::document()
```
