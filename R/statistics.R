#' Confidence interval for a Pearson correlation, adjusted for AR(1) autocorrelation
#'
#' Computes a confidence interval (and optional zero-correlation test) for the
#' Pearson correlation between two serially correlated series. It reduces the
#' nominal sample size N to an effective sample size N_eff based on each series'
#' lag-1 autocorrelation, then applies Fisher's z-transformation using N_eff.
#'
#' @param x,y       Numeric vectors of equal length (paired observations).
#' @param conf      Confidence level (default 0.95).
#' @param bias_correct  Logical; accepted for interface compatibility with
#'   earlier versions of this helper. The lag-1 autocorrelation estimate
#'   used here is the standard (uncorrected) sample ACF in all cases.
#'
#' @return A list with the correlation, effective sample size, CI bounds,
#'         and a t-test of H0: rho = 0 using N_eff.
#' @export
corr_ci_autocorr <- function(x, y, conf = 0.95, bias_correct = TRUE) {

  # Drop pairs with missing values, then require a usable sample size
  ok <- stats::complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 4L) stop("Need at least 4 complete paired observations.")

  r <- stats::cor(x, y)

  # Lag-1 autocorrelation of each series
  rho_x <- stats::acf(x, lag.max = 1, plot = FALSE, demean = TRUE)$acf[2]
  rho_y <- stats::acf(y, lag.max = 1, plot = FALSE, demean = TRUE)$acf[2]

  # Effective sample size for two AR(1) processes; cap at nominal N
  prod_rho <- rho_x * rho_y
  n_eff <- n * (1 - prod_rho) / (1 + prod_rho)
  n_eff <- min(n_eff, n)
  if (n_eff <= 3) {
    warning("Effective sample size <= 3; interval is unreliable.")
    return(list(
      r          = r,
      conf       = conf,
      ci_lower   = NA_real_,
      ci_upper   = NA_real_,
      n          = n,
      n_eff      = n_eff,
      rho_x      = rho_x,
      rho_y      = rho_y,
      t_stat     = NA_real_,
      df         = n_eff - 2,
      p_value    = NA_real_
    ))
  }

  # Back-transform the Fisher-z interval to r-space
  z      <- atanh(r)
  se_z   <- 1 / sqrt(n_eff - 3)
  zcrit  <- stats::qnorm(1 - (1 - conf) / 2)
  ci     <- tanh(z + c(-1, 1) * zcrit * se_z)

  # Significance test of H0: rho = 0, using N_eff
  t_stat <- r * sqrt((n_eff - 2) / (1 - r^2))
  df     <- n_eff - 2
  p_val  <- 2 * stats::pt(-abs(t_stat), df = df)

  list(
    r          = r,
    conf       = conf,
    ci_lower   = ci[1],
    ci_upper   = ci[2],
    n          = n,
    n_eff      = n_eff,
    rho_x      = rho_x,
    rho_y      = rho_y,
    t_stat     = t_stat,
    df         = df,
    p_value    = p_val
  )
}

#' Little simulation, testing if no-multicollinearity == unbiased estimates
#'
#' Simulates `n` observations where `y = a*x1 + b*x2 + e`, with `x1` and `x2`
#' drawn independently (so no multicollinearity), then returns the estimated
#' coefficient of `x1` from a regression of `y` on `x1` alone.
#'
#' @param n Number of observations to simulate.
#' @param a True coefficient of `x1`.
#' @param b True coefficient of `x2`.
#'
#' @return The estimated coefficient of `x1` from `lm(y ~ x1)`.
#' @export
#'
#' @examples
#' hist(replicate(1000, lm_sim(100, -10, 10)))
lm_sim <- function(n, a, b) {

  x1 <- sample(1:100, size = n, replace = TRUE)
  x2 <- sample(1:100, size = n, replace = TRUE)

  e <- stats::rnorm(n)

  y <- a * x1 + b * x2 + e

  stats::coefficients(stats::lm(y ~ x1))[2]
}


#' Simulation to test how multi-collinearity affects estimates
#'
#' Simulates bivariate normal predictors `X1` and `X2` with covariance matrix
#' `covmat`, generates `Y = b0 + b1*X1 + b2*X2 + epsilon`, then fits three
#' regressions: the joint model (`Y ~ X2 + X1`) and each single-predictor
#' model. Returns the coefficient estimates from all three models stacked
#' together.
#'
#' @param n Number of observations to simulate.
#' @param mu1 Mean of `X1`.
#' @param mu2 Mean of `X2`.
#' @param b0 Intercept of the true model.
#' @param b1 True coefficient of `X1`.
#' @param b2 True coefficient of `X2`.
#' @param covmat 2x2 covariance matrix of `X1` and `X2`.
#'
#' @return A data frame with columns `term`, `value`, and `type` (`"joint"`
#'   for the joint model, `"single"` for the single-predictor models).
#' @export
#'
#' @examples
#' multicol_sim(
#'   n = 100, mu1 = 0, mu2 = 0, b0 = 1, b1 = 2, b2 = 3,
#'   covmat = matrix(c(1, 0.9, 0.9, 1), nrow = 2)
#' )
multicol_sim <- function(n, mu1, mu2, b0, b1, b2, covmat) {

  # Generate n random samples from a bivariate normal with covariance covmat
  L <- chol(covmat)
  z <- matrix(stats::rnorm(n * 2), nrow = n, ncol = 2)
  bivariate_data <- sweep(z %*% L, 2, c(mu1, mu2), "+")

  # Convert to a data frame for easier use
  df <- as.data.frame(bivariate_data)
  colnames(df) <- c("X1", "X2")

  # make linear model
  epsilon <- stats::rnorm(n)
  df$Y <- b0 + b1 * df$X1 + b2 * df$X2 + epsilon

  # Stack coefficient estimates from the joint model and each single-predictor
  # model into one data frame
  coef_table <- function(formula, type) {
    coefs <- stats::coefficients(stats::lm(formula, data = df))
    data.frame(term = names(coefs), value = unname(coefs), type = type)
  }

  rbind(
    coef_table(Y ~ X2 + X1, "joint"),
    coef_table(Y ~ X1, "single"),
    coef_table(Y ~ X2, "single")
  )
}


#' Approximate Bayesian computation cross-validation
#'
#' Evaluates an ABC model the way a neural network is validated: after
#' "training" on a set of simulated parameter/summary-statistic pairs, each
#' validation data point is predicted using the training data alone. Each
#' validation simulation is used as the ABC target in turn, and the posterior
#' is summarized by its credible interval, mean, median, and mode.
#'
#' Requires the `abc` package (installed separately; it is only suggested by
#' this package).
#'
#' @param train_sample Data frame of training simulations; columns selected by
#'   `outcomes` and `summary_stats`.
#' @param val_sample Data frame of validation simulations with the same
#'   columns.
#' @param outcomes Character vector naming the outcome (parameter) column(s).
#' @param summary_stats Character vector naming the summary statistic
#'   column(s).
#' @param tol ABC acceptance tolerance, passed to `abc::abc()`.
#' @param method ABC method, passed to `abc::abc()` (`"rejection"`,
#'   `"loclinear"`, `"neuralnet"`, or `"ridge"`).
#' @param verbose If `TRUE`, emit a message for each validation example.
#'
#' @return A data frame with one row per validation simulation and columns
#'   `truth`, `lwr_cred_int`, `pred_mean`, `pred_median`, `pred_mode`, and
#'   `upr_cred_int`.
#' @export
#'
#' @examples
#' \donttest{
#' sim_stats <- data.frame(
#'   Ne = c(100, 200, 500, 1000),
#'   S1 = c(1, 2, 3, 4),
#'   S2 = c(5, 6, 7, 8)
#' )
#' train_abc(
#'   train_sample = sim_stats, val_sample = sim_stats,
#'   outcomes = "Ne", summary_stats = c("S1", "S2"),
#'   tol = 0.5, method = "rejection"
#' )
#' }
train_abc <- function(train_sample, val_sample, outcomes, summary_stats, tol, method, verbose = FALSE) {

  if (!requireNamespace("abc", quietly = TRUE)) {
    stop("The abc package is required for train_abc(). Install it with install.packages('abc').")
  }

  # subset frames to relevant metrics
  train_outcomes <- train_sample[, outcomes]
  train_stats <- train_sample[, summary_stats]
  val_stats <- val_sample[, summary_stats]

  # Preallocate result vectors, one entry per validation simulation
  n_val <- nrow(val_stats)
  lwr_results <- numeric(n_val)
  median_results <- numeric(n_val)
  mean_results <- numeric(n_val)
  mode_results <- numeric(n_val)
  upr_results <- numeric(n_val)

  # loop over validation data
  for (i in seq_len(n_val)) {
    if (verbose) message("Validation example: ", i)

    # get one validation simulation
    val_target <- val_stats[i, ]

    # run abc on validation sim; abc() and summary.abc() write to stdout with
    # cat()/print(), which suppressMessages()/suppressWarnings() cannot touch,
    # so divert stdout with capture.output() and use print = FALSE
    utils::capture.output(
      result <- suppressWarnings(suppressMessages(
        abc::abc(val_target, train_outcomes, train_stats, tol, method)
      )),
      type = "output"
    )
    result_summary <- suppressWarnings(suppressMessages(
      summary(result, print = FALSE)
    ))

    cred_int <- round(c(result_summary[[2]], result_summary[[6]]))

    # save prediction
    lwr_results[i] <- cred_int[1]
    median_results[i] <- round(result_summary[[3]])
    mean_results[i] <- round(result_summary[[4]])
    mode_results[i] <- round(result_summary[[5]])
    upr_results[i] <- cred_int[2]
  }

  return(
    data.frame(
      truth = val_sample[, outcomes],
      lwr_cred_int = lwr_results,
      pred_mean = mean_results,
      pred_median = median_results,
      pred_mode = mode_results,
      upr_cred_int = upr_results
    )
  )
}

# ---------------------------------------------------------------------------
# Block bootstrap
# ---------------------------------------------------------------------------

# Internal: one block-resampled copy of `data`. Rows are never shuffled
# freely; whole blocks are concatenated so serially correlated observations
# (e.g., adjacent genomic windows) move together.
resample_windows <- function(data, block_length = NULL, block_col = NULL,
                             scaffold_col = NULL) {
  n <- nrow(data)

  # Cluster bootstrap: resample whole groups (block ids) with replacement.
  if (!is.null(block_col)) {
    rows_by_id <- split(seq_len(n), data[[block_col]], drop = TRUE)
    picked <- sample(names(rows_by_id), length(rows_by_id), replace = TRUE)
    idx <- unlist(rows_by_id[picked], use.names = FALSE)
    return(data[idx, , drop = FALSE])
  }

  # Circular moving-block bootstrap, one circular series per scaffold.
  if (!is.null(scaffold_col)) {
    keys <- factor(data[[scaffold_col]], levels = unique(data[[scaffold_col]]))
    idx <- unlist(lapply(split(seq_len(n), keys, drop = TRUE), function(rows) {
      ng <- length(rows)
      lg <- min(block_length, ng)
      starts <- sample.int(ng, ceiling(ng / lg), replace = TRUE)
      blocks <- unlist(
        lapply(starts, function(st) ((st:(st + lg - 1) - 1) %% ng) + 1)
      )
      rows[blocks[seq_len(ng)]]
    }), use.names = FALSE)
    return(data[idx, , drop = FALSE])
  }

  # Single circular series: wrap-around blocks over row order.
  starts <- sample.int(n, ceiling(n / block_length), replace = TRUE)
  idx <- unlist(
    lapply(starts, function(st) ((st:(st + block_length - 1) - 1) %% n) + 1)
  )
  data[idx[seq_len(n)], , drop = FALSE]
}

#' Block bootstrap confidence intervals for linear model coefficients
#'
#' Estimates confidence intervals for the coefficients of a linear model via
#' block bootstrapping, designed for data whose rows are serially correlated —
#' e.g., genomic windows, where adjacent windows share recombination history.
#' Resampling whole blocks keeps the correlation structure intact, giving more
#' accurate interval estimates than ordinary case (row-wise) resampling.
#'
#' Two blocking schemes are supported (exactly one of `block_length` or
#' `block_col` must be supplied):
#'
#' - **Positional circular moving blocks** (`block_length`): each replicate
#'   draws `ceiling(n / L)` block starts uniformly and wraps around the end of
#'   the series (row *n* rejoins row 1), truncating to *n* rows so each
#'   replicate matches the original sample size. If `scaffold_col` is given,
#'   each scaffold is treated as its own circular series: block wrap-around
#'   never crosses scaffold boundaries and every replicate keeps each
#'   scaffold's original window count. Scaffolds shorter than `block_length`
#'   are resampled as a single whole-scaffold block. Without `scaffold_col`
#'   the whole data set is assumed to be one ordered sequence — the caller is
#'   responsible for row order (see `order_by`).
#' - **Cluster blocks** (`block_col`): the unique values of this column
#'   (e.g., chromosome chunk ids) are resampled with replacement and all rows
#'   of a sampled id stay together.
#'
#' Coefficients are estimated on each resampled data set with
#' \code{\link[stats]{lm}}; intervals are computed per coefficient with
#' \code{\link[boot]{boot.ci}} (default `type = "basic"`, matching the
#' reverse-percentile interval in the `boot` package's vignette workflows).
#' The raw \code{\link[boot]{boot}} object is attached as attribute
#' `"boot_out"` for further interval types or diagnostics.
#'
#' @param data A data frame of observations (e.g., genome windows).
#' @param formula A linear model formula passed to \code{\link[stats]{lm}}.
#' @param R Number of bootstrap replicates (default 1000).
#' @param block_length Block length in rows (windows) for the positional
#'   circular moving-block bootstrap. Must be a positive integer smaller than
#'   `nrow(data)` when `scaffold_col` is not used.
#' @param block_col Name of the column whose unique values define the
#'   resampling units in cluster mode.
#' @param order_by Optional name of a column used to sort rows (within each
#'   scaffold, when `scaffold_col` is given) before positional blocking, e.g.,
#'   a window midpoint. Ignored in cluster mode.
#' @param scaffold_col Optional name of the column identifying the
#'   scaffold/chromosome each row belongs to; required with `block_length`
#'   when rows span multiple sequences. Ignored in cluster mode.
#' @param conf Confidence level (default 0.95).
#' @param type Interval type passed to \code{\link[boot]{boot.ci}}:
#'   `"basic"` (default), `"perc"`, `"norm"`, or `"bca"`.
#' @param seed Optional integer seed for reproducibility. When `NULL`, the
#'   caller's random state (e.g., a preceding \code{\link[base]{set.seed}})
#'   is used.
#'
#' @return A data frame with one row per model coefficient: `Term`,
#'   `Estimate`, `Lower_CI`, `Upper_CI`. The raw
#'   \code{\link[boot]{boot}} object is attached as attribute `"boot_out"`.
#'
#' @examples
#' set.seed(42)
#' windows <- data.frame(
#'   pooled_theta_pi = rnorm(120),
#'   cds_density     = rnorm(120),
#'   mean_r2_direct  = rnorm(120),
#'   window_midpoint = rep(seq_len(40), times = 3),
#'   scaffold        = rep(paste0("chr", 1:3), each = 40)
#' )
#' \donttest{
#' block_bootstrap_lm(
#'   windows,
#'   pooled_theta_pi ~ cds_density + mean_r2_direct,
#'   R = 500, block_length = 4,
#'   scaffold_col = "scaffold", order_by = "window_midpoint",
#'   seed = 123
#' )
#' }
#' @export
block_bootstrap_lm <- function(data, formula, R = 1000, block_length = NULL,
                               block_col = NULL, order_by = NULL,
                               scaffold_col = NULL, conf = 0.95,
                               type = "basic", seed = NULL) {

  if (!is.data.frame(data) || nrow(data) < 1L) {
    stop("data must be a non-empty data frame.")
  }
  n <- nrow(data)

  if (is.null(block_length) == is.null(block_col)) {
    stop("Exactly one of block_length or block_col must be supplied.")
  }
  if (!is.null(block_length)) {
    if (!is.numeric(block_length) || length(block_length) != 1L ||
        is.na(block_length) || block_length < 1L ||
        block_length != floor(block_length)) {
      stop("block_length must be a positive integer.")
    }
    if (is.null(scaffold_col) && block_length >= n) {
      stop("block_length must be smaller than the number of rows.")
    }
  }
  if (!is.character(type) || length(type) != 1L || is.na(type)) {
    stop("type must be a single character string.")
  }
  if (!is.numeric(R) || length(R) != 1L || is.na(R) || R < 1L ||
      R != floor(R)) {
    stop("R must be a positive integer.")
  }

  if (is.null(block_col)) {
    # Positional mode: put rows into blocking order up front; the resampler
    # relies on that order.
    if (!is.null(scaffold_col) && !scaffold_col %in% names(data)) {
      stop("scaffold_col is not a column of data.")
    }
    if (!is.null(order_by) && !order_by %in% names(data)) {
      stop("order_by is not a column of data.")
    }
    sort_keys <- c(scaffold_col, order_by)
    if (length(sort_keys) > 0L) {
      data <- data[do.call(order, data[sort_keys]), , drop = FALSE]
    }
  } else {
    if (!block_col %in% names(data)) {
      stop("block_col is not a column of data.")
    }
  }

  ran_gen <- function(original, mle) {
    resample_windows(
      original,
      block_length = block_length,
      block_col = block_col,
      scaffold_col = scaffold_col
    )
  }
  statistic <- function(dat, ...) stats::coef(stats::lm(formula, data = dat))

  if (!is.null(seed)) set.seed(seed)
  boot_out <- boot::boot(
    data,
    statistic = statistic,
    R = R,
    sim = "parametric",
    ran.gen = ran_gen,
    mle = NULL
  )

  coefs <- boot_out$t0
  bounds <- t(vapply(seq_along(coefs), function(i) {
    ci <- boot::boot.ci(boot_out, conf = conf, type = type, index = i)
    if (identical(type, "norm")) {
      ci$normal[c(2, 3)]
    } else {
      ci[[type]][c(4, 5)]
    }
  }, numeric(2)))

  out <- data.frame(
    Term = names(coefs),
    Estimate = as.numeric(coefs),
    Lower_CI = as.numeric(bounds[, 1]),
    Upper_CI = as.numeric(bounds[, 2])
  )
  attr(out, "boot_out") <- boot_out
  out
}
