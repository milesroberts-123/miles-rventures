#' Confidence interval for a Pearson correlation, adjusted for AR(1) autocorrelation
#'
#' Computes a confidence interval (and optional zero-correlation test) for the
#' Pearson correlation between two serially correlated series. It reduces the
#' nominal sample size N to an effective sample size N_eff based on each series'
#' lag-1 autocorrelation, then applies Fisher's z-transformation using N_eff.
#'
#' @param x,y       Numeric vectors of equal length (paired observations).
#' @param conf      Confidence level (default 0.95).
#'
#' @return A list with the correlation, effective sample size, CI bounds,
#'         and a t-test of H0: rho = 0 using N_eff.
#' @export
corr_ci_autocorr <- function(x, y, conf = 0.95) {

  # --- input handling: drop pairs with NA, check length ---------------------
  ok <- stats::complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 4L) stop("Need at least 4 complete paired observations.")

  # --- Pearson correlation --------------------------------------------------
  r <- stats::cor(x, y)

  # --- lag-1 autocorrelation of each series ---------------------------------
  rho_x <- stats::acf(x, lag.max = 1, plot = FALSE, demean = TRUE)$acf[2]
  rho_y <- stats::acf(y, lag.max = 1, plot = FALSE, demean = TRUE)$acf[2]

  # --- effective sample size for two AR(1) processes ------------------------
  #   N_eff = N * (1 - rho_x * rho_y) / (1 + rho_x * rho_y)
  prod_rho <- rho_x * rho_y
  n_eff <- n * (1 - prod_rho) / (1 + prod_rho)
  # Guard against pathological values; N_eff must exceed 3 for the SE below.
  n_eff <- min(n_eff, n)            # cap at nominal N (negative prod_rho case)
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

  # --- Fisher z-transform and CI -------------------------------------------
  z      <- atanh(r)
  se_z   <- 1 / sqrt(n_eff - 3)
  zcrit  <- stats::qnorm(1 - (1 - conf) / 2)
  ci     <- tanh(z + c(-1, 1) * zcrit * se_z)   # back-transform to r-space

  # --- significance test of H0: rho = 0, using N_eff ------------------------
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

  # vectors to store results
  lwr_results <- c()
  median_results <- c()
  mean_results <- c()
  mode_results <- c()
  upr_results <- c()

  # subset frames to relevant metrics
  train_outcomes <- train_sample[, outcomes]
  train_stats <- train_sample[, summary_stats]
  val_stats <- val_sample[, summary_stats]

  # loop over validation data
  for (i in 1:nrow(val_stats)) {
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
    median_est <- round(result_summary[[3]])
    mean_est <- round(result_summary[[4]])
    mode_est <- round(result_summary[[5]])

    # save prediction
    lwr_results <- c(lwr_results, cred_int[1])
    median_results <- c(median_results, median_est)
    mean_results <- c(mean_results, mean_est)
    mode_results <- c(mode_results, mode_est)
    upr_results <- c(upr_results, cred_int[2])
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
