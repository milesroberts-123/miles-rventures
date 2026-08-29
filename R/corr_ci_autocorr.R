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
