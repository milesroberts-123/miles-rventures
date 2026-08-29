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

  # estimate parameters
  modf <- data.frame(
    term = names(stats::coefficients(stats::lm(Y ~ X2 + X1, data = df))),
    value = unname(stats::coefficients(stats::lm(Y ~ X2 + X1, data = df)))
  )
  modf$type <- "joint"

  mod1 <- data.frame(
    term = names(stats::coefficients(stats::lm(Y ~ X1, data = df))),
    value = unname(stats::coefficients(stats::lm(Y ~ X1, data = df)))
  )
  mod1$type <- "single"

  mod2 <- data.frame(
    term = names(stats::coefficients(stats::lm(Y ~ X2, data = df))),
    value = unname(stats::coefficients(stats::lm(Y ~ X2, data = df)))
  )
  mod2$type <- "single"

  all_mod <- rbind(modf, mod1, mod2)

  return(all_mod)
}
