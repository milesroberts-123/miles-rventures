test_that("corr_ci_autocorr returns a list with expected names", {
  set.seed(1)
  x <- rnorm(100)
  y <- rnorm(100)
  out <- corr_ci_autocorr(x, y)
  expect_named(
    out,
    c("r", "conf", "ci_lower", "ci_upper", "n", "n_eff",
      "rho_x", "rho_y", "t_stat", "df", "p_value")
  )
})

test_that("corr_ci_autocorr CI contains 0 for uncorrelated data", {
  set.seed(1)
  x <- rnorm(200)
  y <- rnorm(200)
  out <- corr_ci_autocorr(x, y)
  expect_true(out$ci_lower <= 0 && out$ci_upper >= 0)
})

test_that("corr_ci_autocorr drops NA pairs", {
  x <- c(1, 2, 3, 4, 5, NA)
  y <- c(2, 4, 6, 8, 10, 1)
  out <- corr_ci_autocorr(x, y)
  expect_equal(out$n, 5)
})

test_that("corr_ci_autocorr errors with fewer than 4 observations", {
  expect_error(corr_ci_autocorr(1:3, 1:3), "at least 4")
})

test_that("corr_ci_autocorr accepts bias_correct argument for interface parity", {
  set.seed(1)
  x <- rnorm(100)
  y <- rnorm(100)
  out_default <- corr_ci_autocorr(x, y)
  out_flag <- corr_ci_autocorr(x, y, bias_correct = TRUE)
  out_off <- corr_ci_autocorr(x, y, bias_correct = FALSE)
  expect_equal(out_default$r, out_flag$r)
  expect_equal(out_flag$r, out_off$r)
})

test_that("lm_sim returns a scalar", {
  set.seed(1)
  out <- lm_sim(100, -10, 10)
  expect_length(out, 1)
  expect_true(is.numeric(out))
})

test_that("multicol_sim returns 7 rows with term/value/type", {
  set.seed(1)
  out <- multicol_sim(
    n = 100, mu1 = 0, mu2 = 0, b0 = 1, b1 = 2, b2 = 3,
    covmat = matrix(c(1, 0.9, 0.9, 1), nrow = 2)
  )
  expect_equal(nrow(out), 7)
  expect_named(out, c("term", "value", "type"))
  expect_true(all(out$type %in% c("joint", "single")))
})

test_that("train_abc predicts validation outcomes and is silent", {
  skip_if_not_installed("abc")

  set.seed(42)
  n_sim <- 200
  train_sample <- data.frame(
    Ne = sample(c(100, 200, 500, 1000, 2000), n_sim, replace = TRUE),
    S1 = rnorm(n_sim),
    S2 = rnorm(n_sim)
  )
  # summary statistics correlated with the outcome so predictions are informative
  train_sample$S1 <- log10(train_sample$Ne) + rnorm(n_sim, sd = 0.2)
  train_sample$S2 <- rnorm(n_sim, sd = 0.5)

  val_sample <- data.frame(Ne = c(100, 500, 2000), S1 = NA, S2 = NA)
  val_sample$S1 <- log10(val_sample$Ne) + rnorm(nrow(val_sample), sd = 0.2)
  val_sample$S2 <- rnorm(nrow(val_sample), sd = 0.5)

  expect_silent({
    out <- train_abc(
      train_sample, val_sample,
      outcomes = "Ne", summary_stats = c("S1", "S2"),
      tol = 0.1, method = "rejection"
    )
  })

  expect_s3_class(out, "data.frame")
  expect_named(out, c(
    "truth", "lwr_cred_int", "pred_mean",
    "pred_median", "pred_mode", "upr_cred_int"
  ))
  expect_equal(nrow(out), nrow(val_sample))
  expect_equal(out$truth, val_sample$Ne)

  # posterior summaries must fall inside the credible interval
  expect_true(all(out$lwr_cred_int <= out$pred_mean))
  expect_true(all(out$lwr_cred_int <= out$pred_median))
  expect_true(all(out$lwr_cred_int <= out$pred_mode))
  expect_true(all(out$pred_mean <= out$upr_cred_int))
  expect_true(all(out$pred_median <= out$upr_cred_int))
  expect_true(all(out$pred_mode <= out$upr_cred_int))
})

test_that("train_abc errors when abc is not installed", {
  skip_if(requireNamespace("abc", quietly = TRUE))

  expect_error(
    train_abc(
      data.frame(Ne = 1, S1 = 1), data.frame(Ne = 1, S1 = 1),
      outcomes = "Ne", summary_stats = "S1", tol = 0.5, method = "rejection"
    ),
    "abc package is required"
  )
})

# ---------------------------------------------------------------------------
# Ported: rm_na_after_na, conv_cor_wn_env, replicate_gt extras
# (from grenenet-phase2 resources/R/tests.R)
# ---------------------------------------------------------------------------

test_that("conv_cor_wn_env returns a length-2 vector: numerator and denominator", {
  set.seed(42)
  pdiff <- matrix(rnorm(30), nrow = 5, ncol = 6)
  result <- conv_cor_wn_env(pdiff)
  expect_equal(length(result), 2)
})

test_that("conv_cor_wn_env returns numerator smaller than denominator", {
  set.seed(99)
  pdiff <- matrix(rnorm(60), nrow = 10, ncol = 6)
  result <- conv_cor_wn_env(pdiff)
  expect_true(abs(result[1]) <= result[2])
})

test_that("conv_cor_wn_env works with matrices of varying sizes", {
  set.seed(7)
  pdiff_4x3 <- matrix(rnorm(12), nrow = 4, ncol = 3)
  result <- conv_cor_wn_env(pdiff_4x3)
  expect_equal(length(result), 2)
})

test_that("conv_cor_wn_env returns positive numerator for positively correlated data", {
  set.seed(123)
  a <- rnorm(10, mean = 0, sd = 1)
  b <- a + rnorm(10, mean = 0, sd = 0.1)
  pdiff <- cbind(a, b)
  result <- conv_cor_wn_env(pdiff)
  expect_true(result[1] > 0)
})

test_that("conv_cor_wn_env returns negative numerator for negatively correlated data", {
  set.seed(123)
  a <- rnorm(10, mean = 0, sd = 1)
  b <- -1 * a + rnorm(10, mean = 0, sd = 0.1)
  pdiff <- cbind(a, b)
  result <- conv_cor_wn_env(pdiff)
  expect_true(result[1] < 0)
})

test_that("conv_cor_wn_env requires matrix input with ncol > 1 and nrow > 1", {
  expect_error(conv_cor_wn_env(c(1, 2, 3)))
})
