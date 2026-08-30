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
