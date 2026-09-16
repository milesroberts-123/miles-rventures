test_that("expr_matrix validates and classifies a GxS matrix", {
  m <- matrix(c(10, 0, 5, 8, 2, 7), nrow = 3, ncol = 2)
  em <- expr_matrix(m)
  expect_s3_class(em, "expr_matrix")
  expect_identical(dim(em), c(3L, 2L))
  expect_true(is.numeric(unclass(em)))
  # NA and negative values allowed
  em2 <- expr_matrix(matrix(c(1, NA, -2, 4), nrow = 2))
  expect_s3_class(em2, "expr_matrix")
  # errors
  expect_error(expr_matrix(data.frame(a = 1:3)), "x must be a numeric matrix.")
  expect_error(expr_matrix(c(1, 2, 3)), "x must be a numeric matrix.")
  expect_error(expr_matrix(matrix("a", 2, 2)), "x must be a numeric matrix.")
  expect_error(expr_matrix(matrix(numeric(0), nrow = 0, ncol = 2)),
               "at least one row and one column.")
  expect_error(expr_matrix(matrix(numeric(0), nrow = 2, ncol = 0)),
               "at least one row and one column.")
})

test_that("expr_matrix prints a summary header", {
  m <- matrix(c(10, 0, 5, 8, 2, 7), nrow = 3)
  colnames(m) <- c("s1", "s2")
  out <- capture.output(print(expr_matrix(m)))
  expect_match(out[1], "expr_matrix: 3 genes x 2 samples", fixed = TRUE)
  expect_match(out[2], "Samples: s1, s2", fixed = TRUE)
})

test_that("expr_sample_meta validates and classifies an Sx? table", {
  d <- data.frame(
    cell_type = c("neuron", "glia", "neuron"),
    condition = c("ctrl", "ctrl", "treat"),
    batch = c(1, 1, 2)
  )
  em <- expr_sample_meta(d)
  expect_s3_class(em, "expr_sample_meta")
  expect_identical(dim(em), c(3L, 3L))
  # NA values allowed
  d2 <- data.frame(cell_type = c("neuron", NA))
  expect_s3_class(expr_sample_meta(d2), "expr_sample_meta")
  # one-column table is fine
  expect_s3_class(expr_sample_meta(data.frame(tissue = "liver")),
                  "expr_sample_meta")
  # errors
  expect_error(expr_sample_meta(matrix(1:4, nrow = 2)),
               "x must be a data.frame.")
  expect_error(expr_sample_meta(data.frame()),
               "x must have at least one row.")
})

test_that("expr_sample_meta prints a summary header", {
  d <- data.frame(cell_type = c("neuron", "glia"))
  out <- capture.output(print(expr_sample_meta(d)))
  expect_match(out[1], "expr_sample_meta: 2 samples", fixed = TRUE)
})

test_that("tau scores ubiquitous, specific, and intermediate genes", {
  # ubiquitous: 0
  expect_equal(tau(c(1, 1, 1, 1)), 0)
  expect_equal(tau(c(100, 100, 100)), 0)
  # single cell type: 1
  expect_equal(tau(c(10, 0, 0, 0)), 1)
  # hand-computed intermediate: x = (5, 2, 8, 1), max = 8
  # scaled = (0.625, 0.25, 1, 0.125); tau = (0.375 + 0.75 + 0 + 0.875)/3
  expect_equal(tau(c(5, 2, 8, 1)), (0.375 + 0.75 + 0 + 0.875) / 3)
  # scale-invariance
  expect_equal(tau(c(5, 2, 8, 1)), tau(c(5000, 2000, 8000, 1000)))
  expect_true(tau(c(5, 2, 8, 1)) > 0 && tau(c(5, 2, 8, 1)) < 1)
})

test_that("tau returns NA for degenerate input", {
  # fewer than 2 values
  expect_identical(tau(1), NA_real_)
  expect_identical(tau(numeric(0)), NA_real_)
  # all zeros
  expect_identical(tau(c(0, 0, 0)), NA_real_)
})

test_that("tau handles NA values by na_action", {
  # default propagate: silent NA
  expect_identical(tau(c(1, NA, 3)), NA_real_)
  expect_no_warning(tau(c(1, NA, 3)))
  # remove: compute on observed values only
  # observed = (1, 3); scaled = (1/3, 1); tau = (2/3 + 0)/1
  expect_equal(tau(c(1, NA, 3), na_action = "remove"), 2 / 3)
  # remove with fewer than 2 observed values: NA
  expect_identical(tau(c(NA, 1), na_action = "remove"), NA_real_)
  expect_identical(tau(c(NA, NA), na_action = "remove"), NA_real_)
  # error
  expect_error(tau(c(1, NA, 3), na_action = "error"),
               "x must not contain NA values.")
  # no NA: all modes agree
  expect_equal(tau(c(5, 2, 8, 1), na_action = "error"),
               tau(c(5, 2, 8, 1), na_action = "remove"))
  expect_equal(tau(c(5, 2, 8, 1), na_action = "error"),
               tau(c(5, 2, 8, 1), na_action = "propagate"))
  # invalid na_action
  expect_error(tau(c(1, 2), na_action = "skip"))
  expect_error(tau(c(1, 2), na_action = c("propagate", "remove")))
})

test_that("tau validates input type", {
  expect_error(tau(matrix(1:4, nrow = 2)), "x must be a numeric vector.")
  expect_error(tau(c("a", "b")), "x must be a numeric vector.")
  expect_error(tau(list(1, 2)), "x must be a numeric vector.")
})

test_that("mean_no_zeros averages values strictly above the threshold", {
  expect_equal(mean_no_zeros(c(0, 2, 4, 8), thresh = 0), (2 + 4 + 8) / 3)
  # strictly greater: value == thresh is excluded
  expect_equal(mean_no_zeros(c(2, 2, 10), thresh = 2), 10)
  # no value exceeds thresh: NA_real_, not NaN, without warning
  expect_identical(mean_no_zeros(c(0, 0, 0), thresh = 0), NA_real_)
  expect_identical(mean_no_zeros(c(1, 2), thresh = 10), NA_real_)
  expect_no_warning(mean_no_zeros(c(0, 0), thresh = 0))
  # NA handling: NAs in x are dropped by na.rm inside mean
  expect_equal(mean_no_zeros(c(NA, 2, 4), thresh = 0), 3)
  expect_identical(mean_no_zeros(c(NA, NA), thresh = 0), NA_real_)
})

test_that("expression module integrates: per-gene tau over a GxS matrix", {
  set.seed(21)
  G <- 5
  S <- 6
  m <- matrix(runif(G * S, 0, 10), nrow = G, ncol = S)
  m[1, ] <- 7 # ubiquitous gene
  m[2, 1] <- 10
  m[2, -1] <- 0 # single-sample gene
  em <- expr_matrix(m)
  meta <- expr_sample_meta(data.frame(
    cell_type = rep(c("A", "B"), times = 3), batch = rep(1:2, each = 3)
  ))
  # sample metadata rows match matrix columns
  expect_identical(nrow(meta), ncol(em))
  # per-gene tau
  tau_vec <- apply(unclass(em), 1, tau)
  expect_equal(length(tau_vec), G)
  expect_equal(tau_vec[1], 0)
  expect_equal(tau_vec[2], 1)
  expect_true(all(tau_vec[3:5] > 0 & tau_vec[3:5] < 1))
  # per-gene mean_no_zeros
  mnz_vec <- apply(unclass(em), 1, mean_no_zeros, thresh = 0)
  expect_equal(length(mnz_vec), G)
  expect_equal(unname(mnz_vec[1]), 7)
  expect_equal(unname(mnz_vec[2]), 10)
})