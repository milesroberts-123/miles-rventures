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

test_that("create_prob rescales a vector to fractions of its total", {
  expect_equal(unname(create_prob(c(4, 1))), c(0.8, 0.2))
  expect_equal(sum(create_prob(c(3, 3, 2, 2))), 1)
  # single value: the vector of length 1
  expect_equal(create_prob(5), 1)
  # scale invariance
  expect_equal(unname(create_prob(c(4, 1) * 1e6)), c(0.8, 0.2))
  # logical input works like tau's
  expect_equal(unname(create_prob(c(TRUE, TRUE, FALSE, FALSE))), c(0.5, 0.5, 0, 0))
})

test_that("create_prob returns NA_real_ on degenerate totals, silently", {
  # all zeros: no NaN
  expect_identical(create_prob(c(0, 0, 0)), NA_real_)
  expect_no_warning(create_prob(c(0, 0, 0)))
  # negative values canceling the positive ones
  expect_identical(create_prob(c(1, -1)), NA_real_)
  expect_no_warning(create_prob(c(1, -1)))
  # empty input
  expect_identical(create_prob(numeric(0)), NA_real_)
})

test_that("create_prob handles NA values by na_rm", {
  # default propagate: same-length all-NA vector
  expect_identical(as.numeric(create_prob(c(4, NA, 1))), rep(NA_real_, 3))
  expect_length(create_prob(c(4, NA, 1)), 3)
  expect_no_warning(create_prob(c(4, NA, 1)))
  # all-NA input: propagates as a same-length all-NA vector
  expect_identical(create_prob(c(NA, NA)), c(NA_real_, NA_real_))
  # na_rm = TRUE renormalizes over observed values
  expect_equal(create_prob(c(NA, 2, 6), na_rm = TRUE), c(0.25, 0.75))
  expect_equal(sum(create_prob(c(NA, 2, 6), na_rm = TRUE)), 1)
  # na_rm = TRUE with only NAs remains NA_real_
  expect_identical(create_prob(c(NA, NA), na_rm = TRUE), NA_real_)
})

test_that("create_prob validates input type", {
  expect_error(create_prob(matrix(1:4, nrow = 2)), "x must be a numeric vector.")
  expect_error(create_prob(c("a", "b")), "x must be a numeric vector.")
  expect_error(create_prob(list(1, 2)), "x must be a numeric vector.")
  expect_error(create_prob(data.frame(a = 1:3)), "x must be a numeric vector.")
})

test_that("expression module integrates: per-sample create_prob over a GxS matrix", {
  set.seed(21)
  G <- 5
  S <- 6
  m <- matrix(runif(G * S, 0, 10), nrow = G, ncol = S)
  m[, 4] <- 0 # all-zero sample column
  em <- expr_matrix(m)
  # column shares: a zero-total column collapses to a length-1 NA under
  # apply's simplify2array, producing a list of mixed lengths
  probs <- apply(unclass(em), 2, create_prob)
  ok_cols <- setdiff(seq_len(S), 4)
  expect_equal(vapply(probs[ok_cols], sum, numeric(1)), rep(1, S - 1))
  expect_identical(probs[[4]], NA_real_)
})

test_that("entropy matches hand-computed Shannon entropy", {
  # concentrated: 0 bits regardless of base
  expect_equal(entropy(create_prob(c(1, 0, 0))), 0)
  expect_equal(entropy(create_prob(c(1, 0, 0)), log_base = exp(1)), 0)
  # two samples (4, 1): p = (0.8, 0.2); H = -(0.8 lb 0.8 + 0.2 lb 0.2)
  h_bits <- -(0.8 * log(0.8, 2) + 0.2 * log(0.2, 2))
  expect_equal(entropy(create_prob(c(4, 1))), h_bits)
  expect_equal(entropy(c(0.8, 0.2)), h_bits) # probabilities direct
  # nats: same p, natural log
  h_nats <- -(0.8 * log(0.8) + 0.2 * log(0.2))
  expect_equal(entropy(c(0.8, 0.2), log_base = exp(1)), h_nats)
  # uniform: maximum
  expect_equal(entropy(c(0.25, 0.25, 0.25, 0.25)), 2) # bits, 4 entries
  expect_equal(entropy(rep(1/5, 5)), log(5, 2))
})

test_that("entropy treats zero entries as zero contribution", {
  # 0 log 0 = 0: no NaN, no warning
  expect_no_warning(h <- entropy(c(0.5, 0.5, 0)))
  expect_equal(h, 1)
  # a vector that is one 1 and all zeros degenerates to H = 0
  expect_equal(entropy(c(1, 0, 0, 0, 0)), 0)
  # sum < 1 is allowed (no renormalization): p = (0.1, 0.1)
  expect_equal(entropy(c(0.1, 0.1)), 2 * -(0.1 * log(0.1, 2)))
})

test_that("entropy normalized lies in [0, 1]", {
  expect_equal(entropy(c(0.25, 0.25, 0.25, 0.25), normalized = TRUE), 1)
  expect_equal(entropy(c(1, 0, 0), normalized = TRUE), 0)
  # nats base: uniform max is ln(3)
  expect_equal(entropy(rep(1/3, 3), log_base = exp(1), normalized = TRUE), 1)
  # intermediate value hand-check: p = (0.8, 0.2), max = 1 bit
  expect_equal(entropy(c(0.8, 0.2), normalized = TRUE),
               -(0.8 * log(0.8, 2) + 0.2 * log(0.2, 2)))
  # fewer than two entries: NA
  expect_identical(entropy(1, normalized = TRUE), NA_real_)
})

test_that("entropy handles NA values by na_rm", {
  # default: propagate
  expect_identical(entropy(c(0.8, NA, 0.2)), NA_real_)
  expect_no_warning(entropy(c(0.8, NA, 0.2)))
  # na_rm: compute on observed values
  expect_equal(entropy(c(0.8, NA, 0.2), na_rm = TRUE),
               entropy(c(0.8, 0.2)))
  # all-NA
  expect_identical(entropy(c(NA, NA), na_rm = TRUE), NA_real_)
  expect_identical(entropy(c(NA, NA)), NA_real_)
})

test_that("entropy validates input", {
  expect_error(entropy(matrix(1:4, nrow = 2)), "p must be a numeric vector.")
  expect_error(entropy(c("a", "b")), "p must be a numeric vector.")
  expect_error(entropy(c(0.5, -0.1)), "probabilities between 0 and 1")
  expect_error(entropy(c(0.5, 1.2)), "probabilities between 0 and 1")
  expect_error(entropy(c(0.5, 0.5), log_base = 1), "log_base")
  expect_error(entropy(c(0.5, 0.5), log_base = -2), "log_base")
  expect_error(entropy(c(0.5, 0.5), log_base = c(2, 2)), "log_base")
  expect_error(entropy(c(0.5, 0.5), normalized = "yes"), "normalized")
})

test_that("expression module integrates: entropy specificity over a GxS matrix", {
  set.seed(21)
  G <- 5
  S <- 6
  m <- matrix(runif(G * S, 0, 10), nrow = G, ncol = S)
  m[1, ] <- 7 # ubiquitous gene
  m[2, ] <- c(10, rep(0, S - 1)) # single-sample gene
  em <- expr_matrix(m)
  # per-gene sample shares, then entropy
  h_bits <- apply(unclass(em), 1, function(col) entropy(create_prob(col)))
  expect_equal(length(h_bits), G)
  # ubiquitous: log2(S) bits, normalized to 1
  expect_equal(h_bits[1], log(S, 2))
  expect_equal(entropy(create_prob(unclass(em)[1, ]), normalized = TRUE), 1)
  # single sample: 0
  expect_equal(h_bits[2], 0)
  # per-gene share vectors sum to 1
  shares <- apply(unclass(em), 1, create_prob)
  expect_equal(colSums(shares), rep(1, G))
})