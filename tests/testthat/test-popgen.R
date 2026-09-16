test_that("fitfreq returns q unchanged when s = 0", {
  expect_equal(fitfreq(q = 0.3, h = 0.5, s = 0), 0.3)
})

test_that("fitfreq returns a value in [0, 1]", {
  q <- runif(1)
  out <- fitfreq(q = q, h = 0.5, s = 0.1)
  expect_true(out >= 0 && out <= 1)
})

test_that("WF_sel returns a length-G vector of frequencies", {
  out <- WF_sel(N = 100, q = 0.1, h = 0.5, s = 0.01, G = 50)
  expect_length(out, 50)
  expect_true(all(out >= 0 & out <= 1))
})

test_that("WF_sel recycles a scalar s", {
  out <- WF_sel(N = 100, q = 0.1, h = 0.5, s = 0.01, G = 10)
  expect_length(out, 10)
})

test_that("WF_sel accepts a vector s of length G - 1", {
  s <- rep(0.01, 9)
  out <- WF_sel(N = 100, q = 0.1, h = 0.5, s = s, G = 10)
  expect_length(out, 10)
})

test_that("WF_sel errors on wrong-length s", {
  expect_error(
    WF_sel(N = 100, q = 0.1, h = 0.5, s = c(0.01, 0.02), G = 10),
    "scalar or a vector of length G - 1"
  )
})

test_that("fc computes the standardized variance of allele frequency change", {
  # single locus: (0.5 - 0.25)^2 / (0.5 * 0.5) = 0.25
  expect_equal(fc(p0 = 0.5, pt = 0.25), 0.25)

  # two loci: (0.25 + 0.25) / (0.25 + 0.21) = 0.3125 / 0.46
  expect_equal(
    fc(p0 = c(0.5, 0.7), pt = c(0.25, 0.2)),
    0.3125 / 0.46
  )

  # identical frequencies give 0
  expect_equal(fc(p0 = c(0.3, 0.6), pt = c(0.3, 0.6)), 0)
})

test_that("fc skips NA pairs and validates inputs", {
  expect_equal(
    fc(p0 = c(0.5, NA), pt = c(0.25, 0.1)),
    0.25
  )
  expect_error(fc(p0 = c(0.5, 0.6), pt = 0.25), "length")
  expect_error(fc(p0 = 0.5, pt = 1.2), "p0 <= 1|pt <= 1")
})

test_that("waples_ne applies the Waples 1989 selfing correction", {
  t <- 2
  ne <- waples_ne(r_low = 2, r_high = 3, p0 = 0.5, pt = 0.25, S0 = 20, St = 20, t = t)
  expect_length(ne, 2)

  # hand-computed with Fc = fc(0.5, 0.25) = 0.25
  Fc <- 0.25
  expect_equal(
    ne,
    (c(2, 3) * t - 2) / (2 * c(2, 3) * (Fc - 1 / 20 - 1 / 20))
  )
  # here Fc - 1/S0 - 1/St = 0.15 > 0, so the (r*t - 2)/(2r) ratio, which
  # decreases in r, drives Ne downward as r grows
  expect_equal(ne, c(10 / 3, 40 / 9))
})

test_that("hill_weir_r2 matches the closed form at d = 0 and decays with distance", {
  # at d = 0, r2 = 5/11 * (1 + 36/(22n)), independent of C
  n <- 30
  expect_equal(
    hill_weir_r2(d = 0, n = n, C = 0.001),
    (10 / 22) * (1 + 36 / (n * 22))
  )
  expect_equal(
    hill_weir_r2(d = 0, n = n, C = 5),
    hill_weir_r2(d = 0, n = n, C = 0.001)
  )

  # LD decays with distance, and the function is vectorized in d
  r2 <- hill_weir_r2(d = c(0, 1000, 10000), n = n, C = 0.001)
  expect_length(r2, 3)
  expect_true(all(diff(r2) < 0))
})

# ---------------------------------------------------------------------------
# Dataset kit
# ---------------------------------------------------------------------------

make_test_dataset <- function(L = 10) {
  set.seed(1)
  p0 <- runif(L, 0.2, 0.8)
  pops <- c("AA", "BB")
  times <- c(0, 5, 10)
  reps <- c("R1", "R2")
  meta <- expand.grid(replicate = reps, time_point = times, population = pops)
  meta$sample_size <- rep(c(30, 25), length.out = nrow(meta))
  meta <- meta[order(meta$population, meta$time_point, meta$replicate), ]
  meta <- meta[, c("population", "time_point", "replicate", "sample_size")]
  S <- nrow(meta)
  fm <- matrix(runif(L * S, 0.05, 0.95), nrow = L)
  colnames(fm) <- paste(meta$population, meta$time_point, meta$replicate, sep = "_")
  coords <- data.frame(chrom = rep(c("1", "2"), length.out = L),
                       pos = seq_len(L) * 1000)
  list(
    freq_mat = freq_matrix(fm),
    coords = snp_coords(coords),
    p0 = p0_vec(p0),
    meta = sample_info(meta)
  )
}

test_that("freq_matrix validates input and keeps the matrix", {
  m <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2)
  fm <- freq_matrix(m)
  expect_s3_class(fm, "freq_matrix")
  expect_identical(unclass(fm), unname(m))
  expect_error(freq_matrix(matrix(c(-0.1, 0.5), nrow = 1)), "between 0 and 1")
  expect_error(freq_matrix(matrix(c(1.5, 0.5), nrow = 1)), "between 0 and 1")
  expect_error(freq_matrix(data.frame(a = 0.5)), "numeric matrix")
  fm_na <- freq_matrix(matrix(c(NA, 0.5), nrow = 1))
  expect_s3_class(fm_na, "freq_matrix")
})

test_that("snp_coords normalizes names and validates", {
  d <- data.frame(chrom = c("1", "1", "2"), pos = c(100, 250, 80))
  sc <- snp_coords(d)
  expect_s3_class(sc, "snp_coords")
  expect_identical(names(sc), c("CHROM", "POS"))
  expect_identical(as.character(sc$CHROM), c("1", "1", "2"))
  expect_error(snp_coords(data.frame(chr = "1", positions = 100)), "chrom")
  expect_error(snp_coords(data.frame(chrom = "1", pos = 0)), "at least 1")
  expect_error(snp_coords(data.frame(chrom = NA, pos = 100)), "NA")
  expect_error(snp_coords(data.frame(chrom = "1", pos = 100)[, 1, drop = FALSE]),
               "exactly 2 columns")
})

test_that("p0_vec coerces single-column matrices and validates", {
  expect_s3_class(p0_vec(c(0.1, 0.5)), "p0_vec")
  expect_s3_class(p0_vec(matrix(0.5, nrow = 2)), "p0_vec")
  expect_error(p0_vec(1.5), "between 0 and 1")
  expect_error(p0_vec("a"), "numeric vector")
})

test_that("switch_tracked_allele flips a shared random subset of variants", {
  set.seed(42)
  d <- make_test_dataset(L = 50)
  fm_orig <- unclass(d$freq_mat)
  p0_orig <- unclass(d$p0)
  out <- switch_tracked_allele(d$freq_mat, d$p0)
  expect_s3_class(out$freq_mat, "freq_matrix")
  expect_s3_class(out$p0, "p0_vec")
  # derive switched rows from p0, then verify freq_mat flips exactly those
  changed_rows <- which(out$p0 != p0_orig)
  expect_true(length(changed_rows) > 0 && length(changed_rows) < 50)
  expect_true(all(out$p0[changed_rows] == 1 - p0_orig[changed_rows]))
  expect_true(all(out$p0[-changed_rows] == p0_orig[-changed_rows]))
  expect_identical(unclass(out$freq_mat)[changed_rows, ],
                   1 - fm_orig[changed_rows, ])
  expect_identical(unclass(out$freq_mat)[-changed_rows, ],
                   fm_orig[-changed_rows, ])
  # colnames and dimnames survive
  expect_identical(colnames(out$freq_mat), colnames(d$freq_mat))
  expect_identical(dim(out$freq_mat), dim(d$freq_mat))
})

test_that("switch_tracked_allele respects a user-supplied idx", {
  d <- make_test_dataset(L = 5)
  fm_orig <- unclass(d$freq_mat)
  p0_orig <- unclass(d$p0)
  idx <- c(TRUE, FALSE, TRUE, FALSE, TRUE)
  out1 <- switch_tracked_allele(d$freq_mat, d$p0, idx = idx)
  out2 <- switch_tracked_allele(d$freq_mat, d$p0, idx = idx)
  expect_identical(unclass(out1$freq_mat), unclass(out2$freq_mat))
  expect_identical(unclass(out1$p0), unclass(out2$p0))
  expect_identical(unclass(out1$freq_mat)[idx, ], 1 - fm_orig[idx, ])
  expect_identical(unclass(out1$freq_mat)[!idx, ], fm_orig[!idx, ])
  expect_identical(unclass(out1$p0)[idx], 1 - p0_orig[idx])
  expect_identical(unclass(out1$p0)[!idx], p0_orig[!idx])
})

test_that("switch_tracked_allele passes NA frequencies through", {
  m <- matrix(c(0.2, 0.8, NA, 0.6), nrow = 2)
  fm <- freq_matrix(m)
  p0 <- p0_vec(c(0.3, 0.7))
  out <- switch_tracked_allele(fm, p0, idx = c(TRUE, FALSE))
  expect_true(is.na(out$freq_mat[1, 2]))
  expect_identical(unclass(out$freq_mat)[2, ], m[2, ])
  expect_identical(out$p0, p0_vec(c(0.7, 0.7)))
})

test_that("switch_tracked_allele validates inputs", {
  d <- make_test_dataset(L = 5)
  expect_error(
    switch_tracked_allele(unclass(d$freq_mat), d$p0),
    "freq_mat must be a freq_matrix object"
  )
  expect_error(
    switch_tracked_allele(d$freq_mat, unclass(d$p0)),
    "p0 must be a p0_vec object"
  )
  expect_error(
    switch_tracked_allele(d$freq_mat, p0_vec(d$p0[1:3])),
    "Variant count mismatch"
  )
  expect_error(
    switch_tracked_allele(d$freq_mat, d$p0, idx = c(1, 0, 1, 0, 1)),
    "logical"
  )
  expect_error(
    switch_tracked_allele(d$freq_mat, d$p0, idx = rep(TRUE, 4)),
    "one entry per variant"
  )
  expect_error(
    switch_tracked_allele(d$freq_mat, d$p0, idx = c(TRUE, NA, TRUE, FALSE, TRUE)),
    "no NA"
  )
})

test_that("filter_fixations drops boundary and NA p0 variants", {
  m <- matrix(c(
    0.0, 0.3,
    0.4, 0.5,
    0.1, 0.2,
    0.6, 0.7,
    0.9, 0.9,
    0.8, 0.9
  ), nrow = 6, byrow = TRUE)
  fm <- freq_matrix(m)
  # exact tol kept (inclusive, as in the reference workflow); above tol dropped
  p0 <- p0_vec(c(0, 0.4, 1 / 231, 0.6, 0.999, NA))
  out <- filter_fixations(fm, p0, tol = 1 / 231)
  expect_s3_class(out$freq_mat, "freq_matrix")
  expect_s3_class(out$p0, "p0_vec")
  expect_identical(length(out$p0), 3L)
  expect_identical(unclass(out$p0), c(0.4, 1 / 231, 0.6))
  expect_identical(unclass(out$freq_mat), m[2:4, , drop = FALSE])
  # snp_coords not supplied -> NULL
  expect_null(out$snp_coords)
})

test_that("filter_fixations marks boundary entries NA in kept rows", {
  m <- matrix(c(
    0.0,   0.4, 0.4,
    1.0,   0.6, 0.6,
    0.001, 0.5, 0.5
  ), nrow = 3, byrow = TRUE)
  fm <- freq_matrix(m)
  p0 <- p0_vec(c(0.5, 0.5, 0.5))
  out <- filter_fixations(fm, p0, tol = 1 / 1000)
  expect_true(is.na(out$freq_mat[1, 1]))
  expect_true(is.na(out$freq_mat[2, 1]))
  expect_identical(unclass(out$freq_mat)[3, 1], 0.001)
  expect_identical(unclass(out$freq_mat)[, 2:3], m[, 2:3])
})

test_that("filter_fixations subsets snp_coords along with rows", {
  m <- matrix(c(0.0, 0.5, 0.9), nrow = 3)
  fm <- freq_matrix(m)
  p0 <- p0_vec(c(0, 0.5, 0.8))
  coords <- snp_coords(data.frame(chrom = c("1", "1", "2"), pos = c(100, 200, 300)))
  out <- filter_fixations(fm, p0, coords, tol = 1 / 231)
  expect_s3_class(out$snp_coords, "snp_coords")
  expect_identical(nrow(out$snp_coords), 2L)
  expect_identical(out$snp_coords$POS, c(200, 300))
})

test_that("filter_fixations validates inputs", {
  d <- make_test_dataset(L = 5)
  expect_error(
    filter_fixations(unclass(d$freq_mat), d$p0, tol = 0.01),
    "freq_mat must be a freq_matrix object"
  )
  expect_error(
    filter_fixations(d$freq_mat, unclass(d$p0), tol = 0.01),
    "p0 must be a p0_vec object"
  )
  expect_error(
    filter_fixations(d$freq_mat, d$p0, unclass(d$coords), tol = 0.01),
    "snp_coords must be a snp_coords object"
  )
  expect_error(
    filter_fixations(d$freq_mat, p0_vec(d$p0[1:3]), tol = 0.01),
    "Variant count mismatch"
  )
  expect_error(
    filter_fixations(d$freq_mat, d$p0, d$coords[1:2, ], tol = 0.01),
    "Variant count mismatch"
  )
  expect_error(
    filter_fixations(d$freq_mat, d$p0, tol = "a"),
    "single non-negative number"
  )
  expect_error(
    filter_fixations(d$freq_mat, d$p0, tol = c(0.01, 0.02)),
    "single non-negative number"
  )
  expect_error(
    filter_fixations(d$freq_mat, d$p0, tol = NA),
    "single non-negative number"
  )
  expect_error(
    filter_fixations(d$freq_mat, d$p0, tol = -0.01),
    "single non-negative number"
  )
  expect_error(
    filter_fixations(d$freq_mat, d$p0, tol = Inf),
    "single non-negative number"
  )
})

test_that("filter_fixations errors when all variants are dropped", {
  m <- matrix(c(0.0, 1.0, 1 / 231, 0.5), nrow = 2)
  fm <- freq_matrix(m)
  p0 <- p0_vec(c(0, 0.0001))
  expect_error(
    filter_fixations(fm, p0, tol = 1 / 231),
    "nothing to keep"
  )
  p0_na <- p0_vec(c(NA_real_, NA_real_))
  expect_error(
    filter_fixations(freq_matrix(matrix(c(0.5, 0.5), nrow = 2)), p0_na, tol = 0.01),
    "nothing to keep"
  )
})

test_that("sample_info coerces column types and validates", {
  d <- data.frame(
    population = c("AA", "BB"),
    time_point = c(0, 10),
    replicate = c("R1", "R1"),
    sample_size = c(30, 25)
  )
  si <- sample_info(d)
  expect_s3_class(si, "sample_info")
  expect_type(si$population, "character")
  expect_type(si$time_point, "integer")
  expect_type(si$replicate, "character")
  expect_type(si$sample_size, "integer")
  # column order is normalized
  si2 <- sample_info(d[, c("sample_size", "replicate", "time_point",
                           "population")])
  expect_identical(names(si2), c("population", "time_point", "replicate",
                                 "sample_size"))
  expect_error(sample_info(d[, 1:3]), "exactly 4 columns")
  expect_error(
    sample_info(data.frame(pop = c("AA", "BB"), time = c(0, 10),
                           rep = c("R1", "R1"), size = c(30, 25))),
    "population.*time_point.*replicate"
  )
  expect_error(
    sample_info(data.frame(population = c("AA", NA), time_point = c(0, 10),
                           replicate = c("R1", "R1"), sample_size = c(30, 25))),
    "NA"
  )
  expect_error(
    sample_info(data.frame(population = c("AA", "BB"), time_point = c(-1, 10),
                           replicate = c("R1", "R1"), sample_size = c(30, 25))),
    "non-negative"
  )
  expect_error(
    sample_info(data.frame(population = c("AA", "BB"), time_point = c(0, 10),
                           replicate = c("R1", "R1"), sample_size = c(0, 25))),
    "positive integer"
  )
  expect_error(
    sample_info(data.frame(population = c("AA", "BB"), time_point = c(0, 10),
                           replicate = c("R1", "R1"), sample_size = c(-5, 25))),
    "positive integer"
  )
  # non-numeric sample_size coerces to NA (with a coercion warning) and is
  # caught by the NA check
  expect_error(
    expect_warning(
      sample_info(data.frame(population = c("AA", "BB"), time_point = c(0, 10),
                             replicate = c("R1", "R1"),
                             sample_size = c("30", "big"))),
      "coercion"
    ),
    "NA"
  )
})

test_that("print methods summarize the four dataset objects", {
  d <- make_test_dataset()
  expect_output(print(d$freq_mat), "10 variants x 12 samples")
  expect_output(print(d$coords), "10 variants on 2 chromosome")
  expect_output(print(d$p0), "10 initial allele frequencies")
  expect_output(print(d$meta), "12 samples")
})

test_that("validate_af_dataset passes consistent objects and flags mismatches", {
  d <- make_test_dataset()
  expect_true(validate_af_dataset(d$freq_mat, d$coords, d$p0, d$meta))

  # variant mismatch
  expect_error(
    validate_af_dataset(d$freq_mat, d$coords[1:5, ], d$p0, d$meta),
    "Variant count mismatch"
  )
  # sample mismatch
  expect_error(
    validate_af_dataset(d$freq_mat, d$coords, d$p0, d$meta[1:6, ]),
    "Sample count mismatch"
  )
  # colname scheme check
  bad <- unname(d$freq_mat)
  colnames(bad) <- paste0("s", seq_len(ncol(bad)))
  expect_error(
    validate_af_dataset(freq_matrix(bad), d$coords, d$p0, d$meta),
    "colnames do not match"
  )
  # colname check skipped when names are absent
  expect_true(validate_af_dataset(freq_matrix(unname(d$freq_mat)),
                                  d$coords, d$p0, d$meta))
  # class enforcement: reclassing through the constructors keeps it passing
  expect_true(validate_af_dataset(d$freq_mat, d$coords, d$p0, d$meta,
                                  check_colnames = FALSE))
  expect_true(validate_af_dataset(d$freq_mat, d$coords, p0_vec(d$p0),
                                  sample_info(d$meta)))
})

test_that("extract_samples subsets columns, rebuilds names, attaches metadata", {
  d <- make_test_dataset()
  sub <- extract_samples(d$freq_mat, d$meta, population = "AA")
  expect_s3_class(sub, "freq_matrix")
  expect_identical(ncol(sub), 6L)
  expect_true(all(grepl("^AA_", colnames(sub))))
  att <- attr(sub, "sample_info")
  expect_identical(
    colnames(sub),
    paste(att$population, att$time_point, att$replicate, sep = "_")
  )
  expect_s3_class(att, "sample_info")
  expect_true(all(att$population == "AA"))

  # combined filters
  sub2 <- extract_samples(d$freq_mat, d$meta, population = "BB", time = c(0, 5),
                          replicate = "R2")
  expect_identical(ncol(sub2), 2L)
  expect_identical(colnames(sub2), c("BB_0_R2", "BB_5_R2"))

  # no match errors
  expect_error(extract_samples(d$freq_mat, d$meta, population = "ZZ"),
               "No samples match")
})

# ---------------------------------------------------------------------------
# Ported user tests (describe/it blocks converted to test_that)
# ---------------------------------------------------------------------------

test_that("freq_increments correctly computes differences between adjacent generations", {
  pmat <- matrix(c(0.1, 0.2, 0.3, 0.4,
                   0.2, 0.3, 0.4, 0.5), nrow = 2, ncol = 4, byrow = TRUE)
  expected <- matrix(c(0.1, 0.1, 0.1,
                       0.1, 0.1, 0.1), nrow = 2, ncol = 3, byrow = TRUE)
  result <- freq_increments(pmat)
  expect_true(all(abs(result - expected) < 1e-9))
})

test_that("freq_increments returns a matrix with one fewer column than input", {
  pmat <- matrix(1:6, nrow = 2, ncol = 3)
  result <- freq_increments(pmat)
  expect_equal(ncol(result), ncol(pmat) - 1)
  expect_equal(nrow(result), nrow(pmat))
})

test_that("freq_increments handles 1-row matrix", {
  pmat <- matrix(c(0.1, 0.3, 0.6), nrow = 1)
  result <- freq_increments(pmat)
  expected <- matrix(c(0.2, 0.3), nrow = 1)
  expect_true(all(abs(result - expected) < 1e-9))
})

test_that("sum_of_het correctly computes sum of 2*x*(1-x) for a vector", {
  expect_equal(sum_of_het(c(0.5)), 0.5)
  expect_equal(sum_of_het(c(0.0, 1.0)), 0)
  expect_equal(sum_of_het(c(1/3, 1/3, 1/3)), 2 * (1/3) * (2/3) * 3)
})

test_that("sum_of_het returns zero for fixation states", {
  expect_equal(sum_of_het(c(0, 0, 0)), 0)
  expect_equal(sum_of_het(c(1, 1, 1)), 0)
})

test_that("geom_pairwise_mean rejects vectors with 1 element", {
  expect_error(geom_pairwise_mean(c(10)))
})

test_that("geom_pairwise_mean rejects non-numeric vector", {
  expect_error(geom_pairwise_mean(c("hello", "world", "blake")))
})

test_that("geom_pairwise_mean returns 1 when input is all 1s", {
  expect_equal(geom_pairwise_mean(rep(1, times = 10)), 1)
})

test_that("geom_pairwise_mean hand-computed values", {
  # n = 2: sqrt(2*8) = 4
  expect_equal(geom_pairwise_mean(c(2, 8)), 4)
  # n = 3: (sqrt(2*4) + sqrt(2*8) + sqrt(4*8)) / 3 = (2.83 + 4 + 5.66) / 3
  expect_equal(geom_pairwise_mean(c(2, 4, 8)),
               (sqrt(2 * 4) + sqrt(2 * 8) + sqrt(4 * 8)) / 3)
})

test_that("rm_na_after_na returns unchanged vector when no NA present", {
  x <- c(1, 2, 3, 4, 5)
  expect_true(all(rm_na_after_na(x) == x))
})

test_that("rm_na_after_na na-truncates from first NA position", {
  x <- c(1, 2, NA, 4, 5)
  expect_equal(rm_na_after_na(x), c(1, 2, NA, NA, NA))
})

test_that("rm_na_after_na handles vector with NA at start", {
  x <- c(NA, 2, 3)
  expect_true(all(is.na(rm_na_after_na(x))))
})

test_that("rm_na_after_na returns vector unchanged when NA only at the end", {
  x <- c(1, 2, NA)
  expect_equal(rm_na_after_na(x), c(1, 2, NA))
})

test_that("arcsin_sqrt returns endpoints correctly", {
  expect_equal(arcsin_sqrt(1), pi / 2)
  expect_equal(arcsin_sqrt(0), 0)
})

test_that("arcsin_sqrt is monotonic", {
  expect_lt(arcsin_sqrt(0.1), arcsin_sqrt(0.2))
})

test_that("arcsin_sqrt catches bad frequencies", {
  expect_error(arcsin_sqrt(50))
  expect_error(arcsin_sqrt(-0.1))
})

test_that("arcsin_sqrt NA input gives NA output", {
  expect_equal(arcsin_sqrt(c(0, NA, 1)), c(0, NA, pi / 2))
})

test_that("sign_permute_increments returns input as output for none procedure", {
  set.seed(42)
  mat <- as.data.frame(matrix(runif(10 * 10, min = -1, max = 1), nrow = 10, ncol = 10))
  result <- sign_permute_increments(mat, procedure = "none")
  expect_equal(result, mat)
})

test_that("sign_permute_increments returns output with same dim as input", {
  set.seed(42)
  mat <- as.data.frame(matrix(runif(10 * 10, min = -1, max = 1), nrow = 10, ncol = 10))
  result_cell <- sign_permute_increments(mat, procedure = "cell")
  result_window <- sign_permute_increments(
    mat, procedure = "window",
    windows = c(rep(1, times = 5), rep(2, times = 5))
  )
  result_genome <- sign_permute_increments(mat, procedure = "genome")
  expect_equal(dim(result_cell), dim(mat))
  expect_equal(dim(result_window), dim(mat))
  expect_equal(dim(result_genome), dim(mat))
})

test_that("sign_permute_increments preserves absolute values in cell and window modes", {
  set.seed(42)
  mat <- as.data.frame(matrix(runif(20, min = -1, max = 1), nrow = 5, ncol = 4))
  result_cell <- sign_permute_increments(mat, procedure = "cell")
  expect_true(all(abs(as.matrix(result_cell)) == abs(as.matrix(mat))))

  # window mode: one sign per (window, column); constant within each window block
  windows <- c(1, 1, 2, 2, 2)
  result_window <- sign_permute_increments(mat, procedure = "window", windows = windows)
  signed <- as.matrix(result_window) / as.matrix(mat)
  for (w in unique(windows)) {
    idx <- which(windows == w)
    block <- signed[idx, , drop = FALSE]
    expect_true(all(apply(block, 2, function(col) length(unique(col)) == 1)))
  }
})

test_that("sign_permute_increments genome mode flips whole columns", {
  set.seed(7)
  mat <- as.data.frame(matrix(runif(20, min = -1, max = 1), nrow = 5, ncol = 4))
  result_genome <- sign_permute_increments(mat, procedure = "genome")

  # each result column is a single +1/-1 multiple of the corresponding input column
  for (j in seq_len(ncol(mat))) {
    col_signs <- as.matrix(result_genome)[, j] / as.matrix(mat)[, j]
    expect_true(all(col_signs == col_signs[1]))
    expect_true(all(abs(col_signs) == 1))
  }
})

test_that("sign_permute_increments errors on invalid procedure", {
  mat <- as.data.frame(matrix(runif(4), nrow = 2))
  expect_error(sign_permute_increments(mat, procedure = "bogus"), "Procedure not valid")
  expect_error(sign_permute_increments(mat, procedure = "window", windows = NULL),
               "windows")
})

test_that("covmat_from_pmat without correction equals cov of increments", {
  set.seed(1)
  pmat <- matrix(runif(15), nrow = 5, ncol = 3)
  pd <- pmat[, -1] - pmat[, -ncol(pmat), drop = FALSE]
  expect_equal(
    unname(covmat_from_pmat(pmat, correct_for_n = FALSE)),
    unname(stats::cov(pd))
  )
})

test_that("covmat_from_pmat raw sample-size correction matches hand computation", {
  set.seed(1)
  pmat <- matrix(runif(15), nrow = 5, ncol = 3)
  n <- c(50, 50, 50)
  pd <- pmat[, -1] - pmat[, -ncol(pmat), drop = FALSE]

  cm <- stats::cov(pd)
  cm[1, 2] <- cm[2, 1] <- mean(pd[, 1] * pd[, 2]) +
    mean(pmat[, 2] * (1 - pmat[, 2]) / (n[2] - 1))
  v1 <- mean(pd[, 1]^2) - mean(pmat[, 1] * (1 - pmat[, 1]) / (n[1] - 1)) -
    mean(pmat[, 2] * (1 - pmat[, 2]) / (n[2] - 1))
  v2 <- mean(pd[, 2]^2) - mean(pmat[, 2] * (1 - pmat[, 2]) / (n[2] - 1)) -
    mean(pmat[, 3] * (1 - pmat[, 3]) / (n[3] - 1))
  cm[1, 1] <- ifelse(v1 < 0, 0, v1)
  cm[2, 2] <- ifelse(v2 < 0, 0, v2)

  expect_equal(
    unname(covmat_from_pmat(pmat, n = n, correct_for_n = TRUE)),
    unname(cm)
  )
})

test_that("covmat_from_pmat warns and zeroes negative variances", {
  pmat <- matrix(c(0.3, 0.4, 0.3, 0.4), nrow = 2)
  expect_warning(
    cm <- covmat_from_pmat(pmat, n = c(5, 5), correct_for_n = TRUE),
    "variance negative"
  )
  expect_equal(cm[1, 1], 0)
})

test_that("covmat_from_pmat warns on untransformed input with asin correction", {
  set.seed(1)
  pmat <- matrix(runif(15), nrow = 5, ncol = 3)
  expect_warning(
    covmat_from_pmat(pmat, n = c(50, 50, 50), input_asin_trans = TRUE),
    "arcsin transformed"
  )
})

test_that("covmat_from_pmat validates sample sizes", {
  pmat <- matrix(c(0.2, 0.4, 0.3, 0.5), nrow = 2)
  expect_error(covmat_from_pmat(pmat, n = c(1, 1)), "All n must be >= 2")
  expect_error(covmat_from_pmat(pmat, n = c(10, 10, 10)), "sample size for every time point")
})

test_that("covmat_from_pmat works with a single interval (2 time points)", {
  pmat <- matrix(c(0.2, 0.4, 0.3, 0.5), nrow = 2)
  cm <- covmat_from_pmat(pmat, correct_for_n = FALSE)
  expect_equal(dim(cm), c(1, 1))
  pd <- pmat[, 2] - pmat[, 1]
  expect_equal(unname(cm[1, 1]), stats::var(pd))
})

test_that("standard_cov_by_het divides by half heterozygosity sums", {
  pmat <- matrix(c(0.5, 0.2, 0.6, 0.3), nrow = 2) # 2 intervals x 2 variants
  covmat <- matrix(c(4, 1, 1, 9), nrow = 2)
  out <- standard_cov_by_het(pmat, covmat)
  h1 <- 0.5 * sum_of_het(pmat[, 1])
  h2 <- 0.5 * sum_of_het(pmat[, 2])
  expect_equal(out[1, 1], 4 / h1)
  expect_equal(out[2, 2], 9 / h2)
  expect_equal(out[1, 2], 1 / h1) # min(1,2) = 1
  expect_equal(out[2, 1], 1 / h1)
})

test_that("rolling_matrix_sum sums successively larger top-left sub-squares", {
  set.seed(456)
  mat <- matrix(runif(10 * 10, min = -1, max = 1), nrow = 10, ncol = 10)
  result <- rolling_matrix_sum(mat)
  expect_equal(result[1], sum(mat[1, 1]))
  expect_equal(result[2], sum(mat[1:2, 1:2]))
  expect_equal(result[3], sum(mat[1:3, 1:3]))
  expect_equal(result[length(result)], sum(mat))
})

test_that("rolling_matrix_sum errors on non-square input", {
  expect_error(rolling_matrix_sum(matrix(1:9, nrow = 3)), NA)
  expect_error(rolling_matrix_sum(matrix(1:6, nrow = 2, ncol = 3)),
               "nrow\\(mat\\) == ncol\\(mat\\)")
})

test_that("gt_from_covmat returns 0 gt for 2x2 matrix with only one off-diagonal", {
  covmat <- matrix(c(1, 0.5, 0.5, 2), nrow = 2)
  result <- gt_from_covmat(covmat)
  expect_equal(nrow(result), 2)
  # First row: sums_cov=0, sums_var=1 => gt=0/(1+0)=0
  expect_equal(result$gt[1], 0)
})

test_that("gt_from_covmat handles symmetric covariance matrix", {
  covmat <- matrix(c(4, 1, 1, 9), nrow = 2)
  result <- gt_from_covmat(covmat)
  expect_equal(result$sum_var[2], 13) # 4+9
  expect_equal(result$sum_cov[2], 2)
  expect_equal(result$sum_abs_cov[2], 2)
})

test_that("gt_from_covmat hand-computed 3x3 values", {
  covmat <- matrix(c(1, 0.5, 0.2, 0.5, 2, 0.3, 0.2, 0.3, 3), nrow = 3)
  res <- gt_from_covmat(covmat)
  expect_equal(res$gen, 1:3)
  expect_equal(res$gt[1], 0)
  # gen 2: sum_cov = 2*0.5 = 1, sum_var = 3
  expect_equal(res$sum_cov[2], 1)
  expect_equal(res$gt[2], 1 / 4)
  # gen 3: sum_cov = 2*(0.5 + 0.2 + 0.3) = 2, sum_var = 6
  expect_equal(res$sum_cov[3], 2)
  expect_equal(res$gt[3], 2 / 8)
  expect_equal(res$sum_abs_cov[3], 2)
  expect_equal(res$sum_pos_cov[3], 2)
  expect_equal(res$sum_neg_cov[3], 0)
})

test_that("replicate_gt correctly scales denominator of expectation with replicates", {
  set.seed(123)
  pdiff <- matrix(runif(30 * 10, min = -1, max = 1), nrow = 30, ncol = 10)
  rep_labels <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)
  time_labels <- c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2)
  result <- replicate_gt(pdiff, rep_labels, time_labels)
  expect_equal(result$n_covar, 10)
  expect_equal(result$n_var, 5)
})

test_that("replicate_gt hand-computed 2x2 partition", {
  pdiff <- matrix(c(0.1, 0.2, -0.1, -0.2), nrow = 2, ncol = 2)
  res <- replicate_gt(pdiff, rep_labels = c(1, 2), time_labels = c(1, 2))
  # cov(pd1, pd2) for these vectors: mean products - cov of means
  expect_equal(res$n_var, 2)
  expect_equal(res$n_covar, 1)
  expect_equal(res$total_var, sum(diag(stats::cov(pdiff))))
  expect_equal(res$total_covar, 2 * stats::cov(pdiff)[1, 2])
})

test_that("replicate_gt validates labels", {
  pdiff <- matrix(rnorm(20), nrow = 10, ncol = 2)
  expect_error(replicate_gt(pdiff, rep_labels = c(1), time_labels = c(1, 2)),
               "length")
  expect_error(replicate_gt(pdiff, rep_labels = c(1, 1), time_labels = c(1, 2)),
               "unique")
  expect_error(replicate_gt(pdiff, rep_labels = c(1, 2), time_labels = c(1, 1)),
               "unique")
})

test_that("conv_cor_wn_env computes numerator and denominator", {
  set.seed(1)
  pdiff <- matrix(rnorm(40, sd = 0.05), nrow = 20, ncol = 2)
  res <- conv_cor_wn_env(pdiff)
  covmat <- stats::cov(pdiff)
  expect_length(res, 2)
  expect_equal(res[1], covmat[1, 2])
  expect_equal(res[2], geom_pairwise_mean(diag(covmat)))
  expect_lte(abs(res[1]), res[2])
})

test_that("conv_cor_wn_env errors when replicates have zero variance", {
  expect_error(conv_cor_wn_env(cbind(c(0.5, 0.5, 0.5), c(0.5, 0.5, 0.5))),
               "denominator")
})

test_that("g_prime matches manual computation", {
  set.seed(1)
  pmat <- matrix(runif(15), nrow = 5, ncol = 3)
  n <- c(50, 50, 50)
  cvm <- covmat_from_pmat(pmat, n)
  ep0 <- mean(pmat[, 1] * (1 - pmat[, 1]))
  manual <- 1 - (2:3 * ep0) / (2 * 100 * rolling_matrix_sum(cvm))
  expect_equal(g_prime(times = 2:3, pmat = pmat, N = 100, n = n), manual)
  # take_abs variant differs on signed covariances
  g_abs <- g_prime(times = 2:3, pmat = pmat, N = 100, n = n, take_abs = TRUE)
  expect_false(isTRUE(all.equal(g_abs, manual)) || all(g_abs == manual))
  expect_error(g_prime(times = c(2), pmat = pmat, N = 100, n = n), "length")
  expect_error(g_prime(times = 2:3, pmat = pmat, N = 1, n = n), "N > 2")
})

test_that("covmat_pop_pair computes mean standardized between-pop covariance", {
  set.seed(1)
  L <- 10
  gen0 <- runif(L, 0.3, 0.7)
  p1 <- cbind(AA_R1 = gen0 + rnorm(L, 0, 0.02), AA_R2 = gen0 + rnorm(L, 0, 0.02))
  p2 <- cbind(BB_R1 = gen0 + rnorm(L, 0, 0.02), BB_R2 = gen0 + rnorm(L, 0, 0.02))
  val <- covmat_pop_pair(p1, p2, pop1 = "AA", pop2 = "BB", gen0 = gen0)

  pd1 <- p1 - gen0
  pd2 <- p2 - gen0
  cm <- stats::cov(cbind(pd1, pd2))
  rv1 <- diag(cm)[1:2]
  rv2 <- diag(cm)[3:4]
  mean_std <- mean(sqrt(rv1 %o% rv2))
  bw <- mean(cm[1:2, 3:4])
  expect_equal(val, bw / mean_std)
})

test_that("covmat_pop_pair requires column names", {
  gen0 <- c(0.5, 0.4)
  p1 <- cbind(gen0 + 0.01, gen0 - 0.01)
  p2 <- cbind(gen0 - 0.01, gen0 + 0.01)
  expect_error(covmat_pop_pair(p1, p2, pop1 = "AA", pop2 = "BB", gen0 = gen0),
               "colnames")
})

# ---------------------------------------------------------------------------
# Integration: dataset kit feeding the temporal analyses
# ---------------------------------------------------------------------------

test_that("dataset kit integrates with temporal functions end to end", {
  d <- make_test_dataset(L = 15)
  expect_true(validate_af_dataset(d$freq_mat, d$coords, d$p0, d$meta))

  # population AA replicates: time 0 -> 5 -> 10, two replicates
  aa <- extract_samples(d$freq_mat, d$meta, population = "AA")
  aa_meta <- attr(aa, "sample_info")

  # per-replicate pmat for covmat_from_pmat
  rep1_cols <- which(aa_meta$replicate == "R1")
  pmat_rep1 <- as.matrix(aa[, rep1_cols, drop = FALSE])
  cm <- covmat_from_pmat(pmat_rep1, correct_for_n = FALSE)
  expect_equal(dim(cm), c(2, 2))
  expect_equal(unname(cm), unname(stats::cov(freq_increments(pmat_rep1))))

  # g_prime consumes the same pmat
  g <- g_prime(times = 2:3, pmat = pmat_rep1, N = 500, n = rep(100, 3))
  expect_length(g, 2)
  expect_true(all(g <= 1))

  # two-population comparison via covmat_pop_pair
  bb <- extract_samples(d$freq_mat, d$meta, population = "BB")
  val <- covmat_pop_pair(as.matrix(aa), as.matrix(bb),
                         pop1 = "AA", pop2 = "BB",
                         gen0 <- as.matrix(aa[, aa_meta$time_point == 0, drop = FALSE])[, 1])
  expect_type(val, "double")
  expect_lte(abs(val), 1 + 1e-9)

  # replicate_gt with metadata-derived labels: columns from time 5 to 10
  late <- which(aa_meta$time_point > 0)
  pdiff <- as.matrix(aa[, late, drop = FALSE]) -
    as.matrix(aa[, which(aa_meta$time_point == 0), drop = FALSE])[, 1]
  res <- replicate_gt(pdiff,
                      rep_labels = aa_meta$replicate[late],
                      time_labels = aa_meta$time_point[late])
  expect_equal(res$n_var, 2)
  expect_equal(res$n_covar, 1)
})

# ---------------------------------------------------------------------------
# feder_t_test
# ---------------------------------------------------------------------------

test_that("feder_t_test matches a hand-computed pooled t-test", {
  meta <- data.frame(
    population = c("AA", "AA", "AA", "AA", "BB", "BB", "BB"),
    time_point = c(0L, 1L, 2L, 0L, 0L, 1L, 2L),
    replicate = c("R1", "R1", "R1", "R2", "R1", "R1", "R1"),
    sample_size = 30L
  )
  set.seed(1)
  L <- 2
  fm <- freq_matrix(matrix(runif(L * 7, 0.1, 0.9), nrow = L,
                           dimnames = list(NULL, paste(meta$population,
                               meta$time_point, meta$replicate, sep = "_"))))
  p0 <- p0_vec(c(0.4, 0.6))

  res <- feder_t_test(fm, p0, sample_info(meta))

  # AA: R1 has t0 -> t1 -> t2; R2 sampled only at t0 contributes nothing.
  # Anchors 0, 1, 2 with dt = 1 between adjacent anchors.
  p <- unclass(fm)
  y1 <- (p[, "AA_1_R1"] - p0) / sqrt(2 * p0 * (1 - p0) * 1)
  y2 <- (p[, "AA_2_R1"] - p[, "AA_1_R1"]) /
    sqrt(2 * p[, "AA_1_R1"] * (1 - p[, "AA_1_R1"]) * 1)
  yall <- cbind(y1, y2)
  aa <- res[res$population == "AA", ]
  expect_equal(aa$ybar, rowMeans(yall))
  expect_equal(aa$s2, apply(yall, 1, var))
  expect_equal(aa$n, rep(2, L))
  expect_equal(aa$df, rep(1, L))
  expect_equal(aa$tvalue, rowMeans(yall) / sqrt(apply(yall, 1, var) / 2))
  expect_equal(aa$pvalue, 2 * stats::pt(-abs(aa$tvalue), df = 1))

  # BB mirrors the single-replicate structure
  expect_equal(res$n[res$population == "BB"], rep(2, L))
  expect_equal(nrow(res), 2 * L)
})

test_that("feder_t_test agrees with stats::t.test on pooled increments", {
  meta <- data.frame(
    population = rep("AA", 6),
    time_point = c(0L, 0L, 1L, 1L, 2L, 2L),
    replicate = c("R1", "R2", "R1", "R2", "R1", "R2"),
    sample_size = 30L
  )
  set.seed(42)
  L <- 3
  fm_raw <- matrix(runif(L * 6, 0.1, 0.9), nrow = L)
  fm <- freq_matrix(`colnames<-`(fm_raw, paste(meta$population,
      meta$time_point, meta$replicate, sep = "_")))
  p0 <- p0_vec(runif(L, 0.2, 0.8))

  res <- feder_t_test(fm, p0, sample_info(meta))

  y <- cbind(
    (fm_raw[, 3:4] - unclass(p0)) /
      sqrt(2 * unclass(p0) * (1 - unclass(p0)) * 1),
    (fm_raw[, 5:6] - fm_raw[, 3:4]) /
      sqrt(2 * fm_raw[, 3:4] * (1 - fm_raw[, 3:4]) * 1)
  )
  for (i in 1:L) {
    tt <- stats::t.test(y[i, ], alternative = "two.sided")
    expect_equal(res$tvalue[i], unname(tt$statistic))
    expect_equal(res$pvalue[i], tt$p.value)
    expect_equal(res$df[i], unname(tt$parameter))
    expect_equal(res$n[i], 4) # 2 reps x 2 intervals
  }
})

test_that("feder_t_test respects one-sided alternatives", {
  meta <- data.frame(
    population = "AA",
    time_point = c(0L, 1L, 2L),
    replicate = "R1",
    sample_size = 30L
  )
  set.seed(3)
  L <- 4
  fm <- freq_matrix(matrix(runif(L * 3, 0.1, 0.9), nrow = L,
                           dimnames = list(NULL, paste(meta$population,
                               meta$time_point, meta$replicate, sep = "_"))))
  p0 <- p0_vec(runif(L, 0.2, 0.8))

  two <- feder_t_test(fm, p0, sample_info(meta))
  less <- feder_t_test(fm, p0, sample_info(meta), alternative = "less")
  greater <- feder_t_test(fm, p0, sample_info(meta), alternative = "greater")
  # one-sided probabilities partition the two-sided one
  expect_equal(two$pvalue, 2 * pmin(less$pvalue, greater$pvalue))
  expect_equal(less$pvalue + greater$pvalue, rep(1, L))
  expect_error(feder_t_test(fm, p0, sample_info(meta), alternative = "bogus"),
               "should be one of")
})

test_that("feder_t_test uses p0 as baseline with dt from generation gaps", {
  # times 0, 5, 10: first increment is p1 - p0 with dt = 5, second p2 - p1 dt = 5
  meta <- data.frame(
    population = "AA",
    time_point = c(0L, 5L, 10L),
    replicate = "R1",
    sample_size = 30L
  )
  L <- 1
  p0 <- 0.5
  p1 <- 0.6
  p2 <- 0.3
  fm <- freq_matrix(matrix(c(0.5, p1, p2), nrow = L,
                           dimnames = list(NULL, c("AA_0_R1", "AA_5_R1",
                                                   "AA_10_R1"))))
  res <- feder_t_test(fm, p0_vec(p0), sample_info(meta))
  y1 <- (p1 - p0) / sqrt(2 * p0 * (1 - p0) * 5)
  y2 <- (p2 - p1) / sqrt(2 * p1 * (1 - p1) * 5)
  expect_equal(res$n, 2)
  expect_equal(res$ybar, mean(c(y1, y2)))
  expect_equal(res$s2, var(c(y1, y2)))
})

test_that("feder_t_test handles NA frequencies and all-NA variants", {
  meta <- data.frame(
    population = rep("AA", 4),
    time_point = c(0L, 0L, 1L, 1L),
    replicate = c("R1", "R2", "R1", "R2"),
    sample_size = 30L
  )
  fm_raw <- matrix(c(0.2, 0.4, 0.6, 0.8,   # variant 1 usable
                     0.5, 0.5, NA, NA,     # variant 2: t1 missing in both reps
                     0.3, 0.3, 0.7, 0.7),  # variant 3 usable
                   nrow = 3, byrow = TRUE,
                   dimnames = list(NULL, paste(meta$population, meta$time_point,
                                               meta$replicate, sep = "_")))
  fm <- freq_matrix(fm_raw)
  p0 <- p0_vec(c(0.1, 0.5, 0.9))

  expect_no_warning(res <- feder_t_test(fm, p0, sample_info(meta)))
  # variant 1: 2 increments (R1, R2 into t1); variant 2: 0; variant 3: 2
  expect_equal(res$n, c(2, 0, 2))
  v2 <- res[2, ]
  expect_true(is.na(v2$ybar))
  expect_true(is.na(v2$s2))
  expect_true(is.na(v2$tvalue))
  expect_true(is.na(v2$pvalue))
})

test_that("feder_t_test validates inputs", {
  meta <- data.frame(
    population = "AA", time_point = c(0L, 1L), replicate = "R1",
    sample_size = 30L
  )
  meta_si <- sample_info(meta)
  p0 <- p0_vec(c(0.5, 0.5))
  fm_ok <- freq_matrix(matrix(c(0.3, 0.6, 0.4, 0.7), nrow = 2,
                              dimnames = list(NULL, c("AA_0_R1", "AA_1_R1"))))
  coords <- snp_coords(data.frame(chrom = c("1", "1"), pos = c(100, 200)))

  fm_one <- freq_matrix(`colnames<-`(matrix(0.3, nrow = 1), "AA_0_R1"))
  expect_error(feder_t_test(fm_one, p0, meta_si), "p0_vec has 2 entries")
  expect_error(feder_t_test(fm_ok, p0_vec(c(0.5, 0.5, 0.5)), meta_si),
               "p0_vec has 3 entries")
  expect_error(feder_t_test(fm_ok, p0,
                            sample_info(meta_si[1, , drop = FALSE])),
               "sample_info has 1 row")
  expect_error(feder_t_test(fm_ok, p0, meta_si,
                            snp_coords = data.frame(chrom = "1", pos = 1)),
               "snp_coords must be a snp_coords object")
  # with snp_coords, output gains CHROM/POS columns in front
  res <- feder_t_test(fm_ok, p0, meta_si, snp_coords = coords)
  expect_named(res, c("population", "CHROM", "POS", "ybar", "s2", "n", "df",
                      "tvalue", "pvalue"))
  expect_equal(res$CHROM, rep("1", 2))
  expect_equal(res$POS, rep(c(100, 200), 1))
  expect_error(feder_t_test(fm_ok, p0, meta_si,
                            snp_coords = snp_coords(coords[1, , drop = FALSE])),
               "snp_coords has 1 rows")
  expect_error(feder_t_test(fm_ok, p0, p0), "sample_meta must be a sample_info object")
  # no positive time points
  meta0 <- data.frame(
    population = "AA", time_point = 0L, replicate = "R1", sample_size = 30L
  )
  fm_zero <- freq_matrix(`colnames<-`(matrix(0.3, nrow = 1), "AA_0_R1"))
  expect_error(feder_t_test(fm_zero, p0_vec(0.5), sample_info(meta0)),
               "no positive time points")
})

test_that("feder_t_test integrates with the dataset kit", {
  d <- make_test_dataset(L = 12)
  expect_true(validate_af_dataset(d$freq_mat, d$coords, d$p0, d$meta))
  res <- feder_t_test(d$freq_mat, d$p0, d$meta)
  expect_equal(nrow(res), 2 * 12) # 2 populations x 12 variants
  expect_setequal(unique(res$population), c("AA", "BB"))
  expect_true(all(res$n > 0))
  # each row pools both replicates x 2 increments = 4 non-NA increments
  expect_true(all(res$n == 4))
  # coords subset aligned with the make_test_dataset ordering
  res_c <- feder_t_test(d$freq_mat, d$p0, d$meta, snp_coords = d$coords)
  expect_named(res_c, c("population", "CHROM", "POS", "ybar", "s2", "n", "df",
                        "tvalue", "pvalue"))
})

test_that("fit_af_glm matches an inline glm on the equivalent data", {
  set.seed(11)
  L <- 4
  meta <- rbind(expand.grid(population = "AA", time_point = c(0, 1, 2),
                            replicate = c("R1", "R2")),
                expand.grid(population = "BB", time_point = c(0, 1, 2, 3, 4),
                            replicate = "R1"))
  meta$sample_size <- 30
  meta <- meta[order(meta$population, meta$time_point, meta$replicate), ]
  fm <- freq_matrix(matrix(runif(L * nrow(meta)), nrow = L,
                           dimnames = list(NULL, paste(meta$population,
                               meta$time_point, meta$replicate, sep = "_"))))
  out <- fit_af_glm(fm, sample_info(meta))
  expect_named(out, c("population", "estimate", "std.error", "statistic",
                      "p.value"))
  expect_identical(out$population, rep(c("AA", "BB"), times = L))
  expect_equal(nrow(out), 2 * L)
  # gold check: variant 1, population AA
  idx <- which(meta$population == "AA" & meta$time_point > 0)
  d <- data.frame(value = unclass(fm)[1, idx],
                  t = as.numeric(meta$time_point[idx]),
                  n = as.numeric(meta$sample_size[idx]))
  cf <- summary(glm(value ~ t, data = d, weights = n,
                    family = quasibinomial(link = "logit")))$coefficients["t", ]
  mine <- unlist(out[out$population == "AA", ][1, c("estimate", "std.error",
                                                    "statistic", "p.value")])
  expect_equal(unname(cf), unname(mine))
})

test_that("fit_af_glm recovers a known logistic slope", {
  mk <- data.frame(population = "AA", time_point = c(1, 2, 3, 4),
                   replicate = "R1", sample_size = 500L)
  slope_true <- 0.4
  p <- plogis(-0.2 + slope_true * mk$time_point)
  fmk <- freq_matrix(`colnames<-`(matrix(p, nrow = 1),
                                  paste(mk$population, mk$time_point,
                                        mk$replicate, sep = "_")))
  outk <- fit_af_glm(fmk, sample_info(mk))
  expect_true(abs(outk$estimate[1] - slope_true) < 0.01)
})

test_that("fit_af_glm adds the p0 baseline row with p0_weight", {
  set.seed(12)
  L <- 3
  meta <- data.frame(
    population = "AA", time_point = c(0, 1, 2, 3), replicate = "R1",
    sample_size = 30L
  )
  fm <- freq_matrix(matrix(runif(L * nrow(meta), 0.1, 0.9), nrow = L,
                           dimnames = list(NULL, paste(meta$population,
                               meta$time_point, meta$replicate, sep = "_"))))
  p0 <- p0_vec(runif(L, 0.2, 0.8))
  idx <- which(meta$time_point > 0)
  # without p0: fit on positive time points only
  out <- fit_af_glm(fm, sample_info(meta))
  for (l in seq_len(L)) {
    d <- data.frame(value = unclass(fm)[l, idx],
                    t = as.numeric(meta$time_point[idx]),
                    n = as.numeric(meta$sample_size[idx]))
    cf <- summary(glm(value ~ t, data = d, weights = n,
                      family = quasibinomial(link = "logit")))$coefficients["t", ]
    expect_equal(unname(cf), unname(unlist(out[l, c("estimate", "std.error",
                                                    "statistic", "p.value")])))
  }
  # with p0: one extra row t=0, weight 1000
  out0 <- fit_af_glm(fm, sample_info(meta), p0 = p0)
  for (l in seq_len(L)) {
    d0 <- rbind(data.frame(value = unclass(p0)[l], t = 0, n = 1000),
                data.frame(value = unclass(fm)[l, idx],
                           t = as.numeric(meta$time_point[idx]),
                           n = as.numeric(meta$sample_size[idx])))
    cf0 <- summary(glm(value ~ t, data = d0, weights = n,
                       family = quasibinomial(link = "logit")))$coefficients["t", ]
    expect_equal(unname(cf0), unname(unlist(out0[l, c("estimate", "std.error",
                                                      "statistic", "p.value")])))
  }
  # p0_weight = 50 respected
  outw <- fit_af_glm(fm, sample_info(meta), p0 = p0, p0_weight = 50)
  for (l in seq_len(L)) {
    dw <- rbind(data.frame(value = unclass(p0)[l], t = 0, n = 50),
                data.frame(value = unclass(fm)[l, idx],
                           t = as.numeric(meta$time_point[idx]),
                           n = as.numeric(meta$sample_size[idx])))
    cfw <- summary(glm(value ~ t, data = dw, weights = n,
                       family = quasibinomial(link = "logit")))$coefficients["t", ]
    expect_equal(unname(cfw), unname(unlist(outw[l, c("estimate", "std.error",
                                                      "statistic", "p.value")])))
  }
})

test_that("fit_af_glm fixed_intercept pins the logit(p0) offset", {
  set.seed(13)
  L <- 2
  meta <- data.frame(
    population = "AA", time_point = c(0, 1, 2, 3), replicate = "R1",
    sample_size = 30L
  )
  fm <- freq_matrix(matrix(runif(L * nrow(meta), 0.1, 0.9), nrow = L,
                           dimnames = list(NULL, paste(meta$population,
                               meta$time_point, meta$replicate, sep = "_"))))
  p0 <- p0_vec(runif(L, 0.2, 0.8))
  outf <- fit_af_glm(fm, sample_info(meta), p0 = p0, fixed_intercept = TRUE)
  idx <- which(meta$time_point > 0)
  for (l in seq_len(L)) {
    d0 <- rbind(data.frame(value = unclass(p0)[l], t = 0, n = 1000),
                data.frame(value = unclass(fm)[l, idx],
                           t = as.numeric(meta$time_point[idx]),
                           n = as.numeric(meta$sample_size[idx])))
    modf <- glm(value ~ t - 1 + offset(rep(qlogis(unclass(p0)[l]), nrow(d0))),
                data = d0, weights = n,
                family = quasibinomial(link = "logit"))
    cff <- summary(modf)$coefficients["t", ]
    expect_equal(unname(cff), unname(unlist(outf[l, c("estimate", "std.error",
                                                      "statistic", "p.value")])))
  }
  # requires p0
  expect_error(fit_af_glm(fm, sample_info(meta), fixed_intercept = TRUE),
               "requires p0")
})

test_that("fit_af_glm handles NA and boundary values", {
  set.seed(14)
  L <- 3
  meta <- rbind(expand.grid(population = "AA", time_point = c(1, 2, 3),
                            replicate = "R1"),
                expand.grid(population = "BB", time_point = c(1, 2, 3),
                            replicate = "R1"))
  meta$sample_size <- 30L
  meta <- meta[order(meta$population, meta$time_point, meta$replicate), ]
  fm <- freq_matrix(matrix(runif(L * nrow(meta), 0.1, 0.9), nrow = L,
                           dimnames = list(NULL, paste(meta$population,
                               meta$time_point, meta$replicate, sep = "_"))))
  # all-NA variant: NA statistics, correct population labels
  fm2 <- unclass(fm)
  fm2[1, ] <- NA
  out2 <- fit_af_glm(freq_matrix(fm2), sample_info(meta))
  expect_true(all(is.na(out2[1, c("estimate", "std.error", "statistic",
                                  "p.value")])))
  expect_true(all(is.na(out2[2, c("estimate", "std.error", "statistic",
                                  "p.value")])))
  expect_identical(out2$population[1:2], c("AA", "BB"))
  expect_true(!all(is.na(out2[3, c("estimate", "std.error", "statistic",
                                   "p.value")])))
  # boundary 0/1 values are kept in the fit
  idx <- which(meta$population == "AA" & meta$time_point > 0)
  fm3 <- unclass(fm)
  fm3[1, idx[1]] <- 1
  fm3[1, idx[2]] <- 0
  out3 <- fit_af_glm(freq_matrix(fm3), sample_info(meta))
  d3 <- data.frame(value = unclass(fm)[1, idx],
                   t = as.numeric(meta$time_point[idx]),
                   n = as.numeric(meta$sample_size[idx]))
  d3$value[1] <- 1
  d3$value[2] <- 0
  cf3 <- summary(glm(value ~ t, data = d3, weights = n,
                     family = quasibinomial(link = "logit")))$coefficients["t", ]
  mine3 <- unlist(out3[out3$population == "AA", ][1, c("estimate", "std.error",
                                                       "statistic", "p.value")])
  expect_equal(unname(cf3), unname(mine3))
})

test_that("fit_af_glm handles multiple populations and snp_coords", {
  set.seed(15)
  L <- 3
  meta <- rbind(expand.grid(population = "AA", time_point = c(1, 2),
                            replicate = c("R1", "R2")),
                expand.grid(population = "BB", time_point = c(1, 2, 3),
                            replicate = "R1"))
  meta$sample_size <- 30L
  meta <- meta[order(meta$population, meta$time_point, meta$replicate), ]
  fm <- freq_matrix(matrix(runif(L * nrow(meta), 0.1, 0.9), nrow = L,
                           dimnames = list(NULL, paste(meta$population,
                               meta$time_point, meta$replicate, sep = "_"))))
  coords <- snp_coords(data.frame(chrom = rep("2L", L),
                                  pos = seq_len(L) * 100))
  outc <- fit_af_glm(fm, sample_info(meta), snp_coords = coords)
  expect_named(outc, c("population", "CHROM", "POS", "estimate", "std.error",
                       "statistic", "p.value"))
  expect_equal(nrow(outc), 2 * L)
  expect_identical(as.character(outc$CHROM), rep("2L", 2 * L))
  expect_identical(outc$POS, rep(seq_len(L) * 100, times = 2))
})

test_that("fit_af_glm validates inputs", {
  set.seed(16)
  L <- 2
  meta <- data.frame(
    population = "AA", time_point = c(0, 1, 2), replicate = "R1",
    sample_size = 30L
  )
  fm_ok <- freq_matrix(matrix(runif(L * 3), nrow = L,
                              dimnames = list(NULL, c("AA_0_R1", "AA_1_R1",
                                                      "AA_2_R1"))))
  si_ok <- sample_info(meta)
  coords <- snp_coords(data.frame(chrom = c("1", "1"), pos = c(100, 200)))
  expect_error(fit_af_glm(matrix(0.5, 2, 3), si_ok),
               "freq_mat must be a freq_matrix object.")
  expect_error(fit_af_glm(fm_ok, meta), "sample_meta must be a sample_info object.")
  expect_error(fit_af_glm(fm_ok, si_ok, snp_coords = data.frame(chrom = "1",
                                                                pos = 1)),
               "snp_coords must be a snp_coords object")
  expect_error(fit_af_glm(fm_ok, si_ok, p0 = c(0.5, 0.5)),
               "p0 must be a p0_vec object or NULL.")
  expect_error(fit_af_glm(fm_ok, si_ok, p0 = p0_vec(c(0.5, 0.5, 0.5))),
               "p0_vec has 3 entries")
  expect_error(fit_af_glm(fm_ok, si_ok, p0_weight = -1),
               "p0_weight must be a single positive number.")
  expect_error(fit_af_glm(fm_ok, si_ok, p0_weight = c(10, 20)),
               "p0_weight must be a single positive number.")
  expect_error(fit_af_glm(fm_ok, si_ok, fixed_intercept = "yes"),
               "fixed_intercept must be a single TRUE or FALSE.")
  expect_error(fit_af_glm(fm_ok, si_ok, snp_coords = coords,
                          p0 = p0_vec(0.5)),
               "p0_vec has 1 entries")
  expect_error(fit_af_glm(fm_ok, sample_info(meta[1, , drop = FALSE])),
               "sample_info has 1 row")
  expect_error(fit_af_glm(fm_ok, si_ok,
                          snp_coords = snp_coords(coords[1, , drop = FALSE])),
               "snp_coords has 1 rows")
  # no positive time points
  meta0 <- data.frame(
    population = "AA", time_point = 0L, replicate = "R1", sample_size = 30L
  )
  fm_zero <- freq_matrix(`colnames<-`(matrix(0.3, nrow = 1), "AA_0_R1"))
  expect_error(fit_af_glm(fm_zero, p0 = NULL, sample_info(meta0)),
               "nothing to fit")
})

test_that("fit_af_glm integrates with the dataset kit", {
  d <- make_test_dataset(L = 12)
  expect_true(validate_af_dataset(d$freq_mat, d$coords, d$p0, d$meta))
  res <- fit_af_glm(d$freq_mat, d$meta)
  expect_equal(nrow(res), 2 * 12) # 2 populations x 12 variants
  expect_setequal(unique(res$population), c("AA", "BB"))
  expect_true(all(is.finite(res$estimate)))
  # with p0 baseline added, all rows still finite
  res0 <- fit_af_glm(d$freq_mat, d$meta, p0 = d$p0)
  expect_equal(nrow(res0), 2 * 12)
  expect_true(all(is.finite(res0$estimate)))
  res_c <- fit_af_glm(d$freq_mat, d$meta, snp_coords = d$coords)
  expect_named(res_c, c("population", "CHROM", "POS", "estimate", "std.error",
                        "statistic", "p.value"))
})
