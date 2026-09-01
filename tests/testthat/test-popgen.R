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
