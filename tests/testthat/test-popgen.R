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
