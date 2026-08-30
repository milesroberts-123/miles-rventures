test_that("check_file_exists deletes an existing file", {
  f <- tempfile()
  writeLines("x", f)
  expect_true(file.exists(f))
  expect_true(check_file_exists(f))
  expect_false(file.exists(f))
})

test_that("check_file_exists returns FALSE for a missing file", {
  f <- tempfile()
  expect_false(check_file_exists(f))
})

test_that("append_table writes header on first call and appends after", {
  f <- tempfile(fileext = ".csv")
  append_table(data.frame(a = 1, b = 2), f)
  append_table(data.frame(a = 3, b = 4), f)
  lines <- readLines(f)
  expect_equal(lines[1], "a,b")
  expect_length(lines, 3)
})
