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

test_that("read_paf parses a PAF file correctly", {
  # Real minimap2 output: every line carries a consistent tag block
  paf_lines <- paste(
    "q1\t1000\t100\t400\t+\tt1\t2000\t200\t500\t290\t300\t60\tcg:Z:290M\ttp:A:P",
    "q1\t1000\t500\t900\t-\tt1\t2000\t600\t1000\t390\t400\t30\tcg:Z:390M\ttp:A:P",
    "q2\t1000\t0\t100\t+\tt2\t2000\t0\t100\t50\t100\t0\tcg:Z:50M\tSB:i:10",
    sep = "\n"
  )
  f <- tempfile(fileext = ".paf")
  writeLines(paf_lines, f)
  on.exit(unlink(f))

  # ragged extra columns must parse without warnings
  expect_no_warning(paf <- read_paf(f))

  expect_named(
    paf,
    c("qname", "qlen", "qstart", "qend", "strand", "tname", "tlen",
      "tstart", "tend", "nmatch", "alen", "mapq", "pident")
  )
  expect_equal(nrow(paf), 3)
  expect_type(paf$qname, "character")
  expect_type(paf$strand, "character")
  expect_equal(paf$qlen, c(1000, 1000, 1000))
  expect_true(is.numeric(paf$qstart))
  expect_true(is.numeric(paf$mapq))
  expect_equal(paf$pident, c(290 / 300, 390 / 400, 50 / 100))
})

test_that("bind_files binds tables with custom column names and source ids", {
  f1 <- tempfile(fileext = ".txt")
  f2 <- tempfile(fileext = ".txt")
  on.exit(unlink(c(f1, f2)))
  writeLines("0.1\t0.9", f1)
  writeLines("0.2\t0.8\n0.3\t0.7", f2)

  bound <- bind_files(c(f1, f2), col_names = c("bray_curtis", "cosine"))

  expect_equal(nrow(bound), 3)
  expect_named(bound, c("source_file", "bray_curtis", "cosine"))
  expect_equal(bound$source_file, c(f1, f2, f2))
  expect_equal(bound$bray_curtis, c(0.1, 0.2, 0.3))
  expect_s3_class(bound, "data.frame")
})

test_that("bind_files keeps file headers when col_names is NULL", {
  f1 <- tempfile(fileext = ".csv")
  f2 <- tempfile(fileext = ".csv")
  on.exit(unlink(c(f1, f2)))
  writeLines("a,b\n1,2", f1)
  writeLines("a,b\n3,4", f2)

  bound <- bind_files(c(f1, f2))

  expect_named(bound, c("source_file", "a", "b"))
  expect_equal(bound$a, c(1, 3))
})

test_that("bind_files omits the id column when id is NULL", {
  f1 <- tempfile(fileext = ".txt")
  on.exit(unlink(f1))
  writeLines("1\t2", f1)

  bound <- bind_files(f1, col_names = c("x", "y"), id = NULL)

  expect_named(bound, c("x", "y"))
})

test_that("bind_files errors on missing files", {
  good <- tempfile(fileext = ".txt")
  on.exit(unlink(good))
  writeLines("1\t2", good)

  expect_error(
    bind_files(c(good, "/nonexistent/path.txt")),
    "File\\(s\\) not found: /nonexistent/path\\.txt"
  )
})

test_that("bind_files autodetects delimiters via fread", {
  f1 <- tempfile(fileext = ".csv")
  on.exit(unlink(f1))
  writeLines("1;2", f1)

  bound <- bind_files(f1, col_names = c("x", "y"))

  expect_equal(bound$x, 1)
  expect_equal(bound$y, 2)
})

test_that("bind_files validates input type", {
  expect_error(bind_files(list("a.txt")), "is.character\\(file_list\\) is not TRUE")
})
