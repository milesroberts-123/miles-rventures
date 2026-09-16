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

test_that("ensure_parent_dir creates nested folders and is idempotent", {
  base <- tempfile()
  target <- file.path(base, "nested", "data_results", "output.csv")
  on.exit(unlink(base, recursive = TRUE))

  expect_false(dir.exists(dirname(target)))
  first <- ensure_parent_dir(target)
  expect_true(dir.exists(dirname(target)))
  expect_equal(class(first), class(TRUE))  # logical, invisible
  expect_true(first)

  second <- ensure_parent_dir(target)
  expect_false(second)
  expect_true(dir.exists(dirname(target)))
})

test_that("ensure_parent_dir handles plain file names in the working dir", {
  # dirname("f.csv") == "." which already exists
  expect_false(ensure_parent_dir("f.csv"))
})

test_that("ensure_parent_dir validates input", {
  expect_error(ensure_parent_dir(123), "is.character\\(file_path\\) is not TRUE")
  expect_error(ensure_parent_dir(c("a.csv", "b.csv")),
               "length\\(file_path\\) == 1 is not TRUE")
})

test_that("csv_to_parquet converts a CSV and returns a lazy dataset", {
  skip_if_not_installed("arrow")
  csv <- tempfile(fileext = ".csv")
  pdir <- tempfile()
  on.exit(unlink(c(csv, pdir), recursive = TRUE))
  writeLines("x,y\n1,2\n3,4\n5,6", csv)

  ds <- csv_to_parquet(csv, pdir, max_rows_per_file = 2L)

  expect_s3_class(ds, "Dataset")
  expect_equal(dplyr::collect(ds), dplyr::tibble(x = c(1, 3, 5), y = c(2, 4, 6)))
  # rows split across files at 2 rows each
  expect_length(list.files(pdir), 2)
})

test_that("csv_to_parquet skips conversion when the parquet dir exists", {
  skip_if_not_installed("arrow")
  csv <- tempfile(fileext = ".csv")
  pdir <- tempfile()
  on.exit(unlink(c(csv, pdir), recursive = TRUE))
  writeLines("x,y\n1,2", csv)
  csv_to_parquet(csv, pdir)

  files_before <- list.files(pdir, full.names = TRUE)
  mtimes_before <- file.info(files_before)$mtime
  Sys.sleep(0.01)

  ds <- csv_to_parquet(csv, pdir)
  expect_s3_class(ds, "Dataset")
  expect_equal(dplyr::collect(ds), dplyr::tibble(x = 1, y = 2))
  # parquet files untouched by the second call
  expect_identical(
    file.info(list.files(pdir, full.names = TRUE))$mtime,
    mtimes_before
  )
})

test_that("csv_to_parquet errors on a missing CSV", {
  skip_if_not_installed("arrow")
  expect_error(
    csv_to_parquet("/nonexistent/data.csv", tempfile()),
    "File not found: /nonexistent/data\\.csv"
  )
})

test_that("csv_to_parquet validates tuning arguments", {
  skip_if_not_installed("arrow")
  csv <- tempfile(fileext = ".csv")
  pdir <- tempfile()
  on.exit(unlink(c(csv, pdir), recursive = TRUE))
  writeLines("x\n1", csv)

  expect_error(csv_to_parquet(csv, pdir, block_size = 0),
               "block_size must be a single positive number\\.")
  expect_error(csv_to_parquet(csv, pdir, block_size = c(1024, 2048)),
               "block_size must be a single positive number\\.")
  expect_error(csv_to_parquet(csv, pdir, max_rows_per_file = -1),
               "max_rows_per_file must be a single positive number\\.")
  expect_error(csv_to_parquet(csv, pdir, max_rows_per_file = NA_real_),
               "max_rows_per_file must be a single positive number\\.")
})
