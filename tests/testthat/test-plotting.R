test_that("save_plot writes PDF, PNG, and symlinks", {
  skip_if_not_installed("ggplot2")
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = wt, y = mpg)) +
    ggplot2::geom_point()

  base <- tempfile()
  dir.create(file.path(base, "by-date", "2026-01-01"), recursive = TRUE)
  dir.create(file.path(base, "by-analysis"), recursive = TRUE)

  save_plot(p, "test_plot", "2026-01-01", base, height = 4, width = 4)

  pdf_path <- file.path(base, "by-date", "2026-01-01", "test_plot.pdf")
  png_path <- file.path(base, "by-date", "2026-01-01", "test_plot.png")
  expect_true(file.exists(pdf_path))
  expect_true(file.exists(png_path))

  link_pdf <- file.path(base, "by-analysis", "test_plot.pdf")
  link_png <- file.path(base, "by-analysis", "test_plot.png")
  expect_true(file.exists(link_pdf))
  expect_true(file.exists(link_png))
})

test_that("plot_groups_pdf writes a PDF", {
  skip_if_not_installed("ggplot2")
  f <- tempfile(fileext = ".pdf")
  plot_groups_pdf(
    data = iris,
    groupings = "Species",
    plot_fn = function(sub_data) {
      ggplot2::ggplot(sub_data, ggplot2::aes(x = Sepal.Length, y = Sepal.Width)) +
        ggplot2::geom_point()
    },
    filepath = f
  )
  expect_true(file.exists(f))
  expect_gt(file.info(f)$size, 0)
})

test_that("compare_predictions_truth_scatter returns a ggplot", {
  skip_if_not_installed("ggplot2")
  set.seed(1)
  truth <- rnorm(50)
  p <- compare_predictions_truth_scatter(
    predictions = truth + rnorm(50, sd = 0.5),
    truth = truth,
    model_label = "Test model",
    truth_label = "Truth",
    pred_label = "Prediction",
    color_var = truth,
    color_label = "True value"
  )
  expect_s3_class(p, "ggplot")
})

test_that("compare_predictions_truth_scatter drops NA/NaN/Inf pairs", {
  skip_if_not_installed("ggplot2")
  truth <- c(1, 2, NA, 4, NaN, 6, Inf, 8)
  predictions <- c(1.1, 2.2, 3.3, 4.4, 5.5, NaN, 7.7, 8.8)
  p <- compare_predictions_truth_scatter(
    predictions = predictions,
    truth = truth,
    model_label = "Test model",
    truth_label = "Truth",
    pred_label = "Prediction",
    color_var = seq_along(truth),
    color_label = "Index"
  )
  expect_s3_class(p, "ggplot")
  expect_identical(nrow(p$data), 4L)
})

make_test_paf <- function() {
  data.frame(
    qname = c("q1", "q1", "q2"),
    qlen = c(10000, 10000, 10000),
    qstart = c(100, 500, 100),
    qend = c(4000, 9000, 4000),
    strand = c("+", "-", "+"),
    tname = c("t1", "t1", "t2"),
    tlen = c(20000, 20000, 20000),
    tstart = c(200, 600, 200),
    tend = c(5000, 10000, 5000),
    nmatch = c(3900, 450, 3900),
    alen = c(3900, 500, 3900),
    mapq = c(60, 60, 0)
  )
}

test_that("plot_paf_dotplot returns a ggplot for matching query and target", {
  skip_if_not_installed("ggplot2")
  p <- plot_paf_dotplot(
    make_test_paf(),
    q_seq = "q1", t_seq = "t1",
    xlab = "Query", ylab = "Target", title = "Test"
  )
  expect_s3_class(p, "ggplot")
  # only the long q1-t1 block survives the default min_alen filter
  expect_identical(nrow(p$data), 1L)
})

test_that("plot_paf_dotplot orients reverse-strand alignments", {
  skip_if_not_installed("ggplot2")
  paf <- data.frame(
    qname = "q1", qlen = 10000, qstart = 100, qend = 4000, strand = "-",
    tname = "t1", tlen = 20000, tstart = 200, tend = 5000,
    nmatch = 3900, alen = 3900, mapq = 60
  )
  p <- plot_paf_dotplot(
    paf, q_seq = "q1", t_seq = "t1",
    xlab = "Query", ylab = "Target", title = "Test"
  )
  # reverse strand: original tstart (200) becomes the segment's yend
  expect_identical(p$data$tstart_plot[1], 5000)
  expect_identical(p$data$tend_plot[1], 200)
})

test_that("plot_paf_dotplot applies min_alen and min_mapq filters", {
  skip_if_not_installed("ggplot2")
  paf <- make_test_paf()
  # lower thresholds keep both q1-t1 blocks
  p <- plot_paf_dotplot(
    paf, q_seq = "q1", t_seq = "t1",
    xlab = "Query", ylab = "Target", title = "Test",
    min_alen = 100, min_mapq = 1
  )
  expect_identical(nrow(p$data), 2L)

  # default min_alen = 1000 drops the short q1-t1 block
  p_strict <- plot_paf_dotplot(
    paf, q_seq = "q1", t_seq = "t1",
    xlab = "Query", ylab = "Target", title = "Test",
    min_alen = 1000
  )
  expect_identical(nrow(p_strict$data), 1L)

  # asking for a pair with no alignments anywhere in the file errors
  expect_error(
    plot_paf_dotplot(
      paf, q_seq = "q2", t_seq = "t1",
      xlab = "Query", ylab = "Target", title = "Test",
      min_alen = 100
    ),
    "No alignments found"
  )
})

make_test_snp_table <- function() {
  data.frame(
    CHROM = rep(c("chrA", "chrB", "chrC"), c(5, 3, 2)),
    POS = c(100, 250, 400, 400, 500, 220, 550, 700, 400, 900),
    pvalue = c(0.9, 0.5, 1e-6, 0.2, 0.8, 0.01, 0.001, 0.3, 3e-8, 0.7)
  )
}

test_that("plot_manhattan computes cumulative positions across chromosomes", {
  skip_if_not_installed("ggplot2")
  p <- plot_manhattan(make_test_snp_table())
  expect_s3_class(p, "ggplot")

  # chrB offsets by chrA's max POS (500); chrC by chrA + chrB (500 + 700 = 1200)
  expect_identical(
    p$data$BPcum,
    c(100, 250, 400, 400, 500, 720, 1050, 1200, 1600, 2100)
  )
  # x-axis breaks are one per chromosome, at chromosome centers
  expect_identical(
    ggplot2::ggplot_build(p)$layout$panel_scales_x[[1]]$breaks,
    c(300, 960, 1850)
  )
})

test_that("plot_manhattan adds an orange highlight layer only when the column exists", {
  skip_if_not_installed("ggplot2")
  p <- plot_manhattan(make_test_snp_table())
  expect_identical(length(p$layers), 1L)

  snp_h <- make_test_snp_table()
  snp_h$highlight <- c(rep("no", 9), "yes")
  p_h <- plot_manhattan(snp_h)
  expect_identical(length(p_h$layers), 2L)
  expect_identical(p_h$layers[[2]]$aes_params$colour, "orange")
})

test_that("plot_manhattan errors on empty data or missing columns", {
  skip_if_not_installed("ggplot2")
  expect_error(
    plot_manhattan(make_test_snp_table()[0, ]),
    "No SNPs in 'data'"
  )
  expect_error(
    plot_manhattan(data.frame(CHROM = "chrA", POS = 1)),
    "Missing required column"
  )
})
