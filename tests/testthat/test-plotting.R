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

test_that("plot_predictions_truth_scatter returns a ggplot", {
  skip_if_not_installed("ggplot2")
  set.seed(1)
  truth <- rnorm(50)
  p <- plot_predictions_truth_scatter(
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

test_that("plot_predictions_truth_scatter drops NA/NaN/Inf pairs", {
  skip_if_not_installed("ggplot2")
  truth <- c(1, 2, NA, 4, NaN, 6, Inf, 8)
  predictions <- c(1.1, 2.2, 3.3, 4.4, 5.5, NaN, 7.7, 8.8)
  p <- plot_predictions_truth_scatter(
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

test_that("plot_manhattan point colors are customizable", {
  skip_if_not_installed("ggplot2")
  # default stays grey/black, alternating by chromosome
  p <- plot_manhattan(make_test_snp_table())
  expect_identical(
    sort(unique(ggplot2::ggplot_build(p)$data[[1]]$colour)),
    c("black", "grey")
  )
  chr_colors <- ggplot2::ggplot_build(p)$data[[1]]$colour
  expect_identical(chr_colors[1], "grey")
  expect_identical(chr_colors[6], "black")
  expect_identical(chr_colors[9], "grey")

  # two custom colors cycle across the three chromosomes (chrA/chrC = first
  # color, chrB = second)
  p2 <- plot_manhattan(make_test_snp_table(), point_colors = c("red", "blue"))
  build2 <- ggplot2::ggplot_build(p2)$data[[1]]$colour
  expect_identical(unique(build2[c(1, 9)]), "red")
  expect_identical(unique(build2[6]), "blue")
  scale_values <- ggplot2::ggplot_build(p2)$plot$scales$scales[[1]]$palette(
    nlevels(factor(p2$data$CHROM))
  )
  expect_identical(scale_values, c("red", "blue", "red"))

  # a single color makes all chromosomes uniform
  p3 <- plot_manhattan(make_test_snp_table(), point_colors = "green")
  expect_identical(
    unique(ggplot2::ggplot_build(p3)$data[[1]]$colour),
    "green"
  )
})

test_that("plot_manhattan highlight color is customizable and points use shape 16", {
  skip_if_not_installed("ggplot2")
  snp_h <- make_test_snp_table()
  snp_h$highlight <- c(rep("no", 9), "yes")
  p_h <- plot_manhattan(snp_h, highlight_color = "purple")
  expect_identical(p_h$layers[[2]]$aes_params$colour, "purple")

  # both layers draw solid filled circles (shape 16, no border stroke)
  expect_identical(p_h$layers[[1]]$aes_params$shape, 16)
  expect_identical(p_h$layers[[2]]$aes_params$shape, 16)
})

test_that("plot_manhattan errors on invalid color arguments", {
  skip_if_not_installed("ggplot2")
  expect_error(
    plot_manhattan(make_test_snp_table(), point_colors = character(0)),
    "point_colors must be a non-empty character vector"
  )
  expect_error(
    plot_manhattan(make_test_snp_table(), point_colors = 1),
    "point_colors must be a non-empty character vector"
  )
  expect_error(
    plot_manhattan(make_test_snp_table(), highlight_color = c("red", "blue")),
    "highlight_color must be a single color name"
  )
  expect_error(
    plot_manhattan(make_test_snp_table(), highlight_color = 1),
    "highlight_color must be a single color name"
  )
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

make_test_covmat <- function() {
  matrix(
    c(1.0, 0.3, 0.1,
      0.3, 1.0, 0.5,
      0.1, 0.5, 1.0),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("s1", "s2", "s3"), c("s1", "s2", "s3"))
  )
}

test_that("get_upper_tri blanks the strict lower triangle", {
  upper <- get_upper_tri(make_test_covmat())
  expect_true(all(is.na(upper[lower.tri(upper)])))
  expect_identical(unname(upper[1, ]), c(1, 0.3, 0.1))
})

test_that("plot_var_cov_matrix builds tile/heatmap layers from the upper triangle", {
  skip_if_not_installed("ggplot2")
  p <- plot_var_cov_matrix(make_test_covmat())
  expect_s3_class(p, "ggplot")

  # upper triangle only: 3 + 2 + 1 cells, diagonal included
  expect_identical(nrow(p$data), 6L)
  # base layer uses tiles with a white border
  expect_s3_class(p$layers[[1]]$geom, "GeomTile")
  expect_identical(p$layers[[1]]$aes_params$colour, "white")
  # values are printed on the tiles
  expect_identical(length(p$layers), 2L)

  # diagonal (value 1) is included; strictly-lower cells are not
  expect_false(any(p$data$Var1 == "s3" & p$data$Var2 == "s1"))
  expect_true(any(p$data$Var1 == "s1" & p$data$Var2 == "s3"))
})

test_that("plot_var_cov_matrix switches variants correctly", {
  skip_if_not_installed("ggplot2")
  covmat <- make_test_covmat()

  # raster: base geom changes, layer count unchanged
  p_raster <- plot_var_cov_matrix(covmat, tile_or_raster = "raster")
  expect_s3_class(p_raster$layers[[1]]$geom, "GeomRaster")

  # no values: single layer
  p_noval <- plot_var_cov_matrix(covmat, include_values = FALSE)
  expect_identical(length(p_noval$layers), 1L)

  # continuous mode: gradient2 fill scale
  p_cont <- plot_var_cov_matrix(covmat)
  expect_s3_class(p_cont$scales$scales[[1]], "ScaleContinuous")

  # tick-labels-off branch applies the blank theme
  p_bare <- plot_var_cov_matrix(covmat, include_tick_labels = FALSE)
  expect_s3_class(p_bare$theme$axis.text.x, "element_blank")

  # invalid enum values are rejected
  expect_error(
    plot_var_cov_matrix(covmat, tile_or_raster = "dots"),
    "'arg' should be one of"
  )
  expect_error(
    plot_var_cov_matrix(covmat, discrete_or_continuous = "spectral"),
    "'arg' should be one of"
  )

  # non-square or asymmetrically-named matrices are rejected
  expect_error(plot_var_cov_matrix(covmat[1:2, ]), "square")
  bad_names <- covmat
  rownames(bad_names) <- c("a", "b", "c")
  expect_error(plot_var_cov_matrix(bad_names), "identical row and column names")
})

# ---------------------------------------------------------------------------
# plot_effect_sizes
# ---------------------------------------------------------------------------

test_that("plot_effect_sizes builds an effect-size plot with intercept dropped", {
  ci <- data.frame(
    Term = c("(Intercept)", "cds_density", "mean_r2_direct"),
    Estimate = c(0.5, -1.2, 0.8),
    Lower_CI = c(0.1, -1.8, 0.2),
    Upper_CI = c(0.9, -0.6, 1.4)
  )
  p <- plot_effect_sizes(ci)
  expect_s3_class(p, "ggplot")
  expect_false(any(grepl("[Ii]ntercept", p$data$Term)))
  expect_equal(nrow(p$data), 2)
  expect_equal(p$labels$y, "Block bootstrap estimate of\nscaled coefficient")
  expect_equal(p$labels$x, "Model term")
  expect_null(p$labels$title)
})

test_that("plot_effect_sizes keeps the intercept when asked and accepts a title", {
  ci <- data.frame(
    Term = c("Intercept", "x"),
    Estimate = c(1, 2), Lower_CI = c(0.5, 1), Upper_CI = c(1.5, 3)
  )
  p <- plot_effect_sizes(ci, drop_intercept = FALSE, title = "model")
  expect_equal(nrow(p$data), 2)
  expect_equal(p$labels$title, "model")
  # styling arguments reach the built layers
  build <- ggplot2::ggplot_build(p)
  expect_equal(unique(build$data[[1]]$colour), "red")
  expect_equal(unique(build$data[[1]]$linewidth), 3)
  expect_equal(unique(build$data[[2]]$size), 3)
  expect_equal(unique(build$data[[2]]$colour), "black")
})

test_that("plot_effect_sizes validates its input", {
  ci <- data.frame(
    Term = "x", Estimate = 2, Lower_CI = 1, Upper_CI = 3
  )
  expect_error(plot_effect_sizes(ci[, 1:3]), "missing required columns")
  expect_error(
    plot_effect_sizes(
      data.frame(Term = "(Intercept)", Estimate = 1, Lower_CI = 0, Upper_CI = 2)
    ),
    "no rows remain"
  )
})

# ---------------------------------------------------------------------------
# plot_partial_regression
# ---------------------------------------------------------------------------

test_that("plot_partial_regression plots residualised variables", {
  skip_if_not_installed("ggplot2")
  set.seed(42)
  n <- 100
  dat <- data.frame(
    y = rnorm(n), x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n)
  )
  dat$y <- 1 + 2 * dat$x1 - dat$x2 + 0.5 * dat$x3 + rnorm(n)
  mod <- stats::lm(y ~ x1 + x2 + x3, data = dat)

  p <- plot_partial_regression(mod, "x1")
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$x, "Primary variable")
  expect_equal(p$labels$y, "Residuals")

  # coordinates are the residuals from regressing each variable on the others
  expect_equal(
    p$data$x_part,
    unname(stats::residuals(stats::lm(x1 ~ x2 + x3, data = dat)))
  )
  expect_equal(
    p$data$y_part,
    unname(stats::residuals(stats::lm(y ~ x2 + x3, data = dat)))
  )

  # the drawn line is the partial regression line through the origin
  abline <- ggplot2::ggplot_build(p)$data[[2]]
  expect_equal(abline$slope, unname(stats::coef(mod)[["x1"]]))
  expect_equal(abline$intercept, 0)
})

test_that("plot_partial_regression handles factors, weights, and offsets", {
  skip_if_not_installed("ggplot2")
  set.seed(7)
  n <- 80
  dat <- data.frame(
    x1 = rnorm(n), x2 = rnorm(n),
    f = factor(sample(c("a", "b", "c"), n, replace = TRUE)),
    w = runif(n, 0.5, 2)
  )
  dat$y <- 1 + 2 * dat$x1 + as.numeric(dat$f) + rnorm(n)
  dat$z <- as.numeric(dat$f == "c")

  # factor covariate: residualised on the factor's dummy columns
  mod_f <- stats::lm(y ~ x1 + f, data = dat)
  p_f <- plot_partial_regression(mod_f, "x1")
  expect_equal(
    p_f$data$x_part,
    unname(stats::residuals(stats::lm(x1 ~ f, data = dat)))
  )

  # weights: plotted residuals are sqrt(weight)-scaled
  mod_w <- stats::lm(y ~ x1 + x2, data = dat, weights = w)
  p_w <- plot_partial_regression(mod_w, "x1")
  expect_equal(
    p_w$data$x_part,
    unname(stats::weighted.residuals(
      stats::lm(x1 ~ x2, data = dat, weights = w)
    ))
  )
  expect_equal(
    ggplot2::ggplot_build(p_w)$data[[2]]$slope,
    unname(stats::coef(mod_w)[["x1"]])
  )

  # offset: the response is adjusted before residualising
  mod_o <- stats::lm(y ~ x1 + x2 + offset(z), data = dat)
  p_o <- plot_partial_regression(mod_o, "x1")
  expect_equal(
    p_o$data$y_part,
    unname(stats::residuals(stats::lm(I(y - z) ~ x2, data = dat)))
  )

  # no intercept: there are no other columns to residualise on
  mod_0 <- stats::lm(y ~ 0 + x1, data = dat)
  p_0 <- plot_partial_regression(mod_0, "x1")
  expect_equal(p_0$data$x_part, dat$x1)
  expect_equal(
    ggplot2::ggplot_build(p_0)$data[[2]]$slope,
    unname(stats::coef(mod_0)[["x1"]])
  )

  # rows dropped by na.action stay dropped
  dat$y[c(1, 5)] <- NA
  mod_na <- stats::lm(y ~ x1 + x2, data = dat, na.action = na.exclude)
  p_na <- plot_partial_regression(mod_na, "x1")
  expect_equal(nrow(p_na$data), nrow(mod_na$model))
})

test_that("plot_partial_regression passes styling arguments to the lines", {
  skip_if_not_installed("ggplot2")
  set.seed(3)
  dat <- data.frame(y = rnorm(40), x1 = rnorm(40), x2 = rnorm(40))
  mod <- stats::lm(y ~ x1 + x2, data = dat)

  p <- plot_partial_regression(mod, "x1", colour = "red", linewidth = 2)
  build <- ggplot2::ggplot_build(p)$data
  expect_equal(unique(build[[2]]$colour), "red")
  expect_equal(unique(build[[2]]$linewidth), 2)
  # limit lines inherit the styling and stay dashed
  expect_equal(unique(build[[3]]$colour), "red")
  expect_equal(unique(build[[3]]$linewidth), 2)
  expect_equal(unique(build[[3]]$linetype), "dashed")

  # an explicit linetype overrides the dashed limits
  p_dot <- plot_partial_regression(mod, "x1", linetype = "dotted")
  expect_equal(
    unique(ggplot2::ggplot_build(p_dot)$data[[3]]$linetype),
    "dotted"
  )
})

test_that("plot_partial_regression adds t-based confidence limits", {
  skip_if_not_installed("ggplot2")
  set.seed(5)
  dat <- data.frame(y = rnorm(50), x1 = rnorm(50), x2 = rnorm(50))
  mod <- stats::lm(y ~ x1 + x2, data = dat)

  p <- plot_partial_regression(mod, "x1")
  expect_equal(length(p$layers), 4L)

  slopes <- vapply(
    ggplot2::ggplot_build(p)$data[2:4],
    function(d) unique(d$slope), numeric(1)
  )
  expect_equal(
    sort(slopes),
    sort(unname(c(
      stats::coef(mod)[["x1"]],
      stats::confint(mod, parm = "x1")[1, ]
    )))
  )

  # level is adjustable and uses the t distribution, not a fixed 1.96
  p90 <- plot_partial_regression(mod, "x1", level = 0.9)
  slopes90 <- vapply(
    ggplot2::ggplot_build(p90)$data[2:4],
    function(d) unique(d$slope), numeric(1)
  )
  expect_equal(
    sort(slopes90),
    sort(unname(c(
      stats::coef(mod)[["x1"]],
      stats::confint(mod, parm = "x1", level = 0.9)[1, ]
    )))
  )

  # can be switched off
  expect_equal(length(plot_partial_regression(mod, "x1", ci = FALSE)$layers), 2L)

  expect_error(
    plot_partial_regression(mod, "x1", level = 1.5),
    "between 0 and 1"
  )
})

test_that("plot_partial_regression draws supplied interval bounds", {
  skip_if_not_installed("ggplot2")
  set.seed(11)
  dat <- data.frame(y = rnorm(60), x1 = rnorm(60), x2 = rnorm(60))
  mod <- stats::lm(y ~ x1 + x2, data = dat)

  ci_table <- data.frame(
    Term = c("(Intercept)", "x1", "x2"),
    Estimate = unname(stats::coef(mod)),
    Lower_CI = c(-1, -2, -1),
    Upper_CI = c(1, 0.5, 1)
  )
  p <- plot_partial_regression(mod, "x1", ci = ci_table)
  slopes <- vapply(
    ggplot2::ggplot_build(p)$data[2:4],
    function(d) unique(d$slope), numeric(1)
  )
  expect_equal(
    sort(slopes),
    sort(c(stats::coef(mod)[["x1"]], -2, 0.5))
  )

  # a single-row table needs no Term column
  p1 <- plot_partial_regression(
    mod, "x1",
    ci = data.frame(Lower_CI = -2, Upper_CI = 0.5)
  )
  expect_equal(length(p1$layers), 4L)

  # a table built for a different fit warns but still plots
  expect_warning(
    plot_partial_regression(
      mod, "x1",
      ci = data.frame(
        Term = "x1", Estimate = 99, Lower_CI = -2, Upper_CI = 0.5
      )
    ),
    "does not match the model"
  )

  expect_error(
    plot_partial_regression(
      mod, "x1",
      ci = data.frame(Term = "z", Lower_CI = -1, Upper_CI = 1)
    ),
    "no row for coefficient"
  )
  expect_error(
    plot_partial_regression(mod, "x1", ci = data.frame(Term = "x1")),
    "missing column"
  )
  expect_error(
    plot_partial_regression(
      mod, "x1", ci = data.frame(Lower_CI = 1, Upper_CI = -1)
    ),
    "finite and ordered"
  )
  expect_error(
    plot_partial_regression(mod, "x1", ci = "boot"),
    "must be TRUE, FALSE, or a data frame"
  )
})

test_that("plot_partial_regression accepts block_bootstrap_lm intervals", {
  skip_if_not_installed("ggplot2")
  set.seed(21)
  windows <- data.frame(y = rnorm(60), x1 = rnorm(60), x2 = rnorm(60))
  mod <- stats::lm(y ~ x1 + x2, data = windows)

  ci_table <- block_bootstrap_lm(
    windows, y ~ x1 + x2,
    R = 100, block_length = 4, seed = 3
  )
  p <- plot_partial_regression(mod, "x1", ci = ci_table)

  slopes <- vapply(
    ggplot2::ggplot_build(p)$data[2:4],
    function(d) unique(d$slope), numeric(1)
  )
  row <- ci_table[ci_table$Term == "x1", ]
  expect_equal(
    sort(slopes),
    sort(c(stats::coef(mod)[["x1"]], row$Lower_CI, row$Upper_CI))
  )
})

test_that("plot_partial_regression validates its input", {
  skip_if_not_installed("ggplot2")
  set.seed(1)
  dat <- data.frame(
    y = rnorm(30), x1 = rnorm(30),
    f = factor(rep(c("a", "b", "c"), 10))
  )
  mod <- stats::lm(y ~ x1 + f, data = dat)

  expect_error(plot_partial_regression(list(), "x1"), "must be a linear model")
  expect_error(plot_partial_regression(mod, c("x1", "f")), "single term name")
  expect_error(plot_partial_regression(mod, "nope"), "is not a term")
  expect_error(plot_partial_regression(mod, "f"), "model-matrix columns")

  mod_mv <- stats::lm(cbind(y, x1) ~ x1, data = dat)
  expect_error(plot_partial_regression(mod_mv, "x1"), "Multivariate")
})
