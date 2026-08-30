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
