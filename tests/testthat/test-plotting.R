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
