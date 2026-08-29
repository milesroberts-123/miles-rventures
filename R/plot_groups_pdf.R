#' Plot each group of a data frame to a multi-page PDF
#'
#' Splits `data` by the columns in `groupings`, then writes one plot per
#' subset to a single multi-page PDF. The output directory is created if it
#' does not exist.
#'
#' @param data A data frame to split and plot.
#' @param groupings Character vector of column names in `data` to group by.
#' @param plot_fn A function that takes a subset data frame and returns a
#'   ggplot object.
#' @param filepath Path to the output PDF file.
#' @param width PDF page width in inches.
#' @param height PDF page height in inches.
#'
#' @return Invisibly, `filepath`.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_groups_pdf(
#'   data = iris,
#'   groupings = "Species",
#'   plot_fn = function(sub_data) {
#'     ggplot2::ggplot(sub_data, ggplot2::aes(x = Sepal.Length, y = Sepal.Width)) +
#'       ggplot2::geom_point() +
#'       ggplot2::labs(title = unique(sub_data$Species))
#'   },
#'   filepath = "iris_by_species.pdf"
#' )
#' }
plot_groups_pdf <- function(data, groupings, plot_fn, filepath, width = 8, height = 6) {
  dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)

  data_list <- split(data, data[groupings], drop = TRUE)

  grDevices::pdf(filepath, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (sub_data in data_list) {
    print(plot_fn(sub_data))
  }

  invisible(filepath)
}
