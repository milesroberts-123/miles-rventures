#' Save a plot as PDF and PNG, with symlinks in a by-analysis directory
#'
#' Saves `plot_obj` to `plot_base_dir/by-date/<date>/<plot_name>` as both PDF
#' and PNG, then creates (or refreshes) symlinks in
#' `plot_base_dir/by-analysis/` pointing at the saved files.
#'
#' @param plot_obj A ggplot object to save.
#' @param plot_name Base file name for the plot, without extension.
#' @param date Date string used for the `by-date/<date>/` subdirectory.
#' @param plot_base_dir Base directory containing `by-date/` and
#'   `by-analysis/` subdirectories. Both must already exist.
#' @param height Plot height in inches.
#' @param width Plot width in inches.
#'
#' @return Invisibly, the path to the saved PDF.
#'
#' @export
#' @importFrom ggplot2 ggsave
save_plot <- function(plot_obj, plot_name, date, plot_base_dir, height, width) {
  plot_full_name <- file.path(plot_base_dir, "by-date", date, plot_name)

  ggsave(
    paste0(plot_full_name, ".pdf"),
    plot = plot_obj,
    height = height,
    width = width
  )
  ggsave(
    paste0(plot_full_name, ".png"),
    plot = plot_obj,
    height = height,
    width = width
  )

  symlink_pdf <- file.path(plot_base_dir, "by-analysis", paste0(plot_name, ".pdf"))
  symlink_png <- file.path(plot_base_dir, "by-analysis", paste0(plot_name, ".png"))

  if (file.exists(symlink_pdf)) {
    file.remove(symlink_pdf)
  }

  if (file.exists(symlink_png)) {
    file.remove(symlink_png)
  }

  file.symlink(from = paste0(plot_full_name, ".pdf"), to = symlink_pdf)
  file.symlink(from = paste0(plot_full_name, ".png"), to = symlink_png)

  invisible(paste0(plot_full_name, ".pdf"))
}

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

  grDevices::cairo_pdf(filepath, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (sub_data in data_list) {
    print(plot_fn(sub_data))
  }

  invisible(filepath)
}
