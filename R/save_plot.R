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
