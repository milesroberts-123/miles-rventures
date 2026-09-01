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
save_plot <- function(plot_obj, plot_name, date, plot_base_dir, height, width) {
  plot_full_name <- file.path(plot_base_dir, "by-date", date, plot_name)
  pdf_path <- paste0(plot_full_name, ".pdf")
  png_path <- paste0(plot_full_name, ".png")

  for (path in c(pdf_path, png_path)) {
    ggplot2::ggsave(
      path,
      plot = plot_obj,
      height = height,
      width = width
    )
  }

  symlink_pdf <- file.path(plot_base_dir, "by-analysis", paste0(plot_name, ".pdf"))
  symlink_png <- file.path(plot_base_dir, "by-analysis", paste0(plot_name, ".png"))

  unlink(c(symlink_pdf, symlink_png))
  file.symlink(from = c(pdf_path, png_path), to = c(symlink_pdf, symlink_png))

  invisible(pdf_path)
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

#' Scatterplot comparing predictions to true values
#'
#' Plots predicted values against true values for a continuous variable, with
#' the identity line (predictions = truth), a linear fit, and the Pearson
#' correlation with 95% CI, sample size, and RMSE in the plot title. Pairs
#' where either variable is NA, NaN, or infinite are dropped before analysis.
#'
#' @param predictions Numeric vector of predicted values.
#' @param truth Numeric vector of true values, same length as `predictions`.
#' @param model_label Model name shown in the plot title.
#' @param truth_label Axis label for the true values.
#' @param pred_label Axis label for the predicted values.
#' @param color_var Numeric vector mapped to point color, same length as
#'   `predictions`.
#' @param color_label Legend title for the color scale.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' set.seed(1)
#' truth <- rnorm(100)
#' plot_predictions_truth_scatter(
#'   predictions = truth + rnorm(100, sd = 0.5),
#'   truth = truth,
#'   model_label = "Example model",
#'   truth_label = "Truth",
#'   pred_label = "Prediction",
#'   color_var = truth,
#'   color_label = "True value"
#' )
plot_predictions_truth_scatter <- function(predictions, truth, model_label,
                                              truth_label, pred_label,
                                              color_var, color_label) {
  keep <- !(is.na(predictions) | is.nan(predictions) | is.infinite(predictions) |
    is.na(truth) | is.nan(truth) | is.infinite(truth))
  predictions <- predictions[keep]
  truth <- truth[keep]
  color_var <- color_var[keep]

  # Pearson correlation is invariant to log transforms of either variable
  cor_test <- stats::cor.test(predictions, truth, method = "pearson")
  cor_est <- signif(cor_test$estimate, digits = 3)
  cor_lwb <- signif(cor_test$conf.int[1], digits = 3)
  cor_upb <- signif(cor_test$conf.int[2], digits = 3)
  cor_n <- length(truth)

  rmse <- signif(sqrt(mean((predictions - truth)^2)), digits = 3)

  plotdata <- data.frame(predictions = predictions, truth = truth, color_var = color_var)

  ggplot2::ggplot(ggplot2::aes(y = predictions, x = truth, color = color_var), data = plotdata) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "#8C0172", linewidth = 2) +
    ggplot2::geom_smooth(method = "lm", color = "#9B951B", linewidth = 2) +
    ggplot2::scale_color_continuous(low = "#04050A", high = "#F9CCF9") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 16)) +
    ggplot2::labs(
      x = truth_label,
      y = pred_label,
      title = paste(
        model_label, "\nr = ", cor_est, ", n = ", cor_n,
        "\n95 % CI = [", cor_lwb, ", ", cor_upb, "]",
        "\nRMSE = ", rmse,
        sep = ""
      ),
      color = color_label
    )
}

#' Dotplot of PAF alignments between two sequences
#'
#' Plots alignment blocks from a PAF data frame (as returned by [read_paf()])
#' as line segments: query coordinates on the x axis, target coordinates on
#' the y axis, colored by strand. Reverse-strand blocks are oriented so that
#' collinear alignments lie on a shared diagonal. Low-quality blocks can be
#' dropped via `min_alen`/`min_mapq`.
#'
#' @param paf A PAF data frame from [read_paf()].
#' @param q_seq Query sequence name to plot (`qname`).
#' @param t_seq Target sequence name to plot (`tname`).
#' @param xlab X-axis label (query).
#' @param ylab Y-axis label (target).
#' @param title Plot title.
#' @param min_alen Minimum alignment block length to plot (default 1000).
#' @param min_mapq Minimum mapping quality to plot (default 1).
#'
#' @return A ggplot object. Add highlights (e.g. [ggplot2::geom_hline()]) and
#'   save it with [save_plot()].
#' @export
#'
#' @examples
#' paf <- data.frame(
#'   qname = c("q1", "q1", "q1", "q2"),
#'   qlen = 10000, qstart = c(100, 500, 100, 100), qend = c(4000, 9000, 4000, 4000),
#'   strand = c("+", "+", "-", "+"),
#'   tname = c("t1", "t1", "t1", "t2"),
#'   tlen = 10000, tstart = c(200, 600, 700, 200),
#'   tend = c(5000, 10000, 10000, 5000),
#'   nmatch = c(3900, 4500, 3800, 2900),
#'   alen = c(3900, 4500, 4000, 3900),
#'   mapq = c(60, 60, 60, 0)
#' )
#' p <- plot_paf_dotplot(paf, q_seq = "q1", t_seq = "t1",
#'                       xlab = "Query q1", ylab = "Target t1",
#'                       title = "Demo dotplot")
plot_paf_dotplot <- function(paf, q_seq, t_seq, xlab, ylab, title,
                             min_alen = 1000, min_mapq = 1) {
  paf_f <- dplyr::filter(paf, alen >= min_alen, mapq >= min_mapq)

  # Orient reverse-strand alignments for correct diagonal direction
  paf_plot <- dplyr::mutate(
    paf_f,
    tstart_plot = dplyr::if_else(strand == "-", tend, tstart),
    tend_plot = dplyr::if_else(strand == "-", tstart, tend)
  )

  paf_plot <- dplyr::filter(paf_plot, qname == q_seq, tname == t_seq)
  if (nrow(paf_plot) == 0) {
    stop("No alignments found between '", q_seq, "' and '", t_seq, "'.")
  }

  ggplot2::ggplot(paf_plot) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = qstart, xend = qend,
        y = tstart_plot, yend = tend_plot,
        color = strand
      ),
      linewidth = 0.6
    ) +
    ggplot2::scale_color_manual(
      values = c("+" = "steelblue", "-" = "firebrick")
    ) +
    ggplot2::labs(x = xlab, y = ylab, title = title, color = "Strand") +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

#' Manhattan plot from a snp table of stats
#'
#' Plots a Manhattan plot from a table of SNP associations: cumulative
#' genomic position on the x axis (`BPcum`, computed across chromosomes in
#' the order they appear in `CHROM`), `-log10(pvalue)` on the y axis, and
#' alternating grey/black points per chromosome. Chromosome labels are drawn
#' at each chromosome's center. SNPs with `highlight == "yes"` (if the column
#' is present) are overplotted in orange.
#'
#' @param data A data frame (or tibble) with columns `CHROM`, `POS`, and
#'   `pvalue`. An optional `highlight` column (`"yes"` marks points to
#'   overplot in orange) is used if present.
#' @param xlab X-axis label (default `"Chromosome"`).
#' @param ylab Y-axis label (default `"-log10(p)"`).
#' @param title Plot title.
#'
#' @return A ggplot object. Save it with [save_plot()].
#' @export
#'
#' @examples
#' snp_table <- data.frame(
#'   CHROM = rep(c("chr1", "chr2"), c(900, 600)),
#'   POS = c(1:900 * 1000, 1:600 * 3000),
#'   pvalue = runif(1500, 1e-8, 1)
#' )
#' plot_manhattan(snp_table, title = "Demo Manhattan plot")
plot_manhattan <- function(data, xlab = "Chromosome", ylab = "-log10(p)",
                           title = NULL) {
  required_cols <- c("CHROM", "POS", "pvalue")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "), ".")
  }
  if (nrow(data) == 0) {
    stop("No SNPs in 'data'.")
  }

  # Cumulative base-pair offset per chromosome, in order of first appearance
  chr_bounds <- dplyr::summarise(
    dplyr::group_by(data, CHROM),
    chr_len = max(POS),
    .groups = "drop"
  )
  chr_bounds <- dplyr::mutate(chr_bounds, tot = cumsum(chr_len) - chr_len)
  chr_bounds$chr_len <- NULL

  don <- dplyr::arrange(
    dplyr::left_join(data, chr_bounds, by = "CHROM"),
    CHROM, POS
  )
  don$BPcum <- don$POS + don$tot
  don$tot <- NULL

  axisdf <- dplyr::summarise(
    dplyr::group_by(don, CHROM),
    center = (max(BPcum) + min(BPcum)) / 2,
    .groups = "drop"
  )

  p <- ggplot2::ggplot(don) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = BPcum, y = -log10(pvalue),
        color = as.factor(CHROM)
      ),
      alpha = 0.8, size = 1.3, shape = 16
    ) +
    ggplot2::scale_color_manual(
      values = rep(c("grey", "black"), length.out = length(unique(don$CHROM)))
    ) +
    ggplot2::scale_x_continuous(
      breaks = axisdf$center,
      labels = axisdf$CHROM
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme(
      legend.position = "none",
      panel.border = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = xlab, y = ylab, title = title)

  if ("highlight" %in% names(don)) {
    p <- p + ggplot2::geom_point(
      data = dplyr::filter(don, highlight == "yes"),
      ggplot2::aes(x = BPcum, y = -log10(pvalue)),
      color = "orange", size = 1.5
    )
  }

  p
}

# Keep only the upper triangle of a square matrix: the strict lower triangle
# becomes NA, so it drops out when melting (same as the reshape2::melt recipe
# this function was ported from).
get_upper_tri <- function(mat) {
  mat[lower.tri(mat, diag = FALSE)] <- NA
  mat
}

#' Heatmap of a variance-covariance matrix
#'
#' Plots the upper triangle of a square variance-covariance (or correlation)
#' matrix as a heatmap: one tile per cell pair, filled by the cell value
#' (continuous green-magenta gradient centered at 0, or a discrete scico
#' palette). Cell values and/or tick labels can be suppressed for a
#' publication-style matrix figure.
#'
#' @param covmat A square matrix with symmetric row and column names.
#' @param include_values Print each cell's value (rounded to 4 decimals) on
#'   its tile.
#' @param include_tick_labels Show axis tick labels; if `FALSE`, all axis
#'   text and panel decoration is removed.
#' @param low_col Fill color for the most negative values.
#' @param mid_col Fill color at the midpoint (0).
#' @param high_col Fill color for the most positive values.
#' @param legend_name Legend title.
#' @param tile_or_raster `"tile"` (white cell borders) or `"raster"` (no
#'   borders) geometry.
#' @param discrete_or_continuous `"continuous"` uses a diverging color
#'   gradient (`low_col`/`mid_col`/`high_col`); `"discrete"` uses the scico
#'   `"bam"` palette (requires the scico package, a suggested dependency).
#'
#' @return A ggplot object. Save it with [save_plot()].
#'
#' @note The default legend title contains a Unicode character (Δ). Saving
#'   with the default `pdf()` device cannot render it; use a Unicode-capable
#'   device such as `cairo_pdf` (e.g. `ggplot2::ggsave(..., device =
#'   grDevices::cairo_pdf)`).
#' @export
#'
#' @examples
#' covmat <- matrix(
#'   c(1.0, 0.3, 0.1,
#'     0.3, 1.0, 0.5,
#'     0.1, 0.5, 1.0),
#'   nrow = 3, byrow = TRUE,
#'   dimnames = list(c("s1", "s2", "s3"), c("s1", "s2", "s3"))
#' )
#' # ASCII legend title keeps the example compatible with the default pdf()
#' # device; the default "Cov(\u0394pi, \u0394pj)" needs a Unicode-capable device
#' plot_var_cov_matrix(covmat, include_values = TRUE, legend_name = "Cov")
plot_var_cov_matrix <- function(covmat,
                                include_values = TRUE,
                                include_tick_labels = TRUE,
                                low_col = "#65014B",
                                mid_col = "#F5F0F0",
                                high_col = "#0C4C00",
                                legend_name = "Cov(\u0394pi, \u0394pj)",
                                tile_or_raster = c("tile", "raster"),
                                discrete_or_continuous = c("continuous", "discrete")) {
  tile_or_raster <- match.arg(tile_or_raster)
  discrete_or_continuous <- match.arg(discrete_or_continuous)

  if (!is.matrix(covmat) || nrow(covmat) != ncol(covmat)) {
    stop("'covmat' must be a square matrix.")
  }
  if (!identical(rownames(covmat), colnames(covmat))) {
    stop("'covmat' must have identical row and column names.")
  }

  upper_tri <- get_upper_tri(covmat)

  # Long format in Var1/Var2/value columns, upper-triangle cells only,
  # same semantics as reshape2::melt(..., na.rm = TRUE)
  melted_cormat <- na.omit(as.data.frame.table(
    upper_tri,
    responseName = "value"
  ))

  plot_cor <- ggplot2::ggplot(melted_cormat, ggplot2::aes(Var2, Var1, fill = value))

  if (tile_or_raster == "tile") {
    plot_cor <- plot_cor + ggplot2::geom_tile(color = "white")
  } else {
    plot_cor <- plot_cor + ggplot2::geom_raster()
  }

  if (include_values) {
    plot_cor <- plot_cor + ggplot2::geom_text(
      ggplot2::aes(label = round(value, 4)),
      size = 5,
      color = "black"
    )
  }

  if (discrete_or_continuous == "discrete") {
    if (!requireNamespace("scico", quietly = TRUE)) {
      stop("The scico package is required for discrete_or_continuous = 'discrete'.")
    }
    plot_cor <- plot_cor + scico::scale_fill_scico_d(
      palette = "bam",
      begin = 0.1,
      end = 0.9,
      name = legend_name
    )
  } else {
    plot_cor <- plot_cor + ggplot2::scale_fill_gradient2(
      low = low_col,
      mid = mid_col,
      high = high_col,
      midpoint = 0,
      space = "Lab",
      name = legend_name
    )
  }

  plot_cor <- plot_cor + ggplot2::theme_minimal()

  if (include_tick_labels) {
    plot_cor <- plot_cor + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      ),
      text = ggplot2::element_text(size = 14)
    )
  } else {
    plot_cor <- plot_cor + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    )
  }

  # Remove axis names
  plot_cor + ggplot2::coord_fixed() +
    ggplot2::labs(x = "", y = "")
}
