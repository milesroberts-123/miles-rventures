#' Expression matrix (G genes x S samples)
#'
#' Constructs a validated gene expression matrix with one row per gene and
#' one column per sample. Values are expression measurements (counts, CPM,
#' TPM, log-transformed values, ...); `NA` and negative values are allowed.
#'
#' @param x A numeric matrix with at least one row and one column.
#'
#' @return An `expr_matrix` object: a numeric matrix of expression values.
#' @export
#'
#' @examples
#' m <- matrix(c(10, 0, 5, 8, 2, 7), nrow = 3)
#' expr_matrix(m)
expr_matrix <- function(x) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("x must be a numeric matrix.")
  }
  if (nrow(x) < 1 || ncol(x) < 1) {
    stop("x must have at least one row and one column.")
  }
  class(x) <- c("expr_matrix", class(x))
  x
}

#' Sample metadata for expression data (S samples x any columns)
#'
#' Constructs a validated sample metadata table with one row per sample and
#' any number of annotation columns (cell type, tissue, condition, batch,
#' ...). `NA` values are allowed.
#'
#' @param x A data.frame with S rows (at least one) and any number of
#'   columns.
#'
#' @return An `expr_sample_meta` object: a data.frame of sample metadata.
#' @export
#'
#' @examples
#' d <- data.frame(
#'   cell_type = c("neuron", "glia", "neuron"),
#'   condition = c("ctrl", "ctrl", "treat")
#' )
#' expr_sample_meta(d)
expr_sample_meta <- function(x) {
  if (!is.data.frame(x)) {
    stop("x must be a data.frame.")
  }
  if (nrow(x) < 1) {
    stop("x must have at least one row.")
  }
  class(x) <- c("expr_sample_meta", class(x))
  x
}

#' Print an expr_matrix
#'
#' @param x An `expr_matrix` object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The object, invisibly.
#' @export
#'
#' @examples
#' m <- expr_matrix(matrix(c(10, 0, 5, 8, 2, 7), nrow = 3))
#' print(m)
print.expr_matrix <- function(x, ...) {
  cat(sprintf("expr_matrix: %d genes x %d samples\n", nrow(x), ncol(x)))
  if (!is.null(colnames(x))) {
    cat("Samples:", paste(colnames(x), collapse = ", "), "\n")
  }
  print(utils::head(unclass(round(x, 4))), ...)
  invisible(x)
}

#' Print an expr_sample_meta object
#'
#' @param x An `expr_sample_meta` object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The object, invisibly.
#' @export
#'
#' @examples
#' d <- expr_sample_meta(data.frame(
#'   cell_type = c("neuron", "glia"), condition = c("ctrl", "ctrl")
#' ))
#' print(d)
print.expr_sample_meta <- function(x, ...) {
  cat(sprintf("expr_sample_meta: %d samples\n", nrow(x)))
  print(utils::head(unclass(x)), ...)
  invisible(x)
}

#' Cell-type specificity (tau)
#'
#' Computes the tau specificity index of an expression vector (Kryuchkov
#' et al. 2016): 0 means ubiquitous expression across samples, 1 means
#' expression concentrated in a single sample. The vector is scaled by its
#' maximum before averaging, so tau is scale-invariant.
#'
#' @param x A numeric vector of expression values for one gene.
#' @param na_action How to handle `NA` values: `"propagate"` returns
#'   `NA_real_` if any value is missing; `"remove"` computes tau on the
#'   observed values only (returning `NA_real_` if fewer than two remain);
#'   `"error"` stops on missing values.
#'
#' @return A single numeric value in \code{[0, 1]}, or `NA_real_` for
#'   degenerate input (fewer than 2 values, all zeros, or propagated
#'   missingness).
#' @export
#'
#' @examples
#' tau(c(1, 1, 1, 1))   # ubiquitous: 0
#' tau(c(10, 0, 0, 0))  # single sample: 1
#' tau(c(5, 2, 8, 1))
tau <- function(x, na_action = c("propagate", "remove", "error")) {
  na_action <- match.arg(na_action)
  if ((!is.numeric(x) && !is.logical(x)) || is.matrix(x)) {
    stop("x must be a numeric vector.")
  }
  x <- as.numeric(x)
  if (anyNA(x)) {
    if (na_action == "error") {
      stop("x must not contain NA values.")
    }
    if (na_action == "propagate") {
      return(NA_real_)
    }
    x <- x[!is.na(x)]
  }
  n <- length(x)
  if (n <= 1 || max(x) == 0) {
    return(NA_real_)
  }
  x <- x / max(x)
  sum(1 - x) / (n - 1)
}

#' Mean expression of expressed samples only
#'
#' Computes the mean of an expression vector conditioned on the gene being
#' expressed at all: values at or below `thresh` are dropped before
#' averaging. Returns `NA_real_` when no value exceeds `thresh`.
#'
#' @param x A numeric vector of expression values for one gene.
#' @param thresh Numeric expression threshold; only values strictly greater
#'   than `thresh` contribute to the mean.
#'
#' @return A single numeric value: the mean of the values in `x` that are
#'   strictly greater than `thresh` and not `NA`, or `NA_real_` if no value
#'   qualifies.
#' @export
#'
#' @examples
#' mean_no_zeros(c(0, 2, 4, 8), thresh = 0)
#' mean_no_zeros(c(0, 0, 0), thresh = 0)
mean_no_zeros <- function(x, thresh) {
  y <- x[x > thresh & !is.na(x)]
  if (length(y) == 0) {
    return(NA_real_)
  }
  mean(y, na.rm = TRUE)
}