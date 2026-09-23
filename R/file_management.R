#' Check if a file exists and delete it if so
#'
#' @param file_name Path to the file to check.
#'
#' @return Invisibly, `TRUE` if the file existed (and was deleted), `FALSE`
#'   otherwise.
#' @export
#'
#' @examples
#' \dontrun{
#' check_file_exists("output.csv")
#' }
check_file_exists <- function(file_name) {
  if (file.exists(file_name)) {
    file.remove(file_name)
    message(file_name, " from previous run deleted.")
    invisible(TRUE)
  } else {
    message(file_name, " does not exist.")
    invisible(FALSE)
  }
}

#' Append a row to a CSV table, good for for loops
#'
#' Writes `x` to `output_name` as a comma-separated row. If the file does not
#' exist yet, column names are written first.
#'
#' @param x A data frame or matrix to append.
#' @param output_name Path to the output CSV file.
#'
#' @return Invisibly, the return value of [data.table::fwrite()].
#' @export
#'
#' @examples
#' \dontrun{
#' append_table(data.frame(a = 1, b = 2), "results.csv")
#' }
append_table <- function(x, output_name) {
  if (!is.data.frame(x)) {
    x <- as.data.frame(x) # silence fwrite's matrix-coercion message
  }
  data.table::fwrite(x,
                     output_name,
                     col.names = !file.exists(output_name),
                     append = TRUE,
                     na = "NA",
                     quote = FALSE)
}

#' Read a PAF alignment file
#'
#' Reads a PAF (Pairwise mApping Format) file produced by minimap2, keeps the
#' 12 mandatory columns, coerces the numeric fields, and adds a per-alignment
#' percent-identity column.
#'
#' @param file Path to a PAF file. PAF files written by `minimap2 --paf` have
#'   no header row; extra columns beyond the 12 mandatory ones are ignored.
#'
#' @return A tibble with 13 columns: the 12 mandatory PAF fields
#'   (`qname`, `qlen`, `qstart`, `qend`, `strand`, `tname`, `tlen`,
#'   `tstart`, `tend`, `nmatch`, `alen`, `mapq` — all but `qname`, `strand`,
#'   and `tname` numeric) plus `pident` (`nmatch / alen`).
#' @export
#'
#' @examples
#' \dontrun{
#' paf <- read_paf("hap2_vs_hap1.paf")
#' }
read_paf <- function(file) {
  # 12 mandatory PAF columns
  paf_cols <- c(
    "qname", "qlen", "qstart", "qend", "strand",
    "tname", "tlen", "tstart", "tend",
    "nmatch", "alen", "mapq"
  )

  # Read every field as character (PAF extra columns are ragged), then keep
  # the 12 mandatory columns and coerce the numeric fields.
  paf <- readr::read_tsv(
    file,
    col_names = FALSE,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  )
  paf <- paf[, seq_len(12)]
  names(paf) <- paf_cols

  dplyr::mutate(
    paf,
    dplyr::across(c(qlen, qstart, qend, tlen, tstart, tend, nmatch, alen, mapq), as.numeric),
    pident = nmatch / alen
  )
}

#' Read a list of files and bind their tables together
#'
#' Reads each file in `file_list` with [data.table::fread()] and binds the
#' resulting tables row-wise into a single data frame. A source-file column
#' (full paths, in the order given) is prepended unless `id = NULL`.
#'
#' @param file_list Character vector of file paths, e.g. the output of
#'   [list.files()] with `full.names = TRUE`.
#' @param col_names Optional character vector of column names passed to
#'   `fread`'s `col.names`, overriding the headers in the files. `NULL`
#'   keeps the file headers.
#' @param id Name of the source-file column added to the result, holding
#'   the full path of the file each row came from. `NULL` omits the column.
#' @param ... Additional arguments passed to [data.table::fread()], e.g.
#'   `sep`, `select`, or `header`.
#'
#' @return A data frame with one row per row of every input file, plus the
#'   id column if `id` is not `NULL`.
#' @export
#'
#' @examples
#' f1 <- tempfile(fileext = ".txt")
#' f2 <- tempfile(fileext = ".txt")
#' writeLines("0.1\t0.9", f1)
#' writeLines("0.2\t0.8", f2)
#' bind_files(c(f1, f2), col_names = c("bray_curtis", "cosine"))
bind_files <- function(file_list, col_names = NULL, id = "source_file", ...) {
  stopifnot(is.character(file_list), length(file_list) > 0)
  missing_files <- !file.exists(file_list)
  if (any(missing_files)) {
    stop("File(s) not found: ", paste(file_list[missing_files], collapse = ", "))
  }
  tables <- if (is.null(col_names)) {
    lapply(file_list, data.table::fread, ...)
  } else {
    lapply(file_list, data.table::fread, col.names = col_names, ...)
  }
  names(tables) <- file_list
  dplyr::bind_rows(tables, .id = id)
}

#' Create the parent directories of a target file path, if missing
#'
#' Takes the path to a target file, extracts its folder portion with
#' [base::dirname()], and creates the full folder structure with
#' [base::dir.create()] if it does not yet exist. Useful before writing
#' outputs into a not-yet-existing nested directory.
#'
#' @param file_path Path to the target file. Only the folder portion is
#'   created; the file itself is left untouched.
#'
#' @return Invisibly, `TRUE` if the folder structure was created, `FALSE` if
#'   it already existed.
#' @export
#'
#' @examples
#' target <- file.path(tempdir(), "nested", "data_results", "output.csv")
#' ensure_parent_dir(target)   # TRUE: folders created
#' ensure_parent_dir(target)   # FALSE: folders already exist
ensure_parent_dir <- function(file_path) {
  stopifnot(is.character(file_path), length(file_path) == 1)
  folder_path <- dirname(file_path)
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
    invisible(TRUE)
  } else {
    invisible(FALSE)
  }
}

#' Convert a large CSV file to a parquet dataset for lazy loading
#'
#' One-time conversion of a large CSV file into a directory of parquet files
#' that can be opened with [arrow::open_dataset()] and queried lazily (only
#' the needed chunks are read into memory). If `parquet_dir` does not exist
#' yet, the CSV is read with a large block size and written as parquet in
#' chunks of `max_rows_per_file` rows; otherwise the conversion step is
#' skipped. In both cases the lazy dataset is returned, so the function is
#' safe to place at the top of a script that is re-run repeatedly.
#'
#' Note that an existing but empty or incomplete `parquet_dir` is trusted as
#' a completed conversion, mirroring the `dir.exists()` guard this function
#' wraps.
#'
#' @param csv_path Path to the CSV file to convert.
#' @param parquet_dir Path to the output directory of parquet files.
#' @param block_size Block size in bytes for reading the CSV. Larger blocks
#'   speed up parsing of wide files. Default is 16 MiB.
#' @param max_rows_per_file Maximum number of rows per parquet file in the
#'   output dataset. Default is 50000.
#'
#' @return A lazy [arrow::Dataset] object opened from `parquet_dir`.
#' @export
#'
#' @examples
#' \dontrun{
#' # One-time conversion, then lazy reads on every re-run:
#' hapfire <- csv_to_parquet(
#'   "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/data/SV_SNP_INDEL_allele_frequency_trays.csv",
#'   "/global/scratch/users/milesroberts/moi_lab_projects/grenenet-phase2/data/kmate_parquet"
#' )
#' }
csv_to_parquet <- function(csv_path, parquet_dir,
                           block_size = 16L * 1024L * 1024L,
                           max_rows_per_file = 50000L) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("The arrow package is required: install.packages(\"arrow\").")
  }
  stopifnot(is.character(csv_path), length(csv_path) == 1,
            is.character(parquet_dir), length(parquet_dir) == 1)
  if (!file.exists(csv_path)) {
    stop("File not found: ", csv_path)
  }
  if (!is.numeric(block_size) || length(block_size) != 1 ||
      !is.finite(block_size) || block_size <= 0) {
    stop("block_size must be a single positive number.")
  }
  if (!is.numeric(max_rows_per_file) || length(max_rows_per_file) != 1 ||
      !is.finite(max_rows_per_file) || max_rows_per_file <= 0) {
    stop("max_rows_per_file must be a single positive number.")
  }
  if (!dir.exists(parquet_dir)) {
    arrow::open_dataset(csv_path,
                        format = "csv",
                        read_options = arrow::csv_read_options(block_size = block_size)) |>
      arrow::write_dataset(parquet_dir,
                           format = "parquet",
                           max_rows_per_file = max_rows_per_file)
  }
  arrow::open_dataset(parquet_dir)
}

# ---------------------------------------------------------------------------
# Block construction and window helpers (ported from grenenet-phase2)
# ---------------------------------------------------------------------------

#' Block data by chromosome and position
#'
#' Assigns each row of a sorted variant table to a block, computed either by
#' SNP count or by physical distance, and never across chromosome boundaries.
#'
#' @param data Data frame with chromosome and position columns.
#' @param chrom_col_name Chromosome column (tidy-evaluated).
#' @param pos_col_name Position column (tidy-evaluated).
#' @param window_size Block size: SNPs per block if `block_by = "snp"`, base
#'   pairs per block if `block_by = "base"`.
#' @param sep Separator pasted between chromosome and block index to form
#'   the `window` label.
#' @param block_by Either `"snp"` (fixed number of SNPs per block) or
#'   `"base"` (fixed physical span per block).
#'
#' @return `data` with columns `block` (block index within chromosome) and
#'   `window` (chromosome-sep-block label).
#' @export
#'
#' @examples
#' \dontrun{
#' blocks <- create_blocks(snp_table, chrom, pos, 15000,
#'                         block_by = "snp", sep = "_")
#' }
create_blocks <- function(data,
                          chrom_col_name,
                          pos_col_name,
                          window_size,
                          sep,
                          block_by = "snp") {
  # check data is sorted
  sort_check <- data %>%
    dplyr::group_by({{ chrom_col_name }}) %>%
    dplyr::summarise(is_sorted = all(diff({{ pos_col_name }}) >= 0))

  if (!all(sort_check$is_sorted)) {
    stop("Input data should be sorted by position within chromosome.")
  }

  if (block_by == "snp") {
    blocks <- data %>%
      dplyr::group_by({{ chrom_col_name }}) %>%
      dplyr::mutate(
        block = (dplyr::row_number() - 1) %/% window_size + 1,
        window = paste({{ chrom_col_name }}, block, sep = sep)
      )
  } else if (block_by == "base") {
    blocks <- data %>%
      dplyr::group_by({{ chrom_col_name }}) %>%
      dplyr::mutate(
        bin_start = ({{ pos_col_name }} %/% window_size) * window_size,
        bin_end = bin_start + window_size,
        block = ({{ pos_col_name }} %/% window_size),
        window = paste({{ chrom_col_name }}, block, sep = sep)
      )
  } else {
    stop("block_by not valid. Choose from: snp, base.")
  }

  return(blocks)
}

#' Grab sample sizes for a subset of samples
#'
#' Subsets a sample-size table to the requested samples, prepends the
#' generation-0 size, and sorts by sample label.
#'
#' @param n_data Data frame with a sample-label column and a sample-size
#'   column.
#' @param sample_col Sample-label column (tidy-evaluated).
#' @param samples Character vector of sample labels to keep.
#' @param first_n Sample size to prepend for generation 0.
#' @param n_col Sample-size column (tidy-evaluated).
#'
#' @return A data frame with two columns (sample label and size), the first
#'   row holding the generation-0 entry.
#' @export
grab_sample_sizes <- function(n_data, sample_col, samples, first_n, n_col) {
  rep_n <- n_data %>% dplyr::filter({{ sample_col }} %in% samples)
  rep_n <- rbind(c(0, first_n), rep_n %>% dplyr::select({{ sample_col }}, {{ n_col }}))
  rep_n <- rep_n %>% dplyr::arrange({{ sample_col }})
  return(rep_n)
}

# ---------------------------------------------------------------------------
# Polynomial helpers (ported from grenenet-phase2 resources/R/misc.R)
# ---------------------------------------------------------------------------

#' Fit polynomial to data
#'
#' Fits polynomials of degree 1 through `max_degree` by ordinary least
#' squares and returns the model with the lowest AIC.
#'
#' @param data Data frame holding the predictor and response columns.
#' @param x_col Name of the predictor column, as a string.
#' @param y_col Name of the response column, as a string.
#' @param max_degree Highest polynomial degree to consider (default 10).
#'
#' @return The best-fitting `lm` object.
#' @export
find_best_poly <- function(data, x_col, y_col, max_degree = 10) {
  x_val <- data[[x_col]]
  y_val <- data[[y_col]]

  models <- lapply(1:max_degree, function(d) {
    lm(y_val ~ poly(x_val, d, raw = TRUE))
  })
  aics <- sapply(models, AIC)

  best_idx <- which.min(aics)
  best_model <- models[[best_idx]]

  return(best_model)
}

#' Evaluate derivative of a raw-polynomial model at a point
#'
#' Extracts the coefficients of a raw (`raw = TRUE`) polynomial fit and
#' evaluates its derivative at `x_eval`.
#'
#' @param model An `lm` object fitted with `poly(x, d, raw = TRUE)`.
#' @param x_eval Numeric value at which to evaluate the derivative.
#'
#' @return A single numeric value: the derivative of the fitted polynomial
#'   at `x_eval`.
#' @export
eval_deriv <- function(model, x_eval) {
  beta <- coef(model)
  deriv_val <- 0
  if (length(beta) > 1) {
    for (i in 2:length(beta)) {
      deg <- i - 1
      deriv_val <- deriv_val + (deg * beta[i] * (x_eval ^ (deg - 1)))
    }
  }
  return(as.numeric(deriv_val))
}
