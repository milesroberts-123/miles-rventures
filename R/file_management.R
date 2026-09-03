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
#' @return Invisibly, the return value of [utils::write.table()].
#' @export
#'
#' @examples
#' \dontrun{
#' append_table(data.frame(a = 1, b = 2), "results.csv")
#' }
append_table <- function(x, output_name) {
  utils::write.table(x,
              output_name,
              col.names = !file.exists(output_name),
              append = TRUE,
              row.names = FALSE,
              sep = ",",
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
