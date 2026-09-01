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
