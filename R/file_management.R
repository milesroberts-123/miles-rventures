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
