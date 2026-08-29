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
