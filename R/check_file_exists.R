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
    message(paste(file_name, "from previous run deleted."))
    invisible(TRUE)
  } else {
    message(paste(file_name, "does not exist."))
    invisible(FALSE)
  }
}
