#' Load example dataset for gFlora
#'
#' @return A list with \code{M} (matrix/data.frame) and \code{y} (vector).
#' @export
load_example_data <- function() {

  M_path <- system.file("extdata", "M.csv", package = "gFlora")
  y_path <- system.file("extdata", "fv.csv", package = "gFlora")

  if (M_path == "" || y_path == "") {
    stop("Example data files not found. Make sure inst/extdata/M.csv and inst/extdata/fv.csv exist, then reinstall the package.")
  }

  M <- utils::read.csv(M_path, row.names = 1)
  y <- utils::read.csv(y_path, row.names = 1)[, 1]

  list(M = as.matrix(M), y = y)
}
