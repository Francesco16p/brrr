#' Summary function for an object with class brrr
#'
#' @param obj  Object with class brrr
#' @param round  Significant figures for the output
#'
#' @return
#' @export
#'
#' @examples
summary.brrr <- function(obj, round = 3)
{
  print(round(obj$out,round))
  cat("Nuisance parameter : ", obj$param, "\n")
  cat("Estimator : ", obj$method, "\n")
  invisible(obj)
}
