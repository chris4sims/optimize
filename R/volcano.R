#' volcano
#'
#' A "tilted volcano" objective function for testing optimizers
#'
#' @param `x` a vector of arbitrary dimension n
#'
#' @return Minus the height of the volcano at `x`
#'
#' @details Gradient based routines may tend to go to the volcano rim, then
#' only very slowly crawl along it to the peak.
#' @export
#' @md
volcano <- function(x) {
  rho <- sqrt(sum(x))
  f <- rho^2*exp(-rho)+.02*x[1]/(1+.001*x[1]^2)
  return(-f)
}
