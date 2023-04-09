#' numHess
#'
#' Numerical Hessian.  This is a naive calculation, just
#' applying `numgrad()`'s simple forward difference estimate
#' twice.
#'
#' @param fcn The function to be differentiated.  
#' @param x The point at which the Hessian is to be calculated.
#'
#' @export
#' @md
numHess <- function(fcn, x, ...) {
  f1 <- fcn
  n <- length(x)
  h <- matrix(0, n, n)
  f2 <- function(z, ...) { numgrad(fcn=f1, z, ...)$g}
  h <- numgrad(fcn=f2, x=x, ...)$g
  return(h)
}
    
