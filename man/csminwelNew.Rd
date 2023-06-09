% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/csminwelNew.R
\name{csminwelNew}
\alias{csminwelNew}
\title{Gradient-based minimization}
\usage{
csminwelNew(
  fcn,
  x0,
  H0,
  ...,
  grad = NULL,
  crit = 1e-07,
  nit,
  Verbose = TRUE,
  Long = FALSE
)
}
\arguments{
\item{fcn}{The objective function to be minimized.  If its return value has a "gradient"
attribute (like output of \code{deriv()}), that is used as analytic gradient.
If it has a "hessian" attribute, that is used as the Hessian.}

\item{x0}{Initial value of the parameter vector.}

\item{H0}{Initial value for the inverse Hessian.  Must be positive definite,
if used.  (Not used if attr(fcn,"hessian") exists.)}

\item{...}{A list of optional length of additional parameters that get handed
off to \code{fcn} each time it is called.}

\item{grad}{If this is a numerical vector and \code{attr(fcn,"gradient")} is
not present, then \code{grad} is used as an initial gradient vector.
Useful for restarting if numerical gradient calculation is slow.}

\item{crit}{Convergence criterion.  Iteration will cease when it proves impossible
to improve the function value by more than this.}

\item{nit}{Maximum number of iterations.}
}
\value{
\if{html}{\out{<div class="sourceCode">}}\preformatted{   \\item\{fh\}\{Mimimized function value.\}
   \\item\{xh\}\{Argument value at optimum.\}
   \\item\{gh\}\{Gradient at optimum.\}
   \\item\{H\}\{Approximate inverse Hessian at optimum.\}
   \\item\{itct\}\{Number of iterations.\}
   \\item\{fcount\}\{Number of function evaluations\}
   \\item\{retcodeh\}\{Return code of last iteration (see
          \code{\link{csminit}})\}.
    \\item\{match.call()\}\{Call invoking the function\}
   \\item\{...\}\{Same as input parameter ...\}
}\if{html}{\out{</div>}}
}
\description{
Uses numeric or analytic gradients and a "satisficing" line search, building up
an approximate inverse Hessian with BFGS updates.  Looks in random directions or
tries to detect discontinuities when stuck.
}
\details{
Note that if the program ends abnormally, it is possible to retrieve the current
\code{x}, \code{f}, and \code{H} from the files g1.mat and H.mat that are written
at each iteration and at each hessian update, respectively.  (When the routine
hits certain kinds of difficulty, it writes g2.mat and g3.mat as well.  If all
were written at about the same time, any of them may be a decent starting point.
One can also start from the one with best function value.
}
