#' One iteration of minimization
#'
#' Computes step direction and adjusts step length to achieve target rate of decrease, with
#' various attempts to recover if this doesn't work out.  Used by \code{csminwel()}
#'
#' @param fcn Function that evaluates objective function
#' @param x0 starting point of iteration
#' @param f0 function value at start of iteration
#' @param g0 gradient from previous iteration
#' @param badg indicator for whether gradient is problematic
#' @param H0 inverse Hessian approximation at the start of iteration
#' @param Verbose Print out each iteration of step length change?
#'
#' @return
#'  \item{fhat}{New function value}
#' 
#'         \item{xhat}{New x value}
#' 
#'         \item{fcount}{Function evaluation count}
#' 
#'         \item{retcode}{return code (see Details)}
#'        
#'          
#' @details{ Meaning of \code{retcode} values: 0, normal step.  5, largest step still
#' improves too fast. 4,2 back and forth adjustment of stepsize didn't finish.
#' 3, smallest stepsize still improves too slow.  6, no improvement found.
#' 7, possible bad H matrix.  1, zero gradient.}
#' 
#'@export
#'
csminit <- function(fcn,x0,f0,g0,badg,H0, Verbose,...){
### Fixed 7/17/93 to use inverse-hessian instead of hessian itself in bfgs
### update.
###
### Fixed 7/19/93 to flip eigenvalues of H to get better performance when
### it's not psd.
###
  ## ANGLE <- .0005
  ANGLE <- 1e-7
  THETA <- .01 #(0<THETA<.5) THETA near .5 makes long line searches, possibly fewer iterations.
  FCHANGE <- 1000
  MINLAMB <- 1e-9
### fixed 7/15/94
### MINDX <- .0001;
### MINDX <- 1e-6;
  MINDFAC <- .01
  fcount<-0
  lambda<-1
  xhat<-x0
  f<-f0
  fhat<-f0
  g <- g0
  gnorm <- sqrt(sum(g^2))
###
  if (!badg && (gnorm < 1.e-12)) {      # put !badg 8/4/94
    retcode <- 1
    dxnorm <- 0
    ## gradient convergence
  } else {
    ## with badg true, we don't try to match rate of improvement to directional
    ## derivative.  We're satisfied just to get some improvement in f.
    ##
    ##if(badg)
    ##   dx = -g*FCHANGE/(gnorm*gnorm);
    ##  dxnorm = norm(dx)
    ##  if dxnorm > 1e12
    ##     disp('Bad, small gradient problem.')
    ##     dx = dx*FCHANGE/dxnorm;
    ##   end
    ##else
    ## Gauss-Newton step;
    ##---------- Start of 7/19/93 mod ---------------
    ##[v d] = eig(H0);
    ##toc
    ##d=max(1e-10,abs(diag(d)));
    ##d=abs(diag(d));
    ##dx = -(v.*(ones(size(v,1),1)*d'))*(v'*g);
    dx <- -H0 %*% g
    dxnorm <- sqrt(sum(dx^2))
    if (dxnorm > 1e12){
      cat("Near-singular H problem.\n")
      dx <- dx*FCHANGE/dxnorm
    }
    dfhat <- crossprod(dx,g0)
    if(!badg){
      ## test for alignment of dx with gradient and fix if necessary
      a <- -dfhat/(gnorm*dxnorm)
      ## We don't need'to handle negative a differently.  When
      ## a is close to zero, the adjustment below will make a small
      ## change in dx to get it to point slightly uphill, which is
      ## likely to be much better than switching to pure gradient.
      if(a<ANGLE){
        # if (a < 0) {
        #   dx <- -g / gnorm^2
        #   dfhat <- -1
        #   cat("H unused\n")
        #   ## Don't bother with H.  It's not psd. Step length here appropriate for log LH,
        #   ## where 1.0 is a reasonable scale for changes.
        # } else {
        dx <- dx - as.numeric(ANGLE*dxnorm/gnorm+dfhat/(gnorm*gnorm))*g
        dx <- dx * dxnorm / sqrt(sum(dx^2)) # This line added 2/18/2004
        dfhat <- crossprod(dx,g)
        ## dxnorm <- sqrt(sum(dx^2)) # No longer needed with 2/18/2004 change
        cat("Correct for low angle:" ,a,"\n")
      }
    }
    cat(sprintf("Predicted improvement            = %18.9f\n",-dfhat/2))
    ## cat("Predicted improvement:", sprintf("%18.9f",-dfhat/2),"\n")
    ##
    ## Have OK dx, now adjust length of step (lambda) until min and
    ## max improvement rate criteria are met.
    done <- 0
    factor <- 3
    shrink <- 1
    lambdaMin <- 0
    lambdaMax <- Inf
    lambdaPeak <- 0
    fPeak <- f0
    lambdahat <- 0
    while(!done){
      ## argument of fcn retains its dim (or lack thereof), but g
      ## always emerges as n x 1
      ddx <- dx*lambda
      dim(ddx) <- dim(x0)
      dxtest <- x0+ddx
      f <- fcn(dxtest,...)
      if (Verbose) cat(sprintf("lambda = %10.5g; f = %20.7f",lambda,f ),"\n")
      if(f < fhat){
        fhat <- f
        xhat <- dxtest
        lambdahat <- lambda
      }
      fcount <- fcount+1
      shrinkSignal <- (!badg && (f0-f < max(-THETA*dfhat*lambda,0))) || (badg && ((f0-f) < 0)) 
      growSignal <- (!badg && ( (lambda > 0)  &&  (f0-f > -(1-THETA)*dfhat*lambda) ))
      if(  shrinkSignal  &&   ( (lambda>lambdaPeak) || (lambda<0) )){
        if ((lambda>0) && ((!shrink || (lambda/factor <= lambdaPeak)))){
          shrink <- 1
          factor <- factor^.6
          while(lambda/factor <= lambdaPeak){
            factor <- factor^.6
          }
          if( abs(factor-1)<MINDFAC){
            if( abs(lambda) < 4){
              retcode <- 2
            }else{
              retcode <- 7
            }
            done <- 1
          }
        }
        if(  (lambda<lambdaMax) && (lambda>lambdaPeak) ){
          lambdaMax <- lambda
        }
        lambda <- lambda/factor
        if( abs(lambda) < MINLAMB ){
          if( (lambda > 0) & (f0 <= fhat) ){
            ## try going against gradient, which may be inaccurate
            lambda <- -lambda*factor^6
          }else{
            if( lambda < 0 ){
              retcode <- 6
            }else{
              retcode <- 3
            }
            done <- 1
          }
        }
      }else{
        if(  (growSignal && lambda>0) ||  (shrinkSignal && ((lambda <= lambdaPeak) && (lambda>0)))  ) {
          if( shrink ){
            shrink <- 0
            factor <- factor^.6
            if( abs(factor-1)<MINDFAC ) {
              if( abs(lambda)<4 ) {
                retcode <- 4
              }else{
                retcode <- 7
              }
              done <- 1
            }
          }
          if( ( f<fPeak ) && (lambda>0) ) {
            fPeak <- f
            lambdaPeak <- lambda
            if( lambdaMax<=lambdaPeak ) {
              lambdaMax <- lambdaPeak*factor*factor
            }
          }
          lambda <- lambda*factor
          if( abs(lambda) > 1e20 ) {
            retcode <- 5
            done <-1
          }
        } else {
          done <- 1
          if( factor < 1.2 ){
            retcode <- 7
          } else {
            retcode <- 0
          }
        }
      }
    }
  }
  if(Verbose) cat(sprintf("Norm of dx %10.5g\n", dxnorm))
  return(list(fhat=fhat,xhat=xhat,fcount=fcount,retcode=retcode))
}
