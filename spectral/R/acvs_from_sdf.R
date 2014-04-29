
## ======================================================================
## Calculate the ACVS from lags 0 to lag 'max.lag',
## given a function specifying the SDF, 'the.sdf'.
##
## if use.integrate is TRUE use the 'integrate' funtion to calculate the
## integral.  Otherwise use Gauss-Legendre quadrature from the 'statmod'
## library using the 'gq' object, if available, or creating a new object
## with 'num.weights' weights.
##
## The ACVS at lag $k$ is given by
##
## s_k = \int_{-1/2}^{1/2} e^{i2 pi f k} S(f) df
##     = \int_{-1/2}^{1/2} \cos(2 pi f k) S(f) df
##     = 2 \int_{0}^{1/2} \cos(2 pi f k) S(f) df
##
## Letting u = 2f-1, so that du = 2 df, we get
##         f = (u+1)/2, u=-1..1
##
## s_k = \int_{-1}^{1} \cos(pi (u+1) k) S((u+1)/2) du
## ======================================================================

acvs.from.sdf <- function (max.lag, the.sdf, use.integrate=TRUE,
                           num.weights=200, gq, ...) {

  integrand <- function (f, lag, the.sdf, ...) {
    cos(2 * pi * f * lag) * the.sdf(f, ...)
  }

  if (use.integrate) {
    
    2 * sapply(0:max.lag, function (lag)
               integrate(integrand, 0, 0.5, the.sdf=the.sdf, lag=lag, ...)$value)

    
  } else {

    ## only needed if used outside of the R package.
    ## require(statmod)

    if (missing(gq))
      gq <- gauss.quad(num.weights)
    
    absc <- (gq$nodes+1)*0.5
    mult <- 0.5 * gq$weights * the.sdf(absc, ...)    
    colSums(mult * cos(outer(2 * pi * absc, 0:max.lag)))
  }
}




## Some examples

if (F) {
  
  WN.sdf <- function (f, sigma2=1)
    rep(sigma2, length(f))
  
  acvs.from.sdf(5, WN.sdf)
  
  acvs.from.sdf(5, WN.sdf, use.integrate=FALSE)
}
