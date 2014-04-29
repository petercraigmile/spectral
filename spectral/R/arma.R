
## ======================================================================
## File     : arma.R
## Purpose  : Functions to calculate the spectral density functions (sdfs)
##            of ARMA processes, definitions of the
##            coefficients of two AR processes used in P&W (1993),
##            and a function to calculate the variance of an ARMA
##            process from using integrals of the sdf.
## Updated  : pfc@stat.osu.edu, Apr 2007.
## ======================================================================


ar.sdf <- function (freqs, phi, sigma2=1, delta.t=1) {
  ## ======================================================================
  ## Purpose: calculates the spectral density function (sdf) for an AR
  ##          process with coefficients 'phi' evaluated at the
  ##          frequencies 'freqs'.
  ##          Assumes the process is sampled at rate 'delta.t',
  ##          and the innovation variance is 'sigma2'.
  ## Updated: pfc@stat.osu.edu, Apr 2007.
  ## ======================================================================

  ws <- -2.0 * pi * delta.t * freqs
  js <- seq(length(phi))
  reals <- sapply(ws, function (w, js, phi) 1-sum(phi*cos(js*w)),
                  js = js, phi=phi)
  imags <- sapply(ws, function (w, js, phi) sum(phi*sin(js*w)),
                  js = js, phi=phi)
  (sigma2 * delta.t) / (reals*reals + imags*imags)
}




ma.sdf <- function (freqs, theta, sigma2=1, delta.t=1) {
  ## ======================================================================
  ## Purpose: calculates the spectral density function (sdf) for an MA
  ##          process with coefficients 'theta' evaluated at the
  ##          frequencies 'freqs'.
  ##          Assumes the process is sampled at rate 'delta.t',
  ##          and the innovation variance is 'sigma2'.
  ## Updated: pfc@stat.osu.edu, Apr 2007.
  ## ======================================================================
  
  ws <- -2.0 * pi * delta.t * freqs
  js <- seq(length(theta))
  reals <- sapply(ws, function (w, js, theta) 1+sum(theta*cos(js*w)),
                  js = js, theta=theta)
  imags <- sapply(ws, function (w, js, theta) sum(theta*sin(js*w)),
                  js = js, theta=theta)
  (sigma2 * delta.t) * (reals*reals + imags*imags)
}





arma.sdf <- function (freqs, phi, theta, sigma2=1, delta.t=1)
  ## ======================================================================
  ## Purpose: calculates the spectral density function (sdf) for an ARMA
  ##          process with AR coefficients 'phi' and MA coefficient 'theta'
  ##          evaluated at the frequencies 'freqs'.
  ##          Assumes the process is sampled at rate 'delta.t', and the
  ##          innovation variance is 'sigma2'.
  ## Updated: pfc@stat.osu.edu, Apr 2007.
  ## ======================================================================
{
  ws <- -2.0 * pi * delta.t * freqs

  ## calculate the AR part
  if (!missing(phi)) {
    js <- seq(length(phi))
    reals.ar <- sapply(ws, function (w, js, phi) 1-sum(phi*cos(js*w)),
                       js = js, phi=phi)
    imags.ar <- sapply(ws, function (w, js, phi) sum(phi*sin(js*w)),
                       js = js, phi=phi)
    denom <- (reals.ar^2 + imags.ar^2)
  }
  else denom <- rep(1, length(ws))

  ## calculate the MA part
  if (!missing(theta)) {
    ks <- seq(length(theta))
    reals.ma <- sapply(ws, function (w, ks, theta) 1+sum(theta*cos(ks*w)),
                       ks = ks, theta=theta)
    imags.ma <- sapply(ws, function (w, ks, theta) sum(theta*sin(ks*w)),
                       ks = ks, theta=theta)
    numer <- sigma2 * delta.t * (reals.ma^2 + imags.ma^2)
  }
  else numer <- sigma2 * delta.t

  ## now calculate the sdf
  numer/denom
}


arma.var.from.sdf <- function (phi, theta, sigma2=1, delta.t=1) {
  ## ======================================================================
  ## Purpose: calculates the variance for an ARMA process with AR
  ##          coefficients 'phi' and MA coefficient 'theta' and
  ##          innovation variance 'sigma2'.
  ## Updated: pfc@stat.osu.edu, Apr 2007.
  ## ======================================================================

  if (missing(phi) && missing(theta))
    2 * integrate(arma.sdf, 0, 0.5, sigma2=sigma2, delta.t=delta.t)$value
  else if (missing(phi))
    2 * integrate(arma.sdf, 0, 0.5, theta=theta, sigma2=sigma2, delta.t=delta.t)$value
  else if (missing(theta))
    2 * integrate(arma.sdf, 0, 0.5, phi=phi, sigma2=sigma2, delta.t=delta.t)$value
  else
    2 * integrate(arma.sdf, 0, 0.5, phi=phi, theta=theta,
                  sigma2=sigma2, delta.t=delta.t)$value
}




## ======================================================================
## Two example AR processes
## ======================================================================

##  Equation 45 of P&W (1993)
ar2.example <- c(0.75, -0.5)

##  Equation 46a of P&W (1993)
ar4.example <- c(2.7607, -3.8106, 2.6535, -0.9238)


