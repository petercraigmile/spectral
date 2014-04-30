

## Calculate covariances between DFT ordinates

## ======================================================================
## Calculating the covariance between discrete Fourier transform (DFT)
## ordinates.
## pfc@stat.osu.edu
##
## Reference: L. Wei and P. F. Craigmile (2010). Global and local
## spectral-based tests for periodicities. Biometrika, 97, 223-230.
## ======================================================================




DFT.cov <- function (eta1, eta2, N, acvs, taper, delta.t=1) {

  ws <- -2i * pi * seq(N) * delta.t
  delta.t * sum(outer(taper * exp(ws * eta1), taper * exp(-ws * eta2)) *
                toeplitz(acvs))
}

DFT.cov.matrix.slow <- function (etas, N, acvs, taper, delta.t=1,
                                 debug=FALSE) {

  DFT.cov.one <- function (eta1, eta2, Sig, taper)
    sum(outer(taper * exp(ws * eta1), taper * exp(-ws * eta2)) * Sig)
  
  ws <- -2i * pi * seq(N) * delta.t
  Sig <- toeplitz(acvs)

  delta.t * sapply(etas, function (f1)
                   sapply(etas, function (f2)
                          DFT.cov.one(f1, f2, Sig, taper)) })
}

DFT.cov.matrix <- function (N, acvs, taper, delta.t=1,
                            skip.zero=TRUE, every=1) {
  
  ts <- seq(N)
  ms <- seq(1, floor(N/2)+1, every)
  if (skip.zero) ms <- ms[-1]
  
  U <- sapply(ts, function (t)
              fft(taper * acvs[abs(t-ts)+1], inverse=T))
  delta.t * apply(U, 1, function (x) fft(x*taper)[ms])[,ms]
}
