

## ======================================================================
## Calculating the covariance between discrete Fourier transform (DFT)
## ordinates.
## pfc@stat.osu.edu
##
## Reference: L. Wei and P. F. Craigmile (2010). Global and local
## spectral-based tests for periodicities. Biometrika, 97, 223-230.
## ======================================================================

arma.acvs <- function (ar=numeric(0), ma=numeric(0), sigma2=1, lag.max) {

  arma.var.from.sdf(ar, ma, sigma2) * ARMAacf(ar, ma, lag.max)
}



library(spectral)
source("DFT_covariance.R")


ar2.example.sigma2 <- 1/arma.var.from.sdf(ar2.example)
ar4.example.sigma2 <- 1/arma.var.from.sdf(ar4.example)

N <- 256

the.taper <- cosine.taper(N, 0.1)
acvs <- arma.acvs(ar=ar4.example, sigma2=ar4.example.sigma2, lag.max=N-1)

fs <- 1:(floor(N/2))/N

AA <- DFT.cov.matrix.slow(fs, N, acvs, the.taper, debug=TRUE)
BB <- DFT.cov.matrix(N, acvs, the.taper, skip.zero=TRUE)

## Compare the slow and fast versions to calculate the covariance.
summary(as.numeric(Re(AA)-Re(BB)))
summary(as.numeric(Im(AA)-Im(BB)))

plot(Re(diag(AA)))
points(Re(diag(BB)), col="red")
