
## ======================================================================
## Functions to carry our power simulations for spectral-based tests
## of periodicity.
## pfc@stat.osu.edu
##
## Reference: L. Wei and P. F. Craigmile (2010). Global and local
## spectral-based tests for periodicities. Biometrika, 97, 223-230.
## ======================================================================

require(spectral)

ar2.example.sigma2 <- 1/arma.var.from.sdf(ar2.example)
ar4.example.sigma2 <- 1/arma.var.from.sdf(ar4.example)

ar2.example.sigma <- sqrt(ar2.example.sigma2)
ar4.example.sigma <- sqrt(ar4.example.sigma2)


gen.sim <- function (process, N, amp, freq) {

  if (process=="WN")
    x <- rnorm(N)
  else if (process=="AR2")
    x <- arima.sim(n=N, model=list(ar=ar2.example), sd=ar2.example.sigma)
  else ## if (process=="AR4")
    x <- arima.sim(n=N, model=list(ar=ar4.example), sd=ar4.example.sigma)
  
  x + periodicity(N, amp, freq)
}


power.sim <- function (amp, freq, process, N, nbins, REPS, test="localF",
                       taper=rectangular.taper(N), mtapers,
                       amp2, freq2) {

  if (missing(nbins))
    cat(N, " ", freq, "\n")
  else
    cat(N, " ", freq, " ", nbins, "\n")
  
  rejects <- NULL

  if (test=="hearing") {
    hcrit <- hcrits[nbins]
    #hcrit <- hearing.stat.crit(nbins)
  }
  else if (test=="Thomson")
    tcrit <-  Thomson.crit(ncol(mtapers))

  if (missing(amp2))
    nuisance <- 0
  else
    nuisance <- periodicity(N, amp2, freq2)
  
  for (k in 1:REPS) {
    x <- gen.sim(process, N, amp, freq) + nuisance
    sp <- direct.spectral(x, taper=taper)
    
    if (test=="global")
      rejects[k] <- global.Ftest(sp, freq=freq)$reject
    else if (test=="localF")
      rejects[k] <- local.Ftest(sp, freq=freq, nbins=nbins)$reject
    else if (test=="localF.linear")
      rejects[k] <- local.Ftest(sp, freq=freq, nbins=nbins,
                                local.model="lm")$reject
    else if (test=="hearing")
      rejects[k] <- hearing.stat(sp, freq=freq, nbins=nbins) >= hcrit
    else if (test=="Thomson")
      rejects[k] <- Thomson.stat(x, freq=freq, tapers=mtapers) >= tcrit
  }
  mean(rejects)
}


## An example

REPS <- 500

amps <- seq(0, 0.5, length=11)

pows <- sapply(amps, function (a) power.sim(a, freq=0.25, process="AR2",
                                            N=200, nbins=4, REPS=REPS, test="localF"))

plot(amps, pows, xlab="Amplitude", ylab="Power of local F test", type="b")

## Calculate pointwise 95% CIs (not the correct thing to do at the higher power values)
pm <- sqrt(pows * (1-pows) / REPS)

segments(amps, pows+pm, amps, pows-pm)

abline(h=0.05, lty=2)
