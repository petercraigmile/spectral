
hearing.stat.pvalue <- function (x, nbins) {
  
  if (nbins==1)
    2/((x+2)*(x+1))
  else if (nbins==2)
    32 * (3*x+8)/((x+4)^3*(x+2)^2)
  else if (nbins==3)
    1458 * (144 + 75*x + 10*x^2) / ((3 + x)^3 * (6 + x)^5)
  else if (nbins==4)
    131072*(504*x^2+2464*x+35*x^3+4096)/((8+x)^7*(x+4)^4)
  else if (nbins==5)
    19531250*(104625*x+160000+2940*x^3+126*x^4+26100*x^2)/((10+x)^9*(x+5)^5)
  else if (nbins==6)
    26121388032 * (1327104 + 917136*x + 36630*x  + 2640*x  + 77*x  + 257400*x)/((12 + x)^11*(6 + x)^6)
  else if (nbins==7)
    1356446145698 * (346526726*x + 17265248*x^3  + 1716*x^6  + 1611610*x^4  + 81081*x^5 + 105250236*x^2  + 481890304) / ((x + 7)^7 * (x + 14)^13 )
  else
    stop("The pvalue is only defined for nbins=1,...,7.")
}


hearing.stat.crit <- function (nbins=4, alpha=0.05) {
## assumes that the solution lies in the interval (0, 500)

  uniroot(function (x) (hearing.stat.pvalue(x, nbins)-alpha), c(0, 500))$root
}


hearing.stat <- function (sp, freq, nbins = 4, offset = 1) {
  freq.index <- closest.Fourier.index(freq, sp$N, sp$delta.t)
  test.index <- freq.index+1

  bs <- (1:nbins) + offset
  sp$spec[test.index] /
    max(mean(sp$spec[test.index + 1 - bs]), mean(sp$spec[test.index + 1 + bs]))
}



## ======================================================================
## Calculate the hearing test critical values to 4 decimal places
## for B = 1, ..., 10 bins via Monte Carlo simulation
## ======================================================================

## sps <- lapply(1:50000, function (j) periodogram(rnorm(N)))
## hcrits <- sapply(1:10, function (B)
##                  quantile(sapply(sps, function (sp)
##                                  hearing.stat(sp, 0.25, B)), 0.95))

hcrits <- c(4.8494, 3.3527, 3.0627, 2.9387, 2.8788, 2.8512, 2.8146, 2.8060, 2.8017, 2.8015) 

