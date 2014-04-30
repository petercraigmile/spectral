## ======================================================================
## File    : Thomson.R
## Purpose : R/Splus code to carry out Thomson's F test
##           (See Chapter 10 of Percival and Walden (1993) -- SAPA).
## Note:     The C routine 'R_Thomson_F' needs to be dynamically loaded.
## Updated  : pfc@stat.osu.edu, June 2004.
## 
## Copyright 2002--2006, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
## ======================================================================


Thomson.stat <- function (x, freq, index = floor(freq * length(x) * deltat(x)),
                          tapers, n.tapers = dim(tapers)[2]) {
  ## ======================================================================
  ## Purpose: Calculate Thomson's F statistic.
  ## History: pfc@stat.ohio-state.edu, June 2002.
  ## ======================================================================

  z <- as.double(colSums(tapers))
  .C("R_Thomson_F",
     stat=as.double(0), as.integer(index),
     as.double(x), as.integer(length(x)),
     as.double(tapers), as.integer(n.tapers), z, z, z,
     PACKAGE="spectral")$stat
}


Thomson.crit <- function (n.tapers, alpha=0.05) {
  ## ======================================================================
  ## Purpose: Calculate the critical value.
  ## History: pfc@stat.osu.edu, June 2002.
  ## ======================================================================
  
  qf(1-alpha, 2, 2*(n.tapers-1))
}



Thomson.test <- function (x, freq, index = floor(freq * length(x) * deltat(x)),
                          tapers, n.tapers = dim(tapers)[2], alpha=0.05) {
  
  Fstat  <- Thomson.stat(x, freq, index, tapers, n.tapers)

  df2    <- 2*(n.tapers-1)
  
  crit   <- qf(1 - alpha, 2, df2)

  pvalue <- 1 - pf(Fstat, 2, df2)
  
  list(test   = "Thomson",
       Fstat  = Fstat,
       crit   = crit,
       reject = (Fstat >= crit),
       pvalue = pvalue) 
}
