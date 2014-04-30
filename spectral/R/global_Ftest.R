
global.Ftest <- function (sp, freq, index = floor(freq * sp$N * sp$delta.t),
                          nuisance.freqs,
                          nuisance.indexes = floor(nuisance.freqs * sp$N * sp$delta.t),
                          alpha=0.05) {

  m  <- length(sp$spec)
  sp.nz <- sp$spec[-c(1,m)]

  if (missing(nuisance.freqs))
    denom.omit <- index
  else
    denom.omit <- c(index, nuisance.indexes)
  
  Fstat <- sp.nz[index] / mean(sp.nz[-denom.omit])

  ss <- sp.nz[-denom.omit]
  ff <- (sp$freq[-c(1,m)])[-denom.omit]

  require(gam)
  model <- gam(log(ss) ~ s(ff, 3))
  
  Sat <- Ftest.Satterthwaite(freq, model)
  
  crit   <- qf(1 - alpha, Sat$df1, Sat$df2) * Sat$mult
  pvalue <- 1 - pf(Fstat / Sat$mult, Sat$df1, Sat$df2)
    
  list(test   = "GlobalF",
       stat   = Fstat,
       crit   = crit,
       reject = (Fstat>=crit),
       pvalue = pvalue)
}

