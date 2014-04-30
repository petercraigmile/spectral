
local.Ftest <- function (sp, freq, nbins, mult=1, conf.level=0.95, local.model=NULL) {
  
  freq.index <- closest.Fourier.index(freq, sp$N, sp$delta.t)
  test.index <- freq.index+1

  lower1 <- test.index-mult*(nbins+1)
  lower2 <- test.index-mult
  upper1 <- test.index+mult
  upper2 <- test.index+mult*(nbins+1)
  
  fit.subset <- c(lower1:lower2, upper1:upper2)
  ss <- sp$spec[fit.subset]
  ff <- sp$freq[fit.subset]

  if (is.null(local.model)) {
    model <- NULL
  } else if (local.model=="linear") {
    model <- lm(log(ss) ~ ff)
  } else if (local.model=="quadratic") {
    model <- lm(log(ss) ~ ff + I(ff^2))
#  } else if (local.model=="gam") {
#    require(gam)
#    model <- gam(log(ss) ~ s(ff, 3))
  }

  denom.indexes <- c(seq(lower1, lower2, mult), seq(upper1, upper2, mult))
  Fstat <- sp$spec[test.index] / mean(sp$spec[denom.indexes])

  Sat    <- local.Ftest.Satterthwaite(freq, model, nbins)
  crit   <- qf(conf.level, Sat$df1, Sat$df2) * Sat$mult 
  pvalue <- 1 - pf(Fstat/Sat$mult, Sat$df1, Sat$df2) 

  
  list(test   = paste("LocalF", local.model, sep=""),
       Fstat  = Fstat,
       crit   = crit,
       reject = (Fstat >= crit),
       pvalue = pvalue)
}
  


local.Ftest.Satterthwaite <- function (freq, model=NULL, nbins) {

  if (is.null(model)) {

    list(mult=1, df1=2, df2=4*nbins)
  }
  else {
    
    est.sdf <- exp(fitted(model))
    est.sdf.fk <- exp(predict(model, data.frame(ff = freq)))
    
    sum.est.sdf <- sum(est.sdf)
    sum.est.sdf.sq <- sum(est.sdf^2)
    
    a.hat <- sum.est.sdf.sq/sum.est.sdf
    b.hat <- 2 * sum(est.sdf)^2 / sum.est.sdf.sq
    
    mult <- 2 * length(est.sdf) * est.sdf.fk / (a.hat * b.hat)
    
    list(mult=mult, df1=2, df2=b.hat)
  }
}
