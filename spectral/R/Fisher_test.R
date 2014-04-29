
Fisher.test <- function (sp, alpha=.05)
{
  x       <- sp$spec[2:(sp$m - 1 + sp$N%%2)]
  df      <- length(x)
  reject  <- 1 - (alpha/df)^(1/(df-1))
  stat    <- max(x)/sum(x)
  p.value <- df * (1 - stat)^(df-1)
  where   <- order(-sp$spec)[1]
  
  structure(list(statistic = stat,
                 parameter = df,
                 where     = where,
                 p.value   = p.value,
                 reject    = stat>reject,
                 method    = "Fisher's test for a single periodicity"),
            class="h.test")
}
