## ======================================================================
## File     : Deni_Wald.R
## Purpose  : R code to carry out the Denison-Walden test for periodicity.
## Updated  : pfc@stat.osu.edu, Aug 2004.
## 
## Copyright 2002--2004, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.ohio-state.edu.
##
## ======================================================================

Deni.Wald.crit <- function (n.tapers, fs, is.zero=rep(T,length(fs)),
                            alpha=0.05)
{
  m <- length(fs)
  r <- sum(!is.zero)
  qf(1 - alpha, 2 * (m-r), 2 * m * (n.tapers - 1))
}



Deni.Wald.setup <- function (N, deltat, fs, tapers, is.zero=rep(T,length(fs)))
{
  ## calculate Sigma and V
  ws <- -2i * pi * (0:(N-1)) * deltat
  m <- length(fs)
  V <- Sigma <- NULL
  for (k in 1:m)
  {
    W <- A <- NULL
    for (l in 1:m)
     {
      B <- t(tapers) * exp(ws*(fs[k]-fs[l]))
      A <- cbind(A, tapers %*% B)
      W <- cbind(W, colSums(sqrt(deltat) * B))
    }
    Sigma <- rbind(Sigma, A)
    V <- rbind(V, W)
  }

  ## calculate L, LV, LV^H for the full model.
  L <- solve(complex.chol(Sigma, upper=FALSE))
  LV <- L %*% V
  hLV <- t(Conj(LV))
  I <- diag(rep(1,dim(L)[1]))
  H <- L %*% (I - V %*% solve(hLV %*% LV) %*% hLV %*% L)

  ## calculate L, LV, LV^H for the reduced model.
  if (sum(!is.zero)==0)
    H2 <- L
  else
  {
    D <- as.matrix(diag(rep(1,m))[!is.zero,])
    LV2 <- LV %*% D
    hLV2 <- t(Conj(LV2))
    H2 <- L %*% (I - V %*% D %*% solve(hLV2 %*% LV2) %*% hLV2 %*% L)
  }
  
  list(H=H, H2=H2, ws=ws)
}


Deni.Wald.stat <- function (x, fs, tapers, is.zero = rep(T, length(fs)),
                            setup = Deni.Wald.setup(length(x),  deltat(x), fs, tapers, is.zero)) 
{
    K <- dim(tapers)[1]
    m <- length(fs)
    r <- sum(!is.zero)
    Y <- sqrt(deltat(x)) * as.vector(sapply(fs, function(f, z, 
        ws) colSums(z * exp(ws * f)), z = t(tapers) * as.vector(x), 
        setup$ws))
    SSEn <- sum(Mod(setup$H %*% Y)^2)
    SSEr <- sum(Mod(setup$H2 %*% Y)^2)
    m * (K - 1) * (SSEr - SSEn)/((m - r) * SSEn)
}
