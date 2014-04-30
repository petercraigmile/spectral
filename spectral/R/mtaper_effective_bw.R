
multitaper.effective.bw <- function (tapers)
  ## Calculate the effective bandwidth of a set of multitapers
  ## includes some optimizations to speed up calculation
  ## (The commented code is the original version.)
{

  ##calc.tau <- function (tau, h) {
  ##  ks <- 1:(N-abs(tau))
  ##  sum(h[ks] * h[ks+abs(tau)])
  ##}

  if (is.vector(tapers))
    tapers <- cbind(tapers)
  
  N <- nrow(tapers)

  AA <- apply(tapers, 2, function (x)
              N * acf(x, demean=FALSE, type="cov", plot=F, lag.max=(N-1))$acf)
  ##  AA <<- apply(tapers, 2, function (h)
  ##              sapply(0:(N-1), function (tau) calc.tau(tau, h)))
  if (is.vector(AA))
    AA <- cbind(AA)

  1/(1+2*sum(rowMeans(AA)[-1]^2))
}


