
periodicity <- function (N, amps, freqs,
                         phases=runif(length(amps), 0, 2*pi),
                         delta.t=1) {
  ## ======================================================================
  ##
  ## ======================================================================

  ws <- 2 * pi * seq(N) * delta.t
  drop(rowSums(sapply(1:length(amps), function (k)
                      amps[k] * cos(freqs[k] * ws + phases[k]))))
}

