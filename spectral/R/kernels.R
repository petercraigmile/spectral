
Dirichlet.kernel <- function (freqs, N)
  ## ======================================================================
  ## Exercise 1.3c of Percival and Walden (1993)
  ## ======================================================================
{
  ws <- pi * freqs
  sin(N * ws) / (N * sin(ws))
}


Fejer.kernel <- function (freqs, N, delta.t=1)
  ## ======================================================================
  ## Equation 198b of Percival and Walden (1993)
  ## ======================================================================
{
  ws <- pi * freqs * delta.t
  (delta.t * sin(N * ws)^2) / (N * sin(ws)^2)
}
