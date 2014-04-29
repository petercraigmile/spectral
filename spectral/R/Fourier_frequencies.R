
Fourier.frequencies <- function (N, delta.t=1) {
  ## ======================================================================
  ## Purpose : Calculates all the Fourier frequencies obtained when
  ##           performing spectral analysis on a time series of length 'N',
  ##           with sampling interval 'delta.t'.
  ## Updated : pfc@stat.osu.edu, April 2007.
  ## ======================================================================
  
  (seq(floor(N/2+1))-1) / (N * delta.t)
}



closest.Fourier.frequency <- function (freq, N, delta.t=1) {
  
  ## ======================================================================
  ## Purpose : Calculates the closest Fourier frequency to 'freq'
  ##           obtained when performing spectral analysis on a time series
  ##           of length 'N',  recorded at the sampling interval 'delta.t'.
  ## Updated : pfc@stat.osu.edu, April 2007.
  ## ======================================================================

  round(freq * N * delta.t) / (N * delta.t)
}


closest.Fourier.index <- function (freq, N, delta.t=1) {
  
  ## ======================================================================
  ## Purpose : Calculates the closest Fourier index to 'freq'
  ##           obtained when performing spectral analysis on a time series
  ##           of length 'N',  recorded at the sampling interval 'delta.t'.
  ## Updated : pfc@stat.osu.edu, April 2007.
  ## ======================================================================

  round(freq * N * delta.t)
}
