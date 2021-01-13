## Last updated 2007-05-09, pfc@stat.osu.edu
## Created 2002-06-01, pfc@stat.osu.edu



spect <- function (spec, N, delta.t=1)
  ## ======================================================================
  ## Purpose : Create a class of type "spect" from a spectrum 'spec' of a 
  ##           time series of length 'N' regularly sampled at interval 'delta.t'.
  ## ======================================================================
{
  freq <- Fourier.frequencies(N, delta.t)
  n.Fourier <- length(freq)

  if (is.null(dim(spec)))
    actual.spec <- spec[1:n.Fourier]
  else
    actual.spec <- spec[1:n.Fourier,] 
  
  structure(list(freq      = freq,
                 spec      = actual.spec,
                 delta.t   = delta.t,
                 N         = N,
                 n.Fourier = n.Fourier),
            class = "spect")
}




mean.spect <- function (x, ...) 
  ## ======================================================================
  ## Purpose : Average the different spectra in the spect class 'x'
  ##           (normally used for multitapering).
  ## ======================================================================
{
  if (!is.vector(x$sp))
    x$spec <- drop(rowMeans(x$spec, ...))
  else
    x$spec <- mean(x$spec, ...)
  x  
}

plot.spect <- function (x, taper.num, xlab="Frequency", ylab="Spectrum",
                        log.freq=FALSE, add=FALSE,
                        type="l", ...)
  ## ======================================================================
  ## Purpose : Plot the object 'sp' of class 'spect' in decibels,
  ##           showing frequencies on the log scale if log.freq=T.
  ## ======================================================================
  
{
    if (missing(taper.num))
        spec <- x$spec
    else
        spec <- x$spec[,taper.num]

    if (log.freq) freqs <- log(x$freq)/log(2)
    else freqs <- x$freq
    
    if (add) lines(freqs, dB(spec), ...)
    else plot(freqs, dB(spec), type=type, xlab=xlab, ylab=ylab, ...)
}


spect.max <- function (sp, f1, f2)
  ## ==========================================================================
  ## Purpose: Returns the frequency of the maximum value of the spectra
  ##          sp of class 'spect', between frequencies 'f1' and 'f2' inclusive.

  ## ==========================================================================
{
  sel <- ((f1 <= sp$freq) & (sp$freq <= f2))
  sp$freq[sel][order(-sp$spec[sel])[1]]
}





