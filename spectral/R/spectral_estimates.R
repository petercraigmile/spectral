


periodogram <- pgram <- function (x, delta.t=deltat(x))
{
  N <- length(x)
  spect(Mod(fft(x))^2/N, N, delta.t)
}


direct.spectral <-  function (x, taper.name, delta.t=deltat(x),
                              taper = create.taper(taper.name, length(x), ...), ...)
{
  if (is.null(dim(taper)))
    taper <- matrix(taper, ncol=1)
  
  spect(apply(taper, 2, function (h) Mod(fft(x * h)^2)), length(x), delta.t)
}





