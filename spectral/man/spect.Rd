
\name{spect}
\alias{spect}
\title{Create a object of class 'spect'}
\description{}
\usage{spect(spec, N, delta.t=1)}
\arguments{
  \item{spec}{ a numerical vector }
  \item{N}{ length of the original time series }
  \item{delta.t}{ sampling interval }
  }
\value{ An object of class 'spect' which contains the frequencies
  'freq', spectral estimates 'spec', sampling interval 'delta.t', sample
  size 'N' and number of Fourier frequencies, 'n.Fourier'. }
\details{Create a class of type 'spect' from a spectrum 'spec'of a time
  series of length 'N' regularly sampled at rate 'delta.t'.}
\references{}
\examples{}
\author{Peter F. Craigmile}
\keyword{}

