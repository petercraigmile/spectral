
\name{direct.spectral}
\alias{direct.spectral}
\title{Direct spectral estimate}
\description{Calculate the tapered direct spectral estimate}
\usage{direct.spectral(x, taper.name, delta.t=deltat(x), taper=create.taper(taper.name, length(x), ...), ...)}
\arguments{
  \item{x}{Time series object (the data)}
  \item{taper.name}{The name of the taper}
  \item{delta.t}{Sampling interval}
  \item{taper}{The taper}
  \item{...}{Optional arguments used to create the taper (if necessary)}
}
\value{ An object of class 'spect' containing the direct spectral
  estimate of 'x' using the taper, 'taper'. }
\details{Calculate the tapered direct spectral estimate of a univariate time series of 'x'.}
\references{}
\seealso{}
\examples{}
\author{Peter F. Craigmile}
\keyword{}
