
\name{sine.taper}
\alias{sine.taper}
\title{Calculate a set of sine tapers}
\description{}
\usage{ sine.taper(N, n.tapers) }
\arguments{
  \item{N}{ The length of the time series }
  \item{n.tapers}{ The number of tapers to calculate }
  }
\value{ If n.tapers is one then the single sine taper is returned as a
  vector.  If n.tapers is greater than one, this function returns a
  matrix of dimension N by n.tapers.  Each column of the matrix contains
  a different sine taper. }
\references{ Percival and Walden (1993) \emph{Spectral Analysis for
    Physical Applications: Multitaper and Conventional
    Univariate Techniques}. Cambridge: Cambridge University Press. }
\details{}
\author{Peter F. Craigmile}
\keyword{}
