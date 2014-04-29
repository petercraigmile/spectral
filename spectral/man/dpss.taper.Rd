
\name{dpss.taper}
\alias{dpss.taper}
\title{Calculate a set of dpss tapers}
\description{}
\usage{ dpss.taper(N, n.tapers, NW = 4, W = NW/N) }
\arguments{
  \item{N}{ The length of the time series }
  \item{n.tapers}{ The number of tapers to calculate }
  \item{NW}{ N times by the bandwidth, W }
  \item{W}{ The bandwidth of the tapers, W }
  }
\value{ If n.tapers is one then the first dpss taper of bandwidth W is
  returned as a vector.  If n.tapers is greater than one, this function
  returns a matrix of dimension N by n.tapers.  Each column of the
  matrix contains a different dpss taper. }
\references{ Percival and Walden (1993) \emph{Spectral Analysis for
    Physical Applications: Multitaper and Conventional Univariate
    Techniques}. Cambridge: Cambridge University Press. }
\details{}
\author{Peter F. Craigmile}
\keyword{}


