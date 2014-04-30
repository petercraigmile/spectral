
\name{Thomson.stat}
\alias{Thomson.stat}
\title{Thomson's F statistic}
\description{}
\usage{Thomson.stat(x, index, tapers, n.tapers)}
\arguments{
\item{x}{ The time series of interest }
\item{index}{ The Fourier index to test at }
\item{tapers}{ The matrix of tapers to use for the test }
\item{n.tapers}{ The number of tapers to evaluate }
}
\value{ The F statistic of Thomson's F test }
\details{ Perform Thomson's F test to test for a non-zero amplitude at
          Fourier index 'index' in the time series 'x' using 'n.tapers'
          tapers from the matrix of tapers 'tapers'.}
\references{ Percival and Walden (1993) \emph{Spectral Analysis for Physical Applications: Multitaper and Conventional Univariate Techniques}. Cambridge: Cambridge University Press, Section 10.11. }
\examples{}
\author{Peter F. Craigmile}
\keyword{}
