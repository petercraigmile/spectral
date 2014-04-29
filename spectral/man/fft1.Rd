
\name{fft1}
\alias{fft1}
\title{Calculate fft of a vector using a starting time index of 1}
\description{  }
\usage{ fft1(x) }
\arguments{
\item{x}{ a numerical vector }
}
\value{ Calculates
\deqn{ J( k/n ) = \sum_{t=1}^n x_t \; e^{-i 2 \pi (k/n) t} }
for \eqn{k = 0, ..., (n-1)}, where n=length(x)}
\examples{}
\author{Peter F. Craigmile}
\keyword{}
