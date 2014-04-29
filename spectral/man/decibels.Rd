
\name{decibels}
\alias{decibels}
\alias{dB}
\title{Convert to decibels}
\description{ Converts a numerical vector to the decibels value }
\usage{ decibels(x) }
\arguments{
\item{x}{ a numerical vector }
}
\value{ The vector of x values converted to decibels }
\examples{
decibels(1)  # return 0
decibels(10) # returns 10
decibels(0)  # return -Inf
}
\author{Peter F. Craigmile}
\keyword{}
