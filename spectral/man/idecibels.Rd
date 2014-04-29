
\name{idecibels}
\alias{idecibels}
\alias{idB}
\title{Convert from decibels}
\description{ Converts a numerical vector from the decibels value }
\usage{ idecibels(x) }
\arguments{
\item{x}{ a numerical vector }
}
\value{ The vector of x values converted from decibels to the original scale}
\examples{
idecibels(0)  # return 1
idecibels(10) # returns 10
idecibels(-Inf)  # return 0
}
\author{Peter F. Craigmile}
\keyword{}
