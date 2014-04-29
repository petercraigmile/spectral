\name{acvs.from.sdf}
\alias{acvs.from.sdf}
\title{Calculate the acvs via a Fourier transform of the SDF}
\description{
Calculate the ACVS from lags 0 to lag 'max.lag',
 given a function specifying the SDF, 'the.sdf'.
}
\usage{
acvs.from.sdf(max.lag, the.sdf, use.integrate = TRUE, num.weights = 200, gq, ...)
}
\arguments{
  \item{max.lag}{
calculate the acvs from 0 to the lag 'max.lag'.
}
  \item{the.sdf}{
an R function defining the SDF
}
  \item{use.integrate}{
if use.integrate is TRUE use the 'integrate' function to calculate the
 integral.  Otherwise use Gauss-Legendre quadrature from the 'statmod'
 library using the 'gq' object, if available, or creating a new object
 with 'num.weights' weights.
}
  \item{num.weights}{
the number of weights to use for the Gauss-Legendre quadrature
}
  \item{gq}{
a Gauss-Legendre quadrature from the 'statmod'  library
}
  \item{\dots}{ Other arguments are passed to 'the.sdf' function }
}
\details{
}
\value{
a vector of max.lag+1 acvs values.
}
\references{
}
\author{Peter F. Craigmile}
\note{}
\seealso{}
\examples{}
