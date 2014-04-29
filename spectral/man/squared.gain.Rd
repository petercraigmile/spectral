
\name{squared.gain}
\alias{squared.gain}
\title{Calculate the squared gain function}
\description{}
\usage{squared.gain(fs, filter, multiplier=1, use.C=TRUE)}
\arguments{
  \item{fs}{ a numerical vector of frequencies }
  \item{filter}{ a vector of filter coefficients }
  \item{multiplier}{ a multipler for the frequencies }
  \item{use.C}{ if use.C is TRUE use C code to calculate; otherwise using R code }
  }
\value{ A numerical vector of squared gain function values, the same length as 'fs'. }
\details{Calculate the squared gain function for the 'filter' at the frequencies 'multiplier * fs'.}
\references{}
\examples{}
\author{Peter F. Craigmile}
\keyword{}

