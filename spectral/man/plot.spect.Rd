\name{plot.spect}
\alias{plot.spect}         
\title{plot.spect}         
\description{}
\usage{
	\method{plot}{spect}(x, taper.num, xlab="Frequency", ylab="Spectrum",   
   log.freq=FALSE, add=FALSE, type="l", \ldots)}      
\arguments{
 \item{x}{ an object of class 'spect'} 
 \item{taper.num}{if 'x' contains many tapers, which taper to show.}   
 \item{xlab}{x label for the figure}       
 \item{ylab}{y label for the figure}       
 \item{log.freq}{if TRUE present the figure on the log frequency scale}  
 \item{add}{if TRUE add to current plot, otherwise create new plot}   
 \item{type}{line type is used by default in this figure}    
 \item{...}{further arguments that can be passed to the plot}    
}
\value{ No value is returned }       
\details{Plot the object 'sp' of class 'spect' in decibels,    
 showing frequencies on the log scale if log.freq=T.    
}
\references{}
\seealso{}
\examples{}
\author{Peter F. Craigmile}        
\keyword{} 
