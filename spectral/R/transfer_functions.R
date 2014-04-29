


transfer <- function (fs, filter) {
  ## ======================================================================
  ## Calculate the transfer function for the 'filter' at the
  ## frequencies 'fs'.
  ## ======================================================================
  
  ws   <- 2 * pi * 0:(length(filter)-1)
  term <- outer(ws, fs)

  colSums(filter * cos(term)) + 1i * colSums(filter * sin(term))
}


squared.gain <- function (fs, filter, multiplier=1, use.C=TRUE) {
  ## ======================================================================
  ## Calculate the squared gain function for the 'filter' at the
  ## frequencies 'multiplier * fs'.
  ## ======================================================================

  if (use.C) {

    n <- length(fs)
    
    .C("squared_gain",
       sq.gain=double(n), as.double(fs), as.integer(n),
       as.double(multiplier),
       as.double(filter), as.integer(length(filter)),
       PACKAGE="spectral")$sq.gain
    
  } else {
    
    ws   <- 2 * pi * multiplier * 0:(length(filter)-1)
    term <- outer(ws, fs)
    
    colSums(filter * cos(term))^2 + colSums(filter * sin(term))^2
  }
}




## some tests

if (F) {
  
library(dwt)

a <- c(1,2,3,2,1)
transfer(0.2, a)

g <- dwt.filter("LA8")$scaling
fs <-  seq(0, 0.5, len=1000)
#plot(fs, squared.gain(fs, a))

print(system.time(for (i in 1:1000) a<-squared.gain(fs, g, use.C=T)))

print(system.time(for (i in 1:1000) a<-squared.gain(fs, g, use.C=F)))


#  ff <- sapply(0:4, function (j)
#             squared.gain(2^j*fs, dwt.filter("LA8")$scaling))


## if (F) {
## dwt.squared.gains <- function (fs, scaling, J, wavelet=T) {
##   if (wavelet) 
## #    out <- squared.gain(0.5-2fs, dwt.filter("LA8")$scaling)
##   if (J>1)
##     out <- cbind(sapply(0:(J-2), function (j)
##                         squared.gain(2^j*fs, dwt.filter("LA8")$scaling))
## }
}
