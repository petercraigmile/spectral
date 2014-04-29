## Updated   : pfc@stat.osu.edu, Jun 2009.


create.taper <- function (taper.name, N, ...) {
  ## ======================================================================
  ## Purpose: Create a taper called of the form 'taper.name', of length 'N'.
  ##          Extra parameters can be passed to the taper function via ...
  ## ======================================================================

  all.tapers <- c("rectangular", "cosine", "Hanning", "dpss", "sine")
                   
  pos <- pmatch(tolower(taper.name), tolower(all.tapers))

  if (is.na(pos))
    stop(paste("Taper '", taper.name,"' unknown.", sep=""))
  else
    do.call(paste(all.tapers[pos], ".taper", sep=""), list(N, ...))
}



renormalize.taper <- function (taper) {
  ## ======================================================================
  ## Purpose:  Renormalize the taper so that it has unit norm
  ## ======================================================================
  
  taper / sqrt(sum(taper^2))
}



rectangular.taper <- function (N)
  ## ======================================================================
  ## Purpose   : Calculates the rectangular taper of length 'N'.
  ## Reference : Page 209 of Percival and Walden, 1993.
  ## ======================================================================
{
  if (N<=0)
    stop ("N must be positive")
  
  rep(1/sqrt(N), N)
}



cosine.taper <- function (N, p)
  ## ======================================================================
  ## Equation (209) of Percival and Walden
  ## ======================================================================
{
  if (N<=0)
    stop ("N must be positive")
  if (missing(p))
    stop("the p parameter is missing")
  if ((p<0) || (p>1))
    stop("p must lie in the interval [0,1]")
  
  K <- floor(p * N)
  halfK <- floor(K/2)
  taper <- rep(1,N)

  if ((halfK > 0) && (N>1)) {
    as <- 0.5 * (1 - cos(2 * pi * seq(halfK) / (K+1)))
    taper[N:(N-halfK+1)] <- taper[1:halfK] <- as
  }

  renormalize.taper(taper)
}



Hanning.taper <- function (N)
  ## ======================================================================
  ## Purpose   : Calculates the Hanning taper of length 'N'.
  ##             (The Hanning taper is the same as a 100% cosine taper).
  ## Reference : Page 210 of Percival and Walden, 1993.
  ## ======================================================================
{
  cosine.taper(N, 1)
}




sine.taper <- function (N, n.tapers) {
  ## ======================================================================
  ## Purpose   : Calculates n.tapers sine tapers of length 'N'.
  ## Reference : Chapters ?? and ?? of Percival and Walden, 1993.
  ## ======================================================================

  if (N<=0)
    stop ("N must be positive")

  tapers <-sqrt(2/(N+1)) * sin( pi * outer(seq(N), 1:(n.tapers)) / (N+1) )
  
  if (n.tapers==1)
    as.vector(tapers)
  else
    tapers
}



dpss.taper <- function (N, n.tapers, NW=4, W=NW/N)
  ## ======================================================================
  ## Purpose   : Calculates n.tapers DPSS tapers of length 'N',
  ##             with window length 'W'.
  ## Reference : Chapters ?? and ?? of Percival and Walden, 1993.
  ## Requires  : Fortran functions 'tridib' and 'tinvit"
  ## ======================================================================
  
{
  
  if (missing(n.tapers))
    stop("The number of tapers, 'n.tapers', must be specified")
  
  if (N<=0)
    stop ("N must be positive")
  else if (N==1)
    matrix(1, 1, n.tapers)
  else {
    tpW  <- cos(2 * pi * W)
    ts   <- 0:(N-1)
    js   <- (N - 1 - 2*ts)
    diag <- (js/2)^2 * tpW
    off  <- ts*(N-ts)/2
    off2 <- off^2
    z    <- double(N)
    
    out <- .Fortran("tridib",
                    as.integer(N), as.double(-1.0),
                    as.double(diag), as.double(off), as.double(off2),
                    as.double(0), as.double(0),
                    as.integer(N-n.tapers+1), as.integer(n.tapers),
                    eigen = double(n.tapers), ind = integer(N),
                    as.integer(0), z, z, PACKAGE="spectral")
    
    out2 <- .Fortran("tinvit",
                     as.integer(N), as.integer(N),
                     as.double(diag), as.double(off),
                     as.double(off2),
                     as.integer(n.tapers),
                     as.double(out$eigen), as.integer(out$ind),
                     tapers = double(N * n.tapers),
                     as.integer(0),
                     z, z, z, z, z,
                     PACKAGE="spectral")$tapers

    if (n.tapers==1)
      tapers <- cbind(out2)
    else
      tapers <- matrix(out2, N)[,n.tapers:1]

    ## Now make sure that each taper follows the
    ## Slepian convention (see P.379 of Percival and Walden (1993)).
    Slepian.signs <- sign(sapply(1:n.tapers, function (k)
                          if (k%%2==0) sum(js*tapers[,k]) else sum(tapers[,k])))
    t(Slepian.signs * t(tapers))
  }
}
