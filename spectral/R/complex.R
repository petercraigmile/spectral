
h <- function (x) {
  ## ======================================================================
  ## Given a matrix or 'data.frame' 'x', 'h' returns the
  ## conjugate transpose of 'x'.
  ## ======================================================================
  
  t(Conj(x))
}



complex.chol <- function (x, upper=TRUE) {
  ## requires zpotrf.f

  if (is.matrix(x)) {
    if (nrow(x) != ncol(x)) 
      stop("non-square matrix in chol")
    n <- nrow(x)
  }
  else {
    if (length(x) != 1) 
      stop("non-matrix argument to chol")
    n <- as.integer(1)
  }

  xx <- matrix(as.complex(x), n)
  if (upper) {
    xx[lower.tri(xx)] <- 0
    type <- as.character('U')
  }
  else {
    xx[upper.tri(xx)] <- 0
    type <- as.character('L')
  }
  
  z <- .Fortran("zpotrf", type, n, x=xx, n, info=as.integer(0),
                PACKAGE="spectral")

  if (z$info!=0)
    stop("Error: matrix is not positive definite!")
  else
    z$x
}
