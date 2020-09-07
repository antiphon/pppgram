#' Bartlett's Periodogram for 2D Point Patterns
#'
#' @param x ppp-object, i.e. spatstat's point pattern, with rectangular window.
#' @param omega The wavenumber grid
#' @param ... ignored
#' @param M max steps in Fourier wavenumber grid when omega missing.
#' @param debias_early Do debiasing before squaring (default: TRUE if omega given)? No difference for square windows and Fourier
#' grids, but great reduction in bias for non-square windows (future) or non-Fourier grids.
#' 
#' @details
#'
#' We estimate the spectral density defined by Bartlett 1964 for a stationary process with 
#' intensity $\lambda$ and a pair correlation $g$ as
#'
#' \deqn{\mathcal{F}(\omega) = \lambda + \lambda^2\int_{R^d}[g(z)-1]e^{-2\pi i\omega^T z}dz}
#'
#' using a "periodogram", viz. a quadratic transformation of the Fourier transform of the data. 
#'
#' Note how we parameterise the Fourier transform using wavenumbers.
#'
#' No scaling of the pattern is conducted, only centering to 0 for computational simplicity
#' (does not affect the estimate as pattern assumed stationary).
#'
#' The wavenumbers that the spectrum is estimated are given by `omega`:
#'
#' 1) If omega is a vector, we expand it to d-dimensional grid using expand.grid. 
#'
#' 2) If omega is a column matrix of dimension m x d, each row vector is
#' interpreted as a wavenumber.
#'
#' If omega is missing, we will expand -M:M/l where l is the minimum sidelength of the data window.
#' 
#' If `debias_early = TRUE`, the systematic bias component due to non-centered stochastic 
#' process is cancelled before squaring the discrete Fourier transform.
#'
#' @seealso \code{\link[=periodogram_mt]{periodogram_mt}} for multi-tapered peridograms; 
#' \code{\link{smoothen.periodogram}} for operations on the estimated periodograms.
#'
#' @import Rcpp spatstat
#' @export

periodogram <- function(x, 
                        omega = NULL,
                        ...,
                        M = 25,
                        debias_early = !is.null(omega)) {
  # check input
  x <- check_pp(x)
  # Check wavenumbers
  sl <- sidelengths(Window(x))
  Vol <- area(x)
  if(is.null(omega)) {
    slm <- min( sl )
    omega <- -M:M / slm
  }
  o <- omega
  if( is.vector(o) ){
    o <- as.matrix( expand.grid(w1 = o, w2 = o, KEEP.OUT.ATTRS = FALSE) )
  }
  if( !is.matrix( o )) stop("Omega can be either a vector or a matrix.")
  #
  # input ok.
  #
  # center
  x <- center_pp(x)
  lambda <- intensity(x)
  tloc <- as.matrix(rbind( x$x, x$y) )
  # Here we only do the 'raw' or 'untapered' periodogram.
  weight <- 1/sqrt(Vol) # raw
  #
  # First, compute DFT. 
  DFT <- rowSums(weight * exp( -1i * 2 * pi * o %*% tloc))
  #
  # Debias before modulus?
  if(debias_early) {
    # Fourier transform of a rect window. Assuming symmetric, shouldn't matter.
    arg <- t(o) * (sl * pi) 
    sinc <- sin(arg)/(arg)
    is1 <- (arg == 0)
    sinc[ is1 ] <- 1
    Hv <- apply(sinc, 2, prod) * sqrt(Vol)
    bias <- Hv * lambda
    DFT <- DFT - bias
  }
  #
  # Periodogram is the modulus^2,
  sdf <- Re(DFT)^2 + Im(DFT)^2
  #
  # Drop (0,0,..) value
  zero <- which( apply(o==0, 1, all) )
  zerov <- sdf[zero]
  sdf[ zero ] <- NA
  #
  stops <- lapply(1:2, function(d) sort(unique(o[, d])))
  #
  # Compile output
  colnames(o) <- paste0("w", 1:2)
  #
  out <- list( sdf_estimate = sdf,
               omega = o,
               stops = stops,
               data = list(window = Window(x),
                           n = x$n,
                           lambda = lambda),
               zeroidx = zero,
               zerov=zerov,
               dim = 2,
               type = "Bartlett's",
               lambda = lambda,
               debias_early = debias_early)
  # have a constructor?
  class(out) <- c("pppgram", is(out))
  
  out
}



