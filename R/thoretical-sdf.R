#' Spectral density function of Thomas process
#' 
#' Isotropic sdf of the (modified) Thomas point process
#' 
#' @param k wavenumbers
#' @param lambda Target process intensity (not the generator intensity)
#' @param mu Offspring intensity
#' @param sigma Cluster stdev
#' 
#' @details The parametrisation takes the target intensity, not the generator intensity (usually `kappa`). 
#' But it does not really matter here, scaling by intensity removes connection to 
#' lambda or kappa, only mu remains.
#' 
#' @return Isotropic spectral density function scaled by intensity
#' 
#' @export
sdf_thomas <- function(k, lambda, mu, sigma) {
  if(is.null(nrow(k))) k <- cbind(k)
  kappa <- lambda/mu
  u <- rbind(k)*2*pi
  d <- rowSums(u^2)
  1 + exp(-sigma^2 * d) * mu 
}

#' Spectral density function of Matern type II process
#' 
#' @param k wavenumbers
#' @param lambda Target intensity
#' @param hc hard core range
#' @param nint numerical integration scheme steps
#' 
#' @details The parametrisation takes the desired intensity `lambda`, not the generator intensity.
#' The generator intensity will be solved from `lambda` and `hc`, 
#' and will fail if the combination is not possible.  
#' 
#' @return spectral density function scaled by intensity
#' 
#' @export
sdf_maternII <- function(k, lambda, hc, nint = 2000) {
  if(is.null(nrow(k))) k <- cbind(k)
  k <- rbind(k)
  r <- seq(0, 3*hc, l = nint)
  w <- sqrt(rowSums(k^2))
  g <- pcf_maternII(r = r, target_intensity = lambda, hc = hc)
  pcf_to_sdf(g, r, w, lambda)
}

#' Matern type II pcf
#'
#' Isotropic pair correlation function of Matern type II
#'
#' @param r spatial scales
#' @param target_intensity intensity of the process (not the generator)
#' @param hc hard-core distance
#' @details Generator, or dominant Poisson to be thinned, intensity will be solved 
#' and an error will be given if target intensity is too high given hc distance.
#' 
#' @export
pcf_maternII <- function(r = seq(0, 2, l = 50), 
                            target_intensity = 50,
                            hc = 0.05) {
  lam <- maternII_generator_intensity(target_intensity, hc)
  V <- pi * hc^2
  rr <- r
  rr[ rr > 2*hc] <- 2*hc
  z <- 2 * hc^2 * acos(rr/(2*hc)) - (rr/2) * sqrt(4 * hc^2 - rr^2)
  G <- 2 * V - z
  a <- 2 * G * (1-exp(-lam * V)) - 2 * V * (1- exp(-lam * G))
  b <- V * G * (G - V)
  e <- (a/b) / target_intensity^2
  e[r < hc] <- 0
  e
}

#' Solve for the MaternII generator intensity
#' 
#' @param target_intensity Desired intensity in the end
#' @param hc the hard core parameter
#' @export

maternII_generator_intensity <- function(target_intensity, hc){
  V <- pi * hc^2
  if(1 < target_intensity * V) {
  mh <- 1/sqrt(pi * target_intensity)
  ml <- 1/V
  stop(paste0("max hc for target_intensity ", target_intensity, " is ", mh,"\nmax target_intensity for hc ",hc, " is ", ml) )
  }
  lam <- -log(1 - target_intensity * V)/V
  lam
}


#' Generic function to integrate isotropic pcf to sdf
#' 
#' Numerical Hankel transform of pcf to sdf
#'
#' @param g pcf values at r
#' @param r spatial arguments
#' @param w wavenumber arguments
#' @param lam intensity
#' @details Not the most robust Hankel transform, so use with care.
#'
pcf_to_sdf <- function(g, r, w, lam) {
  # trapezoidal rule, any grid
  arg <- 2 * pi * w 
  # bessel zero's
  ra <- outer(r, arg)
  j <- besselJ(ra, 0)
  h <- g-1
  fk <-  h * r * j
  n <- length(r)
  dr <- diff(r)
  s <- colSums((fk[-1,]+fk[-n,])*dr )/2
  1 + 2 * pi * lam * s
} 

