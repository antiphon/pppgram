#' Isotropic Periodogram Estimator
#'
#' @param x point pattern ppp-object with a rectangular window
#' @param t The wavenumber magnitude grid ("frequencies")
#' @param ... ignored
#' @param debias Apply debiasing after estimation? Will reduce bias near 0, but might lead to a negative estimates.
#' @param normalise Divide by lambda? Poisson process will then be in theory constant 1.
#' @param taper_a Use a squared-exponential taper with range taper_a>0? Default: 25.
#' 
#' @details Basically, the Hankel-transform of the data. Note that the taper-parameter (taper_a) is related to the inverse of 
#' the taper variance, so a smaller number means less tapering. And taper_a = 0 equates to no tapering.
#' 
#' 
#' @importFrom spatstat im
#' @useDynLib pppgram
#' @import Rcpp
#' @export

periodogram_iso <- function(x, t, ..., debias = TRUE, normalise=TRUE, taper_a = 25
                            ) {
  x <- check_pp(x)
  # Check wavenumbers
  sl <- sidelengths(Window(x))
  Vol <- area(x)
  lambda <- intensity(x)
  #
  if(missing(t)) {
    slm <- min( sl )
    t <-seq(0, 25 / slm, length.out = 51) [-1]
  }
  ##################################################################
  # Compute:
  # periodogram
  value <- c_iso_sum(cbind(x$x, x$y), t, sl, taper_a)
  sdf <- lambda + 2 * value/Vol
  
  if(debias) {
    bias <- ppgram_iso_bias(sl, t, taper_a)
    lambda2 <- lambda^2 - lambda/Vol
    #browser()
    sdf <- sdf - lambda2 * bias
  }
  #
  if(normalise) sdf <- sdf/lambda
  #
  type <- "Isotropic periodogram"
  if(taper_a > 0) type <- sprintf("exp(x^2)-tapered %s", type)
  # compile output
  out <- fv( data.frame(t = t, I = sdf),
             valu = "I" ,
             argu = "t",
             desc = c("wave number", type),
             labl = c("t", "%s(t)"),
             fname = "I")
  attr(out, "is_periodogram") <- TRUE
  attr(out, "normalised") <- normalise
  attr(out, "W") <- Window(x)
  attr(out, "lambda") <- lambda
  attr(out, "taper_a") <- taper_a
  attr(out, "type") <- type
  out
}



#' Isotropic Periodogram Bias
#' 
#' Assuming a rectangular window, 2D
#' 
#' @param sl window sidelenghts
#' @param t wavenumber magnitudes
#' @param taper_a parameter of the squared-exponential taper. if 0, no taper.
#' @param n grid lengths for numerical integration, order (ang, r)
#' 
#' @export
ppgram_iso_bias <- function(sl, t, taper_a, n = c(150, 200)) {
  V <- prod(sl)
  # Sq-exp taper.
  h2d <- function(x) exp( -taper_a/4 * ((x[1,]/sl[1])^2 + (x[2,]/sl[2])^2 )  ) / V
  # choose.
  h <- if(taper_a > 0) h2d else function(...) 1/V 
  #
  ## Numerical integration. A bit slow, need to optimsize.
  #
  angs <- seq(0, pi, length = n[1])[-1] # symmetric, half is enough
  da <- angs[2] - angs[1]
  rmax <- sqrt(sum(sl^2)) # afterwhich setcov = 0
  rgrid <- seq(0, rmax, length = n[2])[-1]
  dr <- rgrid[2] - rgrid[1]
  #
  # taper-weighted isotropized set_cov i.e. average of h*set_cov over directions.
  funh <- function(r) { # rotave.
    u <- rbind(r*cos(angs), r*sin(angs))
    set_cov <- pmax(apply(sl - abs(u), 2, prod), 0)
    2 * sum( h(u) * set_cov) * da / (2*pi)
  }
  funhr <- sapply(rgrid, funh) # need to compute only once
  Han <- sapply(t, function(ti) 
    sum( besselJ(2*pi*ti*rgrid, 0) * rgrid * funhr ) * dr
  )
  Han  * ( 2 * pi )
}

#' 
#' #' Isotropised Set Covariance of a Rectangle
#' #' 
#' #' @param sl sidelengths
#' #' @param r argument
#' #' 
#' #' @export
#' ppgram_iso_set_cov <- function(sl, r) {
#'   if(length(r) > 1) return(sapply(r, ppgram_iso_set_cov, sl=sl))
#'   aa <- sl[1]
#'   bb <- sl[2]
#'   a <- min(aa,bb)
#'   b <- max(aa,bb)
#'   A <- a*b # area
#'   d <- b/a #ratio
#'   x <- r/sqrt(A/d)
#'   V <- A/pi
#'   if(x <= 1) {
#'     V * ( pi - 2*x - 2*x/d + x^2/d )
#'   }
#'   else if(x <= d){
#'     u <- sqrt(x^2-1)
#'     V * ( 2*asin(1/x) - 1/d - 2*(x-u) )
#'   }
#'   else if( x < sqrt(d^2+1)){
#'     u <- sqrt(x^2-1)
#'     v <- sqrt(x^2-d^2)
#'     V * ( 2*asin((d-u*v)/x^2) + 2*u + 2*v/d - d - (1+x^2)/d )
#'   }
#'   else 0
#' }
#' 
#' 
