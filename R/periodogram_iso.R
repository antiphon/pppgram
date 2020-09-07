#' Isotropic Periodogram Estimator
#'
#' @param x point pattern ppp-object with a rectangular window
#' @param t The wavenumber magnitude grid ("frequencies")
#' @param ... ignored
#' @param debias Apply debiasing after estimation? Will reduce bias near 0, but might lead to a negative estimates.
#' @param normalise Divide by lambda? Poisson process will then be approximately constant 1.
#' @param taper_a Use a squared-exponential taper with range taper_a>0? Default: 25.
#' @details Basically, the Hankel-transform of the data. 
#' 
#' @importFrom spatstat im
#' @useDynLib pppgram
#' @import Rcpp
#' @export

periodogram_iso <- function(x, t, ..., debias = TRUE, normalise=TRUE, taper_a = 25) {
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
  # d <- as.matrix(dist(cbind(x$x, x$y)))
  # d <- d[lower.tri(d, diag = FALSE)]
  # rr <- 2*pi*t
  # o <- outer(d, rr)
  # value <- colSums(besselJ(o, 0)) 
  # 
  # periodogram
  value <- c_iso_sum(cbind(x$x, x$y), t, sl, taper_a)
  sdf <- lambda + 2 * value/Vol 
  
  if(debias) {
    bias <- 
  }
  
  
  #
  if(normalise) sdf <- sdf/lambda
  type <- "Isotropic periodogram"
  if(taper_a > 0) type <- sprintf("exp(x^2)-tapered %s", type)
  # compile output
  out <- fv( data.frame(t = t, I = sdf),
             valu = "I" ,
             argu = "t",
             desc = c("wave number", type),
             labl = c("t", "%s(t)"),
             fname = type)
  attr(out, "is_periodogram") <- TRUE
  attr(out, "normalised") <- normalise
  attr(out, "W") <- Window(x)
  attr(out, "lambda") <- lambda
  attr(out, "taper_a") <- taper_a
  out
}


ppgram_iso_bias <- function(sl, taper_a) {
  BIASES <<- list()
  f1d <- function(x, a, l) exp( -a * x^2 / (2*l)^2)
  f2d <- function(x, a, l) f1d(x[,1], a, l[1])*f1d(x[,2], a, l[2]) / prod(l)
  precompute_bias <- function(nvec, tvec, apar = 25) {
    # compute for each window
    fun <- function(r, t, W, sl) iso_set_cov(W, r) * 
      f2d(cbind(0,r), apar, sl) * besselJ(2*pi*t*r, 0) * r
    for(n in nvec){
      L <- sqrt(n/lambda)
      W <- square( L )
      sl <- c(L, L)
      BiasIsoCov <- sapply(tvec, 
                           function(t) 
                             integrate(fun, 0, sqrt(2) * L, t = t, W=W, sl = sl, 
                                       subdivisions = 1000)$value ) # 
      BIASES[[as.character(n)]] <<- BiasIsoCov * 2 * pi
}

