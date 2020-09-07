#' Multi-tapered Periodogram for 2D Point Patterns
#' 
#' Using sine-tapers.
#' 
#' @param x ppp-object, i.e. spatstat's point pattern, with rectangular window.
#' @param omega The wavenumber grid
#' @param ... ignored
#' @param M max steps in Fourier wavenumber grid when omega missing.
#' @param P Number of orthogonal sine-tapers to apply will P-squared. (default: 2)
#' @param debias_early Do debiasing before squaring (default: TRUE if omega given)? Very important as the Fourier grid is non-trivial.
#' 
#' @details The multitaper estimator is the average of P^2 taper estimators I_m, where
#' I_p is a tapered estimator using the product of 1D-tapers of the form sin(p * pi * x).

#' @seealso \code{\link[=periodogram]{periodogram}} for the original Bartlett's non-tapered peridogram; 
#'
#' @import Rcpp spatstat
#' @export

periodogram_mt <- function(x, 
                        omega = NULL,
                        ...,
                        M = 25,
                        P = 2,
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
  loc <- as.matrix(cbind( x$x, x$y) )
  #
  #
  # Here we check if early debiasing is required:
  if(debias_early) {
    # the taper FTs
    H <- multitaper_FT_theo(sl, omega = o, P = P)
    H <- H * sqrt(Vol)/2
  }
  # precompute arguments
  OT <- exp( -1i * 2 * pi * loc %*% t(o) )
  #browser()
  # precompute centering and scaling for tapering
  locs <- pi * (t(loc) + sl/2)/sl
  #
  # The tapering coefficients
  pvec <- as.matrix( expand.grid(1:P, 1:P) )
  # Compute tapered periodograms, one for each taper-parameter p:
  Ip <- sapply(1:nrow(pvec), function(p) {
    weight <- sin(locs[1,] * pvec[p,1]) * sin(locs[2,] * pvec[p,2])
    DFT <- colSums(weight * OT)
    if(debias_early) {
      DFT <- DFT - lambda * H[,p]
    }
    Re(DFT)^2 + Im(DFT)^2
  })
  weight0 <- 2/sqrt(Vol) # normalisation
  Ip <- Ip * weight0^2
  # average to get the multitapered
  sdf <- rowMeans(Ip)
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
               type = paste0(P, "-sine tapered"),
               lambda = lambda,
               debias_early = debias_early)
  # have a constructor?
  class(out) <- c("pppgram", is(out))
  
  out
}





#' Fourier Transform of Sine Multitaper over Rectangular Window
#'
#' @param sl sidelengths of the window
#' @param omega wavenumber grid
#' @param P number of tapers
#'
#' @details
#' Compute prod int_[0,lj]sin(pi*pj*xj/lj)exp(-2*pi*i*xj*wj)dxj where lj=sl[j] is jth sidelength
#' and wj is jth coordinate of omega, pj is taper number.
#' Then prod over dimensions (2 most likely).
#'
#' The result is normalised with 2^(d/2)/sqrt(|win|) so that ||h^2||_2=1.
#'
#'
#' @return
#' Matrix with each column corresponding to H_p(.) where p=(1,1),(2,1),...,(P,P) and . in omega.
#' @export

multitaper_FT_theo <- function(sl, omega, P = 1) {
  d <- 2 # only 2D for now
  # all combinations of taper coefficients
  pg <- as.matrix(expand.grid(1:P, 1:P))
  # precompute arguments
  a <- t(2*pi*omega) * sl
  CO <- exp(1i * a/2) * sl * sqrt(2/sl)
  # for one p-vec
  Hp <- function(p) {
    pik <- pi * p
    px  <- CO *  pik * (exp(-1i*a) * cos(pik) - 1)/(a^2-pik^2)
    px[a ==  pik] <- -.5i * CO[a ==  pik]
    px[a == -pik] <- +.5i * CO[a == -pik]
    #browser()
    apply(px, 2, prod)
  }
  out <- apply(pg, 1, Hp)
  out
}

