#' Ring-pass Filter
#' 
#' @param x point pattern, 2D, such as `ppp`
#' @param a lower annulus radius
#' @param b upper annulus radius
#' @param nx resolution, in x-dimension 
#' @param minus ignored for now
#' @param fft use FFT-based approximation? default: TRUE
#' @details
#' y-dimension resolution computed to have square pixels
#' 
#' The non-fft version is more accurate, but slow so if `fft=FALSE` keep `nx` low.
#' 
#' `fft=TRUE` will buffer the window and then clip back to avoid circularity issues (toroidality).
#' 
#' @return
#' spatstat's `im` object
#' 
#' @import spatstat.geom
#' @export
ringpass <- function(x, # point coordinates
                     a = 0.1, 
                     b = 0.2, # annuli radii
                     nx = 2^4, # resolution of output image
                     minus = FALSE, # drop data withing too close to border?
                     fft = TRUE
) {
  # ignore minus for now
  minus <- FALSE
  # 
  co    <- centroid.owin(  boundingbox.ppp(x) )
  bb    <- with(x$window, cbind(xrange, yrange))
  bb_sl <- apply(bb, 2, diff)
  ux    <- seq(bb[1,1], bb[2,1], l = nx+1)
  ny    <- round( bb_sl[2]/bb_sl[1] * nx )
  uy    <- seq(bb[1,2], bb[2,2], l = ny+1)
  # pairwise differences, only distance matters
  # if(minus > 0) {
  #   pd  <- bdist.points(x)
  #   txy  <- txy[, pd > minus, drop=FALSE]
  # }
  if(!fft) {
    txy   <- rbind(x$x - co$x, x$y - co$y)
    ugrid <- expand.grid(uy - co$y, 
                         ux - co$x )[,2:1] |> as.matrix()
    delta <- apply(ugrid, 1, \(u) colSums((u-txy)^2)   ) |> sqrt() |> t()
    # DFT of the annulus-indicator
    ha_a  <- (a/delta * besselJ(2*pi*delta*a, 1) )
    ha_b  <- (b/delta * besselJ(2*pi*delta*b, 1) )
    Ha    <- ha_b - ha_a
    Ha[is.na(Ha)] <- 0 # zero
    # then the convolution is
    Wu    <- rowSums(Ha)
  }
  else{
    # add buffers to top+right
    dx  <- diff(ux[1:2])
    dy  <- diff(uy[1:2])
    uxb <- seq(bb[1,1], bb[2,1]*2, by = dx) - co$x
    uyb <- seq(bb[1,2], bb[2,2]*2, by = dy) - co$y
    nxb <- nx*2
    nyb <- ny*2
    ugrid <- expand.grid(uyb, 
                         uxb)[,2:1] |> as.matrix()
    ud    <- rowSums(ugrid^2) |> sqrt()
    ha_a  <- (a/ud * besselJ(2*pi*ud*a, 1) )
    ha_b  <- (b/ud * besselJ(2*pi*ud*b, 1) )
    ha    <- (ha_b - ha_a)   |> matrix(ncol = nxb+1)
    Ha    <- abs( fft(ha) ) # ?
    # DFT of data pattern.
    xb    <- x
    xb$window <- as.owin(c( range(uxb), range(uyb)  ))
    Xxy   <- (as.im(xb, dimyx = c(nyb+1, nxb+1)))
    Mxy   <- unclass( Xxy$v )
    Mxy[is.na(Mxy)] <- 0 # for masked pixels etc
    Hx    <- fft(Mxy)
    # sum
    HH    <- Ha * Hx
    # then the convolution is
    Wub    <- Re( fft(HH, TRUE)/nrow(ugrid) )
    # browser()
    # Wu <- Wub
    # nx <- nxb
    # ny <- nyb
    # ux <- uxb
    # uy <- uyb
    Wu     <- Wub[1:(ny+1) + ny/2, 1:(nx+1) + nx/2]
  }
  # 
  Wum   <- matrix(Wu, ncol = nx+1)
  
  #  browser()
  im(Wum, xcol = ux, yrow = uy)
}