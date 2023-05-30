#' Rotation Average of pppgram
#'
#' @param x periodogram object
#' @param t the wave number magnitude grid to interpolate the periodogram on
#' @param t_steps number of steps in the output length scale grid
#' @param adjust smoothing bandwidth adjustment
#' @param h Box-car kernel radius
#' @param ... ignored
#' @param normalise divide by intensity? For Poisson the periodogram will then be constant 1.
#' @details
#' Box-car kernel with radius h = adjust * mean(dx, dy) where dx, dy are the step sizes of the grid used to compute the periodogram.
#'
#' @import spatstat.explore
#' @export

rotmean.pppgram <- function(x, t, t_steps = 50, adjust = 1, h, ..., normalise=TRUE) {
  # Rotation averaging: use a grid upto a circle that fits the wx x wy box
  wx <- x$stops[[1]]
  wy <- x$stops[[2]]
  w <- x$omega
  fw <- x$sdf_estimate
  tmax <- min(max(wx), max(wy))
  #
  # smoothing width: problem: frequencies might not be on a uniform grid, as
  # box side lengths might be different.
  h <- if(missing(h))  mean( diff(wx[1:2]), diff(wy[1:2]) ) * adjust else h
  
  kern <- function(d) 1*(abs(d) < h) #dnorm(d, 0, h)
  
  # amplitude, |w|, grid
  tvec <- if(missing(t)) seq(0, tmax, l = t_steps)[-1] else t # drop zero
  #
  # grid ampl's
  wd <- sqrt(rowSums(w^2))
  #
  # the differences
  diffs <- outer(wd, tvec, "-")
  #
  # kern weights
  kvals <- kern(diffs)
  kvals[ is.na(fw) ] <- NA # to get right averaging
  #
  # Estimator: sum of vals / total weights per amplitude:
  kweights <- colSums(kvals, na.rm=TRUE)
  ksums <- colSums(kvals * fw, na.rm=TRUE)
  ifw <- ksums / kweights
  #
  type <- paste(x$type, collapse=" ")
  if(normalise) {
    ifw <- ifw/x$lambda
    type <- paste(type, "normalised")
  }
  # compile output
  out <- fv( data.frame(t = tvec, I = ifw),
             valu = "I" ,
             argu = "t",
             desc = c("wave number (1/range)", "Rot. averaged periodogram"),
             fname = "bar[I]")
  attr(out, "n") <- kweights
  out
}
