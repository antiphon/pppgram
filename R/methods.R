#' Plot 2D Periodogram Estimate
#' 
#' @param x pppgram-object
#' @param ... passed on to plot.im of spatstat
#' @importFrom spatstat.geom as.im
#' @export

plot.pppgram <- function(x, ...) {
  plot( as.im(x), ... )
}

#' Convert pppgram to im-object
#' 
#' @param x pppgram-object
#' @param ... ignored
#' 
#' @importFrom spatstat.geom as.im im
#' @export
as.im.pppgram <- function(x, ...) {
  stops <- x$stops
  sdf_matrix <- matrix(x$sdf_estimate, nrow = length(stops[[1]]) )
  im( t(sdf_matrix) , stops[[1]], stops[[2]])
}

#' Wrapper for Smoothing Periodogram
#' 
#' @param x pppgram-object
#' @param ... passed on to spatstat::blur
#' @details this is convenience function for doing blur(as.im(x), ...) so that the result is restored as pppgram-object.
#' @importFrom spatstat.core Smooth
#' @export
Smooth.pppgram <- function(x, ...)  {
  y <- as.im(x)
  ys <- Smooth(y, ...)
  x$sdf_estimate <- c( t(ys$v) )
  x$zerov <- x$sdf_estimate[x$zeroidx]
  x$sdf_estimate[x$zeroidx] <- NA
  x
}

#' #' Wrapper for Rotation Averaging pppgram
#' #' 
#' #' @param x pppgram-object
#' #' @param ... passed on to spatstat::rotmean
#' #' 
#' #' @export
#' rotmean.pppgram <- function(x, ...) {
#'   y <- as.im(x)
#'   rotmean(y, ...)
#' }


#' Convert pppgram to Data Frame
#' 
#' @param x pppgram-object
#' @param ... passed on to `data.frame`
#' 
#' @export
as.data.frame.pppgram <- function(x, ...) 
  data.frame(x$omega, sdf_estimate = x$sdf_estimate, type = x$type, ...)



#' Print for pppgram
#' 
#' @param x pppgram-object
#' @param ... not used
#' 
#' @export
print.pppgram <- function(x, ...) {
  message( sprintf("periodogram\ntype: %s\ndata dimension: 2\ndata points: %d\nintensity:%f", x$type, x$data$n, x$data$lambda))
  for(i in 1:2) message( sprintf("wave axis-%d: [%f, %f] (%d stops)",i, 
                                 min(x$stops[[i]]), max(x$stops[[i]]), length(x$stops[[i]])))
  message(sprintf("debiased: %s", as.character( x$debias_early )  ))
}
