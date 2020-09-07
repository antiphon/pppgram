#' Plot 2D Periodogram Estimate
#' 
#' @param x pppgram-object
#' @param ... passed on to plot.im of spatstat
#' 
#' @export

plot.pppgram <- function(x, ...) {
  plot( as.im(x), ... )
}

#' Convert pppgram to im-object
#' 
#' @param x pppgram-object
#' @param ... ignored
#' 
#' @importFrom spatstat plot.im as.im
#' @export
as.im.pppgram <- function(x, ...) {
  stops <- x$stops
  sdf_matrix <- matrix(x$sdf_estimate, nrow = length(stops[[1]]) )
  im( t(sdf_matrix) , stops[[1]], stops[[2]])
}

#' Convert pppgram to Data Frame
#' 
#' @param x pppgram-object
#' @param ... passed on to `data.frame`
#' 
#' @export
as.data.frame.ppgram <- function(x, ...) 
  data.frame(x$omega, estimate = x$sdf_estimate, type = x$type, ...)


#' Print for pppgram
#' 
#' @param x pppgram-object
#' @param ... not used
#' 
#' @export
print.pppgram <- function(x, ...) {
  message( sprintf("pppgram-object\ntype: %s\ndata dimension: 2\ndata points: %d\nintensity:%f", x$type, x$data$n, x$data$lambda))
  for(i in 1:2) message( sprintf("wave axis-%d: [%f, %f] (%d stops)",i, 
                                 min(x$stops[[i]]), max(x$stops[[i]]), length(x$stops[[i]])))
  message(sprintf("debiased: %s", as.character( x$debias_early )  ))
}
