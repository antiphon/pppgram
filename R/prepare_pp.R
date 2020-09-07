#' Check Input Pattern Format
#'
#' @param x pattern candidate
#' @details Only support spatstat::ppp objects with rectangular windows.
#' 
#' @importFrom spatstat is.ppp is.rectangle
#' @export

check_pp <- function(x){
  if(!is.ppp(x)) stop("only ppp-objects supported.")
  window <- x$window
  if(!is.rectangle(window)) stop("only rectangle windows supported.")
  x
}

#' Center Data
#' 
#' @param x data passing check_pp
#' 
#' @importFrom spatstat affine.ppp centroid.owin

center_pp <- function(x) {
  c0 <- -unlist( centroid.owin(x) )
  affine.ppp(x, vec = c0)
}