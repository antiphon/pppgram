# Optimise the isotropic bias term computation
devtools::load_all()
sl <- c(1,1)
t <- seq(0.1, 25, l = 51)
taper_a <- 25
####
# dev code
ppgram_iso_bias2 <- function(sl, t, taper_a) {
  V <- prod(sl)
  # Numerical, slow but works for non-squares.
  # squared exponential
  h2d <- function(x) exp( -taper_a/4 * ((x[1,]/sl[1])^2 + (x[2,]/sl[2])^2 )  ) / V
  h <- if(taper_a > 0) h2d else function(...) 1/V 
  angs <- seq(0, 2*pi, l = 50)[-1]
  da <- diff(angs[1:2])
  #browser()
  rmax <- sqrt(sum(sl^2))
  rgrid <- seq(0, rmax, l = 500)
  dr <- rgrid[2] - rgrid[1]
  
  # set covariance:
  set_cov <- function(u)  pmax(apply(sl - abs(u), 2, prod), 0)
  funh <- function(r) { # rotave.
    u <- rbind(r*cos(angs), r*sin(angs))
    # integrator
    sum(h(u) * set_cov(u) * da) / (2*pi)
  }
  funhr <- sapply(rgrid, funh) # need only once
  Han <- sapply(t, function(t) 
    sum( besselJ(2*pi*t*rgrid, 0) * rgrid * funhr * dr )
  )
  Han
}




###
t0 <- microbenchmark::microbenchmark(v1 <- ppgram_iso_bias(sl, t, taper_a), 
                                     v2 <- ppgram_iso_bias2(sl, t, taper_a),
                                     times = 5)
print(t0)
print(identical(v1, v2))
