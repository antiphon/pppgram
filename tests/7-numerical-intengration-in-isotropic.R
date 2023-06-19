# Isotropic part, integration 

devtools::load_all()
library(spatstat)

x <- rpoint(n=150, win = as.owin(100*c(0,1,0,2)))

i2 <- periodogram_iso(x, int_n = c(150, 200))

i <- periodogram_iso(x, int_n = c(500, 500))
plot(i)
lines(i2, col = 2)
