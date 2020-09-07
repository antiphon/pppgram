# Isotropic periodogram

devtools::load_all()
library(spatstat)
set.seed(1)
x <- rpoispp(500)


t0 <- system.time( p <- periodogram_iso(x) )
print(t0)

(p)

plot(p)
