# Multitapered periodogram

devtools::load_all()
library(spatstat)
x <- bei

p0 <- periodogram(x)

p <- periodogram_mt(x)

plot(spatstat::listof(p0,p))

# Check second.
p2 <- periodogram_mt(x, omega = seq(-1,1, l= 101)/20)
p20 <- periodogram(x, omega = p2$omega)
plot(spatstat::listof(p20, p2))
