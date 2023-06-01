# Multitapered periodogram

devtools::load_all()
x <- data("bei", package = "spatstat.data") |> get()

p0 <- periodogram(x, M = 2^6)

p <- periodogram_mt(x, omega = p0$omega)

plot(spatstat.geom::listof(p0,p), col = co <- pals::brewer.brbg(20))

# Check second.
p2 <- periodogram_mt(x, omega = seq(-1,1, l= 201)/50)
p20 <- periodogram(x, omega = p2$omega)
plot(spatstat.geom::listof(p20, p2), col = co)
