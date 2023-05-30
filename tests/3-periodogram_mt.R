# Multitapered periodogram

devtools::load_all()
x <- data("bei", package = "spatstat.data") |> get()

p0 <- periodogram(x)

p <- periodogram_mt(x)

plot(spatstat.geom::listof(p0,p))

# Check second.
p2 <- periodogram_mt(x, omega = seq(-1,1, l= 101)/20)
p20 <- periodogram(x, omega = p2$omega)
plot(spatstat.geom::listof(p20, p2))
