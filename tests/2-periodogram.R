# Standard Bartlett periodogram

devtools::load_all()

x <- data("bei", package = "spatstat.data")  |> get()

p <- periodogram(x)
plot(p)



# Check second.
p2 <- periodogram(x, omega = seq(-1,1, l= 101)/100, debias_early = TRUE)
plot(p2)


