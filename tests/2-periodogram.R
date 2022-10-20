# Standard Bartlett periodogram

devtools::load_all()

x <- bei

p <- periodogram(x)
plot(p)



# Check second.
p2 <- periodogram(x, omega = seq(-1,1, l= 101)/100)
plot(p2)


