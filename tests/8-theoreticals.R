# 
devtools::load_all()
library(dplyr)

k <- seq(0, 10, l = 51)[-1]


# Thomas
v <- sdf_thomas(k, 100, 10, sigma = 0.05)
x <- spatstat.random::rThomas(kappa = 100/10, mu = 10, 
                              0.05, nsim = 100)
xv <- lapply(x, \(b) periodogram_iso(b, t = k, taper_a = 0, 
                                     debias=TRUE,
                                     int_n = c(300, 300)) |> as.data.frame()) |> 
  do.call(what = rbind)

xvm <- xv |> group_by(t) |> 
  summarise(m = mean(I))


# MaternII
v2 <- sdf_maternII(k, 100, 0.05)
x2 <- spatstat.random::rMaternII(kappa = maternII_generator_intensity(100, 0.05), 
                                0.05, nsim = 100)
x2v <- lapply(x2, \(b) 
              periodogram_iso(b, t = k, taper_a = 0,
                              int_n = c(300, 300)) |> as.data.frame()) |> 
  do.call(what = rbind)

x2vm <- x2v |> group_by(t) |> 
  summarise(m = mean(I))




# Plot

par(mfrow=c(2,2))

plot(k, v, type = "l", ylim = range(v, xvm$m))
points(xvm$t, xvm$m)
abline(v = 1)

plot(k, v2, type = "l", ylim = c(0,1))
points(x2vm$t, x2vm$m, col = 3)
abline(v = 1)

# the bias?
bias <- pppgram:::pppgram_iso_bias(c(1,1), k, 0, n = c(500, 1000))
plot(k, bias, log="y")
abline(v = 1)
# 
# 1/window size is the failure point

