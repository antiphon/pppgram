# Ring-pass filter

devtools::load_all()

x <- spatstat.random::rpoispp(lambda = 300/200^2, win = as.owin(c(-1,1,-1,1)*100 - 1000))
#x <- data("bei", package = "spatstat.data")  |> get()


t0 <- system.time( p  <- ringpass(x, nx = 2^6, a = 0, b = 0.01, fft=FALSE) )
t1 <- system.time( pf <- ringpass(x, nx = 2^6, a = 0, b = 0.01, fft=TRUE) )


par(mfrow=c(1,2))
plot(p, col = co <- pals::brewer.rdbu(20))
points(x, cex = .1)

plot(pf, col = co)
points(x, cex = .1)



print(rbind(t0,t1))