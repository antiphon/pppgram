# Ring-pass filter

devtools::load_all()

x <- data("bei", package = "spatstat.data")  |> get()

t0 <- system.time( p  <- ringpass(x, nx = 2^8, a = 0, b = 0.01, fft=FALSE) )
t1 <- system.time( pf <- ringpass(x, nx = 2^8, a = 0, b = 0.01, fft=TRUE) )


par(mfrow=c(2,1))
plot(p, col = co <- pals::brewer.rdbu(20))
points(x, cex = .1)

plot(pf, col = co)
points(x, cex = .1)



print(rbind(t0,t1))