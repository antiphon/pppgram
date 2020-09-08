# Isotropic periodogram

devtools::load_all()
library(spatstat)
set.seed(13)
x <- rpoint(n=150, win = as.owin(100*c(0,1,0,2)))


if(1){
taper_aa <- 25
t0 <- microbenchmark::microbenchmark( p <-  periodogram_iso(x, debias  = TRUE, taper_a = taper_aa),
                                      pf <- periodogram_iso(x, debias = FALSE, taper_a = taper_aa),
                                      times = 3)
print(t0)

plot(pf, ylim = c(-1, 1)*5)
lines(p, col = 2)

p0 <- ppspectral::iso_periodogram(x, t = p$t)
lines(p0, lty=2)
# a bit slow
}

# 
if(0) {
  tv <- seq(0.1, 1, l = 11)
  en00 <- envelope(x, periodogram_iso, debias = FALSE, taper_a = 0, t = tv)
  en0 <- envelope(x, periodogram_iso, debias = FALSE, t = tv)
  en1 <- envelope(x, periodogram_iso, debias = TRUE, t = tv)
  
  par(mfrow= c(3,1))

  library(ggplot2)
  library(dplyr)
  re<-as.data.frame(en00) %>% mutate(type="notaper biased") %>% 
    bind_rows(as.data.frame(en0) %>% mutate(type="taper biased"),
              as.data.frame(en1) %>% mutate(type="taper") ) %>%
                ggplot() +
                geom_ribbon(aes(r, ymax=hi, ymin=lo, fill = type), alpha=.2)+
    geom_line(aes(r, mmean, col = type)) + ylim(c(0,3))
  print(re)
}