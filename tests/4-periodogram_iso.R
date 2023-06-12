# Isotropic periodogram

devtools::load_all()
library(spatstat)
set.seed(13)
x <- rpoint(n=150, win = as.owin(100*c(0,1,0,2)))


if(0){
taper_aa <- 25
tvec <- seq(0, .2, l = 50)[-1]
t0 <- microbenchmark::microbenchmark( p <-  periodogram_iso(x, debias  = TRUE, taper_a = taper_aa, t = tvec),
                                      pf <- periodogram_iso(x, debias = FALSE, taper_a = taper_aa, t = tvec),
                                      times = 3)
print(t0)

plot(pf, ylim = c(-1, 1)*5, col = "red")
lines(p, col = "green")

p0 <- ppspectral::iso_periodogram(x, t = p$t)
lines(p0, lty=2)
# a bit slow
}

# 
if(1) {
  tv <- seq(0, .1, l = 31)
  en00 <- envelope(x, periodogram_iso, debias = FALSE, taper_a = 0, t = tv)
  en0 <- envelope(x, periodogram_iso, debias = FALSE, t = tv)
  en11 <- envelope(x, periodogram_iso, debias = TRUE, taper_a = 0, t = tv)
  en1 <- envelope(x, periodogram_iso, debias = TRUE, t = tv)
  
  par(mfrow= c(3,1))

  library(ggplot2)
  library(dplyr)
  
  re<-as.data.frame(en00) %>% mutate(type="notaper biased") %>% 
    bind_rows(as.data.frame(en0) %>% mutate(type="taper biased"),
              as.data.frame(en1) %>% mutate(type="taper debiased"),
              as.data.frame(en11) %>% mutate(type="notaper debiased") ) %>%
                ggplot() +
                geom_ribbon(aes(r, ymax=hi, ymin=lo, fill = type), alpha=.2)+
    geom_line(aes(r, mmean, col = type))  +
    coord_cartesian(ylim = c(0, 50))
  print(re)
}