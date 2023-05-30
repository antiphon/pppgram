# pppgram

Periodograms for Spatial Point Patterns. This package is for stable implementations and for CRAN. 

## Supported data

At the moment the package only supports data objects of class `spatstat::ppp` with rectangular observation window via the check


```
# is my 'x' ok? The following needs to be TRUE:
is.ppp(x) & is.rectangle(x$window)
```

Done: first versions of

* [x] Bartlett's periodogram 2D
* [x] Multitapered 2D
* [x] isotropic 2D
* [x] basic methods (plot, print, as.im, rotmean, as.data.frame)

 