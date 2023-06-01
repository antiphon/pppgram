# pppgram

Periodograms for Spatial Point Patterns. This package is for stable implementations and for CRAN. 

The package implements the main ideas of the paper

> T Rajala, S Olhede, J Grainger and DJ Murrell (2023): What is the Fourier Transform of a Spatial Point Process?, *IEEE Transactions on Information Theory*, https://doi.org/10.1109/TIT.2023.3269514

## Supported data

At the moment the package only supports data objects of class `spatstat::ppp` with rectangular observation window. 

Check this first:

```
# is my data point pattern 'x' ok? The following needs to be TRUE:
is.ppp(x) & is.rectangle(x$window)
```

## Functionality

The package implements

* [x] Bartlett's periodogram 2D
* [x] Multitapered 2D
* [x] isotropic 2D
* [x] basic methods (plot, print, as.im, rotmean, as.data.frame)


 