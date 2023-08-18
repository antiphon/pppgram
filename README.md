# pppgram

Spectral statistics, particularly periodograms as in estimators of spectral density function, for Spatial Point Pattern data. This package is for stable implementations. 

The package implements the main ideas of the papers

> T Rajala, S Olhede, J Grainger and DJ Murrell (2023): What is the Fourier Transform of a Spatial Point Process?, *IEEE Transactions on Information Theory*, https://doi.org/10.1109/TIT.2023.3269514

and 
> J Granger, T Rajala, DJ Murrell, S Olhede (2023): Visualizing the Wavenumber Content of a Point Pattern,https://arxiv.org/abs/2306.04198

## Supported data

At the moment the package only supports data objects of class `spatstat::ppp` with rectangular observation window. 

Check this first:

```
# is my data point pattern 'x' ok? The following needs to be TRUE:
is.ppp(x) & is.rectangle(x$window)
```

## Functionality

The package implements

* [x] Bartlett's periodogram 2D (`periodogram`)
* [x] Multitapered 2D (`periodogram_mt`)
* [x] isotropic 2D (`periodogram_iso`)
* [x] ring-pass filter (`ringpass`)
* [x] basic methods (plot, print, as.im, rotmean, as.data.frame)

See the examples in the vignette.
 