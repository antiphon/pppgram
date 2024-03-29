% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/periodogram_mt.R
\name{periodogram_mt}
\alias{periodogram_mt}
\title{Multi-tapered Periodogram for 2D Point Patterns}
\usage{
periodogram_mt(
  x,
  omega = NULL,
  ...,
  M = 25,
  P = 2,
  debias_early = !is.null(omega)
)
}
\arguments{
\item{x}{ppp-object, i.e. spatstat's point pattern, with rectangular window.}

\item{omega}{The wavenumber grid}

\item{...}{ignored}

\item{M}{max steps in Fourier wavenumber grid when omega missing.}

\item{P}{Number of orthogonal sine-tapers to apply will P-squared. (default: 2)}

\item{debias_early}{Do debiasing before squaring (default: TRUE if omega given)? Very important as the Fourier grid is non-trivial.}
}
\description{
Using sine-tapers.
}
\details{
The multitaper estimator is the average of P^2 taper estimators I_m, where
I_p is a tapered estimator using the product of 1D-tapers of the form sin(p * pi * x).
}
\seealso{
\code{\link[=periodogram]{periodogram}} for the original Bartlett's non-tapered peridogram;
}
