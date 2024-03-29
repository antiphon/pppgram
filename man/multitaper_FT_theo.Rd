% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/periodogram_mt.R
\name{multitaper_FT_theo}
\alias{multitaper_FT_theo}
\title{Fourier Transform of Sine Multitaper over Rectangular Window}
\usage{
multitaper_FT_theo(sl, omega, P = 1)
}
\arguments{
\item{sl}{sidelengths of the window}

\item{omega}{wavenumber grid}

\item{P}{number of tapers}
}
\value{
Matrix with each column corresponding to H_p(.) where p=(1,1),(2,1),...,(P,P) and . in omega.
}
\description{
Fourier Transform of Sine Multitaper over Rectangular Window
}
\details{
Compute prod int_[0,lj]sin(pi*pj*xj/lj)exp(-2*pi*i*xj*wj)dxj where lj=sl[j] is jth sidelength
and wj is jth coordinate of omega, pj is taper number.
Then prod over dimensions (2 most likely).

The result is normalised with 2^(d/2)/sqrt(|win|) so that ||h^2||_2=1.
}
