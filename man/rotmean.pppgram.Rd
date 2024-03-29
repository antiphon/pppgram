% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rotmean.R
\name{rotmean.pppgram}
\alias{rotmean.pppgram}
\title{Rotation Average of pppgram}
\usage{
rotmean.pppgram(x, t, t_steps = 50, adjust = 1, h, ..., normalise = TRUE)
}
\arguments{
\item{x}{periodogram object}

\item{t}{the wave number magnitude grid to interpolate the periodogram on}

\item{t_steps}{number of steps in the output length scale grid}

\item{adjust}{smoothing bandwidth adjustment}

\item{h}{Box-car kernel radius}

\item{...}{ignored}

\item{normalise}{divide by intensity? For Poisson the periodogram will then be constant 1.}
}
\description{
Rotation Average of pppgram
}
\details{
Box-car kernel with radius h = adjust * mean(dx, dy) where dx, dy are the step sizes of the grid used to compute the periodogram.
}
