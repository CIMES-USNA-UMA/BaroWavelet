% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CWT TransferFun.R
\name{TransferFunCWT}
\alias{TransferFunCWT}
\title{Continuous Wavelet Transform Transfer Function}
\usage{
TransferFunCWT(
  data,
  HF = 0.4,
  LF = 0.15,
  VLF = 0.04,
  chosen.dj = 1/20,
  dt = 0.25,
  demean = TRUE
)
}
\arguments{
\item{data}{A matrix with 3 columns containing time values (first column), RR and SBP values (second
and third column).}

\item{HF}{Maximum limit of the HF band. Default is 0.4 Hz}

\item{LF}{Maximum limit of the LF band. Default is 0.15 Hz}

\item{VLF}{Maximum limit of the VLF band. Default is 0.04 Hz}

\item{chosen.dj}{Frequency resolution. Default is 1/20}

\item{dt}{Time resolution (inverse of the sample rate). Default is 0.25}

\item{demean}{Boolean, should the data be demeaned before analysis? Default is TRUE}
}
\value{
A list with the estimated components of the baroreflex transfer function in the wavelet domain:
\item{TransferFun}{The computed baroreflex transfer function}
\item{Coherence}{The computed coherence between the two variables}
\item{Freqs}{A vector of frequencies for which the transfer function has been computed}
\item{Cone}{The computed cone of influence for the Continuous Wavelet Transform}
\item{Time}{The original vector of time values}
\item{HF}{The chosen maximum limit of the HF band}
\item{LF}{The chosen maximum limit of the LF band}
\item{VLF}{The chosen maximum limit of the VLF band}
\item{type}{A character string specifying which type of transfer function this is}
}
\description{
Computes barorreflex sensitivity as a transfer function suing the Continuous Wavelet Transform
}
\examples{
# ADD EXAMPLE!
}
\author{
Alvaro Chao-Ecija
}
