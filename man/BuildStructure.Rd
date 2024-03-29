% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AnalysesFunctions.R
\name{BuildStructure}
\alias{BuildStructure}
\title{Build analysis structure}
\usage{
BuildStructure(
  name = NULL,
  HF = 0.4,
  LF = 0.15,
  VLF = 0.04,
  wavelet = "d8",
  coh = 0.5,
  dj = 1/20,
  cwt.type = c("tf", "alpha"),
  use.coherence = TRUE,
  error = 5e-04,
  error.type = c("absolute", "relative"),
  index.method = "median"
)
}
\arguments{
\item{name}{A name for the analysis environment. Default is NULL}

\item{HF}{Maximum limit of the HF band to be used in the analyses. Default is 0.4 Hz}

\item{LF}{Maximum limit of the LF band to be used in the analyses. Default is 0.15 Hz}

\item{VLF}{Maximum limit of the VLF band to be used in the analyses. Default is 0.04 Hz}

\item{wavelet}{A wavelet to be used for the analyses. Default is d4}

\item{coh}{A coherence threshold to be used for the analyses. Default is 0.5}

\item{dj}{A frequency resolution to be used for the analyses. Default is 1/20}

\item{cwt.type}{Type of BRS estimation (either "tf" for transfer function or "alpha" for alpha index)}

\item{use.coherence}{Boolean, should a coherence threshold be used in the analyses? Default is TRUE}

\item{error}{Tolerance for the PMODWPT. Default is 0.0005}

\item{error.type}{Type of tolerance (absolute or relative) for the PMODWPT}

\item{index.method}{Method to compute individual BRS indices for subjects. Default is "median"}
}
\value{
A list containing the generated environment for the analyses
}
\description{
Generates a structure in which analyses can be performed and results can be stored
}
\details{
This function allows users to create environments to organize complex analyses and store
their results.
}
\examples{
Study <- BuildStructure()

}
\author{
Alvaro Chao-Ecija
}
