% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotHRV.R
\name{PlotHRV}
\alias{PlotHRV}
\title{Plot individual HRV}
\usage{
PlotHRV(
  hrv,
  time,
  time_flags = NULL,
  col = "brown",
  use.xlim = FALSE,
  tem = FALSE,
  newPlot = TRUE,
  title = "data",
  plotHF = TRUE,
  plotLF = TRUE,
  ratio = TRUE,
  normalize = FALSE,
  ylim = NULL,
  use.ggplot = TRUE,
  fHF = 0.4,
  fLF = 0.15,
  fVLF = 0.04,
  size.axis = 100,
  size.labels = 100,
  size.title = 100,
  d.labels = 100
)
}
\arguments{
\item{hrv}{HRV obtained either by RHRV function \link[RHRV]{CalculatePowerBand} or from \link[BaroWavelet]{AlphaIndexDWT}. If
the former method is used, it must be a list with variables named HF, LF and LFHF. If the latter method is
used, argument hrv must be set to TRUE.}

\item{time}{Time values for the HRV data}

\item{time_flags}{A vector containing the minimum and maximum limits of a time interval, in minutes.
Default is NULL.}

\item{col}{Color used to highlight a specific time interval. Default is brown}

\item{use.xlim}{Boolean. Should the argument time_flags be used as limits of the x-axis? Default is FALSE.}

\item{tem}{Boolean, creates a temporal file for the plot. Default is FALSE}

\item{newPlot}{Boolean, generates a new plot without overwriting a previous plot. Default is TRUE}

\item{title}{Part of the title for the plot. Used to indicate the data source on the title. Default is "data"}

\item{plotHF}{Boolean, plot results form the HF band. Default is TRUE}

\item{plotLF}{Boolean, plot results from the LF band. Default is TRUE}

\item{ratio}{Boolean. Should the LF/HF ratio be plotted? Default is FALSE. Arguments plotHF and plotLF must also be
set to TRUE.}

\item{normalize}{Boolean. Should normalized HRV be plotted? Default is FALSE}

\item{ylim}{Maximum y axis limit. Default is NULL}

\item{use.ggplot}{Boolean, use methods from \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2} package to plot the results. Default is TRUE}

\item{fHF}{Maximum limit of the HF band, shown at the plot title. Default is 0.4 Hz}

\item{fLF}{Maximum limit of the LF band, shown at the plot title. Default is 0.15 Hz}

\item{fVLF}{Maximum limit of the VLF band, shown at the plot title. Default is 0.04 Hz}

\item{size.axis}{Percentage of scaling of axis values. Default is 100}

\item{size.labels}{Percentage of scaling of axis labels. Default is 100}

\item{size.title}{Percentage of scaling of plot titles. Default is 100}

\item{d.labels}{Percentage of scaling of distance between axis and labels. Default is 100}
}
\value{
None
}
\description{
Plots HRV data computed by function \link[RHRV]{CalculatePowerBand} from package \href{https://CRAN.R-project.org/package=RHRV}{RHRV}
}
\examples{
Data <- InterpolateData(DataSimulation(), f = 1)
AlphaIndex <- AlphaIndexDWT(Data, wv = "d8", error = 0.0005, hrv = TRUE)


PlotHRV(AlphaIndex$HRV, AlphaIndex$Time, plotHF = FALSE, plotLF = TRUE,
newPlot = FALSE)

}
\author{
Alvaro Chao-Ecija
}
