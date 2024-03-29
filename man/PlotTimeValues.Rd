% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AnalysesFunctions.R
\name{PlotTimeValues}
\alias{PlotTimeValues}
\title{Plot heart rate and blood pressure levels}
\usage{
PlotTimeValues(
  framework,
  locator,
  use.names = TRUE,
  HRcol = "grey20",
  SBPcol = "grey80",
  restrict = NULL,
  output = FALSE,
  ymax = NULL,
  tem = FALSE,
  newPlot = TRUE
)
}
\arguments{
\item{framework}{An analysis environment generated by \link[BaroWavelet]{BuildStructure}}

\item{locator}{A locator index indicating which analysis slot needs to be used}

\item{use.names}{Boolean, should the names of the subjects be used?}

\item{HRcol}{Color to be used to represent the HR}

\item{SBPcol}{Color to be used to represent the SBP}

\item{restrict}{Vector used to restrict the plot to specific intervals. Default is NULL}

\item{output}{Boolean. Should the indices for each interval be returned? Default is FALSE}

\item{ymax}{Maximum level for the y axis. Default is NULL}

\item{tem}{Boolean, creates a temporal file for the plot. Default is FALSE}

\item{newPlot}{Boolean, generates a new plot without overwriting a previous plot. Default is TRUE}
}
\value{
None
}
\description{
Plots HR and SBP levels associated to a specific subject at each interval
}
\examples{
Data <- InterpolateData(DataSimulation(), f = 1)
Study <- BuildStructure()
Study <- AddAnalysis(Study, name = "Simulation")
Study <- AddDataToAnalysis(Study, 1, Data$RR, Data$SBP, Data$Time)
Study <- AnalyzeBRS(Study, 1)
Study <- AddAvgCwtData(Study, 1)

Study <- AddTimeInterval(Study, "First interval")
Study <- AddTimeInterval(Study, "Second interval")
Study <- AddTimeInterval(Study, "Third interval")

Study <- AnalyzeBRSIndices(Study, 1, 1, c(1.7, 3.3))
Study <- AnalyzeBRSIndices(Study, 1, 2, c(10, 11.7))
Study <- AnalyzeBRSIndices(Study, 1, 3, c(35, 37))

PlotTimeValues(Study, 1, newPlot = FALSE)
}
\author{
Alvaro Chao-Ecija
}
