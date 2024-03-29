% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AnalysesFunctions.R
\name{AddTimeInterval}
\alias{AddTimeInterval}
\title{Add time interval}
\usage{
AddTimeInterval(
  framework,
  name = NULL,
  use.default = TRUE,
  use.coherence = TRUE
)
}
\arguments{
\item{framework}{An analysis environment generated by \link[BaroWavelet]{BuildStructure}}

\item{name}{An optional name for the interval. The default is NULL}

\item{use.default}{Boolean. Use the default variables from "General Data". Default is TRUE}

\item{use.coherence}{Boolean, should a coherence threshold be used in the analyses? Default is TRUE}
}
\value{
A framework with the new additional interval.
}
\description{
Adds a certain time interval to the structure, so that estimates can be obtained from it.
}
\details{
This function allows users to add new time segments to their study environments.For example, it
         allows users to define a control interval and a stimulation interval for their study.
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
}
\author{
Alvaro Chao-Ecija
}
