% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AnalysesFunctions.R
\name{AddAvgCwtData}
\alias{AddAvgCwtData}
\title{Add averaged transfer function calculated from the CWT}
\usage{
AddAvgCwtData(framework, locator)
}
\arguments{
\item{framework}{An analysis environment generated by \link[BaroWavelet]{BuildStructure}}

\item{locator}{A locator index indicating which analysis slot needs to be used}
}
\value{
The original analysis environment, in which the averaged CWT transfer function has been added
        to a specific analysis slot.
}
\description{
Computes a vector of the frequency averaged BRS estimates previously obtained from the CWT.
}
\examples{
Data <- InterpolateData(DataSimulation(), f = 1)
Study <- BuildStructure()
Study <- AddAnalysis(Study, name = "Simulation")
Study <- AddDataToAnalysis(Study, 1, Data$RR, Data$SBP, Data$Time)
Study <- AnalyzeBRS(Study, 1)
Study <- AddAvgCwtData(Study, 1)
}
\author{
Alvaro Chao-Ecija
}
