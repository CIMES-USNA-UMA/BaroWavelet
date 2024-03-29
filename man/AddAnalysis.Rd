% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AnalysesFunctions.R
\name{AddAnalysis}
\alias{AddAnalysis}
\title{Add analysis}
\usage{
AddAnalysis(framework, name = NULL)
}
\arguments{
\item{framework}{An analysis environment generated by \link[BaroWavelet]{BuildStructure}}

\item{name}{A name for the analysis to be added. Default is NULL}
}
\value{
The original analysis environment with a new analysis slot
}
\description{
Add a new analysis to an analysis environment.
}
\examples{
Study <- BuildStructure()
Study <- AddAnalysis(Study, name = "Subject A")
}
\author{
Alvaro Chao-Ecija
}
