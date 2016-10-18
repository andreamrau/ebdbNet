\name{dataFormat}
\alias{dataFormat}
\title{Change the Format of Longitudinal Data to be Compatible with EBDBN}
\description{
This function changes the format of longitudinal data to be compatible
with the format required by the EBDBN, namely a list (of length R) of
PxT matrices, where R, P, and T are the number of replicates, genes,
and time points, respectively.
}
\usage{
dataFormat(longitudinal.data)
}
\arguments{
  \item{longitudinal.data}{Data in the longitudinal format}
}
\details{
The argument refers to the general data structure of the 'longitudinal' package.
}
\value{
  List of length R of PxT matrices, suitable to be used in the EBDBN algorithm.
}
\author{Andrea Rau}
\examples{
library(ebdbNet)
library(GeneNet) ## Load GeneNet package to use T-cell activation data

data(tcell)  ## Load T-cell activation data
tc44 <- combine.longitudinal(tcell.10, tcell.34)

## Put data into correct format for algorithm
tcell.dat <- dataFormat(tc44)
}
\keyword{methods}
