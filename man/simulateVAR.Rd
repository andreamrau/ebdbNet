\name{simulateVAR}
\alias{simulateVAR}
\title{Simulate Simple Autoregressive Process}
\description{
Function to simulate a simple autoregressive process based on a network
adjacency matrix with a given percentage of non-zero elements.
}
\usage{
simulateVAR(R, T, P, v, perc)
}
\arguments{
  \item{R}{Number of replicates}
  \item{T}{Number of time points}
  \item{P}{Number of observations (genes)}
  \item{v}{(Px1) vector of gene precisions}
  \item{perc}{Percent of non-zero edges in adjacency matrix}
}
\details{
Data are simulated with R replicates, T time points, and P genes, based on
a first-order autoregressive process with Gaussian noise. The user can
specify the percentage of non-zero edges to be randomly selected in the
adjacency matrix.
}
\value{
  \item{Dtrue }{Adjacency matrix used to generate data (i.e., the true network)}
  \item{y }{Simulated data}
}
\author{Andrea Rau}
\seealso{\code{\link{ebdbn}}}
\examples{
library(ebdbNet)
tmp <- runif(1) ## Initialize random number generator
set.seed(125214) ## Save seed

## Simulate data
simData <- simulateVAR(R = 5, T = 10, P = 10, v = rep(10, 10), perc = 0.10)
Dtrue <- simData$Dtrue
y <- simData$y
}
\keyword{methods}
