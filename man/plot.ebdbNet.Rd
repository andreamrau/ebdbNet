\name{plot.ebdbNet}
\alias{plot.ebdbNet}
\title{Visualize EBDBN network}
\description{
A function to visualize graph estimated using the Empirical Bayes Dynamic Bayesian Network (EBDBN) 
algorithm.
}
\usage{
\method{plot}{ebdbNet}(x, sig.level, interactive = FALSE, clarify = "TRUE",
layout = layout.fruchterman.reingold,  ...)
}
\arguments{
  \item{x}{ An object of class \code{"ebdbNet"}}
  \item{sig.level}{Desired significance level (between 0 and 1) for edges in network}
  \item{interactive}{If TRUE, interactive plotting through tkplot}
  \item{clarify}{If TRUE, unconnected nodes should be removed from the plot}
  \item{layout}{Layout parameter for graphing network using igraph0}
  \item{...}{Additional arguments (mainly useful for plotting)}
}
\details{

For input networks, the default colors for nodes representing inputs and genes are green and blue, respectively. For feedback networks, the default color for all nodes is blue.

The interactive plotting option should only be used for relatively small networks (less than 
about 100 nodes).

}

\references{
Andrea Rau, Florence Jaffrezic, Jean-Louis Foulley, and R. W. Doerge (2010). An Empirical
Bayesian Method for Estimating Biological Networks from Temporal Microarray Data. \emph{Statistical
Applications in Genetics and Molecular Biology} 9. Article 9.
}
\author{Andrea Rau}

\seealso{\code{\link{ebdbn}}}
\examples{
library(ebdbNet)
tmp <- runif(1) ## Initialize random number generator
set.seed(125214) ## Save seed

## Simulate data
R <- 5
T <- 10
P <- 10
simData <- simulateVAR(R, T, P, v = rep(10, P), perc = 0.10)
Dtrue <- simData$Dtrue
y <- simData$y

## Simulate 8 inputs
u <- vector("list", R)
M <- 8
for(r in 1:R) {
	u[[r]] <- matrix(rnorm(M*T), nrow = M, ncol = T)
}

####################################################
## Run EB-DBN without hidden states
####################################################
## Choose alternative value of K using hankel if hidden states are to be estimated
## K <- hankel(y)$dim

## Run algorithm	
## net <- ebdbn(y = y, K = 0, input = u, conv.1 = 0.15, conv.2 = 0.10, conv.3 = 0.10,
##	verbose = TRUE)

## Visualize results
## plot(net, sig.level = 0.95)

}
\keyword{methods}
\keyword{models}
