\name{ebdbNet-package}
\alias{ebdbNet-package}
\alias{ebdbNet}
\docType{package}
\title{
Empirical Bayes Dynamic Bayesian Network (EBDBN) Inference
}
\description{
This package is used to infer the adjacency matrix of a network from time course data using an
empirical Bayes estimation procedure based on Dynamic Bayesian Networks.
}
\details{
\tabular{ll}{
Package: \tab ebdbNet\cr
Type: \tab Package\cr
Version: \tab 1.2.3\cr
Date: \tab 2013-07-26\cr
License: \tab GPL (version 3 or later)\cr
LazyLoad: \tab yes\cr
}
Posterior distributions (mean and variance) of network parameters are estimated using time-course
data based on a linear feedback state space model that allows for a set of hidden states to be incorporated.
The algorithm is composed of three principal parts: choice of hidden state dimension
(see \code{\link{hankel}}), estimation of hidden states via the Kalman filter and smoother, and calculation of
posterior distributions based on the empirical Bayes estimation of hyperparameters in a hierarchical
Bayesian framework (see \code{\link{ebdbn}}). 
}
\author{
Andrea Rau

Maintainer: Andrea Rau <andrea.rau AT jouy.inra.fr>
}
\references{
Andrea Rau, Florence Jaffrezic, Jean-Louis Foulley, and R. W. Doerge (2010). An Empirical
Bayesian Method for Estimating Biological Networks from Temporal Microarray Data. \emph{Statistical
Applications in Genetics and Molecular Biology} 9. Article 9.
}
\keyword{ package }
\examples{
library(ebdbNet)
library(GeneNet)  ## Load GeneNet package to use T-cell activation data
tmp <- runif(1) ## Initialize random number generator
set.seed(4568818) ## Set seed

## Load T-cell activation data
data(tcell) 
tc44 <- combine.longitudinal(tcell.10, tcell.34)

## Put data into correct format for algorithm
## (List, with one matrix per replicate (P rows and T columns)
tcell.dat <- dataFormat(tc44)

## Use only subset of T-cell data for faster example
R <- 20 ## 20 replicates
P <- 10 ## 10 genes
tcell.sub.dat <- vector("list", R)
rep.sample <- sample(1:44, R)
for(r in 1:R) {
	tcell.sub.dat[[r]] <- tcell.dat[[rep.sample[r]]][sample(1:58, P),]
}

####################################################
# Run EBDBN (no hidden states) with feedback loops
####################################################
## Choose alternative value of K using hankel if hidden states to be estimated
## K <- hankel(tcell.sub.dat, lag = 1)$dim

## Run algorithm (feedback network, no hidden states)
net <- ebdbn(y = tcell.sub.dat, K = 0, input = "feedback", conv.1 = 0.01,
	conv.2 = 0.01, conv.3 = 0.001, verbose = TRUE)

## Visualize results: in this example, mostly feedback loops
## plot(net, sig.level = 0.5)
}
