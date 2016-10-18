\name{hankel}
\alias{hankel}
\title{Perform Singular Value Decomposition of Block-Hankel Matrix}
\description{
This function constructs a block-Hankel matrix based on time-course data, performs the subsequent
singular value decomposition (SVD) on this matrix, and returns the number of large singular values
as defined by a user-supplied cutoff criterion.
}
\usage{
hankel(y, lag, cutoff, type)
}
\arguments{
  \item{y}{A list of R (PxT) matrices of observed time course profiles}
  \item{lag}{Maximum relevant time lag to be used in constructing the block-Hankel matrix}
  \item{cutoff}{Cutoff to be used, determined by desired percent of total variance explained}
  \item{type}{Method to combine results across replicates ("median" or "mean")}
}
\details{
Constructs the block-Hankel matrix \eqn{H} of autocovariances of time series observations is constructed
(see references for additional information), where the maximum relevant time lag must be specified
as \code{lag}. In the context of gene networks, this corresponds to the maximum relevant biological time
lag between a gene and its regulators. This quantity is experiment-specific, but will generally be
small for gene expression studies (on the order of 1, 2, or 3).

The singular value decomposition of \eqn{H} is performed, and the singular values are ordered by size
and scaled by the largest singular value. Note that if there are T time points in the data, only the
first (T - 1) singular values will be non-zero.

To choose the number of large singular values, we wish to find the point at which the inclusion of
an additional singular value does not increase the amount of explained variation enough to justify
its inclusion (similar to choosing the number of components in a Principal Components Analysis).
The user-supplied value of \code{cutoff} gives the desired percent of variance explained by the first set
of K principal components. The algorithm returns the value of K, which may subsequently be used
as the dimension of the hidden state in \code{\link{ebdbn}}.

The argument 'type' takes the value of "median" or "mean", and is used to determine how results
from replicated experiments are combined (i.e., median or mean of the per-replicate final hidden
state dimension).
}
\value{
  \item{svs }{Vector of singular values of the block-Hankel matrix \eqn{H}}
  \item{dim }{Number of large singular values, as determined by the user-supplied cutoff}
}
\references{Masanao Aoki and Arthur Havenner (1991). State space modeling of multiple time series. \emph{Econometric Reviews} 10(1), 1-59.

Martina Bremer (2006). \emph{Identifying regulated genes through the correlation structure of time dependent microarray data}. Ph. D. thesis, Purdue University.

Andrea Rau, Florence Jaffrezic, Jean-Louis Foulley, and R. W. Doerge (2010). An Empirical Bayesian Method for Estimating Biological Networks from Temporal Microarray Data. 
\emph{Statistical Applications in Genetics and Molecular Biology} 9. Article 9.
}
\author{Andrea Rau}
\examples{
library(ebdbNet)
tmp <- runif(1) ## Initialize random number generator
set.seed(125214) ## Save seed

## Simulate data
y <- simulateVAR(R = 5, T = 10, P = 10, v = rep(10, 10), perc = 0.10)$y

## Determine the number of hidden states to be estimated (with lag = 1)
K <- hankel(y, lag = 1, cutoff = 0.90, type = "median")$dim
## K = 5

}
\keyword{methods}
