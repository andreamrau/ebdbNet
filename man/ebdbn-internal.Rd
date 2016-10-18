\name{ebdbn-internal}
\alias{ebdbn-internal}
\alias{fdbkFunc}
\alias{sumFunc}
\title{Internal functions for Empirical Bayes Dynamic Bayesian Network (EBDBN) Estimation}
\description{
Internal functions for the ebdbNet package.
}
\usage{
sumFunc(x, cutoff) 
fdbkFunc(y)

}
\arguments{
  \item{x}{Vector of singular values from singular value decomposition of block-Hankel matrix}
  \item{cutoff}{Value to determine cutoff to be considered for singular values (e.g., 0.90)}
  \item{y}{A list of R (PxT) matrices of observed time course profiles}
}
\author{Andrea Rau}
\seealso{\code{\link{hankel}}, \code{\link{ebdbn}}}
\keyword{methods}

