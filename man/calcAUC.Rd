\name{calcAUC}
\alias{calcAUC}
\title{Calculate the Approximate Area Under the Curve (AUC)}
\description{
Returns the approximate Area Under the Curve (AUC) of a Receiver Operating Characteristic (ROC) curve.
}
\usage{
calcAUC(sens, cspec)
}
\arguments{
  \item{sens}{Vector of sensitivity values, calculated for varying thresholds}
  \item{cspec}{Vector of complementary specificity values, calculated for the same varying thresholds as \code{sens}}
}
\details{
Let TP, FP, TN, and FN represent the number of true positives, false positives, true negatives and false
negatives of inferred network edges, respectively.  Sensitivity is defined as \deqn{\frac{TP}{TP + FN}}{TP/(TP + FN)} and
complementary specificity is defined as \deqn{\frac{TN}{TN + FP}}{TN/(TN + FP)}

Note that \code{sens} and \code{cspc} should be in the same order with respect to the threshold value 
so that their elements correspond.  That is, if the first element of \code{sens} was calculated at a
threshold value of 0.01 and the second at a threshold value of 0.02, then the first element of \code{cpsec}
should be also be calculated at a threshold value of 0.01 and the second at a threshold value of 0.02, and
so on.  The AUC is approximated using the trapezoid method, and can take real values between 0 and 1.  An
AUC of 0.5 indicates a classifier with random performance, and an AUC of 1 indicates a classifer with
perfect performance.}

\value{
AUC of the ROC curve
}
\author{Andrea Rau}
\examples{
library(ebdbNet)
tmp <- runif(1) ## Initialize random number generator

## Generate artificial values for sensitivity and complementary specificity.
fn <- function(x) {return(-1/(x^7)+1)}
set.seed(1459)
sens <- c(fn(seq(1, 2.7, length = 100)),1) ## Sensitivity
cspec <- seq(0, 1, by = 0.01) ## Complementary specificity

## Calculate the AUC of the ROC curve
AUC <- calcAUC(sens, cspec) ## AUC of this ROC curve is 0.9030868
}
\keyword{methods}
