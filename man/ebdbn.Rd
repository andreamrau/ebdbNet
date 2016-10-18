\name{ebdbn}
\alias{ebdbn}
\title{Empirical Bayes Dynamic Bayesian Network (EBDBN) Estimation}
\description{
A function to infer the posterior mean and variance of network parameters
using an empirical Bayes estimation procedure for a Dynamic Bayesian Network (DBN).
}
\usage{
ebdbn(y, K, input = "feedback", conv.1 = .15, conv.2 = .05, conv.3 = .01, verbose = TRUE,
max.iter = 100, max.subiter = 200)
}
\arguments{
  \item{y}{A list of R (PxT) matrices of observed time course profiles (P genes, T time points)}
  \item{K}{Number of hidden states}
  \item{input}{"feedback" for feedback loop networks, or a list of R (MxT) matrices of input profiles}
  \item{conv.1}{Value of convergence criterion 1}
  \item{conv.2}{Value of convergence criterion 2}
  \item{conv.3}{Value of convergence criterion 3}
  \item{verbose}{Verbose output}
  \item{max.iter}{Maximum overall iterations (default value is 100)}
  \item{max.subiter}{Maximum iterations for hyperparameter updates (default value is 200)}
}
\details{
An object of class \code{ebdbNet}.

This function infers the parameters of a network, based on the state space model
\deqn{x_t = Ax_{t-1} + Bu_t + w_t}{x(t) = Ax(t-1) + Bu(t) + w(t)}
\deqn{y_t = Cx_t + Du_t + z_t}{y(t) = Cx(t) + Du(t) + z(t)}
where \eqn{x_t}{x(t)} represents the expression of K hidden states at time \eqn{t},
\eqn{y_t}{y(t)} represents the expression of P observed states (e.g., genes) at time
\eqn{t}, \eqn{u_t}{u(t)} represents the values of M inputs at time \eqn{t},
\eqn{w_t \sim MVN(0,I)}{w(t) ~ MVN(0,I)}, and \eqn{z_t \sim MVN(0,V^{-1})}{z(t) ~ MVN(0,V^(-1))},
with \eqn{V = diag(v_1, \ldots, v_P)}{V = diag(v_1, ..., v_P)}. Note that the
dimensions of the matrices \eqn{A}, \eqn{B}, \eqn{C}, and \eqn{D} are (KxK),
(KxM), (PxK), and (PxM), respectively. When a network is estimated with
feedback rather than inputs (\code{input} = "feedback"), the state
space model is
\deqn{x_t = Ax_{t-1} + By_{t-1} + w_t}{x(t) = Ax(t-1) + By(t-1) + w(t)}
\deqn{y_t = Cx_t + Dy_{t-1} + z_t}{y(t) = Cx(t) + Dy(t-1) + z(t)}

The parameters of greatest interest are typically
contained in the matrix \eqn{D}, which encodes the direct interactions among
observed variables from one time to the next (in the case of feedback loops),
or the direct interactions between inputs and observed variables at each time point
(in the case of inputs).

The value of K is chosen prior to running the algorithm by using \code{\link{hankel}}.
The hidden states are estimated using the classic Kalman filter. Posterior distributions
of \eqn{A}, \eqn{B}, \eqn{C}, and \eqn{D} are estimated using an empirical
Bayes procedure based on a hierarchical Bayesian structure defined over the
parameter set. Namely, if \eqn{a_{(j)}}{a(j)}, \eqn{b_{(j)}}{b(j)},
\eqn{c_{(j)}}{c(j)}, \eqn{d_{(j)}}{d(j)}, denote vectors made up of the
rows of matrices \eqn{A}, \eqn{B}, \eqn{C}, and \eqn{D} respectively, then
\deqn{a_{(j)} \vert \alpha \sim N(0, diag(\alpha)^{-1})}{a(j)|alpha ~ N(0, diag(alpha)^(-1))}
\deqn{b_{(j)} \vert \beta \sim N(0, diag(\beta)^{-1})}{b(j)|beta ~ N(0, diag(beta)^(-1))}
\deqn{c_{(j)} \vert \gamma \sim N(0, diag(\gamma)^{-1})}{c(j)|gamma ~ N(0, diag(gamma)^(-1))}
\deqn{d_{(j)} \vert \delta \sim N(0, diag(\delta)^{-1})}{d(j)|delta ~ N(0, diag(delta)^(-1))}
where \eqn{\alpha = (\alpha_1, ..., \alpha_K)}{alpha = (alpha_1, ..., alpha_K)},
\eqn{\beta = (\beta_1, ..., \beta_M)}{beta = (beta_1, ..., beta_M)},
\eqn{\gamma = (\gamma_1, ..., \gamma_K)}{gamma = (gamma_1, ..., gamma_K)},
and \eqn{\delta = (\delta_1, ..., \delta_M)}{delta = (delta_1, ..., delta_M)}. An EM-like algorithm 
is used to estimate the hyperparameters in an iterative
procedure conditioned on current estimates of the hidden states.

\code{conv.1}, \code{conv.2}, and \code{conv.3} correspond to convergence
criteria \eqn{\Delta_1}{Delta_1}, \eqn{\Delta_2}{Delta_2}, and
\eqn{\Delta_3}{Delta_3} in the reference below, respectively.  After
terminating the algorithm, the z-scores of the \eqn{D} matrix is
calculated, which in turn determines the presence or absence of edges in the network.

See the reference below for additional details about the implementation of the
algorithm. 
}
\value{
  \item{APost }{Posterior mean of matrix \eqn{A}}
  \item{BPost }{Posterior mean of matrix \eqn{B}}
  \item{CPost }{Posterior mean of matrix \eqn{C}}
  \item{DPost }{Posterior mean of matrix \eqn{D}}
  \item{CvarPost }{Posterior variance of matrix C}
  \item{DvarPost }{Posterior variance of matrix D}
  \item{xPost }{Posterior mean of hidden states x}
  \item{alphaEst }{Estimated value of \eqn{\alpha}{alpha}}
  \item{betaEst }{Estimated value of \eqn{\beta}{beta}}
  \item{gammaEst }{Estimated value of \eqn{\gamma}{gamma}}
  \item{deltaEst }{Estimated value of \eqn{\delta}{delta}}
  \item{vEst }{Estimated value of precisions \eqn{v}}
  \item{muEst }{Estimated value of \eqn{\mu}{mu}}
  \item{sigmaEst }{Estimated value of \eqn{\Sigma}{Sigma}}
  \item{alliterations }{Total number of iterations run}
  \item{z }{Z-statistics calculated from the posterior distribution of matrix D}
  \item{type }{Either "input" or "feedback", as specified by the user}
}
\references{
Andrea Rau, Florence Jaffrezic, Jean-Louis Foulley, and R. W. Doerge (2010). An Empirical
Bayesian Method for Estimating Biological Networks from Temporal Microarray Data. \emph{Statistical
Applications in Genetics and Molecular Biology} 9. Article 9.
}
\author{Andrea Rau}

\seealso{\code{\link{hankel}}, \code{\link{calcSensSpec}}, \code{\link{plot.ebdbNet}}}
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
net <- ebdbn(y = y, K = 0, input = u, conv.1 = 0.15, conv.2 = 0.10, conv.3 = 0.10,
	verbose = TRUE)

## Visualize results
## Note: no edges here, which is unsurprising as inputs were randomly simulated
## (in input networks, edges only go from inputs to genes)
## plot(net, sig.level = 0.95)

}
\keyword{methods}
\keyword{models}
