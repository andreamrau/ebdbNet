# ebdbNet: Empirical Bayes Estimation of Dynamic Bayesian Networks

Author: Andrea Rau 

This package is used to infer the adjacency matrix of a network from time course data using an empirical Bayes 
estimation procedure based on Dynamic Bayesian Networks.

Posterior distributions (mean and variance) of network parameters are estimated using time-course
data based on a linear feedback state space model that allows for a set of hidden states to be in-
corporated.  The algorithm is composed of three principal parts:  choice of hidden state dimension
(see `hankel`), estimation of hidden states via the Kalman filter and smoother, and calculation of
posterior distributions based on the empirical Bayes estimation of hyperparameters in a hierarchical
Bayesian framework (see `ebdbn`).

Plot functionalities are provided via the `igraph` package.

### Reference

A. Rau, F. Jaffrezic, J.-L. Foulley, R. W. Doerge (2010). An empirical Bayesian method for estimating biological networks from 
temporal microarray data. Statistical Applications in Genetics and Molecular Biology, vol. 9, iss. 1, article 9. 
