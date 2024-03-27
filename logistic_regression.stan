//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> n;      // Number of samples
  int<lower=0> q;      // Number of covariates
  matrix[n,q] X;       // Predictor matrix
  int<lower=0> y[n];   // Binomial outcome
  int priormu;         // Mean of normal prior
  int priorsigma;      // sd of normal prior
  int m;               // number of trials for binomial 
}

parameters {
  real beta0;
  vector[q] beta;
}

model {
  //Prior
  beta0 ~ normal(priormu,priorsigma);
  beta ~ normal(priormu,priorsigma);
  
  //likelihood
  y ~ binomial_logit(m,beta0 + X*beta);
}

