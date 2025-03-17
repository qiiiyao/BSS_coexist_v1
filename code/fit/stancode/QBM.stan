data{
  int<lower=0> N; // observations
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  vector<lower=0,upper=1>[N] Y; // observation vector
  vector[N] X; // size vector
}

parameters{
  real a_mu;
  vector[Yrs] a;
  real b1_mu;
  vector[Yrs] b1;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  vector<lower=0>[N] sigmaSq;
}

transformed parameters{
  vector<lower=0>[N] tau;
  vector[N] mu;
  tau = sqrt(sigmaSq);
  for(n in 1:N)
    mu[n] = a[yid[n]] + b1[yid[n]]*X[n];
}

model{
  // Priors
  a_mu ~ normal(0,10);
  b1_mu ~ normal(0,10);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,5);
  a ~ normal(a_mu, sig_a);
  b1 ~ normal(b1_mu, sig_b1);
  for(n in 1:N)
  sigmaSq[n] ~ inv_gamma(1, 1);
 
  //Likelihood 
  for(n in 1:N)
  Y[n] ~ lognormal(mu[n], tau[n]);
}