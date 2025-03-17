data{
  int<lower=0> N; // observations
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> S; // species number
  vector<lower=0,upper=1>[N] Y; // observation vector
  matrix[N,S] X; // size vector
}

parameters{
  real a_mu;
  vector[Yrs] a;
  vector[S] b;
  vector[S] b_mu;
  real<lower=0> sig_a;
  vector<lower=0> [S] sig_b;
  vector<lower=0>[N] sigmaSq;
}

transformed parameters{
  vector<lower=0>[N] tau;
  vector[N] mu;
  tau = sqrt(sigmaSq);
  for(n in 1:N)
   for(s in 1:S)
    mu[n] = a[yid[n]] + sum(b[s]*X[n,]);
}

model{
  // Priors
  a_mu ~ normal(0,10);
  
  for(s in 1:S)
  b_mu[s] ~ normal(0,10);
  
  sig_a ~ cauchy(0,5);
  
  for(s in 1:S)
  sig_b[s] ~ cauchy(0,5);
  
  a ~ normal(a_mu, sig_a);
  
  for(s in 1:S)
  b[s] ~ normal(b_mu, sig_b);
  
  for(n in 1:N)
  sigmaSq[n] ~ inv_gamma(1, 1);
 
  //Likelihood 
  Y ~ lognormal(mu, tau);
}