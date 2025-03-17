data {
  int<lower=0> N;   // number of observations
  int<lower=0> P;  // number of predicators
  int<lower=0> ftime;
  matrix[N,P] X;   
  vector[N] G;     
  int<lower=1, upper=ftime> time[N];
}

parameters {
  real r;
  vector[P] a;
  vector[ftime] ran_t;
  real<lower=0> sigma_t;
  real<lower=0> sigma_e;
}
transformed parameters {
 // Rciker model 
  vector[N] yhat;
  for (i in 1:N)
    yhat[i] = r+X[i]*a+ran_t[time[i]];
  
}
model {
  // priors 
  r ~ normal(0, 1);
  for (i in 1:P){
      a[i] ~ normal(0, 1);
  }  
  sigma_t ~ normal(0, 1);
  ran_t ~ normal(0, sigma_t);
  sigma_e ~ normal(0, 1);
  // likelihood
  G ~ normal(yhat, sigma_e);
}

