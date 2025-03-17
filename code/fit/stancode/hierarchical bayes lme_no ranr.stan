data {
  int<lower=0> N;   // number of observations
  int<lower=0> P;  // number of predicators
  int<lower=0> ftime;
  matrix[N,P] X;  // 
  vector[N] G;     
  int<lower=1,upper=ftime> time[N];
  real<lower=0> sdscal;
}
parameters {
  vector<upper=0>[P] a;
  real<lower=0> r;
  real<lower=0> sigma_t;
  real<lower=0> sigmaeps;
  
  vector[ftime] eta_t;
}
transformed parameters {
  vector[ftime] ran_t;
  vector[N] yhat;
  ran_t = sigma_t * eta_t;
  
  for (i in 1:N)
    yhat[i] = r+X[i]*a+ran_t[time[i]];
  
}
model {
  eta_t ~ normal(0, 1);
  sigma_t ~ cauchy(0, 2.5*sdscal);
  sigmaeps ~ cauchy(0, 2.5*sdscal);
  G ~ normal(yhat, sigmaeps);
}

