data {
  int<lower=0> N;   // number of observations
  int<lower=0> P;  // number of predicators
  int<lower=0> time_diff[N]; // year difference in our data
  int<lower=0> ftime;
  int<lower=0> ffield;
  matrix[N,P-1] X_t_inter; 
  vector[N] X_t_intra;
  vector<lower=0, upper=1>[N] Y_t1;     
  int<lower=1, upper=ftime> time[N];
  int<lower=1> field[N];
}

parameters {
  
  real<lower=0> r; // intrinsic growth rate
  real<upper=0> a_intra; // species interaction coefficients no permission
  vector[P-1] a_inter;
  
  vector[ftime] ran_t_raw;
  real<lower=0> sigma_ran_t;
  vector[ffield] ran_f_raw;
  real<lower=0> sigma_ran_f;
  real<lower=0> sigma_e; // sigma_e is observational error
}

transformed parameters {
  vector[ftime] ran_t;
  vector[ffield] ran_f; 

  //non-centred parameterization to improve efficiency and convergence of hierachical model
  // followed Yu_et_al_EL_2024
  ran_t = sigma_ran_t * ran_t_raw;
  ran_f = sigma_ran_f * ran_f_raw;
  
 // Rciker model 
  vector[N] yhat;
  for (i in 1:N) {
    yhat[i] = X_t_intra[i]*exp((r+X_t_intra[i]*a_intra+X_t_inter[i]*a_inter+ran_t[time[i]]+ran_f[field[i]])*time_diff[i]);
  }
}

model {
  // priors 
  r ~ normal(0, 1);
  for (i in 1:(P-1)){
      a_inter[i] ~ normal(0, 1);
  }  
  a_intra ~ normal(-1, 1);
  sigma_ran_t ~ normal(0, 0.1);
  ran_t_raw ~ normal(0, 1);
  sigma_ran_f ~ normal(0, 0.1);
  ran_f_raw ~ normal(0, 1);
  sigma_e ~ normal(0, 1);
  
  // likelihood
  Y_t1 ~ lognormal(yhat, sigma_e); // robust regression using student distribution following McElreath_2020
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = lognormal_lpdf(Y_t1[n] | yhat[n], sigma_e);
  }
}

