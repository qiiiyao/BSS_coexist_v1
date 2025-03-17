data{
  int<lower=0> n;              // length of time-series
  int<lower=0> nplot;           // number of mixes (i.e. communities)
  int<lower=0> Fieldnum[nplot]; //name of Field
  int<lower=0> time_diff[nplot]; //name of Field
  real g[n,nplot]; // data: dim1=time, dim2=mixcode
  real<lower=0> N[n,50,nplot];
  real<lower=0> time[n];       // time-series
}


parameters {
  real<lower=0> r;
  vector<lower=0>[50] a;
  real<lower=0> sigma_e;
  real<lower=0> sigma_field;
  vector[n] rand_t;
  vector[10] rand_f;
}

transformed parameters {
 // Rciker model 
  vector[n] gsim;
  matrix[n,50] Nsim;

 // Ricker competition model
  for (x in 1:nplot){ // loop for mixes
    Nsim = to_matrix(N[,,x]);
    // fitting
      for (t in 1:n){
          gsim[t] = r-Nsim[t]*a+ rand_t[t] + rand_f[Fieldnum[x]];
        } // t
      } // x
}


model {
  
  // priors
  r ~ normal(0,1);
  for (i in 1:50){
      a[i] ~ normal(0,1);
  }
  
  sigma_e ~ gamma(0.001, 0.001);
  sigma_field ~ normal(0, 0.001);
  
  for(i in 1:10){
    rand_f[i] ~ normal(0, sigma_field);
  }
  
  for (i in 1:n){
    rand_t[i] ~ normal(0, 0.001);
  }
  
  // likelihood
  for (x in 1:nplot){ // loop for mixes
    // fitting
      for (t in 1:n){
          g[t,x] ~ normal(gsim[t], sigma_e);
        } // t
      } // x
  }
  
generated quantities {

  matrix[nplot, n] log_lik_g;	// log-likelihood for the RIM
  
  for (x in 1:nplot){ // loop for mixes
    // fitting
      for (t in 1:n){
          log_lik_g[t,x] = normal_lpdf(g[t,x] | gsim[t], sigma_e);
      }
  }
  
}





  
