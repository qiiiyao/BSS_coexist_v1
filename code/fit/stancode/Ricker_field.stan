data{
  int<lower=0> n;              // length of time-series
  int<lower=0> nplot;           // number of mixes (i.e. communities)
  int<lower=0> Fieldnum[nplot]; //name of Field
  real<lower=0> N[n,50,nplot]; // data: dim1=time, dim2=species, dim3=replicate, dim4=mixcode
  real<lower=0> time[n];       // time-series
  int sp[nplot,50];              // species id for all mixes
}

transformed data {
  real x_r[0];
  int x_i[0];
  #real logN[n,50,nplot]; // log transformed data
  #logN=log(N);
}

parameters {
  vector<lower=0>[50] r;
  matrix<lower=0>[50,50] a;
  real<lower=0> sigma_e;
  matrix[n-1,50] rand_t;
  matrix[10,50] rand_f;
  real<lower=0> sigma_t;
}

model {

  // priors
  for (i in 1:50){
    r[i] ~ normal(0,1);
    for (j in 1:50){
      a[i,j] ~ normal(0,1);
    }
  }
  sigma_e ~ normal(0,1);
  sigma_t ~ normal(0,1);
  
  // Ricker competition model
  for (x in 1:nplot){ // loop for mixes
    
    // intermediate matrix
    matrix[n,50] Nsim; // simulated data: dim1=time, dim2=dim_ODE=2
    
    // fitting
    //for (k in 1:2){ // loop for replicates
      
      // simulate trajectory
      for (t in 1:(n-1)){
        for (m in 1:50){
          Nsim[t,m]=N[t,m,x]; // t Nsim is N from the first year
          Nsim[t+1,m]=N[t+1,m,x]; // t+1 Nsim is N from the first year
          Nsim[t+1,m] = Nsim[t,m]*exp(r[sp[x,m]]-sum(a[sp[x,1:50],sp[x,1:50]][m,].*Nsim[t,1:50])+ rand_t[t,m] + rand_f[Fieldnum[x], m]);
        }
      }
      
      // optimization
      for (j in 1:50){ // sr[x] species
        for (t in 1:n){
          N[t,j,x] ~ normal(Nsim[t,j], sigma_e);
        } // t
      } // j
    //} // k
  } // x
  
}
