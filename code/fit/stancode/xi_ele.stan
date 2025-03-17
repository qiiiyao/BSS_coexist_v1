
data {
  int<lower=0> s; // number of plots
  int<lower=0> timestep; //number of timesteps = ?
  int Y_O[s,timestep]; //species occurence data
  real Y_C[s,timestep]; // species cover data
  real Y_G[s,timestep-1]; // growth rate data
}


parameters {
  vector[timestep] A_t;
  vector[s] C_f;
  real b1; //bio interaction effects on occurence
  vector[timestep] B2_t_1;
  real b3; // bio-interaction effects on survivor
  real b4; // bio-interaction effects on colonization
  real b5; // bio-interaction effects on growth
  real<lower=0> tau;
  real sigmaplot;
  real psi;
}


model {
  //prior
  for(t in 1:timestep-1){
  A_t[t]~normal(0,10^2); // random time effects
  }
  
  sigmaplot~uniform(0, 10);
 
  for(i in 1:s){
  C_f[i]~normal(0,tauplot); // random field effects
  }
  
 
  for(t in 1:timestep-1){
  B2_t_1[t]~normal(0,10^2); // random time effects on survivor
}

  b3~normal(0,10^2); // bio-interaction effects on survivor
  b4~normal(0,10^2); // bio-interaction effects on colonization
  b5~normal(0,10^2);
  tau~inv_gamma(1e-3, 1e-3);
  psi~uniform(0,1)
  matrix[s,timestep] MU; // growth rates
  
  //Occurence
  psi_o~uniform(0,10)
  b1~normal(0,10^2); //bio interaction effects on occurence
  sigmaplot_o~uniform(0,.1)
  tauplot = sqrt(1/(sigmaplot*sigmaplot))
  for(i in 1:10){
        c[i] ~ dnorm(0,tauplot)
        
      }
      
      for (i in 1:(nyear-1)){
        a[i] ~ dnorm(0,.1)
      }
  for(t in 1:timestep-1){
     for(i in 1:s){
  if(Y_O[i,t]==0){
    Y_O[i,t+1]~bernoulli_logit(A_t[t]+b1*Y_C[i,t]+C_f[i]);
  }}}
  
  //Survivor_colonization
    for(i in 1:s){
      Y_0[i,1]~bernoulli(psi)
    for(t in 2:timestep){
      Y_O[i,t+1]~bernoulli_logit(A_t[t]+Y_O[i,t]*B2_t_1[t]+Y_O[i,t]*b3*Y_C[i,t]+(1-Y_O[i,t])*b4*Y_C[i,t]+C_f[i]);
  }}
  
  //Growth
   for(t in 1:timestep-1){
     for(i in 1:s){
  MU[i,t] = A_t[t]+b5*Y_C[i,t]+C_f[i];
  Y_G[i,t]~normal(MU[i,t], tau);
  }}
 

}

