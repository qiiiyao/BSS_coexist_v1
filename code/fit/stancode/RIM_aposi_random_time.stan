/* NDDM & RIM joint model

code to run on all focal species at once

Bimler 2022
*/ 
  
  
data {
  int<lower=1> S;          // total number of focal species (s in the manuscript) 
  int<lower=1> N;          // total number of observations (rows in model matrix) (n in the manuscript)
  int<lower=0> T;          // total number of interaction partners (columns in model matrix) (t in the manuscript)

  int<lower=0> species_ID[N];   // index matching observations to focal species (d in the manuscript)
  real perform[N];      // response variable (p in the manuscript)
    
  int<lower=0> ftime;
  int<lower=0> time_diff[N];
  matrix[N,T] X;         // neighbour abundances (the model matrix)
  int<lower=1, upper=ftime> time[N];
} 

transformed data {
  // Define transformed data variables here
  vector[N] log_perform; // Example transformation: Log-transform of observations

  // Compute the transformations
  for (n in 1:N) {
    log_perform[n] = log(perform[n]);
  }
}

parameters {
  vector<lower=0>[S] r;    // species-specific intercept 
    
  vector<lower=0> [S] sigma; // species-specific dispersion deviation parameter, 
  // defined for the negative binomial distribution used to reflect seed production (perform)
  // disp_dev = 1/sqrt(phi)
  
  vector[ftime] ran_t;
  //real<lower=0> weight;    // weighting value controlling the average strength of interactions (must be positive)
  vector[S] raw_response; // species-specific effect parameter (r in the manuscript but without weight)
  vector[T] raw_effect;   // species-specific effect parameter (e in the manuscript)

} 

transformed parameters {
  
  // transformed parameters constructed from parameters above
  //vector[S] response;
  //{
   // real norm = sqrt(dot_self(raw_response)); // 计算向量的欧几里得范数
   // response = raw_response / norm; // 规范化向量
  //}
  
  //vector[T] effect;
 // {
  //  real norm = sqrt(dot_self(raw_effect)); // 计算向量的欧几里得范数
  //  effect = raw_effect / norm; // 规范化向量
 // }
 
  vector[N] mu;              // the RIM linear predictor for perform (here seed production)

  matrix[S, T] a_rim;   // interaction matrix for response - impact estimates (re in the manuscript)

  // get RIM interactions
  //weight * 
  a_rim = raw_response*raw_effect';  // the apostrophe transposes the effect vector
   
  // Estimate response-impact interactions
  for(n in 1:N) {

       mu[n] =(r[species_ID[n]] - dot_product(X[n], a_rim[species_ID[n], ]) + ran_t[time[n]])*time_diff[n];  
  }
   
} 

model {

  // priors
  r ~ cauchy(0, 10);   // prior for the intercept following Gelman 2008
  sigma ~ normal(0, 10);  // safer to place prior on disp_dev than on phi
  //weight ~ normal(0, 10); // constrained by parameter definition to be positive
  ran_t ~ normal(0, 1); 
  // no prior needed for response or effect as we can use the default prior for the unit_vector

  // maximising the likelihood for the joint model 
  for(n in 1:N) {
  
    // maximise the likelihood of the RIM for all observations
    log_perform[n] ~ normal(mu[n], sigma);
    
    // NB: in our case study, the response variable (seed production) shows a better fit to 
    // negative binomial (modify as necessary)
  }
  
}

generated quantities {

  vector[N] log_lik;	// log-likelihood for the RIM
  for(n in 1:N) {	
    log_lik[n] = normal_lpdf(log_perform[n] | mu[n], sigma);

  }
  
}

