data{
  int<lower = 1> T; // Number of Timestep
  int<lower = 1> S; // Number of species
  real Fecundity[T]; // Fecundity of the focal species in each plot
  //int Fecundity[T]; // Fecundity of the focal species in each plot
  #int reserve[T];   // Indicator variable for the reserve each plot is located in
  matrix[T,S] SpMatrix; // Matrix of abundances for each species (including abundances of non-focal individuals of the focal species)
  //vector[T] env;   // Environmental values for each plot
  int<lower = 0> Intra[S]; // Indicator boolean variable to identify the focal species (0 for non-focal and 1 for focal). Included for easier calculations
  // The below values define the regularized horseshoe priors used for species-specific parameters
  int Inclusion_ij[S];  // Boolean indicator variables to identify the species x reserve intercept parameters identified for inclusion in the final model
}

parameters{
  real lambdas;
  real alpha_generic_tilde;
  real alpha_intra_tilde;
  vector[S] alpha_hat_ij;
  real<lower = 0> theta;
}

transformed parameters{
  real alpha_generic;
  real alpha_intra;

  // scale the alpha values
  alpha_generic = 1 * alpha_generic_tilde - 0;
  alpha_intra = 3 * alpha_intra_tilde - 6;
  //alpha_generic[2] = 0.5 * alpha_generic_tilde[2];
  //alpha_intra[2] = 0.5 * alpha_intra_tilde[2];
}

model{
  // Declare objects necessary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //     a matrix of the species specific alpha values for each species and plot (interaction_effects), and a matrix
  //     of the the alpha*T values for each species.
  vector[T] F_hat;
  vector[T] interaction_effects;
  row_vector[S] alpha_ij;
  real lambda_i;
  
  // set regular priors
  alpha_generic_tilde ~ normal(0,1);
  alpha_intra_tilde ~ normal(0,1);
  lambdas ~ normal(0, 1);
  
  theta ~ normal(0, 1);
  
  // implement the biological model
  for(i in 1:T){
    lambda_i = exp(lambdas);
    for(s in 1:S){
        alpha_ij[s] = exp((1-Intra[s]) * alpha_generic + Intra[s] * alpha_intra + Inclusion_ij[s] * alpha_hat_ij[s]);
    }
    
    interaction_effects[i] = sum(alpha_ij.* SpMatrix[i,]);
    F_hat[i] = lambda_i-interaction_effects[i];
  }
  Fecundity ~ normal(F_hat, theta);
  //Fecundity ~ poisson(F_hat);
}

