// reciprocal-model-1sp.stan
functions {
  real[] ode_Ricker(real t,       // time
               real[] z_init,   // system state {competitor1, cometitor2}
               real[] theta, // parameters
               real[] x_r,   // unused data
               int[] x_i) {
    
    real dydt[2];
    
    dydt[1] = ( exp(theta[1] - theta[2] * z_init[1] - theta[3] * z_init[2] )) * z_init[1]; // Competitor 1
    dydt[2] = ( exp(theta[4] - theta[5] * z_init[2] - theta[6] * z_init[1] )) * z_init[2]; // Competitor 2
    
    return dydt;
  }
}

data {
  int<lower = 0> N;          // number of measurement times
  real ts[N];                // measurement times > 0
  real y_init[2];            // initial measured populations
  real t_init;               // initial time
  real<lower = 0> y[N, 2];   // measured populations
}

transformed data{
  real times_measured[N];    // N-1 because first time is initial state
  for ( i in 2:N+1 ) times_measured[i-1] = i;
}

parameters {
  real theta[6];   // { alpha, beta}
  real<lower = 0> z_init[2];  // initial population
  real<lower = 0> sigma[2];   // measurement errors
}

transformed parameters {
  real z[N, 2]
  = integrate_ode_rk45(ode_Ricker, z_init, t_init, times_measured, theta,
                         rep_array(0.0, 0), rep_array(0, 0),
                         1e-5, 1e-3, 5e2);
}

model {
  theta[{1, 4}] ~ normal(0, 1);
  theta[{2, 3, 5, 6}] ~ normal(0, 1);
  sigma ~ lognormal(-1, 5);
  z_init ~ lognormal(log(1), 1);
  for (k in 1:2) {
    y_init[k] ~ lognormal(log(z_init[k]), sigma[k]);
    y[ , k] ~ lognormal(log(z[, k]), sigma[k]);
  }
}

generated quantities {
  real y_init_rep[2];
  real y_rep[N, 2];
  for (k in 1:2) {
    y_init_rep[k] = lognormal_rng(log(z_init[k]), sigma[k]);
    for (n in 1:N)
      y_rep[n, k] = lognormal_rng(log(z[n, k]), sigma[k]);
  }
}
