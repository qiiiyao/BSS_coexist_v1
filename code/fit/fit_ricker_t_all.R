# rm(list=ls())
# load packages
library('cmdstanr')
library("coda")
library("deSolve")
library("bayesplot")
library("parallel")
library("posterior")
library("dplyr")

# load('code/data preparation/invasion_stage.Rdata')
mod_ricker_2 = cmdstan_model(stan_file = 'code/fit/Ricker_random_time.stan')
numer_cores = parallel::detectCores()

# Fit
load('code/data preparation/fit_stage_fp_t_all.rdata')

samplingfunction = function(x) {
  res = mod_ricker_2$sample(
    data = datalist[[x]],seed = 123, 
    chains = 3, parallel_chains = 3,iter_warmup = 1000,
    iter_sampling = 2000,refresh = 10, # print update every 10 iters
    max_treedepth = 15, thin = 10,#init = inits,
    adapt_delta = 0.99,
    show_messages = F
  )
  #res$summary()
}

other_t = setdiff(c(1:length(inva_stage_fp_t)), c(1,7,14))
for (z in other_t) {
  ### year 
  data_t = inva_stage_fp_t[[z]]
  for (i in 1:length(data_t)) {
    stopCluster(cl)
    cl = makeCluster(numer_cores-2)
    #i = 67
  datalist = data_t[[i]][[5]]
  growD = data_t[[i]][[4]]
  #sp = colnames(growD)[10:ncol(growD)]
  #posterior_l = list()
  #summary_l = list()
  clusterExport(cl,c("mod_ricker_2","datalist"))
  out = parLapply(cl, c(1:length(datalist)),samplingfunction)
  save(out, file = paste0('results/fit_results/equal_moving_win/','equal_t', 
                           z, '/posterior/',
                           unique(growD$f_p), '.rdata'))
  }
}
