### Dave Armitage (david.armitage@oist.jp)
### Code adapted from https://doi.org/10.5281/zenodo.7083314 by M. Van Dyke
### Please download the author's original data from the link above prior to running the analysis
### Last edit: 20 Dec 2022

# This script fits and compares Bayesian models of competition to the data provided in Van Dyke et al.
# Where possible, it is an exact replica of the bootstrap analysis performed in the original paper
rm(list = ls())
library(dplyr)
library(cmdstanr)
set_cmdstan_path('C:/Users/Qi Yao/Documents/.cmdstan/cmdstan-2.30.0')
library(tidyr)
library(loo)
library(parallel)

options(mc.cores = detectCores(logical = F))
set.seed(1234)
setwd("D:/R projects/BSS")

### load work dictionary
file_null = file.path("code/fit/stancode/model_comprison/Exp_r_random_time.stan")
mod_null = cmdstan_model(file_null
                         , cpp_options = list(stan_threads = TRUE)
)
file_lv = file.path("code/fit/stancode/model_comprison/Lotka_volterra_random_time_loo.stan") 
mod_lv = cmdstan_model(file_lv, cpp_options = list(stan_threads = TRUE))
file_ricker = file.path("code/fit/stancode/model_comprison/Ricker_random_time_loo.stan")
mod_ricker = cmdstan_model(file_ricker, cpp_options = list(stan_threads = TRUE))
file_ricker_logn1 = file.path("code/fit/stancode/model_comprison/Ricker_log(n+1)_random_time_loo.stan")
mod_ricker_logn1 = cmdstan_model(file_ricker_logn1, cpp_options = list(stan_threads = TRUE))
file_bh = file.path("code/fit/stancode/model_comprison/BH_random_time_loo.stan")
mod_bh = cmdstan_model(file_bh, cpp_options = list(stan_threads = TRUE))
file_bh_exp = file.path("code/fit/stancode/model_comprison/BH_n^a_random_time.stan")
mod_bh_exp = cmdstan_model(file_bh_exp, cpp_options = list(stan_threads = TRUE))
file_bh_allexp = file.path("code/fit/stancode/model_comprison/BH_allb_random_time_loo.stan")
mod_bh_allexp = cmdstan_model(file_bh_allexp, cpp_options = list(stan_threads = TRUE))


# Load data
load("code/data preparation/transformed data/fit_fp_top50_ages1_35.RData")

model_compare = function(x) {
  #x = dat_fit[[1]]
  data = list(N = x$N,
              P = x$P,
              X = x$X,
              ftime = x$ftime,
              time = x$time,
              time_diff = x$time_diff,
              G = x$G)
  
  focal_sp = gsub("1", "", names(x$G)[1])
  
  # Next, we fit the bayesian models for comparison purposes. 
  # This code will fit a bunch of models so it can take some time to run
  # It takes around 1 hour on an M1 chip with 8 cores in parallel.
  fit_null = mod_null$sample(data = data, seed = 123,
                             chains = chains,
                             #  init = 50,
                             parallel_chains = chains,
                             threads_per_chain = 5,
                             iter_warmup = warmup,
                             iter_sampling = sampling,
                             refresh = 10, # print update every 10 iters
                             thin = thin,
                             max_treedepth = 20,
                             adapt_delta = 0.99)
  fit_null_summary = as.data.frame(fit_null$summary(variables = c('r')))
  fit_null_summary$model = 'null'
  
  fit_lv = mod_lv$sample(data = data, seed = 123,
                         chains = chains,
                         #  init = 50,
                         parallel_chains = chains,
                         threads_per_chain = 5,
                         iter_warmup = warmup,
                         iter_sampling = sampling,
                         refresh = 10, # print update every 10 iters
                         thin = thin,
                         max_treedepth = 20,
                         adapt_delta = 0.99)
  fit_lv_summary = as.data.frame(fit_lv$summary()) %>% filter(variable == 'r' | 
                                                                grepl('a\\[.', variable))
  fit_lv_summary$model = 'lv'
  
  fit_ricker = mod_ricker$sample(data = data, seed = 123,
                                 chains = chains,
                                 #  init = 50,
                                 parallel_chains = chains,
                                 threads_per_chain = 5,
                                 iter_warmup = warmup,
                                 iter_sampling = sampling,
                                 refresh = 10, # print update every 10 iters
                                 thin = thin,
                                 max_treedepth = 20,
                                 adapt_delta = 0.99)
  fit_ricker_summary = as.data.frame(fit_ricker$summary()) %>% filter(variable == 'r' | 
                                                                        grepl('a\\[.', variable))
  fit_ricker_summary$model = 'ricker'
  
  fit_ricker_n1 = mod_ricker_logn1$sample(data = data, seed = 123,
                                          chains = chains,
                                          #  init = 50,
                                          parallel_chains = chains,
                                          threads_per_chain = 5,
                                          iter_warmup = warmup,
                                          iter_sampling = sampling,
                                          refresh = 10, # print update every 10 iters
                                          thin = thin,
                                          max_treedepth = 20,
                                          adapt_delta = 0.99)
  fit_ricker_n1_summary = as.data.frame(fit_ricker_n1$summary()) %>% filter(variable == 'r' | 
                                                                              grepl('a\\[.', variable))
  fit_ricker_n1_summary$model = 'ricker_n1'
  
  #fit_bh = mod_bh$sample(data = data, seed = 123,
  #   chains = chains,
  #  init = 50,
  #   parallel_chains = chains,
  #   threads_per_chain = 5,
  #   iter_warmup = warmup,
  #   iter_sampling = sampling,
  #  refresh = 10, # print update every 10 iters
  #  thin = thin,
  #  max_treedepth = 20,
  # adapt_delta = 0.99)
  #fit_bh_exp = mod_bh_exp$sample(data = data, seed = 123,
  #  chains = chains,
  #  init = 50,
  #   parallel_chains = chains,
  #   threads_per_chain = 5,
  #    iter_warmup = warmup,
  #    iter_sampling = sampling,
  #    refresh = 10, # print update every 10 iters
  #   thin = thin,
  #   max_treedepth = 20,
  #   adapt_delta = 0.99)
  #fit_bh_allexp = mod_bh_allexp$sample(data = data, seed = 123,
  #  chains = chains,
  #  init = 50,
  #   parallel_chains = chains,
  #  threads_per_chain = 5,
  #  iter_warmup = warmup,
  #  iter_sampling = sampling,
  #  refresh = 10, # print update every 10 iters
  #  thin = thin,
  #  max_treedepth = 20,
  #  adapt_delta = 0.99)
  
  fit_null$save_object(file = paste0('results/fit_results/model_comparison/null', '/posterior/',
                                     f_p,'_', focal_sp, '.RDS'))
  fit_lv$save_object(file = paste0('results/fit_results/model_comparison/lv', '/posterior/',
                                   f_p,'_', focal_sp, '.RDS'))
  fit_ricker$save_object(file = paste0('results/fit_results/model_comparison/ricker',
                                       '/posterior/',
                                       f_p,'_', focal_sp, '.RDS'))
  fit_ricker_n1$save_object(file = paste0('results/fit_results/model_comparison/ricker_n1',
                                          '/posterior/',
                                          f_p,'_', focal_sp, '.RDS'))
  #fit_bh$save_object(file = paste0('results/fit_results/model_comparison/bh', '/posterior/',
  #                           f_p,'_', focal_sp, '.RDS'))
  #fit_bh_exp$save_object(file = paste0('results/fit_results/model_comparison/bh_exp',
  #               '/posterior/',
  #       f_p,'_', focal_sp, '.RDS'))
  #fit_bh_allexp$save_object(file = paste0('results/fit_results/model_comparison/bh_allexp',
  #              '/posterior/',
  #             f_p,'_', focal_sp, '.RDS'))
  
  df.loo = as.data.frame(loo_compare(
    list(
      null = fit_null$loo(),
      lv = fit_lv$loo(),
      ricker = fit_ricker$loo(),
      ricker_n1 = fit_ricker_n1$loo()
      #,m5 = fit_bh$loo(),
      #m6 = fit_bh_exp$loo(),
      #m7 = fit_bh_allexp$loo()
    )
  )
  )
  df.loo$models = rownames(df.loo)
  df.loo$trust = 1
  
  ### according to https://mc-stan.org/rstan/reference/Rhat.html
  if(mean(fit_null_summary$rhat) > 1.05 | mean(fit_null_summary$rhat) < 0.95 |
     mean(fit_null_summary$ess_bulk) < 400) {
    df.loo['null', 'trust'] = 0
  }
  
  if(mean(fit_lv_summary$rhat) > 1.05 | mean(fit_lv_summary$rhat) < 0.95 |
     mean(fit_lv_summary$ess_bulk) < 400) {
    df.loo['lv', 'trust'] = 0
  }
  
  if(mean(fit_ricker_summary$rhat) > 1.05 | mean(fit_ricker_summary$rhat) < 0.95 |
     mean(fit_ricker_summary$ess_bulk) < 400) {
    df.loo['ricker', 'trust'] = 0
  }
  
  if(mean(fit_ricker_n1_summary$rhat) > 1.05 | mean(fit_ricker_n1_summary$rhat) < 0.95 |
     mean(fit_ricker_n1_summary$ess_bulk) < 400) {
    df.loo['ricker_n1', 'trust'] = 0
  }
  
  df.summary = rbind(fit_null_summary, fit_lv_summary,
                     fit_ricker_summary, fit_ricker_n1_summary)
  df.loo.summary = list(loo = df.loo, summary = df.summary)
  
  return(df.loo.summary)
  
}

#Create a data frame with all combinations of treatment and species 
#spp_list = sort(na.omit( unique(seed_data$focal)))
#treat_list = sort( na.omit( unique(seed_data$Tr)))
#spp_combos = expand.grid(species = spp_list, treatment = treat_list)
# loop through the species list, fit the model, return the dataframe 
#out_bayes_compare.waic = list() #list for for loop output
out_bayes_compare.loo = list() #list for for loop output
numer_cores = parallel::detectCores(logical = F)
chains = 4 # 4 chains in parallel
sampling = 2500
warmup = 1000
thin = 10

plot_list = sapply(seq(1, 480, 48), function(x){y = c(x:(x+47))})

for(i in plot_list[,2]){ 
  #i = 1
  dat_fit = fit_fp_top50_ages1_35[[i]][[3]]
  out_bayes_compare.waic_1 = list()
  out_bayes_compare.loo_1 = list()
  f_p = unique(fit_fp_top50_ages1_35[[i]]$re_cover_ab$f_p)
  sps = colnames(dat_fit[[1]]$X)
  
  out_bayes_compare.loo_1 = mclapply(dat_fit, model_compare)
  
  names(out_bayes_compare.loo_1) = sps
  save(out_bayes_compare.loo_1,
       file = paste0("results/fit_results/model_comparison/out_bayes_compare_",
                     f_p, ".loo.rdata"))
  
}   




