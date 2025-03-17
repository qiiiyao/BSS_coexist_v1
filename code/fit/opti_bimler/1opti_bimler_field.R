# set up R environment
rm(list = ls())
library(cmdstanr)
options(mc.cores = parallel::detectCores(logical = F)) 
library(dplyr)
library(data.table)
library(parallel)
library(doParallel)


# load required functions
setwd("D:/R projects/BSS")
source('code/fit/functions/data_prep.R')
load("code/data preparation/transformed data/fit_field_top50_ages1_35.RData")

file = file.path("code/fit/stancode/joint2_model2_aposi.stan")
mod_cmd = cmdstan_model(file)
file_2 = file.path("code/fit/stancode/joint2_model2_r&aposi.stan")
mod_cmd_2 = cmdstan_model(file_2)
file_3 = file.path("code/fit/stancode/joint2_model2_norant_aposi.stan")
mod_cmd_3 = cmdstan_model(file_3)
file_4 = file.path("code/fit/stancode/joint2_model2_ar_aposi.stan")
mod_cmd_4 = cmdstan_model(file_3)

stan.seed = 1234
numer_cores = parallel::detectCores(logical = F)
chains = 4 # 4 chains in parallel
sampling = 2500
warmup = 1000
thin = 10

# identify focal and neighbouring species to be matched to parameter estimates
plot_list = sapply(seq(1, 480, 48), function(x){y = c(x:(x+47))})
setwd("D:/R projects/BSS")

posi_sp_start = 3
for (i in plot_list[,1]) {
  i = 1
  
  #df = fit_fp_top50_ages1_35[[i]][[2]]
  df = fit_field_top50_ages1_35[[i]][[2]]
  sp_cover = df[,posi_sp_start:ncol(df)]
  sp_01 = sp_cover
  sp_01[sp_01 != 1e-06] = 1
  sp_order = names(sort(colSums(sp_01), decreasing = T))
  sp_topf_5 = sp_order[1:50]
  
  sp = colnames(df)[posi_sp_start:ncol(df)]
  sp_topf_5_posi = which(sp %in% sp_topf_5)
  
  field = unique(df$Field)
  df_l = lapply(fit_field_top50_ages1_35[[i]][[3]], function(x){
    #x = fit_fp_top50_ages1_35[[i]][[3]][[1]]
    y = cbind(
      focal = gsub('[0-9]+',
                   '',
                   names(x$G)),
      time = x$time,
      time_diff = x$time_diff,
      G = x$G, 
      x$X[,sp_topf_5_posi])})
  
  df = as.data.frame(rbindlist(df_l))
  focalID = unique(df$focal)  # this should return the names of unique focal groups in the order
  # in which they are encountered in the dataframe
  neighbourID = colnames(df[,-c(1:4)])
  
  
  # ensure neighbours are linearly independent across the whole dataset (see S1.2)
  N_all = apply(df[,c(5:ncol(df))], c(1,2), as.numeric)
  X_all = cbind(model.matrix(~as.factor(df$focal)), N_all)
  R_all = pracma::rref(X_all)
  Z_all = t(R_all) %*% R_all
  indep = sapply(seq(1, dim(Z_all)[1], 1), function(k){ 
    ifelse(Z_all[k, k] == 1 & sum(Z_all[k, -k]) == 0, 1, 0)
  }) #
  all(indep == 1) # if TRUE then neighbours are linearly independent and we can continue
  if(!all(indep == 1)) warning('WARNING neighbours are not linearly independent') 
  rm(N_all, X_all, R_all, Z_all)
  
  # prepare the data into the format required by STAN and the model code
  stan.data = data_prep(perform = 'G', 
                        focal = 'focal', 
                        nonNcols = 4, # number of columns that aren't neighbour abundances
                        df = df)
  
  message(paste0('Data dimensions = ', dim(df)[1], ', ', dim(df)[2]))
  message(paste0('Number of focal groups = ', length(focalID)))
  message(paste0('Number of neighbour groups = ', length(neighbourID)))
  message(paste0('Proportion of inferrable interactions = ',
                 sum(stan.data$Q)/(stan.data$S*stan.data$'T')))
  
  # Run the model using rstan! 
  #fit = sampling(mod,
  #          data = stan.data, seed = stan.seed,  
  #           chains = 4,
  #        warmup = 5000,          # number of warmup iterations per chain
  #          iter = 7000,            # total number of iterations per chain
  #           refresh = 100,         # show progress every 'refresh' iterations
  #          thin = 10,
  #          control = list(max_treedepth = 20,
  #                        adapt_delta = 0.99))
  #saveRDS(fit, file = paste0('results/fit_results/joint_plot_ages1_35_top50_aposi', '/posterior/',
  #                                 field, '.RDS'))
  
  
  
  # Save results 
  #summary = summary(fit)
  #head(fit$summary())
  #save(summary, file=paste0('results/fit_results/joint_plot_ages1_35_top50_aposi',
  #                   '/summary/','summary_',
  #                   field, '.rdata'))
  # Run the model using cmdstanr
  fit = mod_cmd_4$sample(data = stan.data,
                       seed = stan.seed, 
                       chains = chains, 
                       parallel_chains = chains,
                       iter_warmup = warmup,
                       iter_sampling = sampling,
                       refresh = 100, # print update every 100 iters
                       thin = thin,
                       max_treedepth = 20,
                       adapt_delta = 0.99)
  head(fit$summary())
  summary = fit$summary()
  summary_r_a = summary %>% filter(grepl("^gamma", variable) | grepl("^response", variable))
  fit$save_object(file = paste0('fit_results/opti_bimler',
                                '/posterior/','fit_',
                                field, '.RDS'))
}

files = list.files("/data/home/shpli3/R_projects/BSS_exclude_tree_raw/fit_results/opti_bimler/posterior", full.names = T)
trans_summary_c_all = data.frame()
for (i in 1:length(files)) {
  
  fit = readRDS(files[i])
  trans_summary = fit$summary()
  trans_summary_c = trans_summary %>% filter(grepl("^gamma", variable) | grepl("^beta", variable))  %>% filter(rhat < 0.95 | rhat > 1.05) 
  trans_summary_c_all = rbind(trans_summary_c_all, trans_summary_c)
  
}