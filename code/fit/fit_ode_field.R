# set up R environment
rm(list = ls())
library(cmdstanr)
set_cmdstan_path('C:/Users/Qi Yao/Documents/.cmdstan/cmdstan-2.30.0')
options(mc.cores = parallel::detectCores(logical = F)) 
library(dplyr)
library(data.table)
library(parallel)
library(doParallel)


# load required functions
setwd("D:/R projects/BSS")
source('code/fit/functions/data_prep.R')
load("code/data preparation/transformed data/fit_field_top50_ages1_35.RData")

#file = file.path("D:/R projects/BSS/code/fit/stancode/RIM_random_time.stan")
#mod_cmd = cmdstan_model(file)
sps = seq(1, 2550, 51)
forms_ricker = c()
for (i in 1:length(sps)) {
  #i = 5
  r = paste0('theta[', sps[i], ']')
  sps1 = c(sps[i]:(sps[i]+50))
  a = c()
  for (j in 1:(length(sps1)-1)) {
    #j = 1
    a1 = paste0('theta[', sps1[j+1], ']*', 'z_init[', j, ']')
    a = paste0(a, '-', a1)
  }
  full_form = paste0(r, a)
  forms_ricker = c(forms_ricker, full_form)
}

file_ode_ricker = file.path("D:/R projects/BSS/code/fit/stancode/ode_Ricker.stan")
mod_ode_ricker = cmdstan_model(file_ode_ricker)

file_ode_ricker_multisp = file.path("D:/R projects/BSS/code/fit/stancode/ode_Ricker_multisp.stan")
mod_ode_ricker_multisp = cmdstan_model(file_ode_ricker_multisp)

stan.seed = 1234
numer_cores = parallel::detectCores(logical = F)
chains = 4 # 4 chains in parallel
sampling = 2500
warmup = 1000
thin = 10

# identify focal and neighbouring species to be matched to parameter estimates
data_t = fit_field_top50_ages1_35

posi_sp_start = 3

for (i in 1:10) {
  i = 2
  
  #df = fit_fp_top50_fake_ages1_35[[i]][[2]]
  df = fit_field_top50_ages1_35[[i]][[2]]
  sp_cover = df[,posi_sp_start:ncol(df)]
  sp_cover_prop = as.data.frame(t(apply(sp_cover, 1, function(x){y = x/sum(x)})))
  df = cbind(df[,1:(posi_sp_start-1)], sp_cover_prop)
  sp_01 = sp_cover
  sp_01[sp_01 != 1e-06] = 1
  sp_order = names(sort(colSums(sp_01), decreasing = T))
  sp_topf_5 = sp_order[c(1:50)]
  
  sp = colnames(df)[posi_sp_start:ncol(df)]
  sp_topf_5_posi = which(sp %in% sp_topf_5)
  
  N = length(df$fake_age) - 1
  S = length(sp_topf_5)
  ts = df$fake_age[-1]
  y_init = unlist(df[1,sp_topf_5])   #
  t_init = df$fake_age[1]
  y = as.matrix(df[2:(N+1),sp_topf_5])
  stan.data1 = list(N = N, S = S, ts = ts, y_init = y_init, t_init = t_init, y = y)
  
  plot(1:nrow(df), df[,sp_topf_5][,1], ylim = c(0, 1),
       xlab = 'age',
       ylab = 'cover', xaxt = 'n', type = 'l', lwd = 1.5)
  at = c(1,18,35)
  axis( 1,at=at,labels=df$fake_age[at])
  lines( 1:nrow(df),df[,sp_topf_5][,2],lwd=1.5,col=rangi2)
  points( 1:29,df[,sp_topf_5][,1],bg="black",col="white",pch=21,cex=1.4)
  points( 1:29,df[,sp_topf_5][,2],bg=rangi2,col="white",pch=21,cex=1.4)
  
  
  # Run the model using cmdstanr
  fit = mod_ode_ricker_multisp$sample(data = stan.data1,
                              seed = stan.seed, 
                              chains = chains, 
                              parallel_chains = chains,
                              iter_warmup = warmup,
                              iter_sampling = sampling,
                              refresh = 100, # print update every 100 iters
                              thin = thin,
                              max_treedepth = 20,
                              adapt_delta = 0.99)
  fit$summary()
  
  fit = mod_ricker$sample(data = stan.data2,
                          seed = stan.seed, 
                          chains = chains, 
                          parallel_chains = chains,
                          iter_warmup = warmup,
                          iter_sampling = sampling,
                          refresh = 100, # print update every 100 iters
                          thin = thin,
                          max_treedepth = 20,
                          adapt_delta = 0.99)
  fit$summary()
  
  
  fit$save_object(file = paste0('fit_results/opti_bimler',
                                '/posterior/','fit_',
                                f_p, '.RDS'))
}

files = list.files("/data/home/shpli3/R_projects/BSS_exclude_tree_raw/fit_results/opti_bimler/posterior", full.names = T)
trans_summary_c_all = data.frame()
for (i in 1:length(files)) {
  
  fit = readRDS(files[i])
  trans_summary = fit$summary()
  trans_summary_c = trans_summary %>% filter(grepl("^gamma", variable) | grepl("^beta", variable))  %>% filter(rhat < 0.95 | rhat > 1.05) 
  trans_summary_c_all = rbind(trans_summary_c_all, trans_summary_c)
  
}