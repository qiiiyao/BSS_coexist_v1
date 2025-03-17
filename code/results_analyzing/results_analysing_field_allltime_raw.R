# Load packages
rm(list = ls())
library(dplyr)
#library(betareg)
library(rstan)
library(loo)
library(stringr)
library(data.table)

load("code/data preparation/transformed data/fit_field_alltime_same_0.80_raw.RData")

setwd('D:/R projects/BSS/results/fit_results/BSS_exculde_trees_raw/')
# Read data
inter_l_2 = list.files('field_all_time_filter_1/parameters',
                       pattern = 'inter',
                       full.names = T)

inter_all = data.frame()
for (i in 1:length(inter_l_2)) {
  load(inter_l_2[i])
  inter_all = rbind(inter_all, inter_all_trans)
}

intra_l_2 = list.files('field_all_time_filter_1/parameters',
                       pattern = 'intra',
                       full.names = T)

# Filter the species fitting whose rhat<0.99/rhat>1.01
### Check the Bayesian R2
## Calculate bayesian R2 for rstan
posterior_l = list.files('field_all_time_filter_1/posterior',
                         full.names = T)

bayes_R2 <- function(mu, sigma) {
  var_mu <- apply(mu, 1, var)
  sigma2 <- sigma^2
  var_mu / (var_mu + sigma2)
}

loo_R2 = function(object, origindata) {
  object = fit
  origindata = sp_data
  require('loo', 'rstan')
  y = origindata[[1]]$G
  
  log_ratios = (-extract_log_lik(object, merge_chains = T))
  psis_object = loo::psis(log_ratios, r_eff = NA)
  
  mu = extract(object)$yhat
  sigma = extract(object)$sigma_e
  # generating posterior predictions
  mu_pred = matrix(nrow = dim(mu)[1], ncol = dim(mu)[2])
  for (i in 1:dim(mu)[1]) {     # for each posterior draw
    for (j in 1:dim(mu)[2]) {    # for each observation 
      # draw from the predicted distribution
      mu_pred[i, j] = rnorm(1, mean = mu[i, j], sd = sigma[i])  
    }
  }
  
  mu_pred_loo = loo::E_loo(mu_pred, psis_object, log_ratios = log_ratios)$value
  
  err_loo = mu_pred_loo - y
  
  S = nrow(mu_pred)
  N = ncol(mu_pred)
  
  # set the random seed as the seed used in the first chain and ensure
  # the old RNG state is restored on exit
  rng_state_old = .Random.seed
  on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  set.seed(object@stan_args[[1]]$seed)
  
  # dirichlet weights 
  exp_draws = matrix(rexp(S * N, rate = 1), nrow = S, ncol = N)
  wts = exp_draws / rowSums(exp_draws)
  
  var_y = (rowSums(sweep(wts, 2, y^2, FUN = "*")) -
             rowSums(sweep(wts, 2, y, FUN = "*"))^2) * (N/(N-1))
  
  var_err_loo = (rowSums(sweep(wts, 2, err_loo^2, FUN = "*")) -
                   rowSums(sweep(wts, 2, err_loo, FUN = "*")^2)) * (N/(N-1))
  
  loo_r_squared = 1 - var_err_loo / var_y
  mean(loo_r_squared)
  loo_r_squared[loo_r_squared < -1] = -1
  loo_r_squared[loo_r_squared > 1] = 1
  return(loo_r_squared)
}

## extract alpha and beta with 'permuted = TRUE' 
bayesian_r2_dat = data.frame()
for (i in 1:length(posterior_l)) {
  #i = 1
  load(posterior_l[[i]])
  plo = paste(str_split(str_split(str_split(posterior_l[[i]], '/')[[1]][3], '.rdata')[[1]][1],
                        '_')[[1]][1], collapse = '_')
  sp = paste(str_split(str_split(str_split(posterior_l[[i]], '/')[[1]][3], '.rdata')[[1]][1],
                       '_')[[1]][2:3], collapse = '_')
  fit_ss = extract(fit, permuted = TRUE) # fit_ss is a li
  plo_data = fit_field_alltime_same_0.80_raw[which(names(fit_field_alltime_same_0.80_raw) == plo)]
  sp_data = plo_data[[1]][[5]][which(sapply(plo_data[[1]][[5]], function(x){y = unique(gsub('[0-9]+','',names(x$G)))}) == sp)]
  bayesian_r2 = bayes_R2(fit_ss$yhat, fit_ss$sigma_e) %>% mean
  loo_r2 = loo_R2(object = fit, origindata = sp_data) %>% mean
  bayesian_r2_dat_1 = data.frame(plot = plo, species = sp, bayesian_r2 = bayesian_r2,
                                 loo_r2 = loo_r2)
  bayesian_r2_dat = rbind(bayesian_r2_dat, bayesian_r2_dat_1)
}
nrow(bayesian_r2_dat)
summary(bayesian_r2_dat)
save(bayesian_r2_dat, file = 'field_all_time_filter_1/bayesian_r2_dat.rdata')


summary_l1 = list.files('field_all_time_filter_1/summary', full.names = T)
bad_fit = data.frame()
for (i in 1:length(summary_l1)) {
  #i = 2
  load(summary_l1[[i]])
  load(intra_l_2[[i]])
  sp = intra_trans$sp
  for (k in 1:length(summary_l)) {
    # k = 1
    trans_summary = summary_l[[k]]
    variable = rownames(trans_summary)
    trans_summary_c = trans_summary %>% filter(grepl("^a", variable) | variable == 'r') %>%
      filter(Rhat < 0.99 | Rhat > 1.01) %>% filter(n_eff < 400)
    if (nrow(trans_summary_c) > 0) {
      field = str_split(str_split(str_split(summary_l1[[i]], '/')[[1]][4], '.rdata')[[1]][1],
                        'summary_')[[1]][2]
      species = sp[k]
      bad_fit_1 = data.frame(field = field, species = species)
      bad_fit = rbind(bad_fit_1, bad_fit)
    } else {print('no bad fit!')}
  }
}
nrow(bad_fit)/length(posterior_l) ### 0.101 bad fit rhat > 1.01

# Clear the bad fit
inter_all_l = split(inter_all, inter_all$field)
inter_all_l_c = list()
for (i in 1:length(inter_all_l)) {
  #i = 3
  trans = inter_all_l[[i]]
  f_field = unique(trans$field)
  sp_c1 = bad_fit %>% filter(field == f_field)
  if (nrow(sp_c1) != 0) {
    sp_c = sp_c1$species
    trans_c = trans %>% filter(!(species_i %in% sp_c) & !(species_j %in% sp_c))
    inter_all_l_c[[i]] = trans_c
  } else {inter_all_l_c[[i]] = trans}
}

inter_all_c = rbindlist(inter_all_l_c)

### Data fixed
inter_all_c[,6:ncol(inter_all_c)] = as.data.frame(apply(inter_all_c[,6:ncol(inter_all_c)], 2,
                                                        function(x){as.numeric(x)}))
inter_all_c$lgfd = log10(inter_all_c$fd)

### Add the invasion stage
inter_all_c$stage_i = 'native'
inter_all_c$stage_j = 'native'

inter_all_c_l = split(inter_all_c, inter_all_c$field)
for (i in 1:length(inter_all_c_l)) {
  #i = 3
  trans = inter_all_c_l[[i]]
  field_name = unique(trans$field)
  i_1 = which(names(fit_field_alltime_same_0.80_raw) == field_name)
  spname_stage = as.data.frame(fit_field_alltime_same_0.80_raw[[i_1]]$stage_sp_name_f)  
  
  sp_intro = strsplit(strsplit(spname_stage$intro_sp,
                               'Plot, f_p, ')[[1]][1], ', ')[[1]]
  sp_estab = strsplit(spname_stage$estab_sp,
                      ', ')[[1]]
  sp_domin = strsplit(spname_stage$domin_sp,
                      ', ')[[1]]
  
  if (length(sp_estab) == 1) {
    inter_all_c[inter_all_c$field == field_name &
                  inter_all_c$species_i == 're_cover_ab_estab_f[, 10:ncol(re_cover_ab_estab_f)]',]$species_i = sp_estab
    inter_all_c[inter_all_c$field == field_name &
                  inter_all_c$species_j == 're_cover_ab_estab_f[, 10:ncol(re_cover_ab_estab_f)]',]$species_j = sp_estab
    
  }
  
  if (length(sp_domin) == 1) {
    inter_all_c[inter_all_c$field == field_name &
                  inter_all_c$species_i == 're_cover_ab_domin_f[, 10:ncol(re_cover_ab_domin_f)]',]$species_i = sp_domin
    inter_all_c[inter_all_c$field == field_name &
                  inter_all_c$species_j == 're_cover_ab_domin_f[, 10:ncol(re_cover_ab_domin_f)]',]$species_j = sp_domin
    
  }
  
  if (length(sp_intro) == 1) {
    inter_all_c[inter_all_c$field == field_name &
                  inter_all_c$species_i == 're_cover_ab_intro_f[, 10:ncol(re_cover_ab_intro_f)]',]$species_i = sp_intro
    inter_all_c[inter_all_c$field == field_name &
                  inter_all_c$species_j == 're_cover_ab_intro_f[, 10:ncol(re_cover_ab_intro_f)]',]$species_j = sp_intro
  }
  
  inter_all_c[inter_all_c$field == field_name & inter_all_c$species_i %in% sp_intro,]$stage_i = 'introduce'
  inter_all_c[inter_all_c$field == field_name & inter_all_c$species_i %in% sp_estab,]$stage_i = 'establish'
  inter_all_c[inter_all_c$field == field_name & inter_all_c$species_i %in% sp_domin,]$stage_i = 'dominant'
  
  inter_all_c[inter_all_c$field == field_name & inter_all_c$species_j %in% sp_intro,]$stage_j = 'introduce'
  inter_all_c[inter_all_c$field == field_name & inter_all_c$species_j %in% sp_estab,]$stage_j = 'establish'
  inter_all_c[inter_all_c$field == field_name & inter_all_c$species_j %in% sp_domin,]$stage_j = 'dominant'
}

unique(inter_all_c$stage_i)
inter_all_c$stage_ij = paste(inter_all_c$stage_i, inter_all_c$stage_j, sep = '_')
inter_all_c$sp_pair = paste(inter_all_c$species_i, inter_all_c$species_j, sep = '_')
unique(inter_all_c$species_i)

## Add the different stage probability
inter_all_c$estab_i = 0
inter_all_c$domin_i = 0
unique(inter_all_c$stage_i)
inter_all_c[inter_all_c$stage_i %in% c('establish', 'dominant'), ]$estab_i = 1
inter_all_c[inter_all_c$stage_i %in% c('dominant'), ]$domin_i = 1

inter_all_c$estab_j = 0
inter_all_c$domin_j = 0
unique(inter_all_c$stage_j)
inter_all_c[inter_all_c$stage_j %in% c('establish', 'dominant'), ]$estab_j = 1
inter_all_c[inter_all_c$stage_j %in% c('dominant'), ]$domin_j = 1

## Add the gain or loss and relative abundance
load("D:/R projects/BSS/code/Phylo_Func/fd_dat_all.rdata")
load("D:/R projects/BSS/code/Phylo_Func/pd.rdata")
load("D:/R projects/BSS/code/data preparation/transformed data/sp_racover_f1_mean_field_alltime.rdata")

inter_all_c$f_p_species_i = paste(inter_all_c$f_p,
                                  inter_all_c$species_i,
                                  sep = '_')
inter_all_c$f_p_species_j = paste(inter_all_c$f_p,
                                  inter_all_c$species_j,
                                  sep = '_')

sp_racover_f1_mean_y_alltime$f_p_species_i = sp_racover_f1_mean_y_alltime$f_p_species
sp_racover_f1_mean_y_alltime$f_p_species_j = sp_racover_f1_mean_y_alltime$f_p_species_i

inter_all_c_test = inter_all_c %>% left_join(sp_racover_f1_mean_y_alltime[,c(3, 5)], 
                                             by = 'f_p_species_i') %>%
  left_join(sp_racover_f1_mean_y_alltime[,c(3, 6)], 
            by = 'f_p_species_j') %>% 
  left_join(fd_dat[,c(3, 4)], 
            by = 'sp_pair') %>% 
  left_join(pd_dat[,c(3, 4)], 
            by = 'sp_pair')

colname_1 = gsub('\\.x', '_i', colnames(inter_all_c_test))          
colname_2 = gsub('\\.y', '_j', colname_1)
colnames(inter_all_c_test) = colname_2
inter_all_c_alltime = inter_all_c_test

save(inter_all_c_alltime, file = 'fit_results/field_all_time_filter_1/inter_all_c_alltime.rdata')

