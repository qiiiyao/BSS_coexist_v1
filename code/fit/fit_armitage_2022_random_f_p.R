### Dave Armitage (david.armitage@oist.jp)
### Code adapted from https://doi.org/10.5281/zenodo.7083314 by M. Van Dyke
### Please download the author's original data from the link above prior to running the analysis
### Last edit: 20 Dec 2022

# This script fits and compares Bayesian models of competition to the data provided in Van Dyke et al.
# Where possible, it is an exact replica of the bootstrap analysis performed in the original paper
rm(list = ls())
library(tidyverse)
library(nlstools)
library(dplyr)
library(rsample)
#library(broomExtra)
library(purrr)
library(ggplot2)
library(ggrepel)
library(HDInterval)
library(ggridges)
library(tidybayes)
library(brms)
library(svMisc)
library(parallel)
library(foreach)
library(doParallel)
library(cowplot)
library(tidyr)
library(ggfortify)
library(vegan)
library(ggpubr)
library(gdata)

options(scipen = 5)
set.seed(1234)
setwd("D:/R projects/BSS")

# Load data
load("code/data preparation/transformed data/fit_fp_same_ages_top50_early_suc.RData")

temp_data = as.data.frame(cbind(G = fit_fp_same_ages_top50_early_suc[[1]][[3]][[1]]$G, 
                                time = fit_fp_same_ages_top50_early_suc[[1]][[3]][[1]]$time,
                                fit_fp_same_ages_top50_early_suc[[1]][[3]][[1]]$X))

focal_sp = gsub("[^A-Za-z]", "", rownames(temp_data)[1])
background_sp = sapply(colnames(temp_data)[3:ncol(temp_data)],
                       function(x){y = gsub("[^A-Za-z]", "", x)})
colnames(temp_data)[3:ncol(temp_data)] = background_sp
a_vector = paste('a',background_sp, sep = '')

### Create 7 models' formula
## 2. Pure Lotka-Volterra model
formula_lv_fixed = "G ~ lambda"

for (i in 1:length(a_vector)) {
  AN = paste(a_vector[i], background_sp[i], sep = "*")
  formula_lv_fixed = paste(formula_lv_fixed, AN, sep = "-")
}

formula_lv_fixed = gsub(" \\+$", "", formula_lv_fixed)
formula_lv_fixed = as.formula(formula_lv_fixed)

## 3. Pure Ricker model
#G~log((lambda)*exp(-aACWR*(N_acwr)-aFEMI*(N_femi)-aHOMU*(N_homu)-aPLER*(N_pler)-aSACO*(N_saco)-aURLI*(N_urli)))
formula_ricker_fixed_1 = "G ~ log(lambda*exp("

for (i in 1:length(a_vector)) {
  AN = paste(a_vector[i], background_sp[i], sep = "*")
  formula_ricker_fixed_1 = paste(formula_ricker_fixed_1, AN, sep = "-")
}

formula_ricker_fixed = paste(formula_ricker_fixed_1, '))', sep = '')
formula_ricker_fixed = as.formula(formula_ricker_fixed)

## 4. Log(N+1) Ricker model
# G~log((lambda)*exp(-aACWR*log(N_acwr+1)-aFEMI*log(N_femi+1)-aHOMU*log(N_homu+1)
# -aPLER*log(N_pler+1)-aSACO*log(N_saco+1)-aURLI*log(N_urli+1)))
formula_logricker_fixed_1 = "G ~ log(lambda*exp("

for (i in 1:length(a_vector)) {
  AN = paste(a_vector[i], '*log(', background_sp[i], '+1)', sep = "")
  formula_logricker_fixed_1 = paste(formula_logricker_fixed_1, AN, sep = "-")
}

formula_logricker_fixed = paste(formula_logricker_fixed_1, '))', sep = '')
formula_logricker_fixed = as.formula(formula_logricker_fixed)

## 5. Pure BH model
#G~log((lambda)/(1+aACWR*N_acwr+aFEMI*N_femi+aHOMU*N_homu+aPLER*N_pler+aSACO*N_saco+aURLI*N_urli))
formula_bh_fixed_1 = "G ~ log(lambda/(1"

for (i in 1:length(a_vector)) {
  AN = paste(a_vector[i], background_sp[i], sep = "*")
  formula_bh_fixed_1 = paste(formula_bh_fixed_1, AN, sep = "+")
}

formula_bh_fixed = paste(formula_bh_fixed_1, '))', sep = '')
formula_bh_fixed = as.formula(formula_bh_fixed)

## 6.  exponential competition BH model
# G~log(lambda/(1+(N_acwr^aACWR)+(N_femi^aFEMI)+(N_homu^aHOMU)+
# (N_pler^aPLER)+(N_saco^aSACO)+(N_urli^aURLI)))
formula_e_bh_fixed_1 = "G ~ log(lambda/(1"

for (i in 1:length(a_vector)) {
  AN = paste('(',a_vector[i], '^', background_sp[i], ')', sep = "")
  formula_e_bh_fixed_1 = paste(formula_e_bh_fixed_1, AN, sep = "+")
}

formula_e_bh_fixed = paste(formula_e_bh_fixed_1, '))', sep = '')
formula_e_bh_fixed = as.formula(formula_e_bh_fixed)

## 7.  exponential all parts BH model
# G~log((lambda/((1+(aACWR*N_acwr+aFEMI*N_femi+aHOMU*N_homu+aPLER*N_pler+
# aSACO*N_saco+aURLI*N_urli)^theta))))
formula_eall_bh_fixed_1 = "G ~ log((lambda/((1+("

for (i in 1:length(a_vector)) {
  AN = paste(a_vector[i], background_sp[i], sep = "*")
  formula_eall_bh_fixed_1 = paste(formula_eall_bh_fixed_1, AN, sep = "+")
}
formula_eall_bh_fixed_1 = gsub("\\(\\+", "(", formula_eall_bh_fixed_1)
formula_eall_bh_fixed = paste(formula_eall_bh_fixed_1, ')^theta))))', sep = '')
formula_eall_bh_fixed = as.formula(formula_eall_bh_fixed)

## Model random effects parts
formula_string_random = "lambda"

for (i in 1:length(a_vector)) {
  formula_string_random = paste(formula_string_random, "+", a_vector[i], sep = "")
}

formula_string_random = gsub(" \\+$", "", formula_string_random)
formula_string_random = paste(formula_string_random, '~1','')
formula_string_random_7 = paste(formula_string_random, '+ theta ~ 1')
formula_string_random = as.formula(formula_string_random)
formula_string_random_7 = as.formula(formula_string_random_7)

# Model specifications following table 1 in the Armitage_2022_bioRix
seed_bayes_fit_formula_1 = bf(G~log((lambda)),
                              lambda ~ 1, nl = TRUE)

seed_bayes_fit_formula_2 = bf(formula_lv_fixed,
                              formula_string_random, nl = TRUE)

seed_bayes_fit_formula_3 = bf(formula_ricker_fixed,
                              formula_string_random, nl = TRUE)

seed_bayes_fit_formula_4 = bf(formula_logricker_fixed,
                              formula_string_random, nl = TRUE)

seed_bayes_fit_formula_5 = bf(formula_bh_fixed,
                              formula_string_random, nl = TRUE)

seed_bayes_fit_formula_6 = bf(formula_e_bh_fixed,
                              formula_string_random, nl = TRUE)

seed_bayes_fit_formula_7 = bf(formula_eall_bh_fixed,
                              formula_string_random_7, nl = TRUE)

## Set priors
prior_list = list()

for (i in 1:length(background_sp)) {
  prior_list[[i]] = list("normal(.1, 5)",
                         nlpar = a_vector[i], lb = .001, ub=2)
  
}

priors = set_prior("normal(500, 10000)", nlpar = "lambda", lb = 0, ub = 5000)
for (i in 1:length(prior_list)) {
  #i = 2
  priors = c(priors, do.call(set_prior, prior_list[[i]]))
}

seedpriors = c(priors,
               set_prior("normal(0, 2)", nlpar = "theta")
)

# Next, we fit the bayesian models for comparison purposes. 
# This code will fit a bunch of models so it can take some time to run
# It takes around 1 hour on an M1 chip with 8 cores in parallel.

#Create a data frame with all combinations of treatment and species 
spp_list = sort(na.omit( unique(seed_data$focal)))
treat_list = sort( na.omit( unique(seed_data$Tr)))
spp_combos = expand.grid(species = spp_list, treatment = treat_list)
# loop through the species list, fit the model, return the dataframe 
out_bayes_compare.waic = list() #list for for loop output
out_bayes_compare.loo = list() #list for for loop output
ncores = 8 # set no. of cores for parallel MCMC runs.

for( i in 1:nrow(spp_combos)){ 
  temp_data = seed_data %>%
    filter(focal == spp_combos[i,1], Tr ==spp_combos[i,2])
  
  seed_bayes_fit_1= brm(seed_bayes_fit_formula_1, family=gaussian(), data = temp_data, prior = seedpriors[1,],  chains = 8, cores = ncores, iter = 10000, backend = 'cmdstanr')  %>%  add_criterion("loo")
  seed_bayes_fit_2= brm(seed_bayes_fit_formula_2, family=gaussian(), data = temp_data, prior = seedpriors[-nrow(seedpriors),],  chains = 8, cores = ncores, iter = 10000, backend = 'cmdstanr') %>% add_criterion("loo")
  seed_bayes_fit_3= brm(seed_bayes_fit_formula_3, family=gaussian(), data = temp_data, prior = seedpriors[-nrow(seedpriors),],  chains = 8, cores = ncores, iter = 10000, backend = 'cmdstanr') %>% add_criterion("loo")
  seed_bayes_fit_4= brm(seed_bayes_fit_formula_4, family=gaussian(), data = temp_data, prior = seedpriors[-nrow(seedpriors),],  chains = 8, cores = ncores, iter = 10000, backend = 'cmdstanr') %>% add_criterion("loo")
  seed_bayes_fit_5= brm(seed_bayes_fit_formula_5, family=gaussian(), data = temp_data, prior = seedpriors[-nrow(seedpriors),],  chains = 8, cores = ncores, iter = 10000, backend = 'cmdstanr') %>% add_criterion("loo")
  seed_bayes_fit_6= brm(seed_bayes_fit_formula_6, family=gaussian(), data = temp_data, prior = seedpriors[-nrow(seedpriors),],  chains = 8, cores = ncores, iter = 10000, backend = 'cmdstanr') %>% add_criterion("loo")
  seed_bayes_fit_7= brm(seed_bayes_fit_formula_7, family=gaussian(), data = temp_data, prior = seedpriors,  chains = 8, cores = ncores, iter = 10000, backend = 'cmdstanr') %>% add_criterion("loo")
  
  df.loo = as.data.frame(loo_compare(
    list(
      m1 = seed_bayes_fit_1$criteria$loo,
      m2 = seed_bayes_fit_2$criteria$loo,
      m3 = seed_bayes_fit_3$criteria$loo,
      m4 = seed_bayes_fit_4$criteria$loo,
      m5 = seed_bayes_fit_5$criteria$loo,
      m6 = seed_bayes_fit_6$criteria$loo,
      m7 = seed_bayes_fit_7$criteria$loo
    )
  )
  )
  df.loo$models = rownames(df.loo)
  
  df.waic = as.data.frame(loo_compare(
    list(
      m1 = waic(seed_bayes_fit_1),
      m2 = waic(seed_bayes_fit_2),
      m3 = waic(seed_bayes_fit_3),
      m4 = waic(seed_bayes_fit_4),
      m5 = waic(seed_bayes_fit_5),
      m6 = waic(seed_bayes_fit_6),
      m7 = waic(seed_bayes_fit_7)
    )
  )
  )
  df.waic$models = rownames(df.waic)
  
  out_bayes_compare.waic[[i]] = df.waic
  out_bayes_compare.loo[[i]] = df.loo
}

# calculate delta waic for model comparisons
for(i in 1:length(out_bayes_compare.waic)){
  out_bayes_compare.waic[[i]]$models = rownames(out_bayes_compare.waic[[i]])
  out_bayes_compare.waic[[i]]$del.waic = out_bayes_compare.waic[[i]]$waic - min(out_bayes_compare.waic[[i]]$waic)
}

comparisons = do.call(rbind.data.frame, out_bayes_compare.loo)
meancomps = comparisons %>% group_by(models) %>% summarise_all(mean)

# plot predictive accuracy for different models with their standard error estimates. Values closest to zero
# indicate better leave-out-out accuracy
ggplot(meancomps) +
  geom_point(aes(x = models, y = elpd_diff)) +
  geom_errorbar(aes(x = models, ymin = elpd_diff -se_diff , ymax = elpd_diff +se_diff), width = 0.2)

#calculate standard deviation of measurements across species x treatment combinations
sdcomps = comparisons %>% group_by(models) %>% summarise_all(sd)



### Next, we use the best-fit model (model seven in the above analysis) to calculate ND and FD
### First, we loop over all species x treatment combos
#Create a data frame with all combinations of treatment and species 
spp_list = sort(na.omit( unique(seed_data$focal)))
treat_list = sort( na.omit( unique(seed_data$Tr)))
spp_combos = expand.grid(species = spp_list, treatment = treat_list)

# loop through the species list, fit the model, return the dataframe 
out_bayes = list() #list for for loop output

# fit models as above
for( i in 1:nrow( spp_combos)){ 
  temp_data = seed_data %>%
    filter(focal == spp_combos[i,1], Tr == spp_combos[i,2])
  mod= brm(seed_bayes_fit_formula_7, family=gaussian(), data = temp_data, prior = seedpriors,  chains = 8, cores = 8, iter = 10000, backend = 'cmdstanr')
  out_bayes[[i]] = mod
}

# check MCMC traces and posterior distributions
plot(out_bayes[[1]])
# perform posterior predictive check on model - compares observed y to 100 simulated y's 
brms::pp_check(out_bayes[[1]], ndraws = 100)


out_dfb = list() #list for for loop output
# extract posteror draws for each parameter of interest and reorganise them
for( i in 1:nrow( spp_combos)){ 
  dfb = out_bayes[[i]] %>%
    spread_draws(b_lambda_Intercept, b_aACWR_Intercept, b_aFEMI_Intercept,
                 b_aHOMU_Intercept, b_aPLER_Intercept, b_aSACO_Intercept, b_aURLI_Intercept)
  dfb = dfb[,3:10]
  names(dfb) = c("draw","lambda", "a_ACWR", "a_FEMI", "a_HOMU", "a_PLER", "a_SACO", "a_URLI")
  dfb$focal = rep(spp_combos[i,1], nrow(dfb))
  dfb$treatment = rep(spp_combos[i,2], nrow(dfb))
  
  dfb = dfb %>% pivot_longer(
    cols = starts_with("a"), 
    names_to = c("competitor"), 
    names_prefix = "a_", 
    values_to ="alpha"
  )
  
  out_dfb[[i]] = dfb 
}

posteriors = do.call(rbind.data.frame, out_dfb)

#Make a data frame for the stabilizing niche and fitness differences calculated from each bootstrap's alphas
# take a random sample of 1000 post-warmup posterior draws
draw_list = sort( na.omit( unique(posteriors$draw)))
draw_list = sort(na.omit(sample(draw_list, 1000, replace = F)))
spp_treat_draw_combos = expand.grid(focal = spp_list, treatment = treat_list, competitor = spp_list, id = draw_list)
posteriors = subset(posteriors, draw  %in% draw_list) # subset posteriors to 1000 draws per focal species x treatment combo

#Stabilizing Niche Difference calculations following Van Dyke et al. 2022
stabilizing_niche_diff_bayes = function(df, species1, species2, treat, draw_id) {
  aij = with(df, alpha[focal == species1 & competitor == species2 & treatment== treat & draw == draw_id])
  aji = with(df, alpha[focal == species2 & competitor == species1 & treatment== treat & draw == draw_id])
  ajj = with(df, alpha[focal == species2 & competitor == species2 & treatment== treat &  draw == draw_id])
  aii = with(df, alpha[focal == species1 & competitor == species1 & treatment== treat & draw == draw_id])
  snd = (1 - sqrt((aij * aji)/(ajj * aii)))
  return(snd)
}

stabilizing_niche_diff_bayes(posteriors, "ACWR", "URLI", treat = 1, draw_id = draw_list[[1]])


#Getting eta (ηi) term following Van Dyke et al. 2022
#ηi describes the seeds produced per seed lost from the seed bank for plant species i 
s_g_data = read.csv("./data/s_g_data.csv") #seed survival and germination data
posteriors = merge(posteriors, s_g_data, by = "focal")
posteriors$X = NULL
write_csv(posteriors, "./output/posteriors_model7.csv")

#ηi equation function
get_ni= function(df, species, treat, draw_id){
  lambda = with(df, lambda[ focal == species & treatment == treat & draw == draw_id])[1]
  #print(lambda)
  gi = with( df, g[focal == species & treatment == treat & draw == draw_id])[1]
  #print(gi)
  si = with( df, s[focal == species & treatment == treat & draw == draw_id])[1]
  #print(si)
  
  ni= ((lambda*gi)/(1-((1-gi)*si)))
  return(ni[1])
  
}

get_ni(posteriors, "ACWR", 1, draw_list[1]) #test

#Add snd column to data frame
spp_treat_draw_combos$snd = 0
spp_treat_draw_combos$ni = 0  # Add ni column to data frame


###### This section speeds up the calculation of ND and FD through parallelization
### Depending on the number of draws sampled it can take a significant amount of time
### On an M1 chip w/ 8 cores each loop takes around 5 minutes.
n.cores = parallel::detectCores()
clus = parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = clus)

x = foreach(
  i = 1:n, 
  .combine = 'rbind', .packages = c("magrittr")
) %dopar% {
  sp1 = spp_treat_draw_combos[i, "focal"] %>% unlist
  sp2 = spp_treat_draw_combos[i, "competitor"] %>% unlist
  trt = spp_treat_draw_combos[i , "treatment"] %>% unlist
  draw = spp_treat_draw_combos[i, "id"] %>% unlist
  snd = stabilizing_niche_diff_bayes(posteriors, sp1, sp2, trt, draw)
  ni = get_ni(posteriors, sp1, trt, draw)
  cbind(sp1, sp2, trt, draw, snd, ni)
}

parallel::stopCluster(cl = clus)

spp_treat_draw_combos$snd = x[,"snd"]
spp_treat_draw_combos$ni = x[,"ni"]


#Get fitness differences ------

fitness_diff = function(df1, df2, species1, species2, treat, draw_id) {
  
  ni = with(df1, ni[focal == species1 & treatment == treat & id == draw_id])[1]
  nj = with(df1, ni[focal == species2 & treatment == treat & id == draw_id])[1]
  aij = with(df2, alpha[focal == species1 & competitor == species2 & treatment== treat & draw == draw_id])
  aji = with(df2, alpha[focal == species2 & competitor == species1 & treatment== treat & draw == draw_id])
  ajj = with(df2, alpha[focal == species2 & competitor == species2 & treatment== treat & draw == draw_id])
  aii = with(df2, alpha[focal == species1 & competitor == species1 & treatment== treat & draw == draw_id])
  
  nn= (nj-1)/(ni-1)
  #print(paste("nn: ", nn))
  aa= sqrt((aij * aii)/(ajj * aji))
  #print(paste("aa: ", aa))
  FDij = nn*aa
  #print(FDij[1])
  return(FDij[1])
}

fitness_diff(spp_treat_draw_combos, posteriors, "ACWR", "FEMI", 2, draw_list[1]) #test

#add fitness difference column to data frame
spp_treat_draw_combos$fd = 0

n.cores = parallel::detectCores()
clus = parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = clus)

x = foreach(
  i = 1:n, 
  .combine = 'rbind', .packages = c("magrittr")
) %dopar% {
  sp1 = spp_treat_draw_combos[i, "focal"] %>% unlist
  sp2 = spp_treat_draw_combos[i, "competitor"] %>% unlist
  trt = spp_treat_draw_combos[i , "treatment"] %>% unlist
  draw_id = spp_treat_draw_combos[i, "id"] %>% unlist
  fd = fitness_diff(spp_treat_draw_combos, posteriors, sp1, sp2, trt, draw_id)
  cbind(sp1, sp2, trt, draw_id, fd)
}

parallel::stopCluster(cl = clus)
spp_treat_draw_combos$fd = x[,"fd"]

write_csv(spp_treat_draw_combos, "./output/spp_treat_draw_combos_1000_model7.csv")
#spp_treat_draw_combos = read.csv("./output/spp_treat_draw_combos_1000_model7_new.csv")
spp_treat_draw_combos$sp_pair = paste(spp_treat_draw_combos$focal,spp_treat_draw_combos$competitor,sep='_')

pairs_abc = c("ACWR_FEMI", "ACWR_HOMU", "ACWR_PLER", "ACWR_SACO", "ACWR_URLI", 
              "FEMI_HOMU", "FEMI_PLER", "FEMI_SACO", "FEMI_URLI", "HOMU_PLER",
              "HOMU_SACO", "HOMU_URLI","PLER_SACO", "PLER_URLI", "SACO_URLI",
              "FEMI_ACWR", "HOMU_ACWR" ,"PLER_ACWR","SACO_ACWR", "URLI_ACWR",
              "HOMU_FEMI", "PLER_FEMI", "SACO_FEMI", "URLI_FEMI", "PLER_HOMU",
              "SACO_HOMU", "URLI_HOMU", "SACO_PLER", "URLI_PLER", "URLI_SACO")


spp_treat_draw_combos$fd_superior = ifelse(spp_treat_draw_combos$fd < 1, 1/spp_treat_draw_combos$fd, spp_treat_draw_combos$fd)
spp_treat_draw_combos$fd_sup_sp = ifelse(spp_treat_draw_combos$fd <= 1, 1, 2)
spp_treat_draw_combos_sup = spp_treat_draw_combos %>%
  filter(fd_sup_sp == 2 )
W_superior = with(spp_treat_draw_combos_sup, sp_pair[treatment==1])

draws_pairs_w_sup = spp_treat_draw_combos %>%
  filter(sp_pair %in% W_superior)
draws_pairs_w_sup$treat = factor(draws_pairs_w_sup$treat, levels = c(1, 2))

draws_pairs_unique =  spp_treat_draw_combos %>%
  filter(sp_pair %in% pairs_abc)
draws_pairs_unique$treatment = factor(draws_pairs_unique$treatment, levels = c(1, 2))

draws_pairs_unique$fd_superior = ifelse(draws_pairs_unique$fd < 1, 1/draws_pairs_unique$fd, draws_pairs_unique$fd)

draws_pairs_unique$coexist = ifelse((draws_pairs_unique$snd > (1-1/draws_pairs_unique$fd_superior)), 1, 0 )
draws_pairs_w_sup$treatment = factor(draws_pairs_w_sup$treatment, levels = c(1, 2))
draws_pairs_w_sup$label = paste0(substr(draws_pairs_w_sup$focal, 1, 2), "-", substr(draws_pairs_w_sup$competitor, 1, 2))
write.csv(draws_pairs_w_sup, "./output/plot_data_mod7.csv")
boots_pairs_w_sup = read.csv("./output/boots_pairs_w_sup.csv")
tokeep = boots_pairs_w_sup$label

lot_dat = subset(draws_pairs_w_sup, label %in% tokeep)
medians = plot_dat %>% group_by(label, treatment) %>% summarize(snd = median(snd), fd = median(fd))

# Create coexistence area for plot - min/max fitness difference that permits coexistence
niche_differentiation = seq(from = -.25, to = 1, by = 0.001)
niche_overlap = 1-niche_differentiation
fitness_ratio_min = niche_overlap
fitness_ratio_max = 1/niche_overlap

df = data.frame(niche_diff = niche_differentiation,
                min_fitness_ratio = fitness_ratio_min,
                max_fitness_ratio = fitness_ratio_max)

#Text size 
geom.text.size = 7*(5/14); geom.text.size2 = 6*(5/14)
label.size = 5*(5/14);element.size = 5; theme.size = 7

NDFD_plot_mod7 = ggplot(plot_dat) + 
  theme_classic()+
  theme(text = element_text(size = theme.size))+
  coord_cartesian(xlim = c(-.1, 1), ylim = c(0.1, 100)) +
  geom_line(data = df, aes(x = niche_diff, y = max_fitness_ratio)) +
  geom_line(data = df,  aes(x = niche_diff, y = min_fitness_ratio)) +
  geom_ribbon(data = df, aes(x = niche_diff, ymin = min_fitness_ratio, ymax = max_fitness_ratio), fill = 'grey80') +
  geom_point(aes(x = snd, y = fd, color = treatment, shape = treatment, group = sp_pair), size = 2, alpha = .1) +
  geom_point(data = medians, aes(x = snd, y = fd, fill = treatment, shape = treatment, group = label), size = 3) +
  geom_point(data = boots_pairs_w_sup, aes(x = snd, y = fd, fill = factor(treatment), shape = factor(treatment), group = sp_pair), size = 2, fill = "black") +
  scale_color_manual(values=c('1'="#4E84C4", '2' = "#D16103"), name = "Treatment:", labels = c("Ambient", "Reduced Rain")) +
  scale_fill_manual(values=c('1'="#4E84C4", '2' = "#D16103"), name = "Treatment:", labels = c("Ambient", "Reduced Rain")) +
  scale_shape_manual(values = c('1' = 22, '2' = 21), name = "Treatment:", labels = c("Ambient", "Reduced Rain")) +
  scale_y_log10(expand = c(0,0)) + #setting cut-off and making y on a log scale
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs( x = "Stabilizing niche difference (1-\u03c1)", y = "Fitness difference (Kj/Ki)") +
  theme(axis.title = element_text(size = theme.size), 
        axis.text.x = element_text(size = (theme.size - 2)),
        axis.text.y = element_text(size = (theme.size - 2)),
        legend.title = element_text(size = theme.size),
        legend.text = element_text(size = theme.size),
        legend.position = c(0.88, .035),
        legend.justification = c("bottom"),
        legend.box.just = "bottom",
        legend.margin = margin(1, 1, 1, 1),
        legend.direction = "vertical", 
        strip.background = element_blank(),
        strip.text.x = element_text(size = 8),
        panel.spacing = unit(0,'lines')) +
  guides(fill = guide_legend(title = "Treatment", title.position = "left", cex = 1), col = guide_legend(nrow = 2) , scale = "none")+
  facet_wrap(vars(label), nrow = 4, scales ="free")

plot(NDFD_plot_mod7)

cowplot::save_plot("NDFD_plot_mod7.pdf", new_NDFD_plot, base_height = 8, base_width = 10)

plot_dat$coexist = ifelse((plot_dat$snd > (1-1/plot_dat$fd_superior)), 1, 0 )
plot_dat$Coex1D = plot_dat$snd - (1-1/plot_dat$fd_superior)


coexpct= plot_dat %>% group_by(label, treatment) %>% summarise(coexpct = sum(coexist)/1000)
switchprobs = c()
for(i in 1:length(unique(coexpct$label))){
  pair = unique(coexpct$label)[i]
  tmp = subset(coexpct, label == pair)
  switchprobs[i] = (tmp[which.max(tmp$coexpct),"coexpct"] * (1-tmp[which.min(tmp$coexpct),"coexpct"]))[1,1]
}
coexpct$switchprobs = rep(switchprobs, each = 2)

# plot posterior coexistence predictions with probability of switch
MedianHDI_mod7 = ggplot(plot_dat, aes(x = Coex1D, y = treat, fill = after_stat(quantile))) + 
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE, quantiles = c(0.025, 0.5, 0.975), color = "black", lwd = 0.25) +
  stat_density_ridges(quantile_lines = T, quantiles = 2, fill = NA, color = "black", lwd = 0.25) +
  scale_fill_manual(
    name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#A0A0A0A0", "#0000FFA0")  ) + 
  cowplot::theme_minimal_vgrid() + 
  scale_y_discrete(labels = c("Ambient", "Reduced Rain")) +
  facet_wrap(vars(label), nrow = 4) + scale_x_continuous(breaks = c(-1,0,1), limits = c(-2,2)) +
  xlab("SND - (1-1/FD) (vals > 0 indicate coexistence)") + ylab("") +
  geom_text(data = coexpct, aes(label=round(switchprobs,3)), 
            x = Inf, y = -Inf, hjust=2, vjust=-0.5,
            inherit.aes = FALSE)

plot(MedianHDI_mod7)

save_plot("MedianHDI_mod7", MedianHDI_bestmod, base_height = 10, base_width = 24)

######## REPLICATE FIGURE 3 
seed_data = read.csv("./data/drought_seed_production_data.csv")
seed_data$Tr = ifelse(seed_data$treat == "W", 1, 2)
seed_data$treatment = seed_data$Tr
posteriors = read.csv("./output/posteriors_model7.csv")

spp_list = sort(na.omit( unique(seed_data$focal)))
treat_list = sort( na.omit( unique(seed_data$treatment)))
spp_combos = expand.grid(species = spp_list, treatment = treat_list)

spp_treat_comp_combos = expand.grid(focal = spp_list, competitor = comp_labels, treatment = treat_list)

draw_pairs = spp_treat_comp_combos %>%
  mutate(alpha = 0, alpha_sd = 0, alpha_low = 0, alpha_high = 0, lambda = 0, lambda_low = 0, lambda_high = 0)
draw_pairs$sp_pair = paste(draw_pairs$focal, draw_pairs$competitor, sep = "_")  
for( i in 1:nrow(spp_treat_comp_combos)) {
  sp1 = spp_treat_comp_combos[i, "focal"] %>% unlist
  sp2 = spp_treat_comp_combos[i, "competitor"] %>% unlist
  treatt = spp_treat_comp_combos[i , "treatment"] %>% unlist
  draw_pairs[i, "alpha"] = median(with(posteriors, 
                                       alpha[treatment == treatt & focal == sp1 & competitor == sp2 
                                       ]), na.rm=TRUE)
  draw_pairs[i, "lambda"] = median(with(posteriors, 
                                        lambda[focal == sp1 & competitor == sp2 & 
                                                 treatment == treatt]), na.rm=TRUE)
}  


final_output = read.csv("./output/spp_treat_draw_combos_1000_model7.csv")


mean(with(final_output, ni[focal == "HOMU" & treatment == 1]), na.rm = T)


igr_ratio = function(foc, comp, trt) { #invasion growth rate ratios
  
  nj = mean(with(final_output, ni[focal == comp & treatment == trt]), na.rm = T)
  ni = mean(with(final_output, ni[focal == foc & treatment == trt]), na.rm = T)
  ajj = with(draw_pairs, alpha[focal == comp & competitor == comp & treatment == trt])
  aij = with(draw_pairs, alpha[focal == foc & competitor == comp & treatment == trt])
  n_ratio = (ni-1)/(nj-1)
  a_ratio = ajj/aij
  return(c( n_ratio, a_ratio, n_ratio*a_ratio))
}
igr_ratio("URLI", "HOMU", 1)

comp_labels = sort( na.omit( unique(seed_data$background) ))

spp_treat_comp_combos = expand.grid(species = spp_list, competitor = comp_labels, treatment = treat_list)

spp_treat_comp_combos$n_ratio = 0
spp_treat_comp_combos$a_ratio = 0
spp_treat_comp_combos$product = 0

for(i in 1:nrow(spp_treat_comp_combos)) {
  ii = spp_treat_comp_combos[i, "species"] %>% unlist
  jj = spp_treat_comp_combos[i, "competitor"] %>% unlist
  tt = spp_treat_comp_combos[i, "treatment"] %>% unlist
  
  spp_treat_comp_combos[i, "n_ratio"] = igr_ratio(ii, jj, tt)[1]
  spp_treat_comp_combos[i, "a_ratio"] = igr_ratio(ii, jj, tt)[2]
  spp_treat_comp_combos[i, "product"] = igr_ratio(ii, jj, tt)[3]
  
}
head(spp_treat_comp_combos)
spp_treat_comp_combos$larger = ifelse(spp_treat_comp_combos$a_ratio > spp_treat_comp_combos$n_ratio, "a", "n")
spp_treat_comp_combos$invade = ifelse(spp_treat_comp_combos$product>1, "yes", "no")
head(spp_treat_comp_combos)

invaders = spp_treat_comp_combos %>%
  filter(invade == "yes")

#Which ratio changes more in invasion growth rate inequality?
igr_change = function(foc, comp) {
  
  nj_d = mean(with(final_output, ni[focal == comp & treatment == 2]), na.rm = T)
  ni_d = mean(with(final_output, ni[focal == foc & treatment == 2]), na.rm = T)
  ajj_d = with(draw_pairs, alpha[focal == comp & competitor == comp & treatment == 2])
  aij_d = with(draw_pairs, alpha[focal == foc & competitor == comp & treatment == 2])
  n_ratio_d = log10((ni_d-1)/(nj_d-1))
  a_ratio_d = log10(ajj_d/aij_d)
  
  nj_w = mean(with(final_output, ni[focal == comp & treatment == 1]), na.rm = T)
  ni_w = mean(with(final_output, ni[focal == foc & treatment == 1]), na.rm = T)
  ajj_w = with(draw_pairs, alpha[focal == comp & competitor == comp & treatment == 1])
  aij_w = with(draw_pairs, alpha[focal == foc & competitor == comp & treatment == 1])
  n_ratio_w = log10((ni_w-1)/(nj_w-1))
  a_ratio_w = log10(ajj_w/aij_w)
  
  nc=abs(n_ratio_w - n_ratio_d)
  ac=abs(a_ratio_w - a_ratio_d)
  return(c(nc, ac))
}
igr_change("HOMU", "URLI")

#Create data frame
spp_list = sort(na.omit( unique(seed_data$focal)))
comp_labels = sort( na.omit( unique(seed_data$background) ))
spp_comp_combos = expand.grid(species = spp_list, competitor = comp_labels)
spp_comp_combos$n_change = 0
spp_comp_combos$a_change = 0


for(i in 1:nrow(spp_comp_combos)) {
  ii = spp_comp_combos[i, "species"] %>% unlist
  jj = spp_comp_combos[i, "competitor"] %>% unlist
  spp_comp_combos[i, "n_change"] = igr_change(ii, jj)[1]
  spp_comp_combos[i, "a_change"] = igr_change(ii, jj)[2]
}

spp_comp_combos$larger = ifelse(abs(spp_comp_combos$a_change)> abs(spp_comp_combos$n_change), "a", 'n')

spp_comp_combos = spp_comp_combos %>%
  filter(a_change != n_change)
total_pairs = c("ACWR_FEMI", "ACWR_HOMU", "ACWR_PLER", "SACO_ACWR", "URLI_ACWR", 
                "HOMU_FEMI", "PLER_FEMI", "SACO_FEMI", "URLI_FEMI", "PLER_HOMU",
                "SACO_HOMU", "URLI_HOMU","SACO_PLER", "URLI_PLER", "URLI_SACO")
spp_comp_combos$sp_pair = paste(spp_comp_combos$species, spp_comp_combos$competitor, sep = "_")

spp_comp_combos =spp_comp_combos %>%
  filter(sp_pair %in% total_pairs)

t.test(spp_comp_combos$n_change, spp_comp_combos$a_change, paired = T)

pairs_coexist_change =c("ACWR_FEMI", "ACWR_HOMU", "ACWR_PLER", "SACO_ACWR", 
                        "URLI_ACWR",  "HOMU_FEMI", "PLER_FEMI", "URLI_FEMI",
                        "SACO_PLER", "URLI_SACO")

spp_comp_combos_coexist =spp_comp_combos %>%
  filter(sp_pair %in% pairs_coexist_change)

t.test(spp_comp_combos_coexist$n_change, spp_comp_combos_coexist$a_change, paired = T)


spp_comp_combos_coexist_pivot = pivot_longer(spp_comp_combos_coexist, cols = c(n_change, a_change), names_to = "type_change", values_to ="change_value")

t.test(spp_comp_combos$n_change, spp_comp_combos$a_change, paired = T)
spp_comp_combos_pivot = pivot_longer(spp_comp_combos, cols = c(n_change, a_change), names_to = "type_change", values_to ="change_value")


ggplot(spp_comp_combos_pivot, aes(x = type_change, y =change_value))+
  theme_classic(base_size = 20) +
  theme(text = element_text(size = 18))+
  geom_boxplot(fill = "light grey", outlier.shape = NA) +
  geom_point(shape = 1, size = 4) +
  ylab("Difference between treatments") +
  xlab (" ") +
  #scale_y_log10() +
  scale_x_discrete(labels=c("a_change"="Competition \n coefficients", "n_change" = "Demographic \n potential")) 

ggsave("./figures/alpha_eta_ratio_box_model7.pdf")
write.csv(spp_comp_combos_pivot,  "./output/n_alph_ratio_output_mod7.csv")


#######



gkt = read.csv("./data/traits_gk.csv")
rownames(gkt) = gkt$species

draw_pairs = read.csv("./output/spp_treat_draw_combos_1000_model7.csv") # or run model script first to get this
final_output_draws = draw_pairs
final_output_draws$fd_superior = ifelse(final_output_draws$fd < 1, 1/final_output_draws$fd, final_output_draws$fd)
gkt_pca = prcomp(gkt[,2:12], scale = T)

gkt_pca$x
gkt_pca$x["ACWR",1]
screeplot(gkt_pca)

#Create PC vectors with the six species in the experiment
pc1= c(gkt_pca$x["ACWR",1], gkt_pca$x["FEMI",1], gkt_pca$x["HOMU",1], gkt_pca$x["PLER",1], gkt_pca$x["SACO",1], gkt_pca$x["URLI",1])
pc2= c(gkt_pca$x["ACWR",2], gkt_pca$x["FEMI",2], gkt_pca$x["HOMU",2], gkt_pca$x["PLER",2], gkt_pca$x["SACO",2], gkt_pca$x["URLI",2])
#Dist matrices for fd and snd
fd_d  = matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))

# The following loop propagates error in FD and ND measurements into the distance matrices used to 
# assess differences in ND and FD among species
dist_fd_d = list()
for( i in 1:100){
  fd_d  = matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))
  for(sp1 in rownames(fd_d)) {  
    for(sp2 in colnames(fd_d)) {    
      current_fd = final_output_draws %>% 
        filter(id == unique(final_output_draws$id)[i] & treatment == 2 & focal == sp1 & competitor == sp2) %>% 
        dplyr::select(fd_superior) %>%
        unlist()
      fd_d[sp1, sp2] = current_fd  
      dist_fd_d[[i]]= as.dist(fd_d)
    }
  }
  print(i)
}

dist_pc1= dist(pc1)

#run a mantel test over all posterior samples to generate distribution of 
# r and p statistics.
mantelr = c()
mantelp = c()
for(i in 1:length(dist_fd_d)){
  test= vegan :: mantel(xdis = dist_pc1, ydis = dist_fd_d[[i]])
  mantelr[i] = test$statistic
  mantelp[i] = test$signif
  print(i)
}

#plot histograms and 95% hd intervals for each metric
hist(mantelr^2); median(mantelr^2); HDInterval::hdi(mantelr^2)
hist(mantelp); median(mantelp); HDInterval::hdi(mantelp)

#check correspondencee between mantel test stats from a single median estimate and the 
# posterior distribution of r calculated above
medians = final_output_draws %>% group_by(focal, competitor, treatment) %>% summarise(fd_median = median(fd_superior))
fd_d  = matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))
for(sp1 in rownames(fd_d)) {  
  for(sp2 in colnames(fd_d)) {    
    current_fd = medians %>% 
      filter(treatment == 2 & focal == sp1 & competitor == sp2) %>% 
      dplyr::select(fd_median)
    fd_d[sp1, sp2] = current_fd$fd_median  }}
dist_fd_d= as.dist(fd_d)

vegan :: mantel(xdis = dist_pc1, ydis = dist_fd_d)
median(mantelr); median(mantelp) #these values should be somewhat close to those from the previous line


#Do traits dissimilarities between pairs correlate with magnitude that their snd and fd change between treatments
diffs = final_output_draws

fd_diff_abs = abs(diffs$fd[diffs$treatment == 1] - diffs$fd[diffs$treatment == 2])

diffs = diffs %>% 
  dplyr :: select(id, focal, competitor, treatment) %>%
  filter(treatment == 1) 

diffs$fd_diff_abs = fd_diff_abs

dist_fd_diffs_abs = list()
for( i in 1:100){
  fd_diffs_abs  = matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))
  for(sp1 in rownames(fd_diffs_abs)) {  
    for(sp2 in colnames(fd_diffs_abs)) {    
      fd_change = diffs %>% 
        filter(id == unique(diffs$id)[i] & focal == sp1 & competitor == sp2) %>% 
        dplyr::select(fd_diff_abs) %>%
        unlist()
      fd_diffs_abs[sp1, sp2] = fd_change  }}
  dist_fd_diffs_abs[[i]]= dist(fd_diffs_abs)
  print(i)
}

for(i in 1:length(dist_fd_diffs_abs)){
  test = vegan :: mantel(ydis = dist_fd_diffs_abs[[i]],xdis = dist_pc1)
  mantelr[i] = test$statistic
  mantelp[i] = test$signif
  print(i)
}

hist(mantelr^2); median(mantelr^2); HDInterval::hdi(mantelr^2)
hist(mantelp); median(mantelp); HDInterval::hdi(mantelp)


medians = final_output_draws %>% group_by(focal, competitor, treatment) %>% summarise(fd_median = median(fd))
diffs = as.data.frame(medians)
fd_diff_abs = abs(diffs$fd_median[diffs$treatment == 1] - diffs$fd_median[diffs$treatment == 2])

diffs = diffs %>% 
  dplyr :: select(focal, competitor, treatment) %>%
  filter(treatment == 1) 

diffs$fd_diff_abs = fd_diff_abs


fd_diffs_abs_all  = matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))
for(sp1 in rownames(fd_diffs_abs_all)) {  
  for(sp2 in colnames(fd_diffs_abs_all)) {    
    fd_change = diffs %>% 
      filter(focal == sp1 & competitor == sp2)  %>%
      dplyr::select(fd_diff_abs) %>%
      unlist()
    fd_diffs_abs_all[sp1, sp2] = fd_change  }}
dist_fd_diffs_abs_all= dist(fd_diffs_abs_all)

vegan ::mantel(xdis = dist_pc1,ydis = dist_fd_diffs_abs_all) 
median(mantelr); median(mantelp) #these values should be somewhat close to those from the previous line

pc1_vec = as.vector(dist_pc1)
fd_diffs_abs_vec = as.vector(dist_fd_diffs_abs_all)
fd_diffs_vec = as.vector(dist_fd_d)

fd_d_df = as.data.frame(gdata :: unmatrix(fd_d)) 
colnames(fd_d_df) = "fd_d"
fd_d_df = tibble::rownames_to_column(fd_d_df, "sp_pair")

fd_diffs_abs_df = as.data.frame(gdata :: unmatrix(fd_diffs_abs))
colnames(fd_diffs_abs_df) = "fd_diff_abs"
fd_diffs_abs_df = tibble::rownames_to_column(fd_diffs_abs_df, "sp_pair")

pc1_mat = as.matrix(dist(pc1, upper = T))
rownames(pc1_mat)= c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")
colnames(pc1_mat)= c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")
pc1_df = as.data.frame(gdata::unmatrix(pc1_mat))
colnames(pc1_df) = "pc1"
pc1_df= tibble::rownames_to_column(pc1_df, "sp_pair")


df_compile = left_join(fd_d_df,fd_diffs_abs_df, by = "sp_pair") %>%
  left_join(., pc1_df, by = "sp_pair")

w_pairs= c("SACO:ACWR", "URLI:ACWR", "ACWR:FEMI", "HOMU:FEMI", 
           "PLER:FEMI", "SACO:FEMI", "URLI:FEMI","ACWR:HOMU", 
           "PLER:HOMU", "SACO:HOMU", "URLI:HOMU", "ACWR:PLER", 
           "SACO:PLER", "URLI:PLER","URLI:SACO")

df_compile= df_compile %>%
  filter(sp_pair %in% w_pairs)
df_compile$label = paste0(substr(df_compile$sp_pair, 1, 2), "-", substr(df_compile$sp_pair, 6, 7))

write_csv(df_compile, "./output/pca_mantel_mod7.csv")

vegan :: mantel(dist_pc1, dist_fd_d)
ggplot(df_compile, aes(x = pc1, y = fd_d)) + geom_point() +geom_smooth(method = "lm")

vegan :: mantel(dist_pc1, dist_fd_diffs_abs_all)
ggplot(df_compile, aes(x = pc1, y = fd_diff_abs)) + geom_point() +geom_smooth(method = "lm")























#######  Replicate the previous analysis with the author's BH model (model 5)

options(scipen = 5)
set.seed(1234)
setwd("Downloads/Sedgwick_public/")
seed_data = read.csv("./data/drought_seed_production_data.csv")
seed_data$Tr = ifelse(seed_data$treat == "W", 1, 2)
seed_data = seed_data %>% replace_na(list(background = "ACWR")) # putting acwr instead of NA for the lambdas
seed_data$N_acwr = ifelse(seed_data$background == "ACWR", seed_data$num_comp, 0)
seed_data$N_femi = ifelse(seed_data$background == "FEMI", seed_data$num_comp, 0)
seed_data$N_homu = ifelse(seed_data$background == "HOMU", seed_data$num_comp, 0)
seed_data$N_pler = ifelse(seed_data$background == "PLER", seed_data$num_comp, 0)
seed_data$N_saco = ifelse(seed_data$background == "SACO", seed_data$num_comp, 0)
seed_data$N_urli = ifelse(seed_data$background == "URLI", seed_data$num_comp, 0)


### Next, we use the best-fit model (model seven in the above analysis) to calculate ND and FD
### First, we loop over all species x treatment combos
#Create a data frame with all combinations of treatment and species 
spp_list = sort(na.omit( unique(seed_data$focal)))
treat_list = sort( na.omit( unique(seed_data$Tr)))
spp_combos = expand.grid(species = spp_list, treatment = treat_list)

# loop through the species list, fit the model, return the dataframe 
out_bayes = list() #list for for loop output

# fit models as above
for( i in 1:nrow( spp_combos)){ 
  temp_data = seed_data %>%
    filter(focal == spp_combos[i,1], Tr == spp_combos[i,2])
  mod= brm(seed_bayes_fit_formula_4, family=gaussian(), data = temp_data, prior = seedpriors[-8,],  chains = 8, cores = 8, iter = 10000, backend = 'cmdstanr')
  out_bayes[[i]] = mod
}

# check MCMC traces and posterior distributions
plot(out_bayes[[1]])
# perform posterior predictive check on model - compares observed y to 100 simulated y's 
brms::pp_check(out_bayes[[12]], ndraws = 100)


out_dfb = list() #list for for loop output
# extract posteror draws for each parameter of interest and reorganise them
for( i in 1:nrow( spp_combos)){ 
  dfb = out_bayes[[i]] %>%
    spread_draws(b_lambda_Intercept, b_aACWR_Intercept, b_aFEMI_Intercept,
                 b_aHOMU_Intercept, b_aPLER_Intercept, b_aSACO_Intercept, b_aURLI_Intercept)
  dfb = dfb[,3:10]
  names(dfb) = c("draw","lambda", "a_ACWR", "a_FEMI", "a_HOMU", "a_PLER", "a_SACO", "a_URLI")
  dfb$focal = rep(spp_combos[i,1], nrow(dfb))
  dfb$treatment = rep(spp_combos[i,2], nrow(dfb))
  
  dfb = dfb %>% pivot_longer(
    cols = starts_with("a"), 
    names_to = c("competitor"), 
    names_prefix = "a_", 
    values_to ="alpha"
  )
  
  out_dfb[[i]] = dfb 
}

posteriors = do.call(rbind.data.frame, out_dfb)

#Make a data frame for the stabilizing niche and fitness differences calculated from each bootstrap's alphas
# take a random sample of 1000 post-warmup posterior draws
draw_list = sort( na.omit( unique(posteriors$draw)))
draw_list = sort(na.omit(sample(draw_list, 1000, replace = F)))
spp_treat_draw_combos = expand.grid(focal = spp_list, treatment = treat_list, competitor = spp_list, id = draw_list)
posteriors = subset(posteriors, draw  %in% draw_list) # subset posteriors to 1000 draws per focal species x treatment combo

#Stabilizing Niche Difference calculations following Van Dyke et al. 2022
stabilizing_niche_diff_bayes = function(df, species1, species2, treat, draw_id) {
  aij = with(df, alpha[focal == species1 & competitor == species2 & treatment== treat & draw == draw_id])
  aji = with(df, alpha[focal == species2 & competitor == species1 & treatment== treat & draw == draw_id])
  ajj = with(df, alpha[focal == species2 & competitor == species2 & treatment== treat &  draw == draw_id])
  aii = with(df, alpha[focal == species1 & competitor == species1 & treatment== treat & draw == draw_id])
  snd = (1 - sqrt((aij * aji)/(ajj * aii)))
  return(snd)
}

stabilizing_niche_diff_bayes(posteriors, "ACWR", "URLI", treat = 1, draw_id = draw_list[[1]])


#Getting eta (ηi) term following Van Dyke et al. 2022
#ηi describes the seeds produced per seed lost from the seed bank for plant species i 
s_g_data = read.csv("./data/s_g_data.csv") #seed survival and germination data
posteriors = merge(posteriors, s_g_data, by = "focal")
posteriors$X = NULL
write_csv(posteriors, "./output/posteriors_model5.csv")

#ηi equation function
get_ni= function(df, species, treat, draw_id){
  lambda = with(df, lambda[ focal == species & treatment == treat & draw == draw_id])[1]
  #print(lambda)
  gi = with( df, g[focal == species & treatment == treat & draw == draw_id])[1]
  #print(gi)
  si = with( df, s[focal == species & treatment == treat & draw == draw_id])[1]
  #print(si)
  
  ni= ((lambda*gi)/(1-((1-gi)*si)))
  return(ni[1])
  
}

get_ni(posteriors, "ACWR", 1, draw_list[1]) #test

#Add snd column to data frame
spp_treat_draw_combos$snd = 0
spp_treat_draw_combos$ni = 0  # Add ni column to data frame


###### This section speeds up the calculation of ND and FD through parallelization
### Depending on the number of draws sampled it can take a significant amount of time
### On an M1 chip w/ 8 cores each loop takes around 5 minutes.
n.cores = parallel::detectCores()
clus = parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = clus)

x = foreach(
  i = 1:n, 
  .combine = 'rbind', .packages = c("magrittr")
) %dopar% {
  sp1 = spp_treat_draw_combos[i, "focal"] %>% unlist
  sp2 = spp_treat_draw_combos[i, "competitor"] %>% unlist
  trt = spp_treat_draw_combos[i , "treatment"] %>% unlist
  draw = spp_treat_draw_combos[i, "id"] %>% unlist
  snd = stabilizing_niche_diff_bayes(posteriors, sp1, sp2, trt, draw)
  ni = get_ni(posteriors, sp1, trt, draw)
  cbind(sp1, sp2, trt, draw, snd, ni)
}

parallel::stopCluster(cl = clus)

spp_treat_draw_combos$snd = x[,"snd"]
spp_treat_draw_combos$ni = x[,"ni"]


#Get fitness differences ------

fitness_diff = function(df1, df2, species1, species2, treat, draw_id) {
  
  ni = with(df1, ni[focal == species1 & treatment == treat & id == draw_id])[1]
  nj = with(df1, ni[focal == species2 & treatment == treat & id == draw_id])[1]
  aij = with(df2, alpha[focal == species1 & competitor == species2 & treatment== treat & draw == draw_id])
  aji = with(df2, alpha[focal == species2 & competitor == species1 & treatment== treat & draw == draw_id])
  ajj = with(df2, alpha[focal == species2 & competitor == species2 & treatment== treat & draw == draw_id])
  aii = with(df2, alpha[focal == species1 & competitor == species1 & treatment== treat & draw == draw_id])
  
  nn= (nj-1)/(ni-1)
  #print(paste("nn: ", nn))
  aa= sqrt((aij * aii)/(ajj * aji))
  #print(paste("aa: ", aa))
  FDij = nn*aa
  #print(FDij[1])
  return(FDij[1])
}

fitness_diff(spp_treat_draw_combos, posteriors, "ACWR", "FEMI", 2, draw_list[1]) #test

#add fitness difference column to data frame
spp_treat_draw_combos$fd = 0

n.cores = parallel::detectCores()
clus = parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = clus)

x = foreach(
  i = 1:n, 
  .combine = 'rbind', .packages = c("magrittr")
) %dopar% {
  sp1 = spp_treat_draw_combos[i, "focal"] %>% unlist
  sp2 = spp_treat_draw_combos[i, "competitor"] %>% unlist
  trt = spp_treat_draw_combos[i , "treatment"] %>% unlist
  draw_id = spp_treat_draw_combos[i, "id"] %>% unlist
  fd = fitness_diff(spp_treat_draw_combos, posteriors, sp1, sp2, trt, draw_id)
  cbind(sp1, sp2, trt, draw_id, fd)
}

parallel::stopCluster(cl = clus)
spp_treat_draw_combos$fd = x[,"fd"]

write_csv(spp_treat_draw_combos, "./output/spp_treat_draw_combos_1000_model5.csv")
#spp_treat_draw_combos = read.csv("./output/spp_treat_draw_combos_1000_model7_new.csv")
spp_treat_draw_combos$sp_pair = paste(spp_treat_draw_combos$focal,spp_treat_draw_combos$competitor,sep='_')

pairs_abc = c("ACWR_FEMI", "ACWR_HOMU", "ACWR_PLER", "ACWR_SACO", "ACWR_URLI", 
              "FEMI_HOMU", "FEMI_PLER", "FEMI_SACO", "FEMI_URLI", "HOMU_PLER",
              "HOMU_SACO", "HOMU_URLI","PLER_SACO", "PLER_URLI", "SACO_URLI",
              "FEMI_ACWR", "HOMU_ACWR" ,"PLER_ACWR","SACO_ACWR", "URLI_ACWR",
              "HOMU_FEMI", "PLER_FEMI", "SACO_FEMI", "URLI_FEMI", "PLER_HOMU",
              "SACO_HOMU", "URLI_HOMU", "SACO_PLER", "URLI_PLER", "URLI_SACO")


spp_treat_draw_combos$fd_superior = ifelse(spp_treat_draw_combos$fd < 1, 1/spp_treat_draw_combos$fd, spp_treat_draw_combos$fd)
spp_treat_draw_combos$fd_sup_sp = ifelse(spp_treat_draw_combos$fd <= 1, 1, 2)
spp_treat_draw_combos_sup = spp_treat_draw_combos %>%
  filter(fd_sup_sp == 2 )
W_superior = with(spp_treat_draw_combos_sup, sp_pair[treatment==1])

draws_pairs_w_sup = spp_treat_draw_combos %>%
  filter(sp_pair %in% W_superior)
draws_pairs_w_sup$treat = factor(draws_pairs_w_sup$treat, levels = c(1, 2))

draws_pairs_unique =  spp_treat_draw_combos %>%
  filter(sp_pair %in% pairs_abc)
draws_pairs_unique$treatment = factor(draws_pairs_unique$treatment, levels = c(1, 2))

draws_pairs_unique$fd_superior = ifelse(draws_pairs_unique$fd < 1, 1/draws_pairs_unique$fd, draws_pairs_unique$fd)

draws_pairs_unique$coexist = ifelse((draws_pairs_unique$snd > (1-1/draws_pairs_unique$fd_superior)), 1, 0 )
draws_pairs_w_sup$treatment = factor(draws_pairs_w_sup$treatment, levels = c(1, 2))
draws_pairs_w_sup$label = paste0(substr(draws_pairs_w_sup$focal, 1, 2), "-", substr(draws_pairs_w_sup$competitor, 1, 2))
write.csv(draws_pairs_w_sup, "./output/plot_data_mod5.csv")
boots_pairs_w_sup = read.csv("./output/boots_pairs_w_sup.csv")
tokeep = boots_pairs_w_sup$label

plot_dat = subset(draws_pairs_w_sup, label %in% tokeep)
medians = plot_dat %>% group_by(label, treatment) %>% summarize(snd = median(snd), fd = median(fd))

# Create coexistence area for plot - min/max fitness difference that permits coexistence
niche_differentiation = seq(from = -.25, to = 1, by = 0.001)
niche_overlap = 1-niche_differentiation
fitness_ratio_min = niche_overlap
fitness_ratio_max = 1/niche_overlap

df = data.frame(niche_diff = niche_differentiation,
                min_fitness_ratio = fitness_ratio_min,
                max_fitness_ratio = fitness_ratio_max)

#Text size 
geom.text.size = 7*(5/14); geom.text.size2 = 6*(5/14)
label.size = 5*(5/14);element.size = 5; theme.size = 7

NDFD_plot_mod5 = ggplot(plot_dat) + 
  theme_classic()+
  theme(text = element_text(size = theme.size))+
  coord_cartesian(xlim = c(-.1, 1), ylim = c(0.1, 100)) +
  geom_line(data = df, aes(x = niche_diff, y = max_fitness_ratio)) +
  geom_line(data = df,  aes(x = niche_diff, y = min_fitness_ratio)) +
  geom_ribbon(data = df, aes(x = niche_diff, ymin = min_fitness_ratio, ymax = max_fitness_ratio), fill = 'grey80') +
  geom_point(aes(x = snd, y = fd, color = treatment, shape = treatment, group = sp_pair), size = 2, alpha = .1) +
  geom_point(data = medians, aes(x = snd, y = fd, fill = treatment, shape = treatment, group = label), size = 3) +
  geom_point(data = boots_pairs_w_sup, aes(x = snd, y = fd, fill = factor(treatment), shape = factor(treatment), group = sp_pair), size = 2, fill = "black") +
  scale_color_manual(values=c('1'="#4E84C4", '2' = "#D16103"), name = "Treatment:", labels = c("Ambient", "Reduced Rain")) +
  scale_fill_manual(values=c('1'="#4E84C4", '2' = "#D16103"), name = "Treatment:", labels = c("Ambient", "Reduced Rain")) +
  scale_shape_manual(values = c('1' = 22, '2' = 21), name = "Treatment:", labels = c("Ambient", "Reduced Rain")) +
  scale_y_log10(expand = c(0,0)) + #setting cut-off and making y on a log scale
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs( x = "Stabilizing niche difference (1-\u03c1)", y = "Fitness difference (Kj/Ki)") +
  theme(axis.title = element_text(size = theme.size), 
        axis.text.x = element_text(size = (theme.size - 2)),
        axis.text.y = element_text(size = (theme.size - 2)),
        legend.title = element_text(size = theme.size),
        legend.text = element_text(size = theme.size),
        legend.position = c(0.88, .035),
        legend.justification = c("bottom"),
        legend.box.just = "bottom",
        legend.margin = margin(1, 1, 1, 1),
        legend.direction = "vertical", 
        strip.background = element_blank(),
        strip.text.x = element_text(size = 8),
        panel.spacing = unit(0,'lines')) +
  guides(fill = guide_legend(title = "Treatment", title.position = "left", cex = 1), col = guide_legend(nrow = 2) , scale = "none")+
  facet_wrap(vars(label), nrow = 4, scales ="free")

plot(NDFD_plot_mod5)

cowplot::save_plot("NDFD_plot_mod5.pdf", NDFD_plot_mod5, base_height = 8, base_width = 10)

plot_dat$coexist = ifelse((plot_dat$snd > (1-1/plot_dat$fd_superior)), 1, 0 )
plot_dat$Coex1D = plot_dat$snd - (1-1/plot_dat$fd_superior)


coexpct= plot_dat %>% group_by(label, treatment) %>% summarise(coexpct = sum(coexist)/1000)
switchprobs = c()
for(i in 1:length(unique(coexpct$label))){
  pair = unique(coexpct$label)[i]
  tmp = subset(coexpct, label == pair)
  switchprobs[i] = (tmp[which.max(tmp$coexpct),"coexpct"] * (1-tmp[which.min(tmp$coexpct),"coexpct"]))[1,1]
}
coexpct$switchprobs = rep(switchprobs, each = 2)

# plot posterior coexistence predictions with probability of switch
MedianHDI_mod5 = ggplot(plot_dat, aes(x = Coex1D, y = treat, fill = after_stat(quantile))) + 
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE, quantiles = c(0.025, 0.5, 0.975), color = "black", lwd = 0.25) +
  stat_density_ridges(quantile_lines = T, quantiles = 2, fill = NA, color = "black", lwd = 0.25) +
  scale_fill_manual(
    name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#A0A0A0A0", "#0000FFA0")  ) + 
  cowplot::theme_minimal_vgrid() + 
  scale_y_discrete(labels = c("Ambient", "Reduced Rain")) +
  facet_wrap(vars(label), nrow = 4) + scale_x_continuous(breaks = c(-1,0,1), limits = c(-2,2)) +
  xlab("SND - (1-1/FD) (vals > 0 indicate coexistence)") + ylab("") +
  geom_text(data = coexpct, aes(label=round(switchprobs,3)), 
            x = Inf, y = -Inf, hjust=2, vjust=-0.5,
            inherit.aes = FALSE)

plot(MedianHDI_mod5)

cowplot::save_plot("MedianHDI_mod5.pdf", MedianHDI_mod5, base_height = 10, base_width = 24)

######## REPLICATE FIGURE 3 
seed_data = read.csv("./data/drought_seed_production_data.csv")
seed_data$Tr = ifelse(seed_data$treat == "W", 1, 2)
seed_data$treatment = seed_data$Tr
posteriors = read.csv("./output/posteriors_model5.csv")

spp_list = sort(na.omit( unique(seed_data$focal)))
treat_list = sort( na.omit( unique(seed_data$treatment)))
spp_combos = expand.grid(species = spp_list, treatment = treat_list)

spp_treat_comp_combos = expand.grid(focal = spp_list, competitor = comp_labels, treatment = treat_list)

draw_pairs = spp_treat_comp_combos %>%
  mutate(alpha = 0, alpha_sd = 0, alpha_low = 0, alpha_high = 0, lambda = 0, lambda_low = 0, lambda_high = 0)
draw_pairs$sp_pair = paste(draw_pairs$focal, draw_pairs$competitor, sep = "_")  
for( i in 1:nrow(spp_treat_comp_combos)) {
  sp1 = spp_treat_comp_combos[i, "focal"] %>% unlist
  sp2 = spp_treat_comp_combos[i, "competitor"] %>% unlist
  treatt = spp_treat_comp_combos[i , "treatment"] %>% unlist
  draw_pairs[i, "alpha"] = median(with(posteriors, 
                                       alpha[treatment == treatt & focal == sp1 & competitor == sp2 
                                       ]), na.rm=TRUE)
  draw_pairs[i, "lambda"] = median(with(posteriors, 
                                        lambda[focal == sp1 & competitor == sp2 & 
                                                 treatment == treatt]), na.rm=TRUE)
}  


final_output = read.csv("./output/spp_treat_draw_combos_1000_model5.csv")


mean(with(final_output, ni[focal == "HOMU" & treatment == 1]), na.rm = T)


igr_ratio = function(foc, comp, trt) { #invasion growth rate ratios
  
  nj = mean(with(final_output, ni[focal == comp & treatment == trt]), na.rm = T)
  ni = mean(with(final_output, ni[focal == foc & treatment == trt]), na.rm = T)
  ajj = with(draw_pairs, alpha[focal == comp & competitor == comp & treatment == trt])
  aij = with(draw_pairs, alpha[focal == foc & competitor == comp & treatment == trt])
  n_ratio = (ni-1)/(nj-1)
  a_ratio = ajj/aij
  return(c( n_ratio, a_ratio, n_ratio*a_ratio))
}
igr_ratio("URLI", "HOMU", 1)

comp_labels = sort( na.omit( unique(seed_data$background) ))

spp_treat_comp_combos = expand.grid(species = spp_list, competitor = comp_labels, treatment = treat_list)

spp_treat_comp_combos$n_ratio = 0
spp_treat_comp_combos$a_ratio = 0
spp_treat_comp_combos$product = 0

for(i in 1:nrow(spp_treat_comp_combos)) {
  ii = spp_treat_comp_combos[i, "species"] %>% unlist
  jj = spp_treat_comp_combos[i, "competitor"] %>% unlist
  tt = spp_treat_comp_combos[i, "treatment"] %>% unlist
  
  spp_treat_comp_combos[i, "n_ratio"] = igr_ratio(ii, jj, tt)[1]
  spp_treat_comp_combos[i, "a_ratio"] = igr_ratio(ii, jj, tt)[2]
  spp_treat_comp_combos[i, "product"] = igr_ratio(ii, jj, tt)[3]
  
}
head(spp_treat_comp_combos)
spp_treat_comp_combos$larger = ifelse(spp_treat_comp_combos$a_ratio > spp_treat_comp_combos$n_ratio, "a", "n")
spp_treat_comp_combos$invade = ifelse(spp_treat_comp_combos$product>1, "yes", "no")
head(spp_treat_comp_combos)

invaders = spp_treat_comp_combos %>%
  filter(invade == "yes")

#Which ratio changes more in invasion growth rate inequality?
igr_change = function(foc, comp) {
  
  nj_d = mean(with(final_output, ni[focal == comp & treatment == 2]), na.rm = T)
  ni_d = mean(with(final_output, ni[focal == foc & treatment == 2]), na.rm = T)
  ajj_d = with(draw_pairs, alpha[focal == comp & competitor == comp & treatment == 2])
  aij_d = with(draw_pairs, alpha[focal == foc & competitor == comp & treatment == 2])
  n_ratio_d = log10((ni_d-1)/(nj_d-1))
  a_ratio_d = log10(ajj_d/aij_d)
  
  nj_w = mean(with(final_output, ni[focal == comp & treatment == 1]), na.rm = T)
  ni_w = mean(with(final_output, ni[focal == foc & treatment == 1]), na.rm = T)
  ajj_w = with(draw_pairs, alpha[focal == comp & competitor == comp & treatment == 1])
  aij_w = with(draw_pairs, alpha[focal == foc & competitor == comp & treatment == 1])
  n_ratio_w = log10((ni_w-1)/(nj_w-1))
  a_ratio_w = log10(ajj_w/aij_w)
  
  nc=abs(n_ratio_w - n_ratio_d)
  ac=abs(a_ratio_w - a_ratio_d)
  return(c(nc, ac))
}
igr_change("HOMU", "URLI")

#Create data frame
spp_list = sort(na.omit( unique(seed_data$focal)))
comp_labels = sort( na.omit( unique(seed_data$background) ))
spp_comp_combos = expand.grid(species = spp_list, competitor = comp_labels)
spp_comp_combos$n_change = 0
spp_comp_combos$a_change = 0


for(i in 1:nrow(spp_comp_combos)) {
  ii = spp_comp_combos[i, "species"] %>% unlist
  jj = spp_comp_combos[i, "competitor"] %>% unlist
  spp_comp_combos[i, "n_change"] = igr_change(ii, jj)[1]
  spp_comp_combos[i, "a_change"] = igr_change(ii, jj)[2]
}

spp_comp_combos$larger = ifelse(abs(spp_comp_combos$a_change)> abs(spp_comp_combos$n_change), "a", 'n')

spp_comp_combos = spp_comp_combos %>%
  filter(a_change != n_change)
total_pairs = c("ACWR_FEMI", "ACWR_HOMU", "ACWR_PLER", "SACO_ACWR", "URLI_ACWR", 
                "HOMU_FEMI", "PLER_FEMI", "SACO_FEMI", "URLI_FEMI", "PLER_HOMU",
                "SACO_HOMU", "URLI_HOMU","SACO_PLER", "URLI_PLER", "URLI_SACO")
spp_comp_combos$sp_pair = paste(spp_comp_combos$species, spp_comp_combos$competitor, sep = "_")

spp_comp_combos =spp_comp_combos %>%
  filter(sp_pair %in% total_pairs)

t.test(spp_comp_combos$n_change, spp_comp_combos$a_change, paired = T)

pairs_coexist_change =c("ACWR_FEMI", "ACWR_HOMU", "ACWR_PLER", "SACO_ACWR", 
                        "URLI_ACWR",  "HOMU_FEMI", "PLER_FEMI", "URLI_FEMI",
                        "SACO_PLER", "URLI_SACO")

spp_comp_combos_coexist =spp_comp_combos %>%
  filter(sp_pair %in% pairs_coexist_change)

t.test(spp_comp_combos_coexist$n_change, spp_comp_combos_coexist$a_change, paired = T)


spp_comp_combos_coexist_pivot = pivot_longer(spp_comp_combos_coexist, cols = c(n_change, a_change), names_to = "type_change", values_to ="change_value")

t.test(spp_comp_combos$n_change, spp_comp_combos$a_change, paired = T)
spp_comp_combos_pivot = pivot_longer(spp_comp_combos, cols = c(n_change, a_change), names_to = "type_change", values_to ="change_value")


ggplot(spp_comp_combos_pivot, aes(x = type_change, y =change_value))+
  theme_classic(base_size = 20) +
  theme(text = element_text(size = 18))+
  geom_boxplot(fill = "light grey", outlier.shape = NA) +
  geom_point(shape = 1, size = 4) +
  ylab("Difference between treatments") +
  xlab (" ") +
  #scale_y_log10() +
  scale_x_discrete(labels=c("a_change"="Competition \n coefficients", "n_change" = "Demographic \n potential")) 

ggsave("./figures/alpha_eta_ratio_box_model5.pdf")
write.csv(spp_comp_combos_pivot,  "./output/n_alph_ratio_output_mod5.csv")


#######  Analysis for Figure 4 in text (trait differences comparisons)

gkt = read.csv("./data/traits_gk.csv")
rownames(gkt) = gkt$species

draw_pairs = read.csv("./output/spp_treat_draw_combos_1000_model5.csv") # or run model script first to get this
final_output_draws = draw_pairs
final_output_draws$fd_superior = ifelse(final_output_draws$fd < 1, 1/final_output_draws$fd, final_output_draws$fd)
gkt_pca = prcomp(gkt[,2:12], scale = T)

gkt_pca$x
gkt_pca$x["ACWR",1]
screeplot(gkt_pca)

#Create PC vectors with the six species in the experiment
pc1= c(gkt_pca$x["ACWR",1], gkt_pca$x["FEMI",1], gkt_pca$x["HOMU",1], gkt_pca$x["PLER",1], gkt_pca$x["SACO",1], gkt_pca$x["URLI",1])
pc2= c(gkt_pca$x["ACWR",2], gkt_pca$x["FEMI",2], gkt_pca$x["HOMU",2], gkt_pca$x["PLER",2], gkt_pca$x["SACO",2], gkt_pca$x["URLI",2])
#Dist matrices for fd and snd
fd_d  = matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))

# The following loop propagates error in FD and ND measurements into the distance matrices used to 
# assess differences in ND and FD among species
dist_fd_d = list()
for( i in 1:100){
  fd_d  = matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))
  for(sp1 in rownames(fd_d)) {  
    for(sp2 in colnames(fd_d)) {    
      current_fd = final_output_draws %>% 
        filter(id == unique(final_output_draws$id)[i] & treatment == 2 & focal == sp1 & competitor == sp2) %>% 
        dplyr::select(fd_superior) %>%
        unlist()
      fd_d[sp1, sp2] = current_fd  
      dist_fd_d[[i]]= as.dist(fd_d)
    }
  }
  print(i)
}

dist_pc1= dist(pc1)

#run a mantel test over all posterior samples to generate distribution of 
# r and p statistics.
mantelr = c()
mantelp = c()
for(i in 1:length(dist_fd_d)){
  test= vegan :: mantel(xdis = dist_pc1, ydis = dist_fd_d[[i]])
  mantelr[i] = test$statistic
  mantelp[i] = test$signif
  print(i)
}

#plot histograms and 95% hd intervals for each metric
hist(mantelr^2); median(mantelr^2); HDInterval::hdi(mantelr^2)
hist(mantelp); median(mantelp); HDInterval::hdi(mantelp)

#check correspondencee between mantel test stats from a single median estimate and the 
# posterior distribution of r calculated above
medians = final_output_draws %>% group_by(focal, competitor, treatment) %>% summarise(fd_median = median(fd_superior))
fd_d  = matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))
for(sp1 in rownames(fd_d)) {  
  for(sp2 in colnames(fd_d)) {    
    current_fd = medians %>% 
      filter(treatment == 2 & focal == sp1 & competitor == sp2) %>% 
      dplyr::select(fd_median)
    fd_d[sp1, sp2] = current_fd$fd_median  }}
dist_fd_d= as.dist(fd_d)

vegan :: mantel(xdis = dist_pc1, ydis = dist_fd_d)
median(mantelr); median(mantelp) #these values should be somewhat close to those from the previous line


#Do traits dissimilarities between pairs correlate with magnitude that their snd and fd change between treatments
diffs = final_output_draws

fd_diff_abs = abs(diffs$fd[diffs$treatment == 1] - diffs$fd[diffs$treatment == 2])

diffs = diffs %>% 
  dplyr :: select(id, focal, competitor, treatment) %>%
  filter(treatment == 1) 

diffs$fd_diff_abs = fd_diff_abs

dist_fd_diffs_abs = list()
for( i in 1:100){
  fd_diffs_abs  = matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))
  for(sp1 in rownames(fd_diffs_abs)) {  
    for(sp2 in colnames(fd_diffs_abs)) {    
      fd_change = diffs %>% 
        filter(id == unique(diffs$id)[i] & focal == sp1 & competitor == sp2) %>% 
        dplyr::select(fd_diff_abs) %>%
        unlist()
      fd_diffs_abs[sp1, sp2] = fd_change  }}
  dist_fd_diffs_abs[[i]]= dist(fd_diffs_abs)
  print(i)
}

for(i in 1:length(dist_fd_diffs_abs)){
  test = vegan :: mantel(ydis = dist_fd_diffs_abs[[i]],xdis = dist_pc1)
  mantelr[i] = test$statistic
  mantelp[i] = test$signif
  print(i)
}

hist(mantelr^2); median(mantelr^2); HDInterval::hdi(mantelr^2)
hist(mantelp); median(mantelp); HDInterval::hdi(mantelp)


medians = final_output_draws %>% group_by(focal, competitor, treatment) %>% summarise(fd_median = median(fd))
diffs = as.data.frame(medians)
fd_diff_abs = abs(diffs$fd_median[diffs$treatment == 1] - diffs$fd_median[diffs$treatment == 2])

diffs = diffs %>% 
  dplyr :: select(focal, competitor, treatment) %>%
  filter(treatment == 1) 

diffs$fd_diff_abs = fd_diff_abs


fd_diffs_abs_all  = matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))
for(sp1 in rownames(fd_diffs_abs_all)) {  
  for(sp2 in colnames(fd_diffs_abs_all)) {    
    fd_change = diffs %>% 
      filter(focal == sp1 & competitor == sp2)  %>%
      dplyr::select(fd_diff_abs) %>%
      unlist()
    fd_diffs_abs_all[sp1, sp2] = fd_change  }}
dist_fd_diffs_abs_all= dist(fd_diffs_abs_all)

vegan ::mantel(xdis = dist_pc1,ydis = dist_fd_diffs_abs_all) 
median(mantelr); median(mantelp) #these values should be somewhat close to those from the previous line

pc1_vec = as.vector(dist_pc1)
fd_diffs_abs_vec = as.vector(dist_fd_diffs_abs_all)
fd_diffs_vec = as.vector(dist_fd_d)

fd_d_df = as.data.frame(gdata :: unmatrix(fd_d)) 
colnames(fd_d_df) = "fd_d"
fd_d_df = tibble::rownames_to_column(fd_d_df, "sp_pair")

fd_diffs_abs_df = as.data.frame(gdata :: unmatrix(fd_diffs_abs))
colnames(fd_diffs_abs_df) = "fd_diff_abs"
fd_diffs_abs_df = tibble::rownames_to_column(fd_diffs_abs_df, "sp_pair")

pc1_mat = as.matrix(dist(pc1, upper = T))
rownames(pc1_mat)= c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")
colnames(pc1_mat)= c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")
pc1_df = as.data.frame(gdata::unmatrix(pc1_mat))
colnames(pc1_df) = "pc1"
pc1_df= tibble::rownames_to_column(pc1_df, "sp_pair")


df_compile = left_join(fd_d_df,fd_diffs_abs_df, by = "sp_pair") %>%
  left_join(., pc1_df, by = "sp_pair")

w_pairs= c("SACO:ACWR", "URLI:ACWR", "ACWR:FEMI", "HOMU:FEMI", 
           "PLER:FEMI", "SACO:FEMI", "URLI:FEMI","ACWR:HOMU", 
           "PLER:HOMU", "SACO:HOMU", "URLI:HOMU", "ACWR:PLER", 
           "SACO:PLER", "URLI:PLER","URLI:SACO")

df_compile= df_compile %>%
  filter(sp_pair %in% w_pairs)
df_compile$label = paste0(substr(df_compile$sp_pair, 1, 2), "-", substr(df_compile$sp_pair, 6, 7))

write_csv(df_compile, "./output/pca_mantel_mod5.csv")

vegan :: mantel(dist_pc1, dist_fd_d)
ggplot(df_compile, aes(x = pc1, y = fd_d)) + geom_point() +geom_smooth(method = "lm")

vegan :: mantel(dist_pc1, dist_fd_diffs_abs_all)
ggplot(df_compile, aes(x = pc1, y = fd_diff_abs)) + geom_point() +geom_smooth(method = "lm")



# set priors for bayesian models.
seedpriors = c(
  set_prior("normal(500, 10000)", nlpar = "lambda", lb = 0, ub = 5000),
  set_prior("normal(.1, 5)", nlpar = "aACWR", lb = .001, ub=2),
  set_prior("normal(.1, 5)", nlpar = "aFEMI", lb = .001, ub=2),
  set_prior("normal(.1, 5)", nlpar = "aHOMU", lb = .001, ub=2),
  set_prior("normal(.1, 5)", nlpar = "aPLER", lb = .001, ub=2),
  set_prior("normal(.1, 5)", nlpar = "aSACO", lb = .001, ub=2),
  set_prior("normal(.1, 5)", nlpar = "aURLI", lb = .001, ub=2),
  set_prior("normal(0, 2)", nlpar = "theta")
)

# Model specifications following table 1
seed_bayes_fit_formula_1 = bf(log(num_seeds)~log((lambda)),
                              lambda ~ 1, nl = TRUE)

seed_bayes_fit_formula_2 = bf(log(num_seeds)~lambda-aACWR*N_acwr-aFEMI*N_femi-aHOMU*N_homu-aPLER*N_pler-aSACO*N_saco-aURLI*N_urli,
                              lambda+ aACWR+ aFEMI+ aHOMU+ aPLER+ aSACO+ aURLI ~ 1, nl = TRUE)

seed_bayes_fit_formula_3 = bf(log(num_seeds)~log((lambda)*exp(-aACWR*(N_acwr)-aFEMI*(N_femi)-aHOMU*(N_homu)-aPLER*(N_pler)-aSACO*(N_saco)-aURLI*(N_urli))),
                              lambda+ aACWR+ aFEMI+ aHOMU+ aPLER+ aSACO+ aURLI ~ 1, nl = TRUE)

seed_bayes_fit_formula_4 = bf(log(num_seeds)~log((lambda)*exp(-aACWR*log(N_acwr+1)-aFEMI*log(N_femi+1)-aHOMU*log(N_homu+1)-aPLER*log(N_pler+1)-aSACO*log(N_saco+1)-aURLI*log(N_urli+1))),
                              lambda+ aACWR+ aFEMI+ aHOMU+ aPLER+ aSACO+ aURLI ~ 1, nl = TRUE)

seed_bayes_fit_formula_5 = bf(log(num_seeds)~log((lambda)/(1+aACWR*N_acwr+aFEMI*N_femi+aHOMU*N_homu+aPLER*N_pler+aSACO*N_saco+aURLI*N_urli)),
                              lambda+ aACWR+ aFEMI+ aHOMU+ aPLER+ aSACO+ aURLI ~ 1, nl = TRUE)

seed_bayes_fit_formula_6 = bf(log(num_seeds)~log(lambda/(1+(N_acwr^aACWR)+(N_femi^aFEMI)+(N_homu^aHOMU)+(N_pler^aPLER)+(N_saco^aSACO)+(N_urli^aURLI))),
                              lambda+ aACWR+ aFEMI+ aHOMU+ aPLER+ aSACO+ aURLI ~ 1, nl = TRUE)

seed_bayes_fit_formula_7 = bf(log(num_seeds)~log((lambda/((1+(aACWR*N_acwr+aFEMI*N_femi+aHOMU*N_homu+aPLER*N_pler+aSACO*N_saco+aURLI*N_urli)^theta)))),
                              lambda+ aACWR+ aFEMI+ aHOMU+ aPLER+ aSACO+ aURLI + theta ~ 1, nl = TRUE)

# Next, we fit the bayesian models for comparison purposes. 
# This code will fit a bunch of models so it can take some time to run
# It takes around 1 hour on an M1 chip with 8 cores in parallel.

#Create a data frame with all combinations of treatment and species 
spp_list = sort(na.omit( unique(seed_data$focal)))
treat_list = sort( na.omit( unique(seed_data$Tr)))
spp_combos = expand.grid(species = spp_list, treatment = treat_list)
# loop through the species list, fit the model, return the dataframe 
out_bayes_compare.waic = list() #list for for loop output
out_bayes_compare.loo = list() #list for for loop output
ncores = 8 # set no. of cores for parallel MCMC runs.

for( i in 1:nrow( spp_combos)){ 
  temp_data = seed_data %>%
    filter(focal == spp_combos[i,1], Tr ==spp_combos[i,2])
  
  seed_bayes_fit_1= brm(seed_bayes_fit_formula_1, family=gaussian(), data = temp_data, prior = seedpriors[1,],  chains = 8, cores = ncores, iter = 10000, backend = 'cmdstanr')  %>%  add_criterion("loo")
  seed_bayes_fit_2= brm(seed_bayes_fit_formula_2, family=gaussian(), data = temp_data, prior = seedpriors[-8,],  chains = 8, cores = ncores, iter = 10000, backend = 'cmdstanr') %>% add_criterion("loo")
  seed_bayes_fit_3= brm(seed_bayes_fit_formula_3, family=gaussian(), data = temp_data, prior = seedpriors[-8,],  chains = 8, cores = ncores, iter = 10000, backend = 'cmdstanr') %>% add_criterion("loo")
  seed_bayes_fit_4= brm(seed_bayes_fit_formula_4, family=gaussian(), data = temp_data, prior = seedpriors[-8,],  chains = 8, cores = ncores, iter = 10000, backend = 'cmdstanr') %>% add_criterion("loo")
  seed_bayes_fit_5= brm(seed_bayes_fit_formula_5, family=gaussian(), data = temp_data, prior = seedpriors[-8,],  chains = 8, cores = ncores, iter = 10000, backend = 'cmdstanr') %>% add_criterion("loo")
  seed_bayes_fit_6= brm(seed_bayes_fit_formula_6, family=gaussian(), data = temp_data, prior = seedpriors[-8,],  chains = 8, cores = ncores, iter = 10000, backend = 'cmdstanr') %>% add_criterion("loo")
  seed_bayes_fit_7= brm(seed_bayes_fit_formula_7, family=gaussian(), data = temp_data, prior = seedpriors,  chains = 8, cores = ncores, iter = 10000, backend = 'cmdstanr') %>% add_criterion("loo")
  
  df.loo = as.data.frame(loo_compare(
    list(
      m1 = seed_bayes_fit_1$criteria$loo,
      m2 = seed_bayes_fit_2$criteria$loo,
      m3 = seed_bayes_fit_3$criteria$loo,
      m4 = seed_bayes_fit_4$criteria$loo,
      m5 = seed_bayes_fit_5$criteria$loo,
      m6 = seed_bayes_fit_6$criteria$loo,
      m7 = seed_bayes_fit_7$criteria$loo
    )
  )
  )
  df.loo$models = rownames(df.loo)
  
  df.waic = as.data.frame(loo_compare(
    list(
      m1 = waic(seed_bayes_fit_1),
      m2 = waic(seed_bayes_fit_2),
      m3 = waic(seed_bayes_fit_3),
      m4 = waic(seed_bayes_fit_4),
      m5 = waic(seed_bayes_fit_5),
      m6 = waic(seed_bayes_fit_6),
      m7 = waic(seed_bayes_fit_7)
    )
  )
  )
  df.waic$models = rownames(df.waic)
  
  out_bayes_compare.waic[[i]] = df.waic
  out_bayes_compare.loo[[i]] = df.loo
}

# calculate delta waic for model comparisons
for(i in 1:length(out_bayes_compare.waic)){
  out_bayes_compare.waic[[i]]$models = rownames(out_bayes_compare.waic[[i]])
  out_bayes_compare.waic[[i]]$del.waic = out_bayes_compare.waic[[i]]$waic - min(out_bayes_compare.waic[[i]]$waic)
}

comparisons = do.call(rbind.data.frame, out_bayes_compare.loo)
meancomps = comparisons %>% group_by(models) %>% summarise_all(mean)

# plot predictive accuracy for different models with their standard error estimates. Values closest to zero
# indicate better leave-out-out accuracy
ggplot(meancomps) +
  geom_point(aes(x = models, y = elpd_diff)) +
  geom_errorbar(aes(x = models, ymin = elpd_diff -se_diff , ymax = elpd_diff +se_diff), width = 0.2)

#calculate standard deviation of measurements across species x treatment combinations
sdcomps = comparisons %>% group_by(models) %>% summarise_all(sd)



### Next, we use the best-fit model (model seven in the above analysis) to calculate ND and FD
### First, we loop over all species x treatment combos
#Create a data frame with all combinations of treatment and species 
spp_list = sort(na.omit( unique(seed_data$focal)))
treat_list = sort( na.omit( unique(seed_data$Tr)))
spp_combos = expand.grid(species = spp_list, treatment = treat_list)

# loop through the species list, fit the model, return the dataframe 
out_bayes = list() #list for for loop output

# fit models as above
for( i in 1:nrow( spp_combos)){ 
  temp_data = seed_data %>%
    filter(focal == spp_combos[i,1], Tr == spp_combos[i,2])
  mod= brm(seed_bayes_fit_formula_7, family=gaussian(), data = temp_data, prior = seedpriors,  chains = 8, cores = 8, iter = 10000, backend = 'cmdstanr')
  out_bayes[[i]] = mod
}

# check MCMC traces and posterior distributions
plot(out_bayes[[1]])
# perform posterior predictive check on model - compares observed y to 100 simulated y's 
brms::pp_check(out_bayes[[1]], ndraws = 100)


out_dfb = list() #list for for loop output
# extract posteror draws for each parameter of interest and reorganise them
for( i in 1:nrow( spp_combos)){ 
  dfb = out_bayes[[i]] %>%
    spread_draws(b_lambda_Intercept, b_aACWR_Intercept, b_aFEMI_Intercept,
                 b_aHOMU_Intercept, b_aPLER_Intercept, b_aSACO_Intercept, b_aURLI_Intercept)
  dfb = dfb[,3:10]
  names(dfb) = c("draw","lambda", "a_ACWR", "a_FEMI", "a_HOMU", "a_PLER", "a_SACO", "a_URLI")
  dfb$focal = rep(spp_combos[i,1], nrow(dfb))
  dfb$treatment = rep(spp_combos[i,2], nrow(dfb))
  
  dfb = dfb %>% pivot_longer(
    cols = starts_with("a"), 
    names_to = c("competitor"), 
    names_prefix = "a_", 
    values_to ="alpha"
  )
  
  out_dfb[[i]] = dfb 
}

posteriors = do.call(rbind.data.frame, out_dfb)

#Make a data frame for the stabilizing niche and fitness differences calculated from each bootstrap's alphas
# take a random sample of 1000 post-warmup posterior draws
draw_list = sort( na.omit( unique(posteriors$draw)))
draw_list = sort(na.omit(sample(draw_list, 1000, replace = F)))
spp_treat_draw_combos = expand.grid(focal = spp_list, treatment = treat_list, competitor = spp_list, id = draw_list)
posteriors = subset(posteriors, draw  %in% draw_list) # subset posteriors to 1000 draws per focal species x treatment combo

#Stabilizing Niche Difference calculations following Van Dyke et al. 2022
stabilizing_niche_diff_bayes = function(df, species1, species2, treat, draw_id) {
  aij = with(df, alpha[focal == species1 & competitor == species2 & treatment== treat & draw == draw_id])
  aji = with(df, alpha[focal == species2 & competitor == species1 & treatment== treat & draw == draw_id])
  ajj = with(df, alpha[focal == species2 & competitor == species2 & treatment== treat &  draw == draw_id])
  aii = with(df, alpha[focal == species1 & competitor == species1 & treatment== treat & draw == draw_id])
  snd = (1 - sqrt((aij * aji)/(ajj * aii)))
  return(snd)
}

stabilizing_niche_diff_bayes(posteriors, "ACWR", "URLI", treat = 1, draw_id = draw_list[[1]])


#Getting eta (ηi) term following Van Dyke et al. 2022
#ηi describes the seeds produced per seed lost from the seed bank for plant species i 
s_g_data = read.csv("./data/s_g_data.csv") #seed survival and germination data
posteriors = merge(posteriors, s_g_data, by = "focal")
posteriors$X = NULL
write_csv(posteriors, "./output/posteriors_model7.csv")

#ηi equation function
get_ni= function(df, species, treat, draw_id){
  lambda = with(df, lambda[ focal == species & treatment == treat & draw == draw_id])[1]
  #print(lambda)
  gi = with( df, g[focal == species & treatment == treat & draw == draw_id])[1]
  #print(gi)
  si = with( df, s[focal == species & treatment == treat & draw == draw_id])[1]
  #print(si)
  
  ni= ((lambda*gi)/(1-((1-gi)*si)))
  return(ni[1])
  
}

get_ni(posteriors, "ACWR", 1, draw_list[1]) #test

#Add snd column to data frame
spp_treat_draw_combos$snd = 0
spp_treat_draw_combos$ni = 0  # Add ni column to data frame


###### This section speeds up the calculation of ND and FD through parallelization
### Depending on the number of draws sampled it can take a significant amount of time
### On an M1 chip w/ 8 cores each loop takes around 5 minutes.
n.cores = parallel::detectCores()
clus = parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = clus)

x = foreach(
  i = 1:n, 
  .combine = 'rbind', .packages = c("magrittr")
) %dopar% {
  sp1 = spp_treat_draw_combos[i, "focal"] %>% unlist
  sp2 = spp_treat_draw_combos[i, "competitor"] %>% unlist
  trt = spp_treat_draw_combos[i , "treatment"] %>% unlist
  draw = spp_treat_draw_combos[i, "id"] %>% unlist
  snd = stabilizing_niche_diff_bayes(posteriors, sp1, sp2, trt, draw)
  ni = get_ni(posteriors, sp1, trt, draw)
  cbind(sp1, sp2, trt, draw, snd, ni)
}

parallel::stopCluster(cl = clus)

spp_treat_draw_combos$snd = x[,"snd"]
spp_treat_draw_combos$ni = x[,"ni"]


#Get fitness differences ------

fitness_diff = function(df1, df2, species1, species2, treat, draw_id) {
  
  ni = with(df1, ni[focal == species1 & treatment == treat & id == draw_id])[1]
  nj = with(df1, ni[focal == species2 & treatment == treat & id == draw_id])[1]
  aij = with(df2, alpha[focal == species1 & competitor == species2 & treatment== treat & draw == draw_id])
  aji = with(df2, alpha[focal == species2 & competitor == species1 & treatment== treat & draw == draw_id])
  ajj = with(df2, alpha[focal == species2 & competitor == species2 & treatment== treat & draw == draw_id])
  aii = with(df2, alpha[focal == species1 & competitor == species1 & treatment== treat & draw == draw_id])
  
  nn= (nj-1)/(ni-1)
  #print(paste("nn: ", nn))
  aa= sqrt((aij * aii)/(ajj * aji))
  #print(paste("aa: ", aa))
  FDij = nn*aa
  #print(FDij[1])
  return(FDij[1])
}

fitness_diff(spp_treat_draw_combos, posteriors, "ACWR", "FEMI", 2, draw_list[1]) #test

#add fitness difference column to data frame
spp_treat_draw_combos$fd = 0

n.cores = parallel::detectCores()
clus = parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = clus)

x = foreach(
  i = 1:n, 
  .combine = 'rbind', .packages = c("magrittr")
) %dopar% {
  sp1 = spp_treat_draw_combos[i, "focal"] %>% unlist
  sp2 = spp_treat_draw_combos[i, "competitor"] %>% unlist
  trt = spp_treat_draw_combos[i , "treatment"] %>% unlist
  draw_id = spp_treat_draw_combos[i, "id"] %>% unlist
  fd = fitness_diff(spp_treat_draw_combos, posteriors, sp1, sp2, trt, draw_id)
  cbind(sp1, sp2, trt, draw_id, fd)
}

parallel::stopCluster(cl = clus)
spp_treat_draw_combos$fd = x[,"fd"]

write_csv(spp_treat_draw_combos, "./output/spp_treat_draw_combos_1000_model7.csv")
#spp_treat_draw_combos = read.csv("./output/spp_treat_draw_combos_1000_model7_new.csv")
spp_treat_draw_combos$sp_pair = paste(spp_treat_draw_combos$focal,spp_treat_draw_combos$competitor,sep='_')

pairs_abc = c("ACWR_FEMI", "ACWR_HOMU", "ACWR_PLER", "ACWR_SACO", "ACWR_URLI", 
              "FEMI_HOMU", "FEMI_PLER", "FEMI_SACO", "FEMI_URLI", "HOMU_PLER",
              "HOMU_SACO", "HOMU_URLI","PLER_SACO", "PLER_URLI", "SACO_URLI",
              "FEMI_ACWR", "HOMU_ACWR" ,"PLER_ACWR","SACO_ACWR", "URLI_ACWR",
              "HOMU_FEMI", "PLER_FEMI", "SACO_FEMI", "URLI_FEMI", "PLER_HOMU",
              "SACO_HOMU", "URLI_HOMU", "SACO_PLER", "URLI_PLER", "URLI_SACO")


spp_treat_draw_combos$fd_superior = ifelse(spp_treat_draw_combos$fd < 1, 1/spp_treat_draw_combos$fd, spp_treat_draw_combos$fd)
spp_treat_draw_combos$fd_sup_sp = ifelse(spp_treat_draw_combos$fd <= 1, 1, 2)
spp_treat_draw_combos_sup = spp_treat_draw_combos %>%
  filter(fd_sup_sp == 2 )
W_superior = with(spp_treat_draw_combos_sup, sp_pair[treatment==1])

draws_pairs_w_sup = spp_treat_draw_combos %>%
  filter(sp_pair %in% W_superior)
draws_pairs_w_sup$treat = factor(draws_pairs_w_sup$treat, levels = c(1, 2))

draws_pairs_unique =  spp_treat_draw_combos %>%
  filter(sp_pair %in% pairs_abc)
draws_pairs_unique$treatment = factor(draws_pairs_unique$treatment, levels = c(1, 2))

draws_pairs_unique$fd_superior = ifelse(draws_pairs_unique$fd < 1, 1/draws_pairs_unique$fd, draws_pairs_unique$fd)

draws_pairs_unique$coexist = ifelse((draws_pairs_unique$snd > (1-1/draws_pairs_unique$fd_superior)), 1, 0 )
draws_pairs_w_sup$treatment = factor(draws_pairs_w_sup$treatment, levels = c(1, 2))
draws_pairs_w_sup$label = paste0(substr(draws_pairs_w_sup$focal, 1, 2), "-", substr(draws_pairs_w_sup$competitor, 1, 2))
write.csv(draws_pairs_w_sup, "./output/plot_data_mod7.csv")
boots_pairs_w_sup = read.csv("./output/boots_pairs_w_sup.csv")
tokeep = boots_pairs_w_sup$label

lot_dat = subset(draws_pairs_w_sup, label %in% tokeep)
medians = plot_dat %>% group_by(label, treatment) %>% summarize(snd = median(snd), fd = median(fd))

# Create coexistence area for plot - min/max fitness difference that permits coexistence
niche_differentiation = seq(from = -.25, to = 1, by = 0.001)
niche_overlap = 1-niche_differentiation
fitness_ratio_min = niche_overlap
fitness_ratio_max = 1/niche_overlap

df = data.frame(niche_diff = niche_differentiation,
                min_fitness_ratio = fitness_ratio_min,
                max_fitness_ratio = fitness_ratio_max)

#Text size 
geom.text.size = 7*(5/14); geom.text.size2 = 6*(5/14)
label.size = 5*(5/14);element.size = 5; theme.size = 7

NDFD_plot_mod7 = ggplot(plot_dat) + 
  theme_classic()+
  theme(text = element_text(size = theme.size))+
  coord_cartesian(xlim = c(-.1, 1), ylim = c(0.1, 100)) +
  geom_line(data = df, aes(x = niche_diff, y = max_fitness_ratio)) +
  geom_line(data = df,  aes(x = niche_diff, y = min_fitness_ratio)) +
  geom_ribbon(data = df, aes(x = niche_diff, ymin = min_fitness_ratio, ymax = max_fitness_ratio), fill = 'grey80') +
  geom_point(aes(x = snd, y = fd, color = treatment, shape = treatment, group = sp_pair), size = 2, alpha = .1) +
  geom_point(data = medians, aes(x = snd, y = fd, fill = treatment, shape = treatment, group = label), size = 3) +
  geom_point(data = boots_pairs_w_sup, aes(x = snd, y = fd, fill = factor(treatment), shape = factor(treatment), group = sp_pair), size = 2, fill = "black") +
  scale_color_manual(values=c('1'="#4E84C4", '2' = "#D16103"), name = "Treatment:", labels = c("Ambient", "Reduced Rain")) +
  scale_fill_manual(values=c('1'="#4E84C4", '2' = "#D16103"), name = "Treatment:", labels = c("Ambient", "Reduced Rain")) +
  scale_shape_manual(values = c('1' = 22, '2' = 21), name = "Treatment:", labels = c("Ambient", "Reduced Rain")) +
  scale_y_log10(expand = c(0,0)) + #setting cut-off and making y on a log scale
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs( x = "Stabilizing niche difference (1-\u03c1)", y = "Fitness difference (Kj/Ki)") +
  theme(axis.title = element_text(size = theme.size), 
        axis.text.x = element_text(size = (theme.size - 2)),
        axis.text.y = element_text(size = (theme.size - 2)),
        legend.title = element_text(size = theme.size),
        legend.text = element_text(size = theme.size),
        legend.position = c(0.88, .035),
        legend.justification = c("bottom"),
        legend.box.just = "bottom",
        legend.margin = margin(1, 1, 1, 1),
        legend.direction = "vertical", 
        strip.background = element_blank(),
        strip.text.x = element_text(size = 8),
        panel.spacing = unit(0,'lines')) +
  guides(fill = guide_legend(title = "Treatment", title.position = "left", cex = 1), col = guide_legend(nrow = 2) , scale = "none")+
  facet_wrap(vars(label), nrow = 4, scales ="free")

plot(NDFD_plot_mod7)

cowplot::save_plot("NDFD_plot_mod7.pdf", new_NDFD_plot, base_height = 8, base_width = 10)

plot_dat$coexist = ifelse((plot_dat$snd > (1-1/plot_dat$fd_superior)), 1, 0 )
plot_dat$Coex1D = plot_dat$snd - (1-1/plot_dat$fd_superior)


coexpct= plot_dat %>% group_by(label, treatment) %>% summarise(coexpct = sum(coexist)/1000)
switchprobs = c()
for(i in 1:length(unique(coexpct$label))){
  pair = unique(coexpct$label)[i]
  tmp = subset(coexpct, label == pair)
  switchprobs[i] = (tmp[which.max(tmp$coexpct),"coexpct"] * (1-tmp[which.min(tmp$coexpct),"coexpct"]))[1,1]
}
coexpct$switchprobs = rep(switchprobs, each = 2)

# plot posterior coexistence predictions with probability of switch
MedianHDI_mod7 = ggplot(plot_dat, aes(x = Coex1D, y = treat, fill = after_stat(quantile))) + 
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE, quantiles = c(0.025, 0.5, 0.975), color = "black", lwd = 0.25) +
  stat_density_ridges(quantile_lines = T, quantiles = 2, fill = NA, color = "black", lwd = 0.25) +
  scale_fill_manual(
    name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#A0A0A0A0", "#0000FFA0")  ) + 
  cowplot::theme_minimal_vgrid() + 
  scale_y_discrete(labels = c("Ambient", "Reduced Rain")) +
  facet_wrap(vars(label), nrow = 4) + scale_x_continuous(breaks = c(-1,0,1), limits = c(-2,2)) +
  xlab("SND - (1-1/FD) (vals > 0 indicate coexistence)") + ylab("") +
  geom_text(data = coexpct, aes(label=round(switchprobs,3)), 
            x = Inf, y = -Inf, hjust=2, vjust=-0.5,
            inherit.aes = FALSE)

plot(MedianHDI_mod7)

save_plot("MedianHDI_mod7", MedianHDI_bestmod, base_height = 10, base_width = 24)

######## REPLICATE FIGURE 3 
seed_data = read.csv("./data/drought_seed_production_data.csv")
seed_data$Tr = ifelse(seed_data$treat == "W", 1, 2)
seed_data$treatment = seed_data$Tr
posteriors = read.csv("./output/posteriors_model7.csv")

spp_list = sort(na.omit( unique(seed_data$focal)))
treat_list = sort( na.omit( unique(seed_data$treatment)))
spp_combos = expand.grid(species = spp_list, treatment = treat_list)

spp_treat_comp_combos = expand.grid(focal = spp_list, competitor = comp_labels, treatment = treat_list)

draw_pairs = spp_treat_comp_combos %>%
  mutate(alpha = 0, alpha_sd = 0, alpha_low = 0, alpha_high = 0, lambda = 0, lambda_low = 0, lambda_high = 0)
draw_pairs$sp_pair = paste(draw_pairs$focal, draw_pairs$competitor, sep = "_")  
for( i in 1:nrow(spp_treat_comp_combos)) {
  sp1 = spp_treat_comp_combos[i, "focal"] %>% unlist
  sp2 = spp_treat_comp_combos[i, "competitor"] %>% unlist
  treatt = spp_treat_comp_combos[i , "treatment"] %>% unlist
  draw_pairs[i, "alpha"] = median(with(posteriors, 
                                       alpha[treatment == treatt & focal == sp1 & competitor == sp2 
                                       ]), na.rm=TRUE)
  draw_pairs[i, "lambda"] = median(with(posteriors, 
                                        lambda[focal == sp1 & competitor == sp2 & 
                                                 treatment == treatt]), na.rm=TRUE)
}  


final_output = read.csv("./output/spp_treat_draw_combos_1000_model7.csv")


mean(with(final_output, ni[focal == "HOMU" & treatment == 1]), na.rm = T)


igr_ratio = function(foc, comp, trt) { #invasion growth rate ratios
  
  nj = mean(with(final_output, ni[focal == comp & treatment == trt]), na.rm = T)
  ni = mean(with(final_output, ni[focal == foc & treatment == trt]), na.rm = T)
  ajj = with(draw_pairs, alpha[focal == comp & competitor == comp & treatment == trt])
  aij = with(draw_pairs, alpha[focal == foc & competitor == comp & treatment == trt])
  n_ratio = (ni-1)/(nj-1)
  a_ratio = ajj/aij
  return(c( n_ratio, a_ratio, n_ratio*a_ratio))
}
igr_ratio("URLI", "HOMU", 1)

comp_labels = sort( na.omit( unique(seed_data$background) ))

spp_treat_comp_combos = expand.grid(species = spp_list, competitor = comp_labels, treatment = treat_list)

spp_treat_comp_combos$n_ratio = 0
spp_treat_comp_combos$a_ratio = 0
spp_treat_comp_combos$product = 0

for(i in 1:nrow(spp_treat_comp_combos)) {
  ii = spp_treat_comp_combos[i, "species"] %>% unlist
  jj = spp_treat_comp_combos[i, "competitor"] %>% unlist
  tt = spp_treat_comp_combos[i, "treatment"] %>% unlist
  
  spp_treat_comp_combos[i, "n_ratio"] = igr_ratio(ii, jj, tt)[1]
  spp_treat_comp_combos[i, "a_ratio"] = igr_ratio(ii, jj, tt)[2]
  spp_treat_comp_combos[i, "product"] = igr_ratio(ii, jj, tt)[3]
  
}
head(spp_treat_comp_combos)
spp_treat_comp_combos$larger = ifelse(spp_treat_comp_combos$a_ratio > spp_treat_comp_combos$n_ratio, "a", "n")
spp_treat_comp_combos$invade = ifelse(spp_treat_comp_combos$product>1, "yes", "no")
head(spp_treat_comp_combos)

invaders = spp_treat_comp_combos %>%
  filter(invade == "yes")

#Which ratio changes more in invasion growth rate inequality?
igr_change = function(foc, comp) {
  
  nj_d = mean(with(final_output, ni[focal == comp & treatment == 2]), na.rm = T)
  ni_d = mean(with(final_output, ni[focal == foc & treatment == 2]), na.rm = T)
  ajj_d = with(draw_pairs, alpha[focal == comp & competitor == comp & treatment == 2])
  aij_d = with(draw_pairs, alpha[focal == foc & competitor == comp & treatment == 2])
  n_ratio_d = log10((ni_d-1)/(nj_d-1))
  a_ratio_d = log10(ajj_d/aij_d)
  
  nj_w = mean(with(final_output, ni[focal == comp & treatment == 1]), na.rm = T)
  ni_w = mean(with(final_output, ni[focal == foc & treatment == 1]), na.rm = T)
  ajj_w = with(draw_pairs, alpha[focal == comp & competitor == comp & treatment == 1])
  aij_w = with(draw_pairs, alpha[focal == foc & competitor == comp & treatment == 1])
  n_ratio_w = log10((ni_w-1)/(nj_w-1))
  a_ratio_w = log10(ajj_w/aij_w)
  
  nc=abs(n_ratio_w - n_ratio_d)
  ac=abs(a_ratio_w - a_ratio_d)
  return(c(nc, ac))
}
igr_change("HOMU", "URLI")

#Create data frame
spp_list = sort(na.omit( unique(seed_data$focal)))
comp_labels = sort( na.omit( unique(seed_data$background) ))
spp_comp_combos = expand.grid(species = spp_list, competitor = comp_labels)
spp_comp_combos$n_change = 0
spp_comp_combos$a_change = 0


for(i in 1:nrow(spp_comp_combos)) {
  ii = spp_comp_combos[i, "species"] %>% unlist
  jj = spp_comp_combos[i, "competitor"] %>% unlist
  spp_comp_combos[i, "n_change"] = igr_change(ii, jj)[1]
  spp_comp_combos[i, "a_change"] = igr_change(ii, jj)[2]
}

spp_comp_combos$larger = ifelse(abs(spp_comp_combos$a_change)> abs(spp_comp_combos$n_change), "a", 'n')

spp_comp_combos = spp_comp_combos %>%
  filter(a_change != n_change)
total_pairs = c("ACWR_FEMI", "ACWR_HOMU", "ACWR_PLER", "SACO_ACWR", "URLI_ACWR", 
                "HOMU_FEMI", "PLER_FEMI", "SACO_FEMI", "URLI_FEMI", "PLER_HOMU",
                "SACO_HOMU", "URLI_HOMU","SACO_PLER", "URLI_PLER", "URLI_SACO")
spp_comp_combos$sp_pair = paste(spp_comp_combos$species, spp_comp_combos$competitor, sep = "_")

spp_comp_combos =spp_comp_combos %>%
  filter(sp_pair %in% total_pairs)

t.test(spp_comp_combos$n_change, spp_comp_combos$a_change, paired = T)

pairs_coexist_change =c("ACWR_FEMI", "ACWR_HOMU", "ACWR_PLER", "SACO_ACWR", 
                        "URLI_ACWR",  "HOMU_FEMI", "PLER_FEMI", "URLI_FEMI",
                        "SACO_PLER", "URLI_SACO")

spp_comp_combos_coexist =spp_comp_combos %>%
  filter(sp_pair %in% pairs_coexist_change)

t.test(spp_comp_combos_coexist$n_change, spp_comp_combos_coexist$a_change, paired = T)


spp_comp_combos_coexist_pivot = pivot_longer(spp_comp_combos_coexist, cols = c(n_change, a_change), names_to = "type_change", values_to ="change_value")

t.test(spp_comp_combos$n_change, spp_comp_combos$a_change, paired = T)
spp_comp_combos_pivot = pivot_longer(spp_comp_combos, cols = c(n_change, a_change), names_to = "type_change", values_to ="change_value")


ggplot(spp_comp_combos_pivot, aes(x = type_change, y =change_value))+
  theme_classic(base_size = 20) +
  theme(text = element_text(size = 18))+
  geom_boxplot(fill = "light grey", outlier.shape = NA) +
  geom_point(shape = 1, size = 4) +
  ylab("Difference between treatments") +
  xlab (" ") +
  #scale_y_log10() +
  scale_x_discrete(labels=c("a_change"="Competition \n coefficients", "n_change" = "Demographic \n potential")) 

ggsave("./figures/alpha_eta_ratio_box_model7.pdf")
write.csv(spp_comp_combos_pivot,  "./output/n_alph_ratio_output_mod7.csv")


#######



gkt = read.csv("./data/traits_gk.csv")
rownames(gkt) = gkt$species

draw_pairs = read.csv("./output/spp_treat_draw_combos_1000_model7.csv") # or run model script first to get this
final_output_draws = draw_pairs
final_output_draws$fd_superior = ifelse(final_output_draws$fd < 1, 1/final_output_draws$fd, final_output_draws$fd)
gkt_pca = prcomp(gkt[,2:12], scale = T)

gkt_pca$x
gkt_pca$x["ACWR",1]
screeplot(gkt_pca)

#Create PC vectors with the six species in the experiment
pc1= c(gkt_pca$x["ACWR",1], gkt_pca$x["FEMI",1], gkt_pca$x["HOMU",1], gkt_pca$x["PLER",1], gkt_pca$x["SACO",1], gkt_pca$x["URLI",1])
pc2= c(gkt_pca$x["ACWR",2], gkt_pca$x["FEMI",2], gkt_pca$x["HOMU",2], gkt_pca$x["PLER",2], gkt_pca$x["SACO",2], gkt_pca$x["URLI",2])
#Dist matrices for fd and snd
fd_d  = matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))

# The following loop propagates error in FD and ND measurements into the distance matrices used to 
# assess differences in ND and FD among species
dist_fd_d = list()
for( i in 1:100){
  fd_d  = matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))
  for(sp1 in rownames(fd_d)) {  
    for(sp2 in colnames(fd_d)) {    
      current_fd = final_output_draws %>% 
        filter(id == unique(final_output_draws$id)[i] & treatment == 2 & focal == sp1 & competitor == sp2) %>% 
        dplyr::select(fd_superior) %>%
        unlist()
      fd_d[sp1, sp2] = current_fd  
      dist_fd_d[[i]]= as.dist(fd_d)
    }
  }
  print(i)
}

dist_pc1= dist(pc1)

#run a mantel test over all posterior samples to generate distribution of 
# r and p statistics.
mantelr = c()
mantelp = c()
for(i in 1:length(dist_fd_d)){
  test= vegan :: mantel(xdis = dist_pc1, ydis = dist_fd_d[[i]])
  mantelr[i] = test$statistic
  mantelp[i] = test$signif
  print(i)
}

#plot histograms and 95% hd intervals for each metric
hist(mantelr^2); median(mantelr^2); HDInterval::hdi(mantelr^2)
hist(mantelp); median(mantelp); HDInterval::hdi(mantelp)

#check correspondencee between mantel test stats from a single median estimate and the 
# posterior distribution of r calculated above
medians = final_output_draws %>% group_by(focal, competitor, treatment) %>% summarise(fd_median = median(fd_superior))
fd_d  = matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))
for(sp1 in rownames(fd_d)) {  
  for(sp2 in colnames(fd_d)) {    
    current_fd = medians %>% 
      filter(treatment == 2 & focal == sp1 & competitor == sp2) %>% 
      dplyr::select(fd_median)
    fd_d[sp1, sp2] = current_fd$fd_median  }}
dist_fd_d= as.dist(fd_d)

vegan :: mantel(xdis = dist_pc1, ydis = dist_fd_d)
median(mantelr); median(mantelp) #these values should be somewhat close to those from the previous line


#Do traits dissimilarities between pairs correlate with magnitude that their snd and fd change between treatments
diffs = final_output_draws

fd_diff_abs = abs(diffs$fd[diffs$treatment == 1] - diffs$fd[diffs$treatment == 2])

diffs = diffs %>% 
  dplyr :: select(id, focal, competitor, treatment) %>%
  filter(treatment == 1) 

diffs$fd_diff_abs = fd_diff_abs

dist_fd_diffs_abs = list()
for( i in 1:100){
  fd_diffs_abs  = matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))
  for(sp1 in rownames(fd_diffs_abs)) {  
    for(sp2 in colnames(fd_diffs_abs)) {    
      fd_change = diffs %>% 
        filter(id == unique(diffs$id)[i] & focal == sp1 & competitor == sp2) %>% 
        dplyr::select(fd_diff_abs) %>%
        unlist()
      fd_diffs_abs[sp1, sp2] = fd_change  }}
  dist_fd_diffs_abs[[i]]= dist(fd_diffs_abs)
  print(i)
}

for(i in 1:length(dist_fd_diffs_abs)){
  test = vegan :: mantel(ydis = dist_fd_diffs_abs[[i]],xdis = dist_pc1)
  mantelr[i] = test$statistic
  mantelp[i] = test$signif
  print(i)
}

hist(mantelr^2); median(mantelr^2); HDInterval::hdi(mantelr^2)
hist(mantelp); median(mantelp); HDInterval::hdi(mantelp)


medians = final_output_draws %>% group_by(focal, competitor, treatment) %>% summarise(fd_median = median(fd))
diffs = as.data.frame(medians)
fd_diff_abs = abs(diffs$fd_median[diffs$treatment == 1] - diffs$fd_median[diffs$treatment == 2])

diffs = diffs %>% 
  dplyr :: select(focal, competitor, treatment) %>%
  filter(treatment == 1) 

diffs$fd_diff_abs = fd_diff_abs


fd_diffs_abs_all  = matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))
for(sp1 in rownames(fd_diffs_abs_all)) {  
  for(sp2 in colnames(fd_diffs_abs_all)) {    
    fd_change = diffs %>% 
      filter(focal == sp1 & competitor == sp2)  %>%
      dplyr::select(fd_diff_abs) %>%
      unlist()
    fd_diffs_abs_all[sp1, sp2] = fd_change  }}
dist_fd_diffs_abs_all= dist(fd_diffs_abs_all)

vegan ::mantel(xdis = dist_pc1,ydis = dist_fd_diffs_abs_all) 
median(mantelr); median(mantelp) #these values should be somewhat close to those from the previous line

pc1_vec = as.vector(dist_pc1)
fd_diffs_abs_vec = as.vector(dist_fd_diffs_abs_all)
fd_diffs_vec = as.vector(dist_fd_d)

fd_d_df = as.data.frame(gdata :: unmatrix(fd_d)) 
colnames(fd_d_df) = "fd_d"
fd_d_df = tibble::rownames_to_column(fd_d_df, "sp_pair")

fd_diffs_abs_df = as.data.frame(gdata :: unmatrix(fd_diffs_abs))
colnames(fd_diffs_abs_df) = "fd_diff_abs"
fd_diffs_abs_df = tibble::rownames_to_column(fd_diffs_abs_df, "sp_pair")

pc1_mat = as.matrix(dist(pc1, upper = T))
rownames(pc1_mat)= c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")
colnames(pc1_mat)= c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")
pc1_df = as.data.frame(gdata::unmatrix(pc1_mat))
colnames(pc1_df) = "pc1"
pc1_df= tibble::rownames_to_column(pc1_df, "sp_pair")


df_compile = left_join(fd_d_df,fd_diffs_abs_df, by = "sp_pair") %>%
  left_join(., pc1_df, by = "sp_pair")

w_pairs= c("SACO:ACWR", "URLI:ACWR", "ACWR:FEMI", "HOMU:FEMI", 
           "PLER:FEMI", "SACO:FEMI", "URLI:FEMI","ACWR:HOMU", 
           "PLER:HOMU", "SACO:HOMU", "URLI:HOMU", "ACWR:PLER", 
           "SACO:PLER", "URLI:PLER","URLI:SACO")

df_compile= df_compile %>%
  filter(sp_pair %in% w_pairs)
df_compile$label = paste0(substr(df_compile$sp_pair, 1, 2), "-", substr(df_compile$sp_pair, 6, 7))

write_csv(df_compile, "./output/pca_mantel_mod7.csv")

vegan :: mantel(dist_pc1, dist_fd_d)
ggplot(df_compile, aes(x = pc1, y = fd_d)) + geom_point() +geom_smooth(method = "lm")

vegan :: mantel(dist_pc1, dist_fd_diffs_abs_all)
ggplot(df_compile, aes(x = pc1, y = fd_diff_abs)) + geom_point() +geom_smooth(method = "lm")























#######  Replicate the previous analysis with the author's BH model (model 5)

options(scipen = 5)
set.seed(1234)
setwd("Downloads/Sedgwick_public/")
seed_data = read.csv("./data/drought_seed_production_data.csv")
seed_data$Tr = ifelse(seed_data$treat == "W", 1, 2)
seed_data = seed_data %>% replace_na(list(background = "ACWR")) # putting acwr instead of NA for the lambdas
seed_data$N_acwr = ifelse(seed_data$background == "ACWR", seed_data$num_comp, 0)
seed_data$N_femi = ifelse(seed_data$background == "FEMI", seed_data$num_comp, 0)
seed_data$N_homu = ifelse(seed_data$background == "HOMU", seed_data$num_comp, 0)
seed_data$N_pler = ifelse(seed_data$background == "PLER", seed_data$num_comp, 0)
seed_data$N_saco = ifelse(seed_data$background == "SACO", seed_data$num_comp, 0)
seed_data$N_urli = ifelse(seed_data$background == "URLI", seed_data$num_comp, 0)


### Next, we use the best-fit model (model seven in the above analysis) to calculate ND and FD
### First, we loop over all species x treatment combos
#Create a data frame with all combinations of treatment and species 
spp_list = sort(na.omit( unique(seed_data$focal)))
treat_list = sort( na.omit( unique(seed_data$Tr)))
spp_combos = expand.grid(species = spp_list, treatment = treat_list)

# loop through the species list, fit the model, return the dataframe 
out_bayes = list() #list for for loop output

# fit models as above
for( i in 1:nrow( spp_combos)){ 
  temp_data = seed_data %>%
    filter(focal == spp_combos[i,1], Tr == spp_combos[i,2])
  mod= brm(seed_bayes_fit_formula_4, family=gaussian(), data = temp_data, prior = seedpriors[-8,],  chains = 8, cores = 8, iter = 10000, backend = 'cmdstanr')
  out_bayes[[i]] = mod
}

# check MCMC traces and posterior distributions
plot(out_bayes[[1]])
# perform posterior predictive check on model - compares observed y to 100 simulated y's 
brms::pp_check(out_bayes[[12]], ndraws = 100)


out_dfb = list() #list for for loop output
# extract posteror draws for each parameter of interest and reorganise them
for( i in 1:nrow( spp_combos)){ 
  dfb = out_bayes[[i]] %>%
    spread_draws(b_lambda_Intercept, b_aACWR_Intercept, b_aFEMI_Intercept,
                 b_aHOMU_Intercept, b_aPLER_Intercept, b_aSACO_Intercept, b_aURLI_Intercept)
  dfb = dfb[,3:10]
  names(dfb) = c("draw","lambda", "a_ACWR", "a_FEMI", "a_HOMU", "a_PLER", "a_SACO", "a_URLI")
  dfb$focal = rep(spp_combos[i,1], nrow(dfb))
  dfb$treatment = rep(spp_combos[i,2], nrow(dfb))
  
  dfb = dfb %>% pivot_longer(
    cols = starts_with("a"), 
    names_to = c("competitor"), 
    names_prefix = "a_", 
    values_to ="alpha"
  )
  
  out_dfb[[i]] = dfb 
}

posteriors = do.call(rbind.data.frame, out_dfb)

#Make a data frame for the stabilizing niche and fitness differences calculated from each bootstrap's alphas
# take a random sample of 1000 post-warmup posterior draws
draw_list = sort( na.omit( unique(posteriors$draw)))
draw_list = sort(na.omit(sample(draw_list, 1000, replace = F)))
spp_treat_draw_combos = expand.grid(focal = spp_list, treatment = treat_list, competitor = spp_list, id = draw_list)
posteriors = subset(posteriors, draw  %in% draw_list) # subset posteriors to 1000 draws per focal species x treatment combo

#Stabilizing Niche Difference calculations following Van Dyke et al. 2022
stabilizing_niche_diff_bayes = function(df, species1, species2, treat, draw_id) {
  aij = with(df, alpha[focal == species1 & competitor == species2 & treatment== treat & draw == draw_id])
  aji = with(df, alpha[focal == species2 & competitor == species1 & treatment== treat & draw == draw_id])
  ajj = with(df, alpha[focal == species2 & competitor == species2 & treatment== treat &  draw == draw_id])
  aii = with(df, alpha[focal == species1 & competitor == species1 & treatment== treat & draw == draw_id])
  snd = (1 - sqrt((aij * aji)/(ajj * aii)))
  return(snd)
}

stabilizing_niche_diff_bayes(posteriors, "ACWR", "URLI", treat = 1, draw_id = draw_list[[1]])


#Getting eta (ηi) term following Van Dyke et al. 2022
#ηi describes the seeds produced per seed lost from the seed bank for plant species i 
s_g_data = read.csv("./data/s_g_data.csv") #seed survival and germination data
posteriors = merge(posteriors, s_g_data, by = "focal")
posteriors$X = NULL
write_csv(posteriors, "./output/posteriors_model5.csv")

#ηi equation function
get_ni= function(df, species, treat, draw_id){
  lambda = with(df, lambda[ focal == species & treatment == treat & draw == draw_id])[1]
  #print(lambda)
  gi = with( df, g[focal == species & treatment == treat & draw == draw_id])[1]
  #print(gi)
  si = with( df, s[focal == species & treatment == treat & draw == draw_id])[1]
  #print(si)
  
  ni= ((lambda*gi)/(1-((1-gi)*si)))
  return(ni[1])
  
}

get_ni(posteriors, "ACWR", 1, draw_list[1]) #test

#Add snd column to data frame
spp_treat_draw_combos$snd = 0
spp_treat_draw_combos$ni = 0  # Add ni column to data frame


###### This section speeds up the calculation of ND and FD through parallelization
### Depending on the number of draws sampled it can take a significant amount of time
### On an M1 chip w/ 8 cores each loop takes around 5 minutes.
n.cores = parallel::detectCores()
clus = parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = clus)

x = foreach(
  i = 1:n, 
  .combine = 'rbind', .packages = c("magrittr")
) %dopar% {
  sp1 = spp_treat_draw_combos[i, "focal"] %>% unlist
  sp2 = spp_treat_draw_combos[i, "competitor"] %>% unlist
  trt = spp_treat_draw_combos[i , "treatment"] %>% unlist
  draw = spp_treat_draw_combos[i, "id"] %>% unlist
  snd = stabilizing_niche_diff_bayes(posteriors, sp1, sp2, trt, draw)
  ni = get_ni(posteriors, sp1, trt, draw)
  cbind(sp1, sp2, trt, draw, snd, ni)
}

parallel::stopCluster(cl = clus)

spp_treat_draw_combos$snd = x[,"snd"]
spp_treat_draw_combos$ni = x[,"ni"]


#Get fitness differences ------

fitness_diff = function(df1, df2, species1, species2, treat, draw_id) {
  
  ni = with(df1, ni[focal == species1 & treatment == treat & id == draw_id])[1]
  nj = with(df1, ni[focal == species2 & treatment == treat & id == draw_id])[1]
  aij = with(df2, alpha[focal == species1 & competitor == species2 & treatment== treat & draw == draw_id])
  aji = with(df2, alpha[focal == species2 & competitor == species1 & treatment== treat & draw == draw_id])
  ajj = with(df2, alpha[focal == species2 & competitor == species2 & treatment== treat & draw == draw_id])
  aii = with(df2, alpha[focal == species1 & competitor == species1 & treatment== treat & draw == draw_id])
  
  nn= (nj-1)/(ni-1)
  #print(paste("nn: ", nn))
  aa= sqrt((aij * aii)/(ajj * aji))
  #print(paste("aa: ", aa))
  FDij = nn*aa
  #print(FDij[1])
  return(FDij[1])
}

fitness_diff(spp_treat_draw_combos, posteriors, "ACWR", "FEMI", 2, draw_list[1]) #test

#add fitness difference column to data frame
spp_treat_draw_combos$fd = 0

n.cores = parallel::detectCores()
clus = parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = clus)

x = foreach(
  i = 1:n, 
  .combine = 'rbind', .packages = c("magrittr")
) %dopar% {
  sp1 = spp_treat_draw_combos[i, "focal"] %>% unlist
  sp2 = spp_treat_draw_combos[i, "competitor"] %>% unlist
  trt = spp_treat_draw_combos[i , "treatment"] %>% unlist
  draw_id = spp_treat_draw_combos[i, "id"] %>% unlist
  fd = fitness_diff(spp_treat_draw_combos, posteriors, sp1, sp2, trt, draw_id)
  cbind(sp1, sp2, trt, draw_id, fd)
}

parallel::stopCluster(cl = clus)
spp_treat_draw_combos$fd = x[,"fd"]

write_csv(spp_treat_draw_combos, "./output/spp_treat_draw_combos_1000_model5.csv")
#spp_treat_draw_combos = read.csv("./output/spp_treat_draw_combos_1000_model7_new.csv")
spp_treat_draw_combos$sp_pair = paste(spp_treat_draw_combos$focal,spp_treat_draw_combos$competitor,sep='_')

pairs_abc = c("ACWR_FEMI", "ACWR_HOMU", "ACWR_PLER", "ACWR_SACO", "ACWR_URLI", 
              "FEMI_HOMU", "FEMI_PLER", "FEMI_SACO", "FEMI_URLI", "HOMU_PLER",
              "HOMU_SACO", "HOMU_URLI","PLER_SACO", "PLER_URLI", "SACO_URLI",
              "FEMI_ACWR", "HOMU_ACWR" ,"PLER_ACWR","SACO_ACWR", "URLI_ACWR",
              "HOMU_FEMI", "PLER_FEMI", "SACO_FEMI", "URLI_FEMI", "PLER_HOMU",
              "SACO_HOMU", "URLI_HOMU", "SACO_PLER", "URLI_PLER", "URLI_SACO")


spp_treat_draw_combos$fd_superior = ifelse(spp_treat_draw_combos$fd < 1, 1/spp_treat_draw_combos$fd, spp_treat_draw_combos$fd)
spp_treat_draw_combos$fd_sup_sp = ifelse(spp_treat_draw_combos$fd <= 1, 1, 2)
spp_treat_draw_combos_sup = spp_treat_draw_combos %>%
  filter(fd_sup_sp == 2 )
W_superior = with(spp_treat_draw_combos_sup, sp_pair[treatment==1])

draws_pairs_w_sup = spp_treat_draw_combos %>%
  filter(sp_pair %in% W_superior)
draws_pairs_w_sup$treat = factor(draws_pairs_w_sup$treat, levels = c(1, 2))

draws_pairs_unique =  spp_treat_draw_combos %>%
  filter(sp_pair %in% pairs_abc)
draws_pairs_unique$treatment = factor(draws_pairs_unique$treatment, levels = c(1, 2))

draws_pairs_unique$fd_superior = ifelse(draws_pairs_unique$fd < 1, 1/draws_pairs_unique$fd, draws_pairs_unique$fd)

draws_pairs_unique$coexist = ifelse((draws_pairs_unique$snd > (1-1/draws_pairs_unique$fd_superior)), 1, 0 )
draws_pairs_w_sup$treatment = factor(draws_pairs_w_sup$treatment, levels = c(1, 2))
draws_pairs_w_sup$label = paste0(substr(draws_pairs_w_sup$focal, 1, 2), "-", substr(draws_pairs_w_sup$competitor, 1, 2))
write.csv(draws_pairs_w_sup, "./output/plot_data_mod5.csv")
boots_pairs_w_sup = read.csv("./output/boots_pairs_w_sup.csv")
tokeep = boots_pairs_w_sup$label

plot_dat = subset(draws_pairs_w_sup, label %in% tokeep)
medians = plot_dat %>% group_by(label, treatment) %>% summarize(snd = median(snd), fd = median(fd))

# Create coexistence area for plot - min/max fitness difference that permits coexistence
niche_differentiation = seq(from = -.25, to = 1, by = 0.001)
niche_overlap = 1-niche_differentiation
fitness_ratio_min = niche_overlap
fitness_ratio_max = 1/niche_overlap

df = data.frame(niche_diff = niche_differentiation,
                min_fitness_ratio = fitness_ratio_min,
                max_fitness_ratio = fitness_ratio_max)

#Text size 
geom.text.size = 7*(5/14); geom.text.size2 = 6*(5/14)
label.size = 5*(5/14);element.size = 5; theme.size = 7

NDFD_plot_mod5 = ggplot(plot_dat) + 
  theme_classic()+
  theme(text = element_text(size = theme.size))+
  coord_cartesian(xlim = c(-.1, 1), ylim = c(0.1, 100)) +
  geom_line(data = df, aes(x = niche_diff, y = max_fitness_ratio)) +
  geom_line(data = df,  aes(x = niche_diff, y = min_fitness_ratio)) +
  geom_ribbon(data = df, aes(x = niche_diff, ymin = min_fitness_ratio, ymax = max_fitness_ratio), fill = 'grey80') +
  geom_point(aes(x = snd, y = fd, color = treatment, shape = treatment, group = sp_pair), size = 2, alpha = .1) +
  geom_point(data = medians, aes(x = snd, y = fd, fill = treatment, shape = treatment, group = label), size = 3) +
  geom_point(data = boots_pairs_w_sup, aes(x = snd, y = fd, fill = factor(treatment), shape = factor(treatment), group = sp_pair), size = 2, fill = "black") +
  scale_color_manual(values=c('1'="#4E84C4", '2' = "#D16103"), name = "Treatment:", labels = c("Ambient", "Reduced Rain")) +
  scale_fill_manual(values=c('1'="#4E84C4", '2' = "#D16103"), name = "Treatment:", labels = c("Ambient", "Reduced Rain")) +
  scale_shape_manual(values = c('1' = 22, '2' = 21), name = "Treatment:", labels = c("Ambient", "Reduced Rain")) +
  scale_y_log10(expand = c(0,0)) + #setting cut-off and making y on a log scale
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs( x = "Stabilizing niche difference (1-\u03c1)", y = "Fitness difference (Kj/Ki)") +
  theme(axis.title = element_text(size = theme.size), 
        axis.text.x = element_text(size = (theme.size - 2)),
        axis.text.y = element_text(size = (theme.size - 2)),
        legend.title = element_text(size = theme.size),
        legend.text = element_text(size = theme.size),
        legend.position = c(0.88, .035),
        legend.justification = c("bottom"),
        legend.box.just = "bottom",
        legend.margin = margin(1, 1, 1, 1),
        legend.direction = "vertical", 
        strip.background = element_blank(),
        strip.text.x = element_text(size = 8),
        panel.spacing = unit(0,'lines')) +
  guides(fill = guide_legend(title = "Treatment", title.position = "left", cex = 1), col = guide_legend(nrow = 2) , scale = "none")+
  facet_wrap(vars(label), nrow = 4, scales ="free")

plot(NDFD_plot_mod5)

cowplot::save_plot("NDFD_plot_mod5.pdf", NDFD_plot_mod5, base_height = 8, base_width = 10)

plot_dat$coexist = ifelse((plot_dat$snd > (1-1/plot_dat$fd_superior)), 1, 0 )
plot_dat$Coex1D = plot_dat$snd - (1-1/plot_dat$fd_superior)


coexpct= plot_dat %>% group_by(label, treatment) %>% summarise(coexpct = sum(coexist)/1000)
switchprobs = c()
for(i in 1:length(unique(coexpct$label))){
  pair = unique(coexpct$label)[i]
  tmp = subset(coexpct, label == pair)
  switchprobs[i] = (tmp[which.max(tmp$coexpct),"coexpct"] * (1-tmp[which.min(tmp$coexpct),"coexpct"]))[1,1]
}
coexpct$switchprobs = rep(switchprobs, each = 2)

# plot posterior coexistence predictions with probability of switch
MedianHDI_mod5 = ggplot(plot_dat, aes(x = Coex1D, y = treat, fill = after_stat(quantile))) + 
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE, quantiles = c(0.025, 0.5, 0.975), color = "black", lwd = 0.25) +
  stat_density_ridges(quantile_lines = T, quantiles = 2, fill = NA, color = "black", lwd = 0.25) +
  scale_fill_manual(
    name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#A0A0A0A0", "#0000FFA0")  ) + 
  cowplot::theme_minimal_vgrid() + 
  scale_y_discrete(labels = c("Ambient", "Reduced Rain")) +
  facet_wrap(vars(label), nrow = 4) + scale_x_continuous(breaks = c(-1,0,1), limits = c(-2,2)) +
  xlab("SND - (1-1/FD) (vals > 0 indicate coexistence)") + ylab("") +
  geom_text(data = coexpct, aes(label=round(switchprobs,3)), 
            x = Inf, y = -Inf, hjust=2, vjust=-0.5,
            inherit.aes = FALSE)

plot(MedianHDI_mod5)

cowplot::save_plot("MedianHDI_mod5.pdf", MedianHDI_mod5, base_height = 10, base_width = 24)

######## REPLICATE FIGURE 3 
seed_data = read.csv("./data/drought_seed_production_data.csv")
seed_data$Tr = ifelse(seed_data$treat == "W", 1, 2)
seed_data$treatment = seed_data$Tr
posteriors = read.csv("./output/posteriors_model5.csv")

spp_list = sort(na.omit( unique(seed_data$focal)))
treat_list = sort( na.omit( unique(seed_data$treatment)))
spp_combos = expand.grid(species = spp_list, treatment = treat_list)

spp_treat_comp_combos = expand.grid(focal = spp_list, competitor = comp_labels, treatment = treat_list)

draw_pairs = spp_treat_comp_combos %>%
  mutate(alpha = 0, alpha_sd = 0, alpha_low = 0, alpha_high = 0, lambda = 0, lambda_low = 0, lambda_high = 0)
draw_pairs$sp_pair = paste(draw_pairs$focal, draw_pairs$competitor, sep = "_")  
for( i in 1:nrow(spp_treat_comp_combos)) {
  sp1 = spp_treat_comp_combos[i, "focal"] %>% unlist
  sp2 = spp_treat_comp_combos[i, "competitor"] %>% unlist
  treatt = spp_treat_comp_combos[i , "treatment"] %>% unlist
  draw_pairs[i, "alpha"] = median(with(posteriors, 
                                       alpha[treatment == treatt & focal == sp1 & competitor == sp2 
                                       ]), na.rm=TRUE)
  draw_pairs[i, "lambda"] = median(with(posteriors, 
                                        lambda[focal == sp1 & competitor == sp2 & 
                                                 treatment == treatt]), na.rm=TRUE)
}  


final_output = read.csv("./output/spp_treat_draw_combos_1000_model5.csv")


mean(with(final_output, ni[focal == "HOMU" & treatment == 1]), na.rm = T)


igr_ratio = function(foc, comp, trt) { #invasion growth rate ratios
  
  nj = mean(with(final_output, ni[focal == comp & treatment == trt]), na.rm = T)
  ni = mean(with(final_output, ni[focal == foc & treatment == trt]), na.rm = T)
  ajj = with(draw_pairs, alpha[focal == comp & competitor == comp & treatment == trt])
  aij = with(draw_pairs, alpha[focal == foc & competitor == comp & treatment == trt])
  n_ratio = (ni-1)/(nj-1)
  a_ratio = ajj/aij
  return(c( n_ratio, a_ratio, n_ratio*a_ratio))
}
igr_ratio("URLI", "HOMU", 1)

comp_labels = sort( na.omit( unique(seed_data$background) ))

spp_treat_comp_combos = expand.grid(species = spp_list, competitor = comp_labels, treatment = treat_list)

spp_treat_comp_combos$n_ratio = 0
spp_treat_comp_combos$a_ratio = 0
spp_treat_comp_combos$product = 0

for(i in 1:nrow(spp_treat_comp_combos)) {
  ii = spp_treat_comp_combos[i, "species"] %>% unlist
  jj = spp_treat_comp_combos[i, "competitor"] %>% unlist
  tt = spp_treat_comp_combos[i, "treatment"] %>% unlist
  
  spp_treat_comp_combos[i, "n_ratio"] = igr_ratio(ii, jj, tt)[1]
  spp_treat_comp_combos[i, "a_ratio"] = igr_ratio(ii, jj, tt)[2]
  spp_treat_comp_combos[i, "product"] = igr_ratio(ii, jj, tt)[3]
  
}
head(spp_treat_comp_combos)
spp_treat_comp_combos$larger = ifelse(spp_treat_comp_combos$a_ratio > spp_treat_comp_combos$n_ratio, "a", "n")
spp_treat_comp_combos$invade = ifelse(spp_treat_comp_combos$product>1, "yes", "no")
head(spp_treat_comp_combos)

invaders = spp_treat_comp_combos %>%
  filter(invade == "yes")

#Which ratio changes more in invasion growth rate inequality?
igr_change = function(foc, comp) {
  
  nj_d = mean(with(final_output, ni[focal == comp & treatment == 2]), na.rm = T)
  ni_d = mean(with(final_output, ni[focal == foc & treatment == 2]), na.rm = T)
  ajj_d = with(draw_pairs, alpha[focal == comp & competitor == comp & treatment == 2])
  aij_d = with(draw_pairs, alpha[focal == foc & competitor == comp & treatment == 2])
  n_ratio_d = log10((ni_d-1)/(nj_d-1))
  a_ratio_d = log10(ajj_d/aij_d)
  
  nj_w = mean(with(final_output, ni[focal == comp & treatment == 1]), na.rm = T)
  ni_w = mean(with(final_output, ni[focal == foc & treatment == 1]), na.rm = T)
  ajj_w = with(draw_pairs, alpha[focal == comp & competitor == comp & treatment == 1])
  aij_w = with(draw_pairs, alpha[focal == foc & competitor == comp & treatment == 1])
  n_ratio_w = log10((ni_w-1)/(nj_w-1))
  a_ratio_w = log10(ajj_w/aij_w)
  
  nc=abs(n_ratio_w - n_ratio_d)
  ac=abs(a_ratio_w - a_ratio_d)
  return(c(nc, ac))
}
igr_change("HOMU", "URLI")

#Create data frame
spp_list = sort(na.omit( unique(seed_data$focal)))
comp_labels = sort( na.omit( unique(seed_data$background) ))
spp_comp_combos = expand.grid(species = spp_list, competitor = comp_labels)
spp_comp_combos$n_change = 0
spp_comp_combos$a_change = 0


for(i in 1:nrow(spp_comp_combos)) {
  ii = spp_comp_combos[i, "species"] %>% unlist
  jj = spp_comp_combos[i, "competitor"] %>% unlist
  spp_comp_combos[i, "n_change"] = igr_change(ii, jj)[1]
  spp_comp_combos[i, "a_change"] = igr_change(ii, jj)[2]
}

spp_comp_combos$larger = ifelse(abs(spp_comp_combos$a_change)> abs(spp_comp_combos$n_change), "a", 'n')

spp_comp_combos = spp_comp_combos %>%
  filter(a_change != n_change)
total_pairs = c("ACWR_FEMI", "ACWR_HOMU", "ACWR_PLER", "SACO_ACWR", "URLI_ACWR", 
                "HOMU_FEMI", "PLER_FEMI", "SACO_FEMI", "URLI_FEMI", "PLER_HOMU",
                "SACO_HOMU", "URLI_HOMU","SACO_PLER", "URLI_PLER", "URLI_SACO")
spp_comp_combos$sp_pair = paste(spp_comp_combos$species, spp_comp_combos$competitor, sep = "_")

spp_comp_combos =spp_comp_combos %>%
  filter(sp_pair %in% total_pairs)

t.test(spp_comp_combos$n_change, spp_comp_combos$a_change, paired = T)

pairs_coexist_change =c("ACWR_FEMI", "ACWR_HOMU", "ACWR_PLER", "SACO_ACWR", 
                        "URLI_ACWR",  "HOMU_FEMI", "PLER_FEMI", "URLI_FEMI",
                        "SACO_PLER", "URLI_SACO")

spp_comp_combos_coexist =spp_comp_combos %>%
  filter(sp_pair %in% pairs_coexist_change)

t.test(spp_comp_combos_coexist$n_change, spp_comp_combos_coexist$a_change, paired = T)


spp_comp_combos_coexist_pivot = pivot_longer(spp_comp_combos_coexist, cols = c(n_change, a_change), names_to = "type_change", values_to ="change_value")

t.test(spp_comp_combos$n_change, spp_comp_combos$a_change, paired = T)
spp_comp_combos_pivot = pivot_longer(spp_comp_combos, cols = c(n_change, a_change), names_to = "type_change", values_to ="change_value")


ggplot(spp_comp_combos_pivot, aes(x = type_change, y =change_value))+
  theme_classic(base_size = 20) +
  theme(text = element_text(size = 18))+
  geom_boxplot(fill = "light grey", outlier.shape = NA) +
  geom_point(shape = 1, size = 4) +
  ylab("Difference between treatments") +
  xlab (" ") +
  #scale_y_log10() +
  scale_x_discrete(labels=c("a_change"="Competition \n coefficients", "n_change" = "Demographic \n potential")) 

ggsave("./figures/alpha_eta_ratio_box_model5.pdf")
write.csv(spp_comp_combos_pivot,  "./output/n_alph_ratio_output_mod5.csv")


#######  Analysis for Figure 4 in text (trait differences comparisons)

gkt = read.csv("./data/traits_gk.csv")
rownames(gkt) = gkt$species

draw_pairs = read.csv("./output/spp_treat_draw_combos_1000_model5.csv") # or run model script first to get this
final_output_draws = draw_pairs
final_output_draws$fd_superior = ifelse(final_output_draws$fd < 1, 1/final_output_draws$fd, final_output_draws$fd)
gkt_pca = prcomp(gkt[,2:12], scale = T)

gkt_pca$x
gkt_pca$x["ACWR",1]
screeplot(gkt_pca)

#Create PC vectors with the six species in the experiment
pc1= c(gkt_pca$x["ACWR",1], gkt_pca$x["FEMI",1], gkt_pca$x["HOMU",1], gkt_pca$x["PLER",1], gkt_pca$x["SACO",1], gkt_pca$x["URLI",1])
pc2= c(gkt_pca$x["ACWR",2], gkt_pca$x["FEMI",2], gkt_pca$x["HOMU",2], gkt_pca$x["PLER",2], gkt_pca$x["SACO",2], gkt_pca$x["URLI",2])
#Dist matrices for fd and snd
fd_d  = matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))

# The following loop propagates error in FD and ND measurements into the distance matrices used to 
# assess differences in ND and FD among species
dist_fd_d = list()
for( i in 1:100){
  fd_d  = matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))
  for(sp1 in rownames(fd_d)) {  
    for(sp2 in colnames(fd_d)) {    
      current_fd = final_output_draws %>% 
        filter(id == unique(final_output_draws$id)[i] & treatment == 2 & focal == sp1 & competitor == sp2) %>% 
        dplyr::select(fd_superior) %>%
        unlist()
      fd_d[sp1, sp2] = current_fd  
      dist_fd_d[[i]]= as.dist(fd_d)
    }
  }
  print(i)
}

dist_pc1= dist(pc1)

#run a mantel test over all posterior samples to generate distribution of 
# r and p statistics.
mantelr = c()
mantelp = c()
for(i in 1:length(dist_fd_d)){
  test= vegan :: mantel(xdis = dist_pc1, ydis = dist_fd_d[[i]])
  mantelr[i] = test$statistic
  mantelp[i] = test$signif
  print(i)
}

#plot histograms and 95% hd intervals for each metric
hist(mantelr^2); median(mantelr^2); HDInterval::hdi(mantelr^2)
hist(mantelp); median(mantelp); HDInterval::hdi(mantelp)

#check correspondencee between mantel test stats from a single median estimate and the 
# posterior distribution of r calculated above
medians = final_output_draws %>% group_by(focal, competitor, treatment) %>% summarise(fd_median = median(fd_superior))
fd_d  = matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))
for(sp1 in rownames(fd_d)) {  
  for(sp2 in colnames(fd_d)) {    
    current_fd = medians %>% 
      filter(treatment == 2 & focal == sp1 & competitor == sp2) %>% 
      dplyr::select(fd_median)
    fd_d[sp1, sp2] = current_fd$fd_median  }}
dist_fd_d= as.dist(fd_d)

vegan :: mantel(xdis = dist_pc1, ydis = dist_fd_d)
median(mantelr); median(mantelp) #these values should be somewhat close to those from the previous line


#Do traits dissimilarities between pairs correlate with magnitude that their snd and fd change between treatments
diffs = final_output_draws

fd_diff_abs = abs(diffs$fd[diffs$treatment == 1] - diffs$fd[diffs$treatment == 2])

diffs = diffs %>% 
  dplyr :: select(id, focal, competitor, treatment) %>%
  filter(treatment == 1) 

diffs$fd_diff_abs = fd_diff_abs

dist_fd_diffs_abs = list()
for( i in 1:100){
  fd_diffs_abs  = matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))
  for(sp1 in rownames(fd_diffs_abs)) {  
    for(sp2 in colnames(fd_diffs_abs)) {    
      fd_change = diffs %>% 
        filter(id == unique(diffs$id)[i] & focal == sp1 & competitor == sp2) %>% 
        dplyr::select(fd_diff_abs) %>%
        unlist()
      fd_diffs_abs[sp1, sp2] = fd_change  }}
  dist_fd_diffs_abs[[i]]= dist(fd_diffs_abs)
  print(i)
}

for(i in 1:length(dist_fd_diffs_abs)){
  test = vegan :: mantel(ydis = dist_fd_diffs_abs[[i]],xdis = dist_pc1)
  mantelr[i] = test$statistic
  mantelp[i] = test$signif
  print(i)
}

hist(mantelr^2); median(mantelr^2); HDInterval::hdi(mantelr^2)
hist(mantelp); median(mantelp); HDInterval::hdi(mantelp)


medians = final_output_draws %>% group_by(focal, competitor, treatment) %>% summarise(fd_median = median(fd))
diffs = as.data.frame(medians)
fd_diff_abs = abs(diffs$fd_median[diffs$treatment == 1] - diffs$fd_median[diffs$treatment == 2])

diffs = diffs %>% 
  dplyr :: select(focal, competitor, treatment) %>%
  filter(treatment == 1) 

diffs$fd_diff_abs = fd_diff_abs


fd_diffs_abs_all  = matrix(NA, nrow = 6, ncol = 6, dimnames = list(c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI"), c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")))
for(sp1 in rownames(fd_diffs_abs_all)) {  
  for(sp2 in colnames(fd_diffs_abs_all)) {    
    fd_change = diffs %>% 
      filter(focal == sp1 & competitor == sp2)  %>%
      dplyr::select(fd_diff_abs) %>%
      unlist()
    fd_diffs_abs_all[sp1, sp2] = fd_change  }}
dist_fd_diffs_abs_all= dist(fd_diffs_abs_all)

vegan ::mantel(xdis = dist_pc1,ydis = dist_fd_diffs_abs_all) 
median(mantelr); median(mantelp) #these values should be somewhat close to those from the previous line

pc1_vec = as.vector(dist_pc1)
fd_diffs_abs_vec = as.vector(dist_fd_diffs_abs_all)
fd_diffs_vec = as.vector(dist_fd_d)

fd_d_df = as.data.frame(gdata :: unmatrix(fd_d)) 
colnames(fd_d_df) = "fd_d"
fd_d_df = tibble::rownames_to_column(fd_d_df, "sp_pair")

fd_diffs_abs_df = as.data.frame(gdata :: unmatrix(fd_diffs_abs))
colnames(fd_diffs_abs_df) = "fd_diff_abs"
fd_diffs_abs_df = tibble::rownames_to_column(fd_diffs_abs_df, "sp_pair")

pc1_mat = as.matrix(dist(pc1, upper = T))
rownames(pc1_mat)= c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")
colnames(pc1_mat)= c("ACWR", "FEMI", "HOMU", "PLER", "SACO", "URLI")
pc1_df = as.data.frame(gdata::unmatrix(pc1_mat))
colnames(pc1_df) = "pc1"
pc1_df= tibble::rownames_to_column(pc1_df, "sp_pair")


df_compile = left_join(fd_d_df,fd_diffs_abs_df, by = "sp_pair") %>%
  left_join(., pc1_df, by = "sp_pair")

w_pairs= c("SACO:ACWR", "URLI:ACWR", "ACWR:FEMI", "HOMU:FEMI", 
           "PLER:FEMI", "SACO:FEMI", "URLI:FEMI","ACWR:HOMU", 
           "PLER:HOMU", "SACO:HOMU", "URLI:HOMU", "ACWR:PLER", 
           "SACO:PLER", "URLI:PLER","URLI:SACO")

df_compile= df_compile %>%
  filter(sp_pair %in% w_pairs)
df_compile$label = paste0(substr(df_compile$sp_pair, 1, 2), "-", substr(df_compile$sp_pair, 6, 7))

write_csv(df_compile, "./output/pca_mantel_mod5.csv")

vegan :: mantel(dist_pc1, dist_fd_d)
ggplot(df_compile, aes(x = pc1, y = fd_d)) + geom_point() +geom_smooth(method = "lm")

vegan :: mantel(dist_pc1, dist_fd_diffs_abs_all)
ggplot(df_compile, aes(x = pc1, y = fd_diff_abs)) + geom_point() +geom_smooth(method = "lm")

