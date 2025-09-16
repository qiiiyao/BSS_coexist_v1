rm(list=ls())

# set up R environment
library(cmdstanr)
options(mc.cores = parallel::detectCores(logical = F)) 
library(dplyr)
library(parallel)
library(doParallel)
setwd("/data/home/shpli3/R_projects/BSS_exclude_tree_raw")

# Load data
load("/data/home/shpli3/my_pc/BSS/code/data preparation/transformed data/fit_fp_top50_ages1_35_equal_interval.RData")

# run model chains separately then combine then to avoid memory issues
stan.seed = 52
file = file.path('/data/home/shpli3/my_pc/BSS/code/fit/stancode/model_comprison/constrained_competition_equal_interval/Exp_r_allposi.stan')
mod = cmdstan_model(file
#, stanc_options = list("allow-undefined")
)
numer_cores = parallel::detectCores(logical = F)
chains = 4 # 4 chains in parallel
sampling = 2500
warmup = 1000
thin = 10


fit_save = function(x) {
  #x = 2
  species = sp[x]
  fit = mod$sample(data = datalist[[x]],
                   seed = stan.seed, 
                   chains = chains, 
                   parallel_chains = chains,
                   iter_warmup = warmup,
                   iter_sampling = sampling,
                   refresh = 100, # print update every 100 iters
                   thin = thin,
                   max_treedepth = 20,
                   adapt_delta = 0.99)
  
  # Raw output
  #------------ 
  fit$save_object(file = paste0('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/null', '/posterior/',
                                unique(growD$f_p),'_', species, '.RDS'))
}

#--------------------------------------------------
# Estimate interactions with a joint NDD*RI model for t1
#--------------------------------------------------
data_t = fit_fp_top50_ages1_35_equal_interval
plot_list = sapply(seq(1, 480, 48), function(x){y = c(x:(x+47))})
names(fit_fp_top50_ages1_35_equal_interval)[plot_list[,5]]

fitted_plots_1 = list.files('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/null/summary')
m = regexec("[0-9]+\\_[0-9]+", fitted_plots_1)
fitted_plots = unlist(regmatches(fitted_plots_1, m))
unfitted = setdiff(names(data_t)[plot_list[,5]], fitted_plots)
unfitted_order = which(names(data_t) %in% unfitted)

for (i in unfitted_order) {
  #i = 1
  datalist = data_t[[i]][[3]]
  growD = data_t[[i]][[2]]
  unique(growD$f_p)
  sp = colnames(growD)[10:ncol(growD)]
  
  {cl = makeCluster(numer_cores-2)      
    registerDoParallel(cl)       
    foreach(x=1:length(datalist), .packages = c('dplyr', 'cmdstanr')) %dopar% {fit_save(x)
      }

      # Save results 
      summary_l = foreach(i=1:length(datalist)) %dopar% {
        #i = 2
      species = sp[i]
      fit = readRDS(paste0('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/null', '/posterior/',
                     unique(growD$f_p),'_', species, '.RDS'))
      summary = as.data.frame(fit$summary())}

      save(summary_l, file=paste0('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/null',
                                  '/summary/','summary_',
                                  unique(growD$f_p), '.rdata'))

  # Extract results
   r_trans_sp = vector()
   r_trans_sp = foreach(i=1:length(datalist), .combine = c, .packages = c("dplyr")) %dopar% {(summary_l[[i]] %>% filter(variable == 'r'))$mean}
   names(r_trans_sp) = sp

   save(r_trans_sp, file=paste0('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/null', '/parameters/','r_',
                               unique(growD$f_p), '.rdata'))
   stopCluster(cl)
   }
}