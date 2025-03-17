fit_save_sp = function(x) {
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
  fit$save_object(file = paste0('fit_results/plot_ages1_35_top50', '/posterior/',
                                unique(growD$f_p),'_', species, '.RDS'))
}

fit_save_plot = function(x) {
  #x = 2
  plot = plots[x]
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
  fit$save_object(file = paste0('fit_results/plot_ages1_35_top50', '/posterior/',
                                 plot, '.RDS'))
}
