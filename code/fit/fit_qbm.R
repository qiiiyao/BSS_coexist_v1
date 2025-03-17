## Script to estimate parameters for the quad-based model using STAN

rm(list=ls())

library('cmdstanr')
library("coda")
library("deSolve")
library("bayesplot")
library("parallel")
library("posterior")
library("dplyr")

# Compile outside of loop
load('D:/BSS study/code/invasion_stage.Rdata')
growD = cover_prop_domin_l[[1]]
growD[growD==0] = 1e-4
Y=growD[2:(nrow(growD)),1:8]
Yrs = nrow(Y)
yid = as.numeric(as.factor(Y$Age))
sp_matrix = log(growD[1:(nrow(growD)-1),8:9])
S = ncol(sp_matrix)
datalist = list(N = Yrs, Yrs=Yrs, yid=yid, S = S,
                 Y=Y$Daucus_carota, X=sp_matrix)

mod_QBM = cmdstan_model(stan_file = 'D:/BSS study/code/QBM_2.stan')

inits <- list()

inits[[1]] <- list(a_mu=0, a=rep(0,Yrs), b_mu=rep(0.01,S), b=rep(0.01,S),
                   sig_b1=rep(0.5,S), sig_a=0.5, sigmaSq=0.5)
inits[[2]] <- list(a_mu=1, a=rep(1,Yrs), b_mu=rep(1,S), b=rep(1,S),
                   sig_b=rep(1,S), sig_a=1, sigmaSq=1)
inits[[3]] <- list(a_mu=0.5, a=rep(0.5,Yrs), b_mu=rep(0.5,S), b=rep(0.5,S),
                    sig_b=rep(0.1,S), sig_a=0.1, sigmaSq=0.1)

inits[[1]] <- list(a_mu=0, a=rep(-5,Yrs), b_mu=rep(0.01,S), b=rep(0.01,S),
                   sig_b1=rep(0.5,S), sig_a=0.5, sigmaSq=0.5)
inits[[2]] <- list(a_mu=1, a=rep(10,Yrs), b_mu=rep(10,S), b=rep(1,S),
                   sig_b=rep(1,S), sig_a=1, sigmaSq=1)
inits[[3]] <- list(a_mu=0.5, a=rep(100,Yrs), b_mu=rep(10000,S), b=rep(0.5,S),
                   sig_b=rep(0.1,S), sig_a=0.1, sigmaSq=0.1)

PrelimFit_l = mod_QBM$sample(data = datalist, seed = 123, 
                                     chains = 4, 
                                     parallel_chains = 4,
                                     iter_warmup = 1000,
                                     iter_sampling = 2000,
                                     refresh = 10, # print update every 10 iters
                                     max_treedepth = 15,
                                     thin = 10,
                                     #init = inits,
                                     adapt_delta = 0.99)

summary = PrelimFit_l$summary()
PrelimFit_l$diagnostic_summary()
posterior = PrelimFit_l$draws(format = 'matrix')

plot_title = ggtitle("Posterior distributions",
                     "with medians and 80% intervals")
mcmc_areas(posterior, pars = c("a_mu", 'b[1]'),
           prob = 0.8) + plot_title

color_scheme_set("mix-blue-pink")
p = mcmc_trace(posterior,  pars = c("a_mu", 'b[1]'), n_warmup = 1000,
               facet_args = list(nrow = 2, labeller = label_parsed))
p + facet_text(size = 15)
               
pars=c("a_mu", "a", "b1_mu",  "b1",
       "tau", "gint", "sig_a", "sig_b1", "sig_G")
mcmc_samples <- stan(file = "qbm_noclimate.stan", data=datalist, pars=pars, chains=0)

##  Loop through species and fit the model
for (do_species in sppList){
  print(paste("fitting model for", do_species))
  
  growD <- subset(growD_all, species==do_species)
  groups <- as.numeric(as.factor(growD$group))
  G <- length(unique(growD$group))
  Yrs <- length(unique(growD$year))
  yid <- as.numeric(as.factor(growD$year))
  
  ## Set reasonable initial values for three chains
  inits <- list()
  inits[[1]] <- list(a_mu=0, a=rep(0,Yrs), b1_mu=0.01, b1=rep(0.01,Yrs),
                     gint=rep(0,G), w=c(0,0), sig_b1=0.5, sig_a=0.5, sigmaSq=0.5,
                     sig_G=0.5)
  inits[[2]] <- list(a_mu=1, a=rep(1,Yrs), b1_mu=1, b1=rep(1,Yrs),
                     gint=rep(1,G), w=c(0.5,0.5), sig_b1=1, sig_a=1, sigmaSq=1,
                     sig_G=1)
  inits[[3]] <- list(a_mu=0.5, a=rep(0.5,Yrs), b1_mu=0.5, b1=rep(0.5,Yrs),
                     gint=rep(0.5,G), w=c(-0.5,-0.5), sig_b1=0.1, sig_a=0.1, sigmaSq=0.1,
                     sig_G=0.1)
  
  datalist <- list(N=nrow(growD), Yrs=Yrs, yid=yid,
                   Y=growD$propCover.t1, X=log(growD$propCover.t0),
                   G=G, gid=groups)
  pars=c("a_mu", "a", "b1_mu",  "b1",
         "tau", "gint", "sig_a", "sig_b1", "sig_G")
  
  rng_seed <- 123
  sflist <-
    mclapply(1:3, mc.cores=3,
             function(i) stan(fit=mcmc_samples, data=datalist, pars=pars,
                              seed=rng_seed, chains=1, chain_id=i, refresh=-1,
                              iter=2000, warmup=1000, init=list(inits[[i]])))
  fit <- sflist2stanfit(sflist)
  long <- ggs(fit)
  saveRDS(long, paste("popgrowth_stanmcmc_noclimate_", do_species, ".RDS", sep=""))
  r_hats <- summary(fit)$summary[,10] 
  write.csv(r_hats, paste("rhat_noclimate_", do_species, ".csv", sep=""))
}

