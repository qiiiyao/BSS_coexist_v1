###pairwise interaction effects
###BSSzmat.dat 1/0 data of 25 native and 25 exotic species
###BSScover.dat cover of each species
library("coda")
library("dplyr")
rm(list= ls())
Fieldnum=c(rep(1,48),rep(2,48),rep(3,48),rep(4,48),rep(5,48),rep(6,48),rep(7,48),rep(8,48),rep(9,48),rep(10,48))
setwd("D:/R projects/BSS")
load("code/Yin_ele_2022/fit_yin_ele.rdata")
load("D:/R projects/BSS/code/data preparation/transformed data/fit_fp_same_ages_top50.RData")
sp_list = unique(sp_cover_f2_top50$Species)

##### fit_r_limit ##### 
library(R2jags)
names(fit_yin_ele) = c('BSSzmat.dat', 'BSScover.dat', 'BSSothercover.dat', 'BSSotherrich.dat',
                       'growth.list')
BSScover.dat = fit_yin_ele$BSScover.dat
BSSothercover.dat = fit_yin_ele$BSSothercover.dat
growth.list = fit_yin_ele$growth.list

`pairgrowth.fn_r_limit` =
  function(){
    # this script fits the model to the BSS invading species by evaluate the invader-invader effects
    growth.dat=growth.dat
    nsite=dim(growth.dat)[1]
    nyear=dim(growth.dat)[2]
    focalcover=focalcover
    extcover=extcover
    othercover = othercover
    Fieldnum=Fieldnum
    
    sim.dat = list()
    sim.dat$growth.dat=growth.dat
    sim.dat$nsite=nsite
    sim.dat$nyear=nyear
    sim.dat$nsp = nsp
    sim.dat$focalcover=focalcover
    sim.dat$extcover=extcover
    sim.dat$othercover=othercover
    sim.dat$Fieldnum=Fieldnum
    
    sim.par=c("a", "b_i", "b_j", "b_o", "r", "c") 
    
    sim.fit=jags(sim.dat, model.file='code/fit/yin_bugs_merge_ralimit.txt',
                 parameters.to.save=sim.par, n.chains=1, n.iter=10000)
    
    return(sim.fit)
  }

intra_i = matrix(NA,50,50)
inter_j = matrix(NA,50,50)
inter_oj = matrix(NA,50,50)
r = matrix(NA,50,50)
sd_intra_i = matrix(NA,50,50)
sd_inter_j = matrix(NA,50,50)
sd_inter_oj = matrix(NA,50,50)
sd_r = matrix(NA,50,50)

summary = list()
for (i in 1:50){
  #i = 1
  growth.dat=growth.list[[i]]
  growth.dat = array(data = unlist(replicate(49, growth.dat, simplify = FALSE)), dim = c(480, 42, 49))
  focalcover=as.matrix(scale(BSScover.dat[[i]],center=TRUE))
  focalcover[is.na(focalcover)]=0
  
  focalcover_1 = as.matrix(BSScover.dat[[i]])
  othercover_1 = BSSothercover.dat[c(1:50)[-i]]
  #for (j in 1:50){
  #j = 2
  #if(i == j)(print(i)) else {
  extcover = BSScover.dat[(c(1:50)[-i])]
  extcover=lapply(extcover, function(x){y = as.matrix(scale(x,center=TRUE))
  y[is.na(y)] = 0
  return(y)})
  extcover = array(unlist(extcover), dim = c(480, 42, 49))
  extcover_1 = BSScover.dat[(c(1:50)[-i])]
  othercover = Map(function(x, y){z = x - y
  z1 = as.matrix(scale(z, center=TRUE))
  z1[is.na(z1)] = 0
  return(z1)}, othercover_1, extcover_1)
  othercover = array(unlist(othercover), dim = c(480, 42, 49))
  nsite=dim(growth.dat)[1]
  nyear=dim(growth.dat)[2]
  nsp = dim(growth.dat)[3]
  Fieldnum = Fieldnum
  
  res_r_limit=pairgrowth.fn_r_limit()
  summary[[i]] = res_r_limit$BUGSoutput$summary
  save(summary[[i]], file = paste0('results/fit_results/yin_ele/merge_r_a_limit/summary/',
                                   sp_list[i], '.rdata'))
  #}   
}

data = list(n=n,
            nrep=nrep,
            N=N,
            time=time,
            sp=sp_com,
            sr=sr)

#### model fitting ####
### stan options
chains = 4 # 4 chains in parallel
sampling = 3500
warmup = 1000
thin = 1

### initial values for sampling
init = rep(list(list(a=matrix(rep(0.02, sr*sr), sr, sr),
                     sdev=2,
                     N0=data$N[1, , ])
),chains)

### run model and print result
fit = mod$sample(data=data, seed = 123,
                 chains = chains, 
                 parallel_chains = chains,
                 iter_warmup = warmup,
                 iter_sampling = sampling,
                 init=init,
                 refresh = 10, # print update every 10 iters
                 max_treedepth = 10,
                 adapt_delta = 0.8)

#save(intra_i,inter_j,inter_oj,r,sd_intra_i,sd_inter_j,sd_inter_oj,sd_r,n.eff, file = "D:/BSS study/code/results_fit_xi_nolimit.rdata")
