#### fit ####
library('cmdstanr')
library('dplyr')

###########################################################
##### fitting Ricker competition models with dom data #####

rm(list=ls())
setwd("D:/R projects/BSS")
file = file.path("code/fit/stancode/Ricker_field_g.stan")
mod = cmdstan_model(file)

#### data preparation ####
load("code/data preparation/transformed data/fit_ranfield_ages1_35_top50.rdata")
load("code/data preparation/transformed data/sp_cover_f1_ages1_35_top50.rdata")

BSScover.dat = array(unlist(fit_ranfield_ages1_35_top50$BSScover.dat), dim = c(480, 29, 50))
BSScover.dat = aperm(BSScover.dat, c(2, 3, 1))
growth.list = array(unlist(fit_ranfield_ages1_35_top50$growth.list), dim = c(480, 28, 50))
growth.list = aperm(growth.list, c(2, 3, 1))

dom = sp_cover_f1_ages1_35_top50
sp_bio = colnames(dom)[10:ncol(dom)]
sp_pa = sapply(sp_bio, function(x){y = gsub('\\_', '.', x)})
sr = length(sp_bio)
nplot = length(unique(dom$f_p))
Fieldnum=c(rep(1,48),rep(2,48),rep(3,48),rep(4,48),rep(5,48),rep(6,48),rep(7,48),rep(8,48),rep(9,48),rep(10,48))
fake_age = sort(unique(dom$fake_age))
time = fake_age[-length(fake_age)]
n = length(time)
sp_com = array(100, dim=c(nplot, length(sp_pa)))
BSScover.dat[BSScover.dat==0] = 1e-6 # change zeros into a small number
N_noend = BSScover.dat[-length(fake_age),,] # data: dim1=time, dim2=species, dim3=PLOT
N = BSScover.dat # data: dim1=time, dim2=species, dim3=PLOT

#### model fitting ####
### stan options
chains = 4 # 4 chains in parallel
sampling = 2500
warmup = 1000
thin = 10

### initial values for sampling
init = rep(list(list(r = 0.1,
                     a=c(rep(0.02, sr)),
                     sigma_e=2,
                     sigma_field = 2)),chains)
sp_list = sapply(seq(1, 50, 5), function(x){y = c(x:(x+4))})

for (i in sp_list[,1]) {
  #i = 1
  focal_sp = sp_bio[i]
  g = growth.list[,i,]
  data = list(n=n,
              nplot=nplot,
              Fieldnum = Fieldnum,
              g = g,
              N = N_noend,
              time=time)
  
  ### run model and print result
  fit = mod$sample(data=data, seed = 123,
                   chains = chains, 
                   parallel_chains = chains,
                   iter_warmup = warmup,
                   iter_sampling = sampling,
                   init=init,
                   thin = thin,
                   refresh = 10, # print update every 10 iters
                   max_treedepth = 20,
                   adapt_delta = 0.99)
  
  saveRDS(fit, file = paste0('/data/home/shpli3/R_projects/BSS_exclude_tree_raw/code/fit/fit_random_field_g/fit_',
                             focal_sp, '.RDS'))
}
