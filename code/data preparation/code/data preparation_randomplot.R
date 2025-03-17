###spatio-temporal data
###data sorting
rm(list = ls())
library(dplyr)
library(lme4)

##### Data loading
sp_cover = read.csv('data/original data/BSS_community_332.csv',
                    header = T)
field = read.csv('data/original data/FIELDS.csv',
                 header = T)
trait = read.csv('data/original data/traits332.csv',
                 header = T)

sp_cover$f_p = paste(sp_cover$Field, sp_cover$Plot, sep = '_')
sp_cover = sp_cover %>% relocate(f_p, .after = Plot)
trait$species_id = paste('sp', trait$species_id, sep = '_')
sp_real = vector()
sp = colnames(sp_cover)[8:ncol(sp_cover)]

for (i in 1:length(sp)) {
  sp_real[i] = trait[trait$species_id == sp[i],]$Species
}

colnames(sp_cover)[8:ncol(sp_cover)] = sp_real # change the species name of data

sp_cover_fp = split(sp_cover, sp_cover$f_p)
sp_cover_year = split(sp_cover, sp_cover$Absolute_year)
sp_cover_age = split(sp_cover, sp_cover$Age)

sapply(sp_cover_fp, function(x) {sort(unique(x$Absolute_year))})
sapply(sp_cover_year, function(x) {sort(unique(x$f_p))})

##### Follow Li et al., (2015), condense the two year (half of field survey) to
##### one year and Follow Yin et al., (2022),, our study restricted to 3-51 age
sp_cover_f1 = sp_cover %>% filter(Absolute_year %in% c(1958:2017))
sp_cover_f1$fake_year = sp_cover_f1$Absolute_year
sp_cover_f1$fake_age = sp_cover_f1$Age
full_f_p = unique(sp_cover_f1$f_p)
for (i in seq(1959, 2017, 2)){
  #i = 1991
  f_p_year_2 = unique(sp_cover_f1[sp_cover_f1$Absolute_year == i-1,]$f_p)
  diff_f_p_year_2 = setdiff(full_f_p, f_p_year_2)
  row_diff = nrow(sp_cover_f1[sp_cover_f1$Absolute_year == i & 
                                sp_cover_f1$f_p %in% diff_f_p_year_2,
  ])
  if (row_diff == 0) {print(paste('cannot make it up', 'for',
                                  i-1, sep = ' '))}
  else {
    sp_cover_f1[sp_cover_f1$Absolute_year == i & sp_cover_f1$f_p %in% 
                  diff_f_p_year_2,]$fake_year = i-1
    age = unique(sp_cover_f1[sp_cover_f1$Absolute_year == i & sp_cover_f1$f_p %in% 
                               diff_f_p_year_2,]$Age)
    sp_cover_f1[sp_cover_f1$Absolute_year == i & sp_cover_f1$f_p %in%
                  diff_f_p_year_2,]$fake_age = rep(age-1,
                                                   each = row_diff/length(age-1))
  }
}

sp_cover_f2 = sp_cover_f1 %>% filter(fake_age %in% seq(3, 52, 2))
sp_cover_f2 = sp_cover_f2 %>% relocate(fake_year, .before = Age)
sp_cover_f2 = sp_cover_f2 %>% relocate(fake_age, .before = Field)
sp_cover_f2_y = split(sp_cover_f2, sp_cover_f2$fake_age)
sapply(sp_cover_f2_y, function(x) {sort(unique(x$f_p))})
unique(sp_cover_f2$fake_age)
sp_cover_f2_f = split(sp_cover_f2, sp_cover_f2$Field)

######## filter top 50 species for 10 fields
pchange.all2_l = list()
for (z in 1:length(sp_cover_f2_f)) {
  #z = 1
sp_cover_f2_f1 = sp_cover_f2_f[[z]]
sp_cover_only = sp_cover_f2_f1[,10:ncol(sp_cover_f2_f1)]

top50 = rev(sort(apply(sp_cover_only,
                       2, mean,na.rm=T)))[1:50] 
top50.short = names(top50)

sp_cover_focal = sp_cover_only %>% select(all_of(top50.short))

Rest = apply(sp_cover_only[, -match(names(top50),
                                    names(sp_cover_only))],
             1, sum,na.rm=T)
sp_cover_focal_all = cbind(sp_cover_focal,
                           Rest)
#sp_cover_focal_all = sp_cover_focal_all/apply(sp_cover_focal_all,
#                     1,sum) # to rscale to no more than 100%
sp_cover_focal_all[sp_cover_focal_all == 0] = 0.01
#plants3 = log(plants3)
pyear = split(sp_cover_focal_all, sp_cover_f2_f1$fake_age)
pchange = list()
pchange.logit = list()
pchange.log = list()
pchange.relative = list()

# This is to prepare the database in order to conduct the analyses where t+1 will be compared to t
for(i in 1:(length(pyear)-1)){
  xx = log(pyear[[i+1]]/pyear[[i]])
  names(xx) = paste(names(xx), "_delta",sep="")
  pchange[[i]] = cbind(xx, pyear[[i]])
}

pchange.all = do.call("rbind", pchange)
Year_change = paste(1:24, 2:25, sep="to")
pchange.all2_l[[z]] = data.frame("f_p" = rep(unique(sp_cover_f2_f1$f_p), 24),
                          "Field" = rep(unlist(lapply(strsplit(unique(sp_cover_f2_f1$f_p),
                                                               '_'),
                                                      function(x){
                                                        x[1]
                                                      })), 24),
                          "Plot" = rep(unlist(lapply(strsplit(unique(sp_cover_f2_f1$f_p),
                                                              '_'),
                                                     function(x){
                                                       x[2]
                                                     })), 24),
                          "Year_change" = rep(Year_change, each = 48),
                          "Yeart" = rep(1:24, each = 48), pchange.all) 

}

yy = names(pchange.all2_l[[z]])[grep("delta", names(pchange.all2))]

##### Let's plot a couple of examples to see how the data looks like 

plot(pchange.all2$Poa_pratensis, pchange.all2$Poa_pratensis_delta,
     xlab = "Poa_pratensis_year_t",
     ylab = "growth(Poa_pratensis_year_t)")

##### Run a brms model and extract its fitting data for formal fit
library(rstan)
library(brms)
rstan_options(auto_write = TRUE)
rstan_options(threads_per_chain = 1)

dd = brm(as.formula(paste(yy[1], " ~ 0 + Intercept + ", 
                          paste(top50.short,
                                collapse="+"), 
                          '+ (1|Plot)',
                          '+ar(time = Yeart,
                          gr = Plot,
                          p = 1)')),
         data = pchange.all2_l[[z]], warmup = 1000, iter = 3000, 
         chains = 4, 
         #core = 2,
         core = parallel::detectCores(),
         seed = 1234)
summary(dd)


dd_noAR = brm(as.formula(paste(yy[1], " ~ 0 + Intercept + ", 
                          paste(top50.short,
                                collapse="+"), 
                          '+ (1|Plot)')),
         data = pchange.all2_l[[z]], warmup = 1000, iter = 3000, 
         chains = 4, 
         #core = 2,
         core = parallel::detectCores(),
         seed = 1234)
summary(dd_noAR)

dd_data = standata(dd)
stancode(dd)
save.image('code/data preparation/fit_ricker_randomplot.rdata')

