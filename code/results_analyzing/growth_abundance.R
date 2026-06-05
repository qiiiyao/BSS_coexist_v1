# Load packages
rm(list = ls())

require(stringr)
require(lme4)
require(lmerTest)
require(performance)
require(dplyr)
require(data.table)
setwd("D:/R projects/BSS_coexist")

load("code/data preparation/transformed data/fit_fp_top50_ages1_35_equal_interval.RData")
load("code/data preparation/transformed data/re_cover_ab_f1_ages1_35_fp.rdata")
load("code/data preparation/transformed data/sp_racover_f1_mean_fp_early_suc_1_35.rdata")
sum(names(fit_fp_top50_ages1_35_equal_interval) == names(re_cover_ab_f1_ages1_35_fp))

growth_abundance_l = lapply(1:length(re_cover_ab_f1_ages1_35_fp), function(x){
  #x = 1
  
  fit_pre = fit_fp_top50_ages1_35_equal_interval[[x]]
  re_cover_ab = re_cover_ab_f1_ages1_35_fp[[x]]
  dat = data.frame()
  f_p_1 = fit_pre$stage_sp_name[[1]]
  racover = sp_racover_f1_mean_fp_early_suc_1_35 %>% filter(f_p == f_p_1)
  
  for(i in 1:length(fit_pre[[3]])){
    growth = fit_pre[[3]][[i]]$G
    sp = gsub('1','',names(growth)[1])
    abundaces = mean(unlist((re_cover_ab %>% select(all_of(sp)))),
                     na.rm = T)
    ra = racover %>% filter(species == sp) %>% pull(ra_m_real_t)
    dat = rbind(dat, data.frame(f_p = f_p_1, 
                     sp = sp,
                     mean_growth = mean(growth),
                     mean_abundaces = abundaces,
                     ra = ra))
    
  }
  return(dat)
})

dat_growth_abundance = as_tibble(rbindlist(growth_abundance_l))

dat_growth_abundance$lg_mean_growth = log(dat_growth_abundance$mean_growth)

plot(dat_growth_abundance$mean_abundaces, dat_growth_abundance$mean_growth)
plot(dat_growth_abundance$mean_abundaces, dat_growth_abundance$lg_mean_growth)
plot(dat_growth_abundance$ra, dat_growth_abundance$lg_mean_growth)

mod1 = lm(lg_mean_growth ~ ra, 
            data = dat_growth_abundance)
summary(mod1)

mod2 = lmer(lg_mean_growth ~ ra + (1|f_p), 
          data = dat_growth_abundance)
summary(mod2)

r2(mod2)
anova(mod2) 

sps = sort(unique(dat_growth_abundance$sp))
for (i in 1:length(sps)){
  i = 16
  dat_growth_abundance_sp1 = dat_growth_abundance %>% filter(sp == sps[i])
  mod_sp = lm(lg_mean_growth ~ mean_abundaces, 
              data = dat_growth_abundance_sp1)
  summary(mod_sp)
  
}



