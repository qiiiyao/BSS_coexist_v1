############### Fast Start ####################
rm(list = ls())
setwd("D:/R projects/BSS")

### load packages
library(dplyr)
library(betareg)
library(stringr)
library(data.table)
library(ape)
library(devEMF)
library(reticulate)


load('results/fit_results/BSS_exclude_trees_raw/plot_ages1_35_top50_equal_interval_model_comparison/dat_all_alltime.rdata')
load("results/fit_results/BSS_exclude_trees_raw/plot_ages1_35_top50_equal_interval_model_comparison/bh_partialb_top40/inter_all_c_alltime_newfd.rdata")
inter_all_c_alltime_top40 = inter_all_c_alltime

inter_all_c = inter_all_c_alltime_top40
rm(inter_all_c_alltime_top40)

## Functions for plot 
source('code/function/plot_func.R')
#legend.position = c(0.11, 0.76)

# get the abs value of lgfd
inter_all_c$ablgfd = abs(inter_all_c$lgfd)

#### Check the correlation between nd and fd
#cor.test(inter_all_c$nd, inter_all_c$ablgfd) 
## Full data but corr coeff just 0.13
## correlation
#plot(inter_all_c$nd, inter_all_c$ablgfd)

## Add stage introduce
inter_all_c$intro_i = 1
inter_all_c$intro_j = 1
inter_all_c[inter_all_c$stage_i == 'native',]$intro_i = 0
inter_all_c[inter_all_c$stage_j == 'native',]$intro_j = 0

inter_all_c$stage_ij_estab = inter_all_c$stage_ij
inter_all_c$stage_ij_domin = inter_all_c$stage_ij

inter_all_c$stage_i_estab = 'native'
inter_all_c$stage_j_estab = 'native'
inter_all_c[inter_all_c$intro_i == 1&
              inter_all_c$estab_i == 0,]$stage_i_estab = 'intro.noestab'
inter_all_c[inter_all_c$intro_i == 1&
              inter_all_c$estab_i == 1,]$stage_i_estab = 'intro.estab'
inter_all_c[inter_all_c$intro_j == 1&
              inter_all_c$estab_j == 0,]$stage_j_estab = 'intro.noestab'
inter_all_c[inter_all_c$intro_j == 1&
              inter_all_c$estab_j == 1,]$stage_j_estab = 'intro.estab'
inter_all_c$stage_ij_estab = paste(inter_all_c$stage_i_estab,
                                   inter_all_c$stage_j_estab,
                                   sep = '_')

inter_all_c$stage_i_domin = 'native'
inter_all_c$stage_j_domin = 'native'
inter_all_c[inter_all_c$intro_i == 1&
              inter_all_c$estab_i == 0&
              inter_all_c$domin_i == 0,]$stage_i_domin = 'intro.nodomin'
inter_all_c[inter_all_c$intro_i == 1&
              inter_all_c$estab_i == 1&
              inter_all_c$domin_i == 0,]$stage_i_domin = 'estab.nodomin'
inter_all_c[inter_all_c$intro_i == 1&
              inter_all_c$estab_i == 1&
              inter_all_c$domin_i == 1,]$stage_i_domin = 'estab.domin'
inter_all_c[inter_all_c$intro_j == 1&
              inter_all_c$estab_j == 0&
              inter_all_c$domin_j == 0,]$stage_j_domin = 'intro.nodomin'
inter_all_c[inter_all_c$intro_j == 1&
              inter_all_c$estab_j == 1&
              inter_all_c$domin_j == 0,]$stage_j_domin = 'estab.nodomin'
inter_all_c[inter_all_c$intro_j == 1&
              inter_all_c$estab_j == 1&
              inter_all_c$domin_j == 1,]$stage_j_domin = 'estab.domin'
inter_all_c$stage_ij_domin = paste(inter_all_c$stage_i_domin,
                                   inter_all_c$stage_j_domin,
                                   sep = '_')

inter_all_c = inter_all_c %>% relocate(intro_i, .after = stage_ij) %>% 
  relocate(intro_j, .after = intro_i) %>% 
  relocate(stage_ij_estab, .after = stage_ij) %>% 
  relocate(stage_i_estab, .after = stage_ij) %>% 
  relocate(stage_j_estab, .after = stage_i_estab) %>% 
  relocate(stage_ij_domin, .after = stage_ij_estab) %>% 
  relocate(stage_i_domin, .after = stage_ij_estab) %>% 
  relocate(stage_j_domin, .after = stage_i_domin)

### Draw the coexistence plot
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)
library(ggpubr)
library(ggeffects)
library(stringr)
library(viridis)
library(showtext)
library(extrafont)
windowsFonts(Arial=windowsFont("Arial"))
library(scatterpie)
#loadfonts(device = "win")

gre = "#55a868ff"
ora = "#dd8452ff"
yel = "#ccb974ff"
blu = "#4c72b0ff"
colors_4d = c(blu,ora,gre,yel)
colors_2d = c(ora, blu)



#### Analyze invasion success ####
inter_all_c$ra_m_real_t_i = as.numeric(inter_all_c$ra_m_real_t_i)
inter_all_c$ra_m_real_t_j = as.numeric(inter_all_c$ra_m_real_t_j)
inter_all_c_l_2 = split(inter_all_c, inter_all_c$f_p)
dat_all_2 = split(dat_all_alltime, dat_all_alltime$f_p)

##### Species level #####
dat_suc_sp = data.frame()

for (i in 1:length(inter_all_c_l_2)) {
  #i = 1
  trans_plot = inter_all_c_l_2[[i]]
  inv_suc = trans_plot %>% filter(stage_i != 'native' & stage_j == 'native')
  inv_sp = unique(inv_suc$species_i)
  inv_suc_all = dat_all_2[[i]]
  
  for (j in 1:length(inv_sp)) {
    #j = 4
    inv_suc_sp = inv_suc %>% filter(species_i == inv_sp[j])
    inv_suc_sp_all = inv_suc_all %>% filter(species_i == inv_sp[j])
    
    mnd = mean(inv_suc_sp$nd)
    mlgfd = mean(inv_suc_sp$lgfd)
    mablgfd = mean(inv_suc_sp$ablgfd)
    mpd_all = mean(inv_suc_sp_all$Phylo_dis)
    mfunc_d_all = mean(inv_suc_sp_all$Multi_traits)
    mconti_func_d_all = mean(inv_suc_sp_all$Multi_conti_traits)
    
    
    mnd.a = sum(inv_suc_sp$nd*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/
      sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mlgfd.a = sum(inv_suc_sp$lgfd*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/
      sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mablgfd.a = sum(inv_suc_sp$ablgfd*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/
      sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mpd.a_all = sum(inv_suc_sp_all$Phylo_dis*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/
      sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mfunc_d.a_all = sum(inv_suc_sp_all$Multi_traits*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/
      sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mconti_func_d.a_all = sum(inv_suc_sp_all$Multi_conti_traits*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/
      sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    
    mnnd = min(inv_suc_sp$nd, na.rm = T)
    mnlgfd = min(inv_suc_sp$lgfd, na.rm = T)
    mnablgfd = min(inv_suc_sp$ablgfd, na.rm = T)
    mntd_all = min(inv_suc_sp_all$Phylo_dis, na.rm = T)
    mnfunc_d_all = min(inv_suc_sp_all$Multi_traits, na.rm = T)
    mnconti_func_d_all = min(inv_suc_sp_all$Multi_conti_traits, na.rm = T)
    
    dat_suc_sp_1 = data.frame(f_p = unique(trans_plot$f_p),
                              plot = unique(trans_plot$f_p),
                              field = unique(trans_plot$field),
                              species = inv_sp[j], 
                              stage = unique(inv_suc_sp$stage_i),
                              estab = unique(inv_suc_sp$estab_i),
                              domin = unique(inv_suc_sp$domin_i),
                              
                              mnd = mnd, mlgfd = mlgfd, mablgfd = mablgfd,
                              mpd_all = mpd_all, mfunc_d_all = mfunc_d_all,
                              mconti_func_d_all = mconti_func_d_all,
                              
                              mnd.a = mnd.a, mlgfd.a = mlgfd.a, mablgfd.a = mablgfd.a,
                              mpd.a_all = mpd.a_all,
                              mfunc_d.a_all = mfunc_d.a_all,
                              mconti_func_d.a_all = mconti_func_d.a_all,
                              
                              mnnd = mnnd, mnlgfd = mnlgfd, mnablgfd = mnablgfd,
                              mntd_all = mntd_all,
                              mnfunc_d_all = mnfunc_d_all, 
                              mnconti_func_d_all = mnconti_func_d_all)
    dat_suc_sp = rbind(dat_suc_sp, dat_suc_sp_1)
  }
}

summary(dat_suc_sp)
dat_suc_sp$stage_level = NA
dat_suc_sp[dat_suc_sp$stage == 'introduce',]$stage_level = 1
dat_suc_sp[dat_suc_sp$stage == 'establish',]$stage_level = 2
dat_suc_sp[dat_suc_sp$stage == 'dominant',]$stage_level = 3
dat_suc_sp$stage_level = ordered(dat_suc_sp$stage_level)
str(dat_suc_sp)
save(dat_suc_sp,
     file = 'code/results_analyzing/analysing_ages1_35_top40_equal_interval_bh_partialb_data/dat_suc_sp.rdata')



#### Single regression results for analyzing invasion success probability ~ mnd+mfd+mpd+mfunc_d for all species ####
##### establishment predictive curves for md.ab #####
setwd("D:/R projects/BSS")
library(phyr)
library(tibble)
library(lme4)
require(ape)
library(scales)
library(ggthemes)
load('code/results_analyzing/analysing_ages1_35_top40_equal_interval_bh_partialb_data/dat_suc_sp.rdata')
numcols = grep("^m",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])

dat_suc_sp %>% group_by(estab) %>% 
  summarise_at(vars(mnd), c(length))

pc_prior = list(prec=list("pc.prec", param=c(0.1,0.01)))
dat_suc_sps$species_1 = as.factor(dat_suc_sps$species)
estab_sp_names = unique(dat_suc_sps$species)

tree = read.tree('data/original data/phylo_tree332.txt')
estab_tree_fit = keep.tip(tree, estab_sp_names)
estab_vcv_tree = ape::vcv(estab_tree_fit, model = "Brownian", corr = FALSE)
estab_vcv_tree_sparse = inla.as.sparse(solve(estab_vcv_tree))

dat_dom_sp = dat_suc_sp %>% filter(stage %in% c('establish', 'dominant'))
dat_dom_sps = dat_suc_sps %>% filter(stage %in% c('establish', 'dominant'))
dat_dom_sp %>% group_by(domin) %>% 
  summarise_at(vars(mnd), c(length))

### get establishment predict data for mnd.a
lincombs.data.estab.mnd.a_single = data.frame(mnd.a=seq(min(dat_suc_sp$mnd.a),
                                                        max(dat_suc_sp$mnd.a),
                                                        length=100))

lincombs.matrix.estab.mnd.a_single=model.matrix(~mnd.a,
                                                data=lincombs.data.estab.mnd.a_single)
lincombs.matrix.estab.mnd.a_single=as.data.frame(lincombs.matrix.estab.mnd.a_single)
lincombs.estab.mnd.a_single=inla.make.lincombs(lincombs.matrix.estab.mnd.a_single)

inla.model_lincombs.estab.mnd.a_single = pglmm(estab ~ mnd.a+#(1|species) + 
                                                 (1|f_p) + (1|field), data = dat_suc_sp,
                                               family = "binomial", cov_ranef = list(species = tree),
                                               bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                           config = TRUE),
                                                                    quantiles=c(0.025,0.5,0.975),
                                                                    lincomb=lincombs.estab.mnd.a_single,
                                                                    control.predictor=list(compute=T)),
                                               bayes = T)

lincombs.posterior.estab.mnd.a_single = inla.model_lincombs.estab.mnd.a_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mnd.a_single$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

# inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mnd.a_single$predicted.value=unlist(lapply(lincombs.posterior.estab.mnd.a_single,
                                                               function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mnd.a_single$lower=unlist(lapply(lincombs.posterior.estab.mnd.a_single,
                                                     function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mnd.a_single$upper=unlist(lapply(lincombs.posterior.estab.mnd.a_single,
                                                     function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mnd.a_single
save(lincombs.data.estab.mnd.a_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top40_equal_interval_bh_partialb_data/lincombs.data.estab.mnd.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mnd.a
load("code/results_analyzing/analysing_ages1_35_top40_equal_interval_bh_partialb_data/lincombs.data.estab.mnd.a_single.rdata")
(estab.mnd.a_single.logistic=ggplot(data=lincombs.data.estab.mnd.a_single,aes(x=mnd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1],alpha=0.2)+
    geom_line(color=colors_4d[1],size=2,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #scale_x_continuous(breaks = seq(-0.35, 0.7, 0.35))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnd.a, y=estab),
               color = colors_4d[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='Establishment probability')+
    annotate(geom="text",
             x=min(dat_suc_sp$mnd.a)+(max(dat_suc_sp$mnd.a)-min(dat_suc_sp$mnd.a))*0.3,
             y=c(0.900,0.80),
             label=c("italic(β)['ND'[ab]] == '7.87'",
                     "'95%CI' == '[6.13, 9.61]'"),
             parse=T,size=4)+
    theme_regular()
)
#mnd.a   7.87       6.13       9.61



### get establishment predict data for mlgfd.a
lincombs.data.estab.mlgfd.a_single = data.frame(mlgfd.a=seq(min(dat_suc_sp$mlgfd.a),
                                                            max(dat_suc_sp$mlgfd.a),
                                                            length=100))

lincombs.matrix.estab.mlgfd.a_single=model.matrix(~mlgfd.a,
                                                  data=lincombs.data.estab.mlgfd.a_single)
lincombs.matrix.estab.mlgfd.a_single=as.data.frame(lincombs.matrix.estab.mlgfd.a_single)
lincombs.estab.mlgfd.a_single=inla.make.lincombs(lincombs.matrix.estab.mlgfd.a_single)

inla.model_lincombs.estab.mlgfd.a_single = pglmm(estab ~ mlgfd.a+#(1|species) + 
                                                   (1|f_p) + (1|field), data = dat_suc_sp,
                                                 family = "binomial", cov_ranef = list(species = tree),
                                                 bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                             config = TRUE),
                                                                      quantiles=c(0.025,0.5,0.975),
                                                                      lincomb=lincombs.estab.mlgfd.a_single,
                                                                      control.predictor=list(compute=T)),
                                                 bayes = T)

lincombs.posterior.estab.mlgfd.a_single = inla.model_lincombs.estab.mlgfd.a_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mlgfd.a_single$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mlgfd.a_single$predicted.value=unlist(lapply(lincombs.posterior.estab.mlgfd.a_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mlgfd.a_single$lower=unlist(lapply(lincombs.posterior.estab.mlgfd.a_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mlgfd.a_single$upper=unlist(lapply(lincombs.posterior.estab.mlgfd.a_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mlgfd.a_single
save(lincombs.data.estab.mlgfd.a_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top40_equal_interval_bh_partialb_data/lincombs.data.estab.mlgfd.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mlgfd.a
load("code/results_analyzing/analysing_ages1_35_top40_equal_interval_bh_partialb_data/lincombs.data.estab.mlgfd.a_single.rdata")
(estab.mlgfd.a_single.logistic=ggplot(data=lincombs.data.estab.mlgfd.a_single,aes(x=mlgfd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=2,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mlgfd.a, y=estab),
               color = colors_4d[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y=' ')+
    annotate(geom="text",
             x=min(dat_suc_sp$mlgfd.a)+(max(dat_suc_sp$mlgfd.a)-min(dat_suc_sp$mlgfd.a))*0.3,
             y=c(0.90,0.80),
             label=c("italic(β)['RFD'[ab]] == '3.93'",
                     "'95%CI' == '[3.53, 4.33]'"),
             parse=T,size=4)+
    theme_regular()
)
#mlgfd.a      3.93       3.53       4.33


##### dominance predictive curves for md.ab #####
### get dominance predict data for mnd.a
lincombs.data.domin.mnd.a_single = data.frame(mnd.a=seq(min(dat_dom_sp$mnd.a),
                                                        max(dat_dom_sp$mnd.a),
                                                        length=100))

lincombs.matrix.domin.mnd.a_single=model.matrix(~mnd.a,
                                                data=lincombs.data.domin.mnd.a_single)
lincombs.matrix.domin.mnd.a_single=as.data.frame(lincombs.matrix.domin.mnd.a_single)
lincombs.domin.mnd.a_single=inla.make.lincombs(lincombs.matrix.domin.mnd.a_single)

inla.model_lincombs.domin.mnd.a_single = pglmm(domin ~ mnd.a+#(1|species) + 
                                                 (1|f_p) + (1|field), data = dat_dom_sp,
                                               family = "binomial", cov_ranef = list(species = tree),
                                               bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                           config = TRUE),
                                                                    quantiles=c(0.025,0.5,0.975),
                                                                    lincomb=lincombs.domin.mnd.a_single,
                                                                    control.predictor=list(compute=T)),
                                               bayes = T)

lincombs.posterior.domin.mnd.a_single = inla.model_lincombs.domin.mnd.a_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mnd.a_single$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mnd.a_single$predicted.value=unlist(lapply(lincombs.posterior.domin.mnd.a_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mnd.a_single$lower=unlist(lapply(lincombs.posterior.domin.mnd.a_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mnd.a_single$upper=unlist(lapply(lincombs.posterior.domin.mnd.a_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mnd.a_single
save(lincombs.data.domin.mnd.a_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top40_equal_interval_bh_partialb_data/lincombs.data.domin.mnd.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mnd.a
load("code/results_analyzing/analysing_ages1_35_top40_equal_interval_bh_partialb_data/lincombs.data.domin.mnd.a_single.rdata")
(domin.mnd.a_single.logistic=ggplot(data=lincombs.data.domin.mnd.a_single,aes(x=mnd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1],alpha=0.2)+
    geom_line(color=colors_4d[1],size=2,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #scale_x_continuous(breaks = seq(-0.35, 0.7, 0.35))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnd.a, y=domin),
               color = colors_4d[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Niche difference (ND)', y='Dominance probability')+
    annotate(geom="text",
             x=min(dat_dom_sp$mnd.a)+(max(dat_dom_sp$mnd.a)-min(dat_dom_sp$mnd.a))*0.3,
             y=c(0.90,0.80),
             label=c("italic(β)['ND'[ab]] == '5.41'",
                     "'95%CI' == '[2.77, 8.04]'"),
             parse=T,size=4)+
    theme_regular()
)
#mnd.a        5.41       2.77       8.04


### get dominance predict data for mlgfd.a
lincombs.data.domin.mlgfd.a_single = data.frame(mlgfd.a=seq(min(dat_dom_sp$mlgfd.a),
                                                            max(dat_dom_sp$mlgfd.a),
                                                            length=100))

lincombs.matrix.domin.mlgfd.a_single=model.matrix(~mlgfd.a,
                                                  data=lincombs.data.domin.mlgfd.a_single)
lincombs.matrix.domin.mlgfd.a_single=as.data.frame(lincombs.matrix.domin.mlgfd.a_single)
lincombs.domin.mlgfd.a_single=inla.make.lincombs(lincombs.matrix.domin.mlgfd.a_single)

inla.model_lincombs.domin.mlgfd.a_single = pglmm(domin ~ mlgfd.a+#(1|species) + 
                                                   (1|f_p) + (1|field), data = dat_dom_sp,
                                                 family = "binomial", cov_ranef = list(species = tree),
                                                 bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                             config = TRUE),
                                                                      quantiles=c(0.025,0.5,0.975),
                                                                      lincomb=lincombs.domin.mlgfd.a_single,
                                                                      control.predictor=list(compute=T)),
                                                 bayes = T)

lincombs.posterior.domin.mlgfd.a_single = inla.model_lincombs.domin.mlgfd.a_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mlgfd.a_single$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mlgfd.a_single$predicted.value=unlist(lapply(lincombs.posterior.domin.mlgfd.a_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mlgfd.a_single$lower=unlist(lapply(lincombs.posterior.domin.mlgfd.a_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mlgfd.a_single$upper=unlist(lapply(lincombs.posterior.domin.mlgfd.a_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mlgfd.a_single
save(lincombs.data.domin.mlgfd.a_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top40_equal_interval_bh_partialb_data/lincombs.data.domin.mlgfd.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mlgfd.a
load("code/results_analyzing/analysing_ages1_35_top40_equal_interval_bh_partialb_data/lincombs.data.domin.mlgfd.a_single.rdata")
(domin.mlgfd.a_single.logistic=ggplot(data=lincombs.data.domin.mlgfd.a_single,aes(x=mlgfd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=2,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mlgfd.a, y=domin),
               color = colors_4d[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Relative fitness difference (RFD)', y='')+
    annotate(geom="text",
             x=min(dat_dom_sp$mlgfd.a)+(max(dat_dom_sp$mlgfd.a)-min(dat_dom_sp$mlgfd.a))*0.3,
             y=c(0.9,0.8),
             label=c("italic(β)['RFD'[ab]] == '2.20'",
                     "'95%CI' == '[1.57, 2.83]'"),
             parse=T,size=4)+
    theme_regular()
)
#mlgfd.a      2.20       1.57       2.83



#### Fast start for analyzing invasion success probability ~ mnd+mfd+mpd+mfunc_d for all species ####
setwd("D:/R projects/BSS")
library(phyr)
library(tibble)
library(lme4)
require(ape)
library(scales)
library(ggthemes)
load('code/results_analyzing/analysing_ages1_35_top40_equal_interval_bh_partialb_data/dat_suc_sp.rdata')
numcols = grep("^m",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])

##### Estab #####
pc_prior = list(prec=list("pc.prec", param=c(0.1,0.01)))
dat_suc_sps$species_1 = as.factor(dat_suc_sps$species)
estab_sp_names = unique(dat_suc_sps$species)

tree = read.tree('data/original data/phylo_tree332.txt')
estab_tree_fit = keep.tip(tree, estab_sp_names)
estab_vcv_tree = ape::vcv(estab_tree_fit, model = "Brownian", corr = FALSE)
estab_vcv_tree_sparse = inla.as.sparse(solve(estab_vcv_tree))

## Check the co-linearity
### only continuous functional trait distance
car::vif(
  glmer(estab ~ mnd.a + mlgfd.a + mpd.a_all + mconti_func_d.a_all + (1|f_p)+ (1|field),
        family = binomial, data = dat_suc_sps)
)


###no co-linearity problem, all VIF < 3

### only continuous trait distance 
estab_model_md.a_conti_func_d_all = pglmm(estab~mnd.a+mlgfd.a+mpd.a_all+mconti_func_d.a_all
                                          #+(1|species) 
                                          + (1|f_p) + (1|field), data = dat_suc_sps,
                                          family = "binomial", cov_ranef = list(species = tree),
                                          bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                      config = TRUE),
                                                               quantiles=c(0.025,0.5,0.975)),
                                          bayes = T)

summary(estab_model_md.a_conti_func_d_all)

# plot 
estab_data.inla.md.a_conti_func_d.all_intercept1 = estab_model_md.a_conti_func_d_all$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

estab_data.inla.md.a_conti_func_d.all_intercept = estab_data.inla.md.a_conti_func_d.all_intercept1%>%
  mutate(rowname=c("MND.ab","MRFD.ab","MPD.ab_all","MFD.ab_all"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

### Effect size plot for abundance weighted mean differences
# point + effect size
estab_data.inla.md.a_conti_func_d.all_intercept_ndrfd = 
  estab_data.inla.md.a_conti_func_d.all_intercept %>% 
  filter(rowname %in% c('MRFD.ab', 'MND.ab'))
x_min = min(estab_data.inla.md.a_conti_func_d.all_intercept_ndrfd$lower)
x_max = max(estab_data.inla.md.a_conti_func_d.all_intercept_ndrfd$upper)
(estab_nd_rfd.a.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=estab_data.inla.md.a_conti_func_d.all_intercept_ndrfd,
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
    geom_point(data=estab_data.inla.md.a_conti_func_d.all_intercept_ndrfd,
               mapping = aes(x=mean,y=rowname),
               size = 4.4)+
    geom_errorbar(data=estab_data.inla.md.a_conti_func_d.all_intercept_ndrfd,
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper), width = 0.12, linewidth = 1.2)+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MRFD.ab', 'MND.ab'),
                     label = c(expression(RFD['ab']), expression(ND['ab'])))+
    scale_x_continuous(limits=c(ifelse(x_min>0,-0.15,x_min-0.15),
                                x_max+0.1
    ))+
    scale_color_manual(values = c(colors_4d[1],
                                  colors_4d[2]),
                       labels = c("MND.ab" = "MND.ab",
                                  "MRFD.ab" = "MRFD.ab"),
                       name = ' ')+
    scale_fill_manual(values = c(colors_4d[1],
                                 colors_4d[2]),
                      labels = c("MND.ab" = "MND.ab",
                                 "MRFD.ab" = "MRFD.ab"),
                      name = ' ')+
    annotate(geom="text", x = x_min+(x_max-x_min)*0.6
             , y=2.45, family = 'Arial',
             label='Establishment',
             size = 6)+
    labs(x = '  ', y = '  ')+
    guides(color="none")+
    theme_regular())

estab_data.inla.md.a_conti_func_d.all_intercept_pdfd = 
  estab_data.inla.md.a_conti_func_d.all_intercept %>% 
  filter(rowname %in% c('MFD.ab_all', 'MPD.ab_all'))
x_min = min(estab_data.inla.md.a_conti_func_d.all_intercept_pdfd$lower)
x_max = max(estab_data.inla.md.a_conti_func_d.all_intercept_pdfd$upper)
(estab_pd_fd.a.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=estab_data.inla.md.a_conti_func_d.all_intercept_pdfd,
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
    geom_point(data=estab_data.inla.md.a_conti_func_d.all_intercept_pdfd,
               mapping = aes(x=mean,y=rowname),
               size = 4.4)+
    geom_errorbar(data=estab_data.inla.md.a_conti_func_d.all_intercept_pdfd,
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper), width = 0.12, linewidth = 1.2)+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD.ab_all','MPD.ab_all'),
                     label=c(expression(MFD['ab']), expression(MPD['ab'])))+
    scale_x_continuous(limits=c(ifelse(x_min>0,-0.15,x_min-0.15),
                                ifelse(x_max<0,0.15,x_max+0.15)
    ))+
    scale_color_manual(values = c(colors_4d[4],
                                  colors_4d[3]),
                       labels = c("MPD.ab_all" = "MPD.ab_all",
                                  "MFD.ab_all" = "MFD.ab_all"),
                       name = ' ')+
    scale_fill_manual(values = c(colors_4d[4],
                                 colors_4d[3]),
                      labels = c("MPD.ab_all" = "MPD.ab_all",
                                 "MFD.ab_all" = "MFD.ab_all"),
                      name = ' ')+
    annotate(geom="text", x = x_min+(x_max-x_min)*0.15,
             y=2.45, family = 'Arial',
             label='Establishment', size = 6)+
    labs(x = '  ', y = ' ')+
    guides(color="none")+
    theme_regular())

##### Domin #####
dat_dom_sp = dat_suc_sp %>% filter(stage %in% c('establish', 'dominant'))
dat_dom_sps = dat_suc_sps %>% filter(stage %in% c('establish', 'dominant'))

## Check the co-linearity
car::vif(
  glmer(domin ~ mnd.a + mlgfd.a + mpd.a_all+mconti_func_d.a_all + 
          (1|f_p)+ (1|field),
        family=binomial,data=dat_dom_sps)
)

### no colinearity problem

### only continuous trait distance 
domin_model_md.a_conti_func_d_all = pglmm(domin~mnd.a+mlgfd.a+mpd.a_all+mconti_func_d.a_all
                                          #+(1|species)
                                          + (1|f_p)
                                          + (1|field), data = dat_dom_sps,
                                          family = "binomial", cov_ranef = list(species = tree),
                                          bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                      config = TRUE),
                                                               quantiles=c(0.025,0.5,0.975)),
                                          bayes = T)
summary(domin_model_md.a_conti_func_d_all)

# plot 
domin_data.inla.md.a_conti_func_d.all_intercept1 = domin_model_md.a_conti_func_d_all$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

domin_data.inla.md.a_conti_func_d.all_intercept = domin_data.inla.md.a_conti_func_d.all_intercept1%>%
  mutate(rowname=c("MND.ab","MRFD.ab","MPD.ab_all","MFD.ab_all"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))%>% 
  mutate(significant = ifelse(lower / abs(lower) == upper / abs(upper), 'yes', 'no'))


### Effect size plot for abundance weighted mean differences
# point + effect size
domin_data.inla.md.a_conti_func_d.all_intercept_ndrfd = 
  domin_data.inla.md.a_conti_func_d.all_intercept %>% 
  filter(rowname %in% c('MRFD.ab', 'MND.ab'))
x_min = min(domin_data.inla.md.a_conti_func_d.all_intercept_ndrfd$lower)
x_max = max(domin_data.inla.md.a_conti_func_d.all_intercept_ndrfd$upper)
(domin_nd_rfd.a.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=domin_data.inla.md.a_conti_func_d.all_intercept_ndrfd,
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
    scale_color_manual(values = c(colors_4d[1],
                                  colors_4d[2]),
                       labels = c("MND.ab" = "MND.ab",
                                  "MRFD.ab" = "MRFD.ab"),
                       name = ' ')+
    scale_fill_manual(values = c(colors_4d[1],
                                 colors_4d[2]),
                      labels = c("MND.ab" = "MND.ab",
                                 "MRFD.ab" = "MRFD.ab"),
                      name = ' ')+
    geom_point(data=domin_data.inla.md.a_conti_func_d.all_intercept_ndrfd,
               mapping = aes(x=mean,y=rowname, color = significant),
               size = 4.4)+
    geom_errorbar(data=domin_data.inla.md.a_conti_func_d.all_intercept_ndrfd,
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper, color = significant),
                  width = 0.12, linewidth = 1.2)+
    scale_color_manual(values = c('grey',
                                  'black'),
                       name = ' ')+
    ggnewscale::new_scale_color()+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MRFD.ab', 'MND.ab'),
                     label =c(c(expression(RFD['ab']),
                                expression(ND['ab']))))+
    scale_x_continuous(limits=c(ifelse(x_min>0,-0.15,x_min-0.15),
                                ifelse(x_max<0,0.15,x_max+0.15)))+
    annotate(geom="text", x = x_min+(x_max-x_min)*0.85,
             y=2.45, family = 'Arial',
             label='Dominance', size = 6)+
    labs(x = 'Standardized effects', y = ' ')+
    guides(color="none")+
    theme_regular())

domin_data.inla.md.a_conti_func_d.all_intercept_pdfd = 
  domin_data.inla.md.a_conti_func_d.all_intercept %>% 
  filter(rowname %in% c('MFD.ab_all', 'MPD.ab_all'))
x_min = min(domin_data.inla.md.a_conti_func_d.all_intercept_pdfd$lower)
x_max = max(domin_data.inla.md.a_conti_func_d.all_intercept_pdfd$upper)
(domin_pd_fd.a.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=domin_data.inla.md.a_conti_func_d.all_intercept_pdfd,
             mapping = aes(x=mean,y=rowname,fill=rowname), 
             width = 0.4)+
    geom_point(data=domin_data.inla.md.a_conti_func_d.all_intercept_pdfd,
               mapping = aes(x=mean,y=rowname), color = 'grey',
               size = 4.4)+
    geom_errorbar(data=domin_data.inla.md.a_conti_func_d.all_intercept_pdfd,
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper), color = 'grey',
                  width = 0.12, linewidth = 1.2)+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD.ab_all','MPD.ab_all'),
                     label=c(expression(MFD['ab']), expression(MPD['ab'])))+
    scale_x_continuous(limits=c(ifelse(x_min>0,-0.15,x_min-0.15),
                                ifelse(x_max<0,0.15,x_max+0.15)
    ))+
    scale_color_manual(values = c(colors_4d[4],
                                  colors_4d[3]),
                       labels = c("MPD.ab_all" = "MPD.ab_all",
                                  "MFD.ab_all" = "MFD.ab_all"),
                       name = ' ')+
    scale_fill_manual(values = c(colors_4d[4],
                                 colors_4d[3]),
                      labels = c("MPD.ab_all" = "MPD.ab_all",
                                 "MFD.ab_all" = "MFD.ab_all"),
                      name = ' ')+
    annotate(geom="text", x = x_min+(x_max-x_min)*0.15,
             y=2.45, family = 'Arial',
             label='Dominance', size = 6)+
    labs(x = '  ', y = ' ')+
    guides(color="none")+
    theme_regular())



#### Fig.S3：BH model results for demonstration of how Weighted Mean ND/RFD affects successful colonization and dominance of species, including presentation of raw data and single & multi-regression results ####
library(ggpubr)
library(export)
library(devEMF)
library(cowplot)

# Adjust grid settings
num_rows  =  2
num_cols  =  3
plot_width  =  1 / num_cols
plot_height  =  1 / num_rows

# Create a ggdraw object with all plots
Fig.S3.all = ggdraw() +
  draw_plot(estab.mnd.a_single.logistic, x = 0, y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(estab.mlgfd.a_single.logistic, x = plot_width, y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(estab_nd_rfd.a.all.varied.intercept.plot, x = 2 * plot_width, y = 0.5,
            width = plot_width,height = plot_height) +
  draw_plot(domin.mnd.a_single.logistic, x = 0, y = 0.05,
            width = plot_width, height = plot_height) +
  draw_plot(domin_nd_rfd.a.all.varied.intercept.plot, x = 2 * plot_width, y = 0.05,
            width = plot_width, height = plot_height) +
  draw_plot(domin.mlgfd.a_single.logistic, x = plot_width, y = 0.05, 
            width = plot_width, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f"),
    x = c(0.01, 1 * plot_width+0.01, (2 * plot_width)+0.01,
          0.01, 1 * plot_width+0.01, (2 * plot_width)+0.01),
    y = c(1, 1, 1, 0.55, 0.55, 0.55),
    hjust = 0, vjust = 1.1, size = 18
  )

Fig.S3.all
emf('results/figures_ages1_35_top50_equal_interval_bh_partialb/Fig.S3.emf',
    width = 25.2*1.1, height = 16.2*1.1, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.S3.all
dev.off() #turn off device and finalize file



