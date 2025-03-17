############### Fast Start ####################
rm(list = ls())
setwd("D:/R projects/BSS")

### load packages
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
library(MuMIn)
library(extrafont)
windowsFonts(Arial=windowsFont("Arial"))
library(scatterpie)
source('code/function/plot_func.R')

load("results/fit_results/BSS_exclude_trees_raw/plot_ages1_35_top50_equal_interval_model_comparison/bh_allb/inter_all_c_alltime_rhat_105.rdata")
inter_all_c_alltime_bh_allb = inter_all_c_alltime
load("results/fit_results/BSS_exclude_trees_raw/plot_ages1_35_top50_equal_interval_model_comparison/bh_partialb/inter_all_c_alltime_rhat_105.rdata")
inter_all_c_alltime_bh_partialb = inter_all_c_alltime
rm(inter_all_c_alltime)

inter_all_c_alltime_comparison = inter_all_c_alltime_bh_partialb %>% 
  left_join(inter_all_c_alltime_bh_allb[,c('f_p_sp_pair', "nd" ,"fd", "lgfd")],
            by = 'f_p_sp_pair') %>% 
  filter(!is.na(nd.y))

plot(inter_all_c_alltime_comparison$nd.x,
     inter_all_c_alltime_comparison$nd.y)



### the corrletion between estimated ND values from model_bh_partialb and from model_bh_allb
nd_model_cor = lmer(nd.y ~ nd.x + (1|f_p) + (1|sp_pair),
                      data = inter_all_c_alltime_comparison,
                      REML = T)
nd_model_cor = lm(nd.y ~ nd.x,
                    data = inter_all_c_alltime_comparison)
summary(nd_model_cor)
r.squaredGLMM(nd_model_cor)

#### Get the model predictive value
v = seq(min(inter_all_c_alltime_comparison$nd.x),
        max(inter_all_c_alltime_comparison$nd.x),
        length.out = 100)

nd_model_cor_pre = ggpredict(nd_model_cor,
                             type = 'fixed', terms = 'nd.x[v]')

range_nd_x = range(nd_model_cor_pre$x)
range_nd_y = range(nd_model_cor_pre$predicted)
fig_nd_model_cor = ggplot()+
    #geom_point(data = inter_all_c_alltime_comparison,
            #  mapping = aes(x = nd.x, y = nd.y),
             #   color = "grey50",
              #  alpha = 0.5)+
    geom_ribbon(data = nd_model_cor_pre,
                aes(x = x, ymin = conf.low, ymax = conf.high),
                fill = 'darkgrey', alpha = 0.5) + # CI ribbon 
    geom_line(data = nd_model_cor_pre,
              mapping = aes(x = x, y = predicted),
              linewidth = 2,
              color = 'black')+
    labs(x = 'Model 1', y = 'Model 2')+
    annotate("text", x=range_nd_x[1]+(range_nd_x[2]-range_nd_x[1])*0.1,
             y=range_nd_y[1]+(range_nd_y[2]-range_nd_y[1])*0.97,
             label = 'ND, R2 = 0.16'
          # , size = 18
    )+
    #scale_y_continuous(limits = c(0.5, 5.2), breaks = c(1, 2, 3, 4, 5)) +
    theme_regular()

### the corrletion between estimated lgfd values from model_bh_partialb and from model_bh_allb
lgfd_model_cor = lmer(lgfd.y ~ lgfd.x + (1|f_p) + (1|sp_pair),
                    data = inter_all_c_alltime_comparison,
                    REML = T)
cor.test(inter_all_c_alltime_comparison$lgfd.x,
         inter_all_c_alltime_comparison$lgfd.y)
summary(lgfd_model_cor)

lgfd_model_cor = lm(lgfd.y ~ lgfd.x,
                  data = inter_all_c_alltime_comparison)
summary(lgfd_model_cor)

#### Get the model predictive value
v = seq(min(inter_all_c_alltime_comparison$lgfd.x),
        max(inter_all_c_alltime_comparison$lgfd.x),
        length.out = 100)
lgfd_model_cor_pre = ggpredict(lgfd_model_cor,
                             type = 'fixed', terms = 'lgfd.x[v]')

range_lgfd_x = range(lgfd_model_cor_pre$x)
range_lgfd_y = range(lgfd_model_cor_pre$predicted)
fig_lgfd_model_cor = ggplot()+
  #geom_point(data = inter_all_c_alltime_comparison,
           #  mapping = aes(x = lgfd.x, y = lgfd.y),
            # color = "grey50",
           #  alpha = 0.5)+
  geom_ribbon(data = lgfd_model_cor_pre,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = 'darkgrey', alpha = 0.5) + # CI ribbon 
  geom_line(data = lgfd_model_cor_pre,
            mapping = aes(x = x, y = predicted),
            linewidth = 2,
            color = 'black')+
  annotate("text", x=range_lgfd_x[1]+(range_lgfd_x[2]-range_lgfd_x[1])*0.1,
           y=range_lgfd_y[1]+(range_lgfd_y[2]-range_lgfd_y[1])*0.97,
           label = 'RFD, R2 = 0.78'
           #, size = 18
  )+
  labs(x = 'Model 1', y = 'Model 2')+
  #scale_y_continuous(limits = c(0.5, 5.2), breaks = c(1, 2, 3, 4, 5)) +
  theme_regular()


#Fig.S1
# Merge and export eight plots
sem_nfd_model_cor = ggarrange(fig_nd_model_cor, fig_lgfd_model_cor,
                                  nrow = 1, ncol = 2,
                                  labels = c('a',
                                             'b')
                              #,vjust = 1.8
                              )
emf('results/figures_ages1_35_top50_equal_interval_bh_partialb/Fig.S1_mod_cor.emf',
    width = 30*0.8, height = 16*0.8, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
sem_nfd_model_cor
dev.off() #turn off device and finalize file
