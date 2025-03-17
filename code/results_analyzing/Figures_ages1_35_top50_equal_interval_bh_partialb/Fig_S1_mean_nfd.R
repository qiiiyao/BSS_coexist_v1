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
library(devEMF)
library(ggeffects)
library(stringr)
library(viridis)
library(showtext)
library(MuMIn)
library(extrafont)
windowsFonts(Arial=windowsFont("Arial"))
library(scatterpie)
source('code/function/plot_func.R')

load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_allb_data/dat_suc_sp.rdata")
dat_suc_sp_bh_allb = dat_suc_sp
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/dat_suc_sp.rdata")
dat_suc_sp_bh_partialb = dat_suc_sp
rm(dat_suc_sp)

dat_suc_sp_bh_allb$f_p_sp = paste(dat_suc_sp_bh_allb$f_p,
                                  dat_suc_sp_bh_allb$species,
                                  sep = '_')

dat_suc_sp_bh_partialb$f_p_sp = paste(dat_suc_sp_bh_partialb$f_p,
                                      dat_suc_sp_bh_partialb$species,
                                      sep = '_')

dat_suc_sp_comparison = dat_suc_sp_bh_partialb %>% 
  left_join(dat_suc_sp_bh_allb[,c('f_p_sp', "mnd.a" ,"mlgfd.a")],
            by = 'f_p_sp') %>% 
  filter(!is.na(mnd.a.y))

plot(dat_suc_sp_comparison$mnd.a.x,
     dat_suc_sp_comparison$mnd.a.y)



### the corrletion between estimated ND values from model_bh_partialb and from model_bh_allb
nd_model_cor = lmer(mnd.a.y ~ mnd.a.x + (1|f_p),
                    data = dat_suc_sp_comparison,
                    REML = T)
r.squaredGLMM(nd_model_cor)

nd_model_cor = lm(mnd.a.y ~ mnd.a.x,
                  data = dat_suc_sp_comparison)
summary(nd_model_cor)

#### Get the model predictive value
v = seq(min(dat_suc_sp_comparison$mnd.a.x),
        max(dat_suc_sp_comparison$mnd.a.x),
        length.out = 100)

nd_model_cor_pre = ggpredict(nd_model_cor,
                             type = 'fixed', terms = 'mnd.a.x[v]')

range_nd_x = range(dat_suc_sp_comparison$mnd.a.x)
range_nd_y = range(dat_suc_sp_comparison$mnd.a.y)
fig_nd_model_cor = ggplot()+
  geom_point(data = dat_suc_sp_comparison,
    mapping = aes(x = mnd.a.x, y = mnd.a.y),
     color = "grey50",
    alpha = 0.5)+
  geom_ribbon(data = nd_model_cor_pre,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = 'darkgrey', alpha = 0.5) + # CI ribbon 
  geom_line(data = nd_model_cor_pre,
            mapping = aes(x = x, y = predicted),
            linewidth = 2,
            color = 'black')+
  labs(x = 'Partial_EBH', y = 'EBH')+
  annotate("text", x=range_nd_x[1]+(range_nd_x[2]-range_nd_x[1])*0.225,
           y=range_nd_y[1]+(range_nd_y[2]-range_nd_y[1])*0.99,
           label = expression("Niche difference"["ND"[ab]])
            , size = 5
  )+
  annotate("text", x=range_nd_x[1]+(range_nd_x[2]-range_nd_x[1])*0.1,
           y=range_nd_y[1]+(range_nd_y[2]-range_nd_y[1])*0.92,
           label = expression(R^2 ~ "= 0.39")
            , size = 5
  )+
  scale_x_continuous(limits = c(range_nd_x[1], 0.82)) +
  theme_regular()

### the corrletion between estimated lgfd values from model_bh_partialb and from model_bh_allb
lgfd_model_cor = lmer(mlgfd.a.y ~ mlgfd.a.x + (1|f_p),
                      data = dat_suc_sp_comparison,
                      REML = T)
cor.test(dat_suc_sp_comparison$mlgfd.a.x,
         dat_suc_sp_comparison$mlgfd.a.y)
summary(lgfd_model_cor)

lgfd_model_cor = lm(mlgfd.a.y ~ mlgfd.a.x,
                    data = dat_suc_sp_comparison)
summary(lgfd_model_cor)

#### Get the model predictive value
v = seq(min(dat_suc_sp_comparison$mlgfd.a.x),
        max(dat_suc_sp_comparison$mlgfd.a.x),
        length.out = 100)
lgfd_model_cor_pre = ggpredict(lgfd_model_cor,
                               type = 'fixed', terms = 'mlgfd.a.x[v]')

range_lgfd_x = range(dat_suc_sp_comparison$mlgfd.a.x)
range_lgfd_y = range(dat_suc_sp_comparison$mlgfd.a.y)
fig_lgfd_model_cor = ggplot()+
  geom_point(data = dat_suc_sp_comparison,
    mapping = aes(x = mlgfd.a.x, y = mlgfd.a.y),
   color = "grey50",
    alpha = 0.5)+
  geom_ribbon(data = lgfd_model_cor_pre,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = 'darkgrey', alpha = 0.5) + # CI ribbon 
  geom_line(data = lgfd_model_cor_pre,
            mapping = aes(x = x, y = predicted),
            linewidth = 2,
            color = 'black')+
  annotate("text", x=range_lgfd_x[1]+(range_lgfd_x[2]-range_lgfd_x[1])*0.345,
           y=range_lgfd_y[1]+(range_lgfd_y[2]-range_lgfd_y[1])*0.99,
           label = expression("Relative fitness difference"["RFD"[ab]])
            , size = 5
  )+
  annotate("text", x=range_lgfd_x[1]+(range_lgfd_x[2]-range_lgfd_x[1])*0.1,
           y=range_lgfd_y[1]+(range_lgfd_y[2]-range_lgfd_y[1])*0.92,
           label = expression(R^2 ~ "= 0.72")
            , size = 5
  )+
  labs(x = 'Partial_EBH', y = ' ')+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1) ) +
  theme_regular()


#Fig.S1
# Merge and export eight plots
sem_nfd_model_cor = ggarrange(fig_nd_model_cor, fig_lgfd_model_cor,
                              nrow = 1, ncol = 2,
                              labels = c('a',
                                         'b'),
                              
                              font.label = list(size = 18)
                              #,vjust = 1.8
)

emf('results/figures_ages1_35_top50_equal_interval_bh_partialb/Fig.S1_mod_cor.emf',
    width = 30*0.8, height = 15*0.8, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
sem_nfd_model_cor
dev.off() #turn off device and finalize file
