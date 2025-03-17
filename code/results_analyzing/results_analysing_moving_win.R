### Load data
library(data.table)
load("D:/R projects/BSS/results/fit_results/inter_moving_c_l.rdata")
inter_moving_c_l_pure = lapply(inter_moving_c_l,
                               function(x){x[[1]]}) 
inter_all_c_movingwin = rbindlist(inter_moving_c_l_pure)
save(inter_all_c_movingwin,
     file = 'results/fit_results/inter_all_c_movingwin.rdata')

############ Fast start #############
## Add stage introduce
### load packages ###
rm(list = ls())
library(dplyr)
library(stringr)
library(data.table)
library(parallel)
library(doParallel)

load("D:/R projects/BSS/results/fit_results/inter_all_c_movingwin.rdata")
inter_all_c_movingwin$initial_age = inter_all_c_movingwin$t_win+(inter_all_c_movingwin$t_win-1)
inter_all_c_movingwin$intro_i = 1
inter_all_c_movingwin$intro_j = 1
inter_all_c_movingwin[inter_all_c_movingwin$stage_i == 'native',]$intro_i = 0
inter_all_c_movingwin[inter_all_c_movingwin$stage_j == 'native',]$intro_j = 0

inter_all_c_movingwin$stage_ij_estab = inter_all_c_movingwin$stage_ij
inter_all_c_movingwin$stage_ij_domin = inter_all_c_movingwin$stage_ij

inter_all_c_movingwin$stage_i_estab = 'native'
inter_all_c_movingwin$stage_j_estab = 'native'
inter_all_c_movingwin[inter_all_c_movingwin$intro_i == 1&
              inter_all_c_movingwin$estab_i == 0,]$stage_i_estab = 'intro.noestab'
inter_all_c_movingwin[inter_all_c_movingwin$intro_i == 1&
              inter_all_c_movingwin$estab_i == 1,]$stage_i_estab = 'intro.estab'
inter_all_c_movingwin[inter_all_c_movingwin$intro_j == 1&
              inter_all_c_movingwin$estab_j == 0,]$stage_j_estab = 'intro.noestab'
inter_all_c_movingwin[inter_all_c_movingwin$intro_j == 1&
              inter_all_c_movingwin$estab_j == 1,]$stage_j_estab = 'intro.estab'
inter_all_c_movingwin$stage_ij_estab = paste(inter_all_c_movingwin$stage_i_estab,
                                   inter_all_c_movingwin$stage_j_estab,
                                   sep = '_')

inter_all_c_movingwin$stage_i_domin = 'native'
inter_all_c_movingwin$stage_j_domin = 'native'
inter_all_c_movingwin[inter_all_c_movingwin$intro_i == 1&
              inter_all_c_movingwin$estab_i == 1&
              inter_all_c_movingwin$domin_i == 0,]$stage_i_domin = 'estab.nodomin'
inter_all_c_movingwin[inter_all_c_movingwin$intro_i == 1&
              inter_all_c_movingwin$estab_i == 1&
              inter_all_c_movingwin$domin_i == 1,]$stage_i_domin = 'estab.domin'
inter_all_c_movingwin[inter_all_c_movingwin$intro_j == 1&
              inter_all_c_movingwin$estab_j == 1&
              inter_all_c_movingwin$domin_j == 0,]$stage_j_domin = 'estab.nodomin'
inter_all_c_movingwin[inter_all_c_movingwin$intro_j == 1&
              inter_all_c_movingwin$estab_j == 1&
              inter_all_c_movingwin$domin_j == 1,]$stage_j_domin = 'estab.domin'
inter_all_c_movingwin$stage_ij_domin = paste(inter_all_c_movingwin$stage_i_domin,
                                   inter_all_c_movingwin$stage_j_domin,
                                   sep = '_')

inter_all_c_movingwin = inter_all_c_movingwin %>% relocate(intro_i, .after = stage_ij) %>% 
  relocate(intro_j, .after = intro_i) %>% 
  relocate(stage_ij_estab, .after = stage_ij) %>% 
  relocate(stage_i_estab, .after = stage_ij) %>% 
  relocate(stage_j_estab, .after = stage_i_estab) %>% 
  relocate(stage_ij_domin, .after = stage_ij_estab) %>% 
  relocate(stage_i_domin, .after = stage_ij_estab) %>% 
  relocate(stage_j_domin, .after = stage_i_domin)


############# Draw the coexistence plot ###############
library(lme4)
library(lmerTest)
library(INLA)
library(ggplot2)
library(ggpubr)
library(ggeffects)
library(viridis)
library(ggbreak)

source('code/function/plot_func.R')

########### Draw standard deviation simply #############
inter_all_movingwin_forestab = inter_all_c_movingwin %>%
  filter(stage_ij_estab %in% c("intro.estab_native",
                               "intro.noestab_native"))
str(inter_all_movingwin_forestab)
### Calculate SD at different stage, time windows
#### Established vs. nonestablished #####
dat_movingwin_forcoexplot_estab = inter_all_movingwin_forestab %>%
                            group_by(stage_ij_estab, t_win, initial_age) %>%
                            summarise_at(vars(nd, lgfd),
                                         list(mean = mean,sd = sd))
spss.f = function(x) log10(1-x)
x = seq(-1, 1, 0.001)
spss.ff = function(x) -log10(1-x)

(Fig.4_estab = ggplot(data = dat_movingwin_forcoexplot_estab,
               aes(x = nd_mean,y = lgfd_mean,
                   color = initial_age, shape = stage_ij_estab)) +
  geom_errorbar(mapping = aes(ymin = (lgfd_mean-lgfd_sd),
                              ymax = (lgfd_mean+lgfd_sd)),
                linewidth = 1, 
                show.legend = T)+
  geom_errorbarh(mapping = aes(xmin = (nd_mean-nd_sd),
                               xmax = (nd_mean+nd_sd)), 
                 linewidth = 1, 
                 show.legend = T)+
  geom_point(mapping = aes(x = nd_mean, y = lgfd_mean),
             size=40) + 
  geom_path(aes(x = nd_mean, y = lgfd_mean),
            color = 'grey', linewidth = 4)+
  scale_shape_manual(values=c(19,15),
                     guide = guide_legend(title = 'Stage_ij',
                                            override.aes = list(size = 40)))+
  scale_color_viridis(direction = -1,option = 'A', 
                       discrete = F,
                       guide = guide_colorbar(title = 'Initial age',
                                              barwidth = 4,
                                              barheight = 40))+
  #scale_color_brewer(palette="Set2")+
  #scale_x_continuous(limits=c(0, 0.53)) +
  #scale_y_continuous(limits=c(-0.2, 0.2)) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.5) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=0.5) +
  xlab("Niche difference") +
  ylab(expression(paste(log[10],
                        "(Fitness difference)",
                        sep=' ')))+
  stat_function(fun=spss.f, colour="black")+
  stat_function(fun=spss.ff, colour="black")+
  #scale_x_continuous(breaks=c(-1, 1)) +
  annotate(geom="text", x=0.4, y=0.15, label=str_wrap("Invasion and Coexistence",
                                                        width = 20), size=40) +
  annotate(geom="text", x=0.4, y=-0.15, label=str_wrap("Invasion and Coexistence",
                                                         width = 20), size=40) +
  annotate(geom="text", x=0.1, y=0.15, label=str_wrap("Invasion and exclude residents",
                                                       width = 20), size=40) +
  annotate(geom="text", x=0.1, y=-0.15, label=str_wrap("Residents repel invasion",
                                                        width = 20), size=40) +
    guides(shape = guide_legend(title='Invasion stage'))+
  theme(panel.background = element_rect(fill = 'white',color = 'black',linewidth = 2),
        panel.grid = element_blank(),
        legend.position = c(0.2, 0.88),
        legend.background = element_blank(),
        legend.text = element_text(size=120),
        legend.title = element_text(size=145),
        plot.margin = margin(10,10,10,10),
        plot.background = element_blank(),
        text = element_text(size = 150),
        axis.ticks = element_line(linewidth = 3),
        axis.ticks.length = unit(-5,'lines'),
        axis.title.y = element_text(margin = margin(0,3,0,0),color = '#000000',face = 'bold'),
        axis.title.x = element_text(margin = margin(5,0,0,0),color = '#000000',face = 'bold'),
        axis.text.y = element_text(margin = margin(0,3,0,0),color = '#000000'),
        axis.text.x = element_text(margin = margin(5,0,0,0),color = '#000000')))

#### Dominant vs. nondominant #####
inter_all_movingwin_fordomin = inter_all_c_movingwin %>%
  filter(stage_ij_domin %in% c("estab.domin_native",
                               "estab.nodomin_native"))

dat_movingwin_forcoexplot_domin = inter_all_movingwin_fordomin %>%
  group_by(stage_ij_domin, t_win, initial_age) %>%
  summarise_at(vars(nd, lgfd),
               list(mean = mean,sd = sd))

spss.f = function(x) log10(1-x)
x = seq(-1, 1, 0.001)
spss.ff = function(x) -log10(1-x)

(Fig.4_domin = ggplot(data = dat_movingwin_forcoexplot_domin,
                     aes(x = nd_mean,y = lgfd_mean,
                         color = initial_age, shape = stage_ij_domin)) +
  geom_errorbar(mapping = aes(ymin = (lgfd_mean-lgfd_sd),
                              ymax = (lgfd_mean+lgfd_sd)),
                linewidth = 1, 
                show.legend = T)+
  geom_errorbarh(mapping = aes(xmin = (nd_mean-nd_sd),
                               xmax = (nd_mean+nd_sd)), 
                 linewidth = 1, 
                 show.legend = T)+
  geom_point(mapping = aes(x = nd_mean, y = lgfd_mean),
             size=40) + 
  geom_line(aes(x = nd_mean, y = lgfd_mean),
            color = 'grey', linewidth = 4)+
  scale_shape_manual(values=c(19,15),
                     guide = guide_legend(title = 'Stage_ij',
                                          override.aes = list(size = 25)))+
  scale_color_viridis(direction = -1,option = 'A', 
                      discrete = F,
                      guide = guide_colorbar(title = 'Initial age',
                                             barwidth = 4,
                                             barheight = 40))+
  #scale_color_brewer(palette="Set2")+
  #scale_x_continuous(limits=c(0, 0.53)) +
  #scale_y_continuous(limits=c(-0.2, 0.2)) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.5) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=0.5) +
  labs(x = "Niche difference", 
       y = NULL) +
  stat_function(fun=spss.f, colour="black")+
  stat_function(fun=spss.ff, colour="black")+
  #scale_x_continuous(breaks=c(-1, 1)) +
  annotate(geom="text", x=0.4, y=0.15, label=str_wrap("Invasion and Coexistence",
                                                       width = 20), size=40) +
  annotate(geom="text", x=0.4, y=-0.15, label=str_wrap("Invasion and Coexistence",
                                                        width = 20), size=40) +
  annotate(geom="text", x=0.1, y=0.15, label=str_wrap("Invasion and exclude residents",
                                                       width = 20), size=40) +
  annotate(geom="text", x=0.1, y=-0.15, label=str_wrap("Residents repel invasion",
                                                        width = 20), size=40) +
    guides(shape = guide_legend(title='Invasion stage'))+
  theme(panel.background = element_rect(fill = 'white',color = 'black',linewidth = 2),
        panel.grid = element_blank(),
        legend.position = c(0.2, 0.88),
        legend.background = element_blank(),
        legend.text = element_text(size=120),
        legend.title = element_text(size=145),
        plot.margin = margin(10,10,10,10),
        plot.background = element_blank(),
        text = element_text(size = 150),
        axis.ticks = element_line(linewidth = 3),
        axis.ticks.length = unit(-5,'lines'),
        axis.title.y = element_text(margin = margin(0,3,0,0),color = '#000000',face = 'bold'),
        axis.title.x = element_text(margin = margin(5,0,0,0),color = '#000000',face = 'bold'),
        axis.text.y = element_text(margin = margin(0,3,0,0),color = '#000000'),
        axis.text.x = element_text(margin = margin(5,0,0,0),color = '#000000')))

#### Merge two sd-coexistence plot
gap = ggplot(NULL)+theme_void()
library(ggpubr)

Fig.4_split_stage_2 = ggarrange(gap,Fig.4_estab,gap,
                                Fig.4_domin, nrow = 4, ncol = 1,
                                labels = c('','a)','','b)'), hjust = -0.5,
                                vjust = 0, heights = c(0.1, 1, 0.1, 1),
                                font.label = list(size = 150))
ggsave(plot = Fig.4_split_stage_2,
       'results/figures_movingwin/Fig.4_split_stage_2.svg',
       width = 210,height = 340, dpi = 300, units = 'cm',
       limitsize = F)

Fig.4_split_stage_2_for_ppt = ggarrange(gap,Fig.4_estab,gap,
                                Fig.4_domin, nrow = 1, ncol = 4,
                                labels = c('','a)','','b)'), hjust = 0.5,
                                vjust = 1, widths = c(0.05, 1, 0.05, 1),
                                font.label = list(size = 150))
ggsave(plot = Fig.4_split_stage_2_for_ppt,
       'results/figures_movingwin/Fig.4_split_stage_2_for_ppt.svg',
       width = 400,height = 200, dpi = 300, units = 'cm',
       limitsize = F)


######### LMER methods ############
####### Just split three stages ########
### nd ~ stage
mod_nd_compare_movingwin_1 = lmer(nd ~ stage_ij*t_win+(1|field/plot),
                                  data = inter_all_c_movingwin_forcoexplot, 
                                  REML = F)
mod_nd_compare_movingwin_lm = lm(nd ~ stage_ij*t_win,
                                 data = inter_all_c_movingwin_forcoexplot)
anova(mod_nd_compare_movingwin_1, mod_nd_compare_movingwin_lm)
mod_nd_compare_movingwin_2 = lmer(nd ~ stage_ij*t_win+(1|field/plot)+(1|sp_pair),
                         data = inter_all_c_movingwin_forcoexplot, 
                         REML = F)
anova(mod_nd_compare_movingwin_2, mod_nd_compare_movingwin_1)
mod_nd_compare_movingwin = lmer(nd ~ stage_ij*t_win+(1|field/plot)+(1|sp_pair),
                                  data = inter_all_c_movingwin_forcoexplot, 
                                  REML = T)

### lgfd ~ stage
mod_lgfd_compare_movingwin_1 = lmer(lgfd ~ stage_ij*t_win+(1|field/plot),
                                  data = inter_all_c_movingwin_forcoexplot, 
                                  REML = F)
mod_lgfd_compare_movingwin_lm = lm(lgfd ~ stage_ij*t_win,
                                   data = inter_all_c_movingwin_forcoexplot)
anova(mod_lgfd_compare_movingwin_1, mod_lgfd_compare_movingwin_lm)
mod_lgfd_compare_movingwin_2 = lmer(lgfd ~ stage_ij*t_win+(1|field/plot)+(1|sp_pair),
                                  data = inter_all_c_movingwin_forcoexplot, 
                                  REML = F)
anova(mod_lgfd_compare_movingwin_2, mod_lgfd_compare_movingwin_1)
mod_lgfd_compare_movingwin = lmer(lgfd ~ stage_ij*t_win+(1|field/plot)+(1|sp_pair),
                                data = inter_all_c_movingwin_forcoexplot, 
                                REML = T)

summary(mod_nd_compare_movingwin)
summary(mod_lgfd_compare_movingwin)

pre_mod_nd_compare_movingwin = as.data.frame(ggpredict(mod_nd_compare_movingwin,
                                                terms = c('stage_ij',
                                                't_win[1,2,3,4,5,6,7,8,9,10,11,12,13,14]')))

pre_mod_lgfd_compare_movingwin = as.data.frame(ggpredict(mod_lgfd_compare_movingwin,
                                                terms = c('stage_ij',
                                                't_win[1,2,3,4,5,6,7,8,9,10,11,12,13,14]')))

pre_ndfd_movingwin = data.frame(nd = pre_mod_nd_compare_movingwin$predicted,
                         lgfd = pre_mod_lgfd_compare_movingwin$predicted,
                         invasion_stage = pre_mod_nd_compare_movingwin$x,
                         nd_sd = pre_mod_nd_compare_movingwin$std.error,
                         nd_low = pre_mod_nd_compare_movingwin$conf.low,
                         nd_high = pre_mod_nd_compare_movingwin$conf.high,
                         fd_sd = pre_mod_lgfd_compare_movingwin$std.error,
                         fd_low = pre_mod_lgfd_compare_movingwin$conf.low,
                         fd_high = pre_mod_lgfd_compare_movingwin$conf.high,
                         time_window = unique(pre_mod_lgfd_compare_movingwin$group))

spss.f = function(x) log10(1-x)
x = seq(-1, 1, 0.001)
spss.ff = function(x) -log10(1-x)

Fig.2 = ggplot(data = pre_ndfd_movingwin,
               aes(x = nd,y = lgfd,
                             color = time_window,
                             shape = invasion_stage)) +
  geom_errorbar(mapping = aes(ymin = fd_low,
                                ymax = fd_high),
                  linewidth = 4, 
                  show.legend = T)+
  geom_errorbarh(mapping = aes(xmin = nd_low,
                                xmax = nd_high), 
                  linewidth = 4, 
                  show.legend = T)+
  geom_point(size=20) + 
  scale_shape_manual(values=c(19,15,17))+
  scale_color_viridis(direction = -1,option = 'A', discrete = T)+
  #scale_color_brewer(palette="Set2")+
  #scale_x_break(c(0, 0.28),#?????????????????????
  #              space = 0.3,#????????????
  #              scales = 1.5)+
  scale_x_continuous(limits = c(0.29, 0.36), expand = c(0, 0)) +
  scale_y_continuous(limits=c(-0.1, 0.1)) +
  #coord_cartesian(xlim = c(0.29, 0.36))+
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.5) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=0.5) +
  xlab("Niche difference") +
  ylab(expression(paste(log[10],
                        "(Fitness difference)",
                        sep=' ')))+
  stat_function(fun=spss.f, colour="black")+
  stat_function(fun=spss.ff, colour="black")+
  #scale_x_continuous(breaks=c(-1, 1)) +
  annotate(geom="text", x=0.55, y=0.092, label=str_wrap("Invasion and Coexistence",
                                                        width = 20), size=20) +
  annotate(geom="text", x=0.55, y=-0.092, label=str_wrap("Invasion and Coexistence",
                                                         width = 20), size=20) +
  annotate(geom="text", x=0.1, y=0.092, label=str_wrap("Invasion and exclude residents",
                                                       width = 20), size=20) +
  annotate(geom="text", x=0.1, y=-0.092, label=str_wrap("Residents repel invasion",
                                                        width = 20), size=20) +
  theme_custom_legend()

ggsave('results/figures_movingwin/Fig.2.svg',
       plot=Fig.2, 
       width=105, height=80, dpi=600, units='cm',
       limitsize=F)

########## Stage split two parts: estab/noestab; domin/nodomin ##########
### Normal methods ####
#### estab noestab ####
unique(inter_all_c_movingwin$stage_ij_estab)
inter_all_movingwin_forestab = inter_all_c_movingwin %>%
  filter(stage_ij_estab %in% c("intro.estab_native",
                               "intro.noestab_native"))
mod_nd_compare_movingwin_estab_1 = lmer(nd ~ stage_ij_estab*t_win+(1|field/plot),
                              data = inter_all_movingwin_forestab, 
                              REML = F)
mod_nd_compare_movingwin_estab_2 = lmer(nd ~ stage_ij_estab*t_win
                                        +(1|field/plot)
                                        +(1|sp_pair) ,
                              data = inter_all_movingwin_forestab, 
                              REML = F)
mod_nd_compare_movingwin_estab = lmer(nd ~ stage_ij_estab*t_win+(1|field/plot)+
                              (1|species_i) +(1|species_j) +
                              (1|sp_pair),
                            data = inter_all_movingwin_forestab, 
                            REML = F)
anova(mod_nd_compare_movingwin_estab_1,
      mod_nd_compare_movingwin_estab) ### select mod_nd_compare_movingwin_estab
anova(mod_nd_compare_movingwin_estab_2,
      mod_nd_compare_movingwin_estab) ### select mod_nd_compare_movingwin_estab
mod_nd_compare_movingwin_estab = lmer(nd ~ stage_ij_estab*t_win+(1|field/plot)+
                                        (1|species_i) +(1|species_j) +
                                        (1|sp_pair),
                                      data = inter_all_movingwin_forestab, 
                                      REML = T)
summary(mod_nd_compare_movingwin_estab)

mod_lgfd_compare_movingwin_estab = lmer(lgfd ~ stage_ij_estab*t_win+(1|field/plot)+
                                        (1|species_i) +(1|species_j) +
                                        (1|sp_pair),
                                      data = inter_all_movingwin_forestab, 
                                      REML = T)

summary(mod_lgfd_compare_movingwin_estab)

### Load inla model results adding the time auto-regression ###  
load("D:/R projects/BSS/results/analysing_results/mod_nd_compare_movingwin_estab.rdata")
load("D:/R projects/BSS/results/analysing_results/mod_fd_compare_movingwin_estab.rdata")

load("D:/R projects/BSS/results/analysing_results/estab.nd.predict.rdata")
load("D:/R projects/BSS/results/analysing_results/estab.fd.predict.rdata")

load("D:/R projects/BSS/results/analysing_results/estab_3.nd.predict.rdata")
load("D:/R projects/BSS/results/analysing_results/estab_3.fd.predict.rdata")

pre_ndfd_compare_estab_movingwin = data.frame(nd = estab_3.nd.predict$mean,
                                lgfd = estab_3.fd.predict$mean,
                                invasion_stage = estab_3.fd.predict$stage_ij_estab,
                                nd_low = estab_3.nd.predict$'0.025quant',
                                nd_high = estab_3.nd.predict$'0.975quant',
                                fd_low = estab_3.fd.predict$'0.025quant',
                                fd_high = estab_3.fd.predict$'0.975quant',
                                time_window = unique(estab_3.fd.predict$t_win))
pre_ndfd_compare_estab_movingwin$time_window = ordered(pre_ndfd_compare_estab_movingwin$time_window)
spss.f = function(x) log10(1-x)
x = seq(-1, 1, 0.001)
spss.ff = function(x) -log10(1-x)
# Define a collection of palettes to alter the default based on number of levels
# to encode Template function for creating densities grouped by a variable
str(pre_ndfd_compare_estab_movingwin)
max(pre_ndfd_compare_estab_movingwin$nd_high)
min(pre_ndfd_compare_estab_movingwin$nd_low)
max(pre_ndfd_compare_estab_movingwin$fd_high)
min(pre_ndfd_compare_estab_movingwin$fd_low)

Fig.4_compare_estab = ggplot(data = pre_ndfd_compare_estab_movingwin,
                             aes(x = nd,y = lgfd,
                                 color = time_window,
                                 shape = invasion_stage)) +
  geom_errorbar(mapping = aes(ymin = fd_low,
                              ymax = fd_high),
                linewidth = 0.5, 
                show.legend = T)+
  geom_errorbarh(mapping = aes(xmin = nd_low,
                               xmax = nd_high), 
                 linewidth = 0.5, 
                 show.legend = T)+
  geom_point(size=20) + 
  scale_shape_manual(values=c(19,15,17))+
  scale_color_viridis(direction = -1, option = 'A', discrete = T)+
  #scale_color_brewer(palette="Set2")+
  scale_x_continuous(limits = c(0, 0.39), expand = c(0, 0)) +
  scale_y_continuous(limits=c(-0.2, 0.2)) +
  #coord_cartesian(xlim = c(0.29, 0.36))+
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.5) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=0.5) +
  xlab("Niche difference") +
  ylab(expression(paste(log[10],
                        "(Fitness difference)",
                        sep=' ')))+
  stat_function(fun=spss.f, colour="black")+
  stat_function(fun=spss.ff, colour="black")+
  #annotate(geom="text", x=0.55, y=0.092, label=str_wrap("Invasion and Coexistence",
  #                                                      width = 20), size=20) +
  #annotate(geom="text", x=0.55, y=-0.092, label=str_wrap("Invasion and Coexistence",
  #                                                       width = 20), size=20) +
  #annotate(geom="text", x=0.1, y=0.092, label=str_wrap("Invasion and exclude residents",
  #                                                     width = 20), size=20) +
  #annotate(geom="text", x=0.1, y=-0.092, label=str_wrap("Residents repel invasion",
  #                                                      width = 20), size=20) +
  guides(colour=guide_legend(title='Age',
                             override.aes = list(size = 10)),
         shape = guide_legend(title='Invasion stages',
                              override.aes = list(size = 10))) +
  theme(panel.background = element_rect(fill = 'white',color = 'black',linewidth = 0.5),
        panel.grid = element_blank(),
        #legend.position = c(0.18, 0.9),
        legend.text = element_text(size=80),
        legend.title = element_text(size=80),
        plot.margin = margin(10,10,10,10),
        plot.background = element_blank(),
        text = element_text(size = 150),
        axis.ticks = element_line(linewidth = 3),
        axis.ticks.length = unit(-5,'lines'),
        axis.title.y = element_text(margin = margin(0,3,0,0),color = '#000000',face = 'bold'),
        axis.title.x = element_text(margin = margin(5,0,0,0),color = '#000000',face = 'bold'),
        axis.text.y = element_text(margin = margin(0,3,0,0),color = '#000000'),
        axis.text.x = element_text(margin = margin(5,0,0,0),color = '#000000'))


#### domin nodomin ####
unique(inter_all_c_movingwin$stage_ij_domin)
inter_all_movingwin_fordomin = inter_all_c_movingwin %>%
  filter(stage_ij_domin %in% c("estab.domin_native",
                               "estab.nodomin_native"))
mod_nd_compare_movingwin_domin_1 = lmer(nd ~ stage_ij_domin*t_win+(1|field/plot),
                                        data = inter_all_movingwin_fordomin, 
                                        REML = F)
mod_nd_compare_movingwin_domin_2 = lmer(nd ~ stage_ij_domin*t_win
                                        +(1|field/plot)
                                        +(1|sp_pair) ,
                                        data = inter_all_movingwin_fordomin, 
                                        REML = F)
mod_nd_compare_movingwin_domin = lmer(nd ~ stage_ij_domin*t_win+(1|field/plot)+
                                        (1|species_i) +(1|species_j) +
                                        (1|sp_pair),
                                      data = inter_all_movingwin_fordomin, 
                                      REML = F)
anova(mod_nd_compare_movingwin_domin_1,
      mod_nd_compare_movingwin_domin) ### select mod_nd_compare_movingwin_domin
anova(mod_nd_compare_movingwin_domin_2,
      mod_nd_compare_movingwin_domin) ### select mod_nd_compare_movingwin_domin
mod_nd_compare_movingwin_domin = lmer(nd ~ stage_ij_domin*t_win+(1|field/plot)+
                                        (1|species_i) +(1|species_j) +
                                        (1|sp_pair),
                                      data = inter_all_movingwin_fordomin, 
                                      REML = T)

summary(mod_nd_compare_movingwin_domin)

mod_lgfd_compare_movingwin_domin = lmer(lgfd ~ stage_ij_domin*t_win+(1|field/plot)+
                                          (1|species_i) +(1|species_j) +
                                          (1|sp_pair),
                                        data = inter_all_movingwin_fordomin, 
                                        REML = T)

summary(mod_lgfd_compare_movingwin_domin)

pre_mod_nd_compare_domin = as.data.frame(ggpredict(mod_nd_compare_movingwin_domin,
                                                   terms = c('stage_ij_domin',
                                                             't_win[1,2,3,4,5,6,7,8,9,10,11,12,13,14]'),
                                                   type = 'fixed'))

pre_mod_fd_compare_domin = as.data.frame(ggpredict(mod_lgfd_compare_movingwin_domin,
                                                   terms = c('stage_ij_domin',
                                                             't_win[1,2,3,4,5,6,7,8,9,10,11,12,13,14]'),
                                                   type = 'fixed'))


load("D:/R projects/BSS/results/analysing_results/mod_nd_compare_movingwin_domin.rdata")
load("D:/R projects/BSS/results/analysing_results/mod_fd_compare_movingwin_domin.rdata")

load("D:/R projects/BSS/results/analysing_results/domin.nd.predict.rdata")
load("D:/R projects/BSS/results/analysing_results/domin.fd.predict.rdata")

load("D:/R projects/BSS/results/analysing_results/domin_3.nd.predict.rdata")
load("D:/R projects/BSS/results/analysing_results/domin_3.fd.predict.rdata")

pre_ndfd_compare_domin_movingwin = data.frame(nd = domin_3.nd.predict$mean,
                                              lgfd = domin_3.fd.predict$mean,
                                              invasion_stage = domin_3.fd.predict$stage_ij_domin,
                                              nd_low = domin_3.nd.predict$'0.025quant',
                                              nd_high = domin_3.nd.predict$'0.975quant',
                                              fd_low = domin_3.fd.predict$'0.025quant',
                                              fd_high = domin_3.fd.predict$'0.975quant',
                                              time_window = unique(domin_3.fd.predict$t_win))
pre_ndfd_compare_domin_movingwin$time_window = ordered(pre_ndfd_compare_domin_movingwin$time_window)
spss.f = function(x) log10(1-x)
x = seq(-1, 1, 0.001)
spss.ff = function(x) -log10(1-x)
# Define a collection of palettes to alter the default based on number of levels
# to encode Template function for creating densities grouped by a variable
str(pre_ndfd_compare_domin_movingwin)
max(pre_ndfd_compare_domin_movingwin$nd_high)
min(pre_ndfd_compare_domin_movingwin$nd_low)
max(pre_ndfd_compare_domin_movingwin$fd_high)
min(pre_ndfd_compare_domin_movingwin$fd_low)

Fig.4_compare_domin = ggplot(data = pre_ndfd_compare_domin_movingwin,
                             aes(x = nd,y = lgfd,
                                 color = time_window,
                                 shape = invasion_stage)) +
  geom_errorbar(mapping = aes(ymin = fd_low,
                              ymax = fd_high),
                linewidth = 0.5, 
                show.legend = T)+
  geom_errorbarh(mapping = aes(xmin = nd_low,
                               xmax = nd_high), 
                 linewidth = 0.5, 
                 show.legend = T)+
  geom_point(size=20) + 
  scale_shape_manual(values=c(19,15,17))+
  scale_color_viridis(direction = -1, option = 'A', discrete = T)+
  #scale_color_brewer(palette="Set2")+
  scale_x_continuous(limits = c(0, 0.39), expand = c(0, 0)) +
  scale_y_continuous(limits=c(-0.2, 0.2)) +
  #coord_cartesian(xlim = c(0.29, 0.36))+
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.5) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=0.5) +
  xlab("Niche difference") +
  ylab(expression(paste(log[10],
                        "(Fitness difference)",
                        sep=' ')))+
  stat_function(fun=spss.f, colour="black")+
  stat_function(fun=spss.ff, colour="black")+
  #annotate(geom="text", x=0.55, y=0.092, label=str_wrap("Invasion and Coexistence",
  #                                                      width = 20), size=20) +
  #annotate(geom="text", x=0.55, y=-0.092, label=str_wrap("Invasion and Coexistence",
  #                                                       width = 20), size=20) +
  #annotate(geom="text", x=0.1, y=0.092, label=str_wrap("Invasion and exclude residents",
  #                                                     width = 20), size=20) +
  #annotate(geom="text", x=0.1, y=-0.092, label=str_wrap("Residents repel invasion",
  #                                                      width = 20), size=20) +
  guides(colour=guide_legend(title='Age',
                             override.aes = list(size = 10)),
         shape = guide_legend(title='Invasion stages',
                              override.aes = list(size = 10))) +
  theme(panel.background = element_rect(fill = 'white',color = 'black',linewidth = 0.5),
        panel.grid = element_blank(),
        #legend.position = c(0.18, 0.9),
        legend.text = element_text(size=80),
        legend.title = element_text(size=80),
        plot.margin = margin(10,10,10,10),
        plot.background = element_blank(),
        text = element_text(size = 150),
        axis.ticks = element_line(linewidth = 3),
        axis.ticks.length = unit(-5,'lines'),
        axis.title.y = element_text(margin = margin(0,3,0,0),color = '#000000',face = 'bold'),
        axis.title.x = element_text(margin = margin(5,0,0,0),color = '#000000',face = 'bold'),
        axis.text.y = element_text(margin = margin(0,3,0,0),color = '#000000'),
        axis.text.x = element_text(margin = margin(5,0,0,0),color = '#000000'))





ggsave('results/figures_movingwin/Fig.2_compare_domin.svg',
       plot=Fig.2_compare_domin, 
       width=105, height=80, dpi=600, units='cm',
       limitsize=F)

gap = ggplot(NULL)+theme_void()
library(ggpubr)

Fig.4_split_stage = ggarrange(gap, Fig.4_compare_estab,gap,
                              Fig.4_compare_domin, nrow = 4, ncol = 1,
                                labels = c('','a)','','b)'), hjust = -0.5,
                                vjust = 0, heights = c(0.1, 1, 0.1, 1),
                                font.label = list(size = 150))

ggsave(plot = Fig.4_split_stage,
       'results/figures_movingwin/Fig.4_split_stage.svg',
       width = 105,height = 170, dpi = 300, units = 'cm',
       limitsize = F)


################ Invasion success ~ nd*t_windows+lgfd*t_windows #########
###### Create dat_suc_sp data, running in linux ######
inter_all_c_movingwin_l = split(inter_all_c_movingwin, inter_all_c_movingwin$t_win)

number_cores = detectCores()
cl = makeCluster(number_cores-2)      
registerDoParallel(cl)
dat_suc_sp_movingwin = foreach(i = 1:14, .packages = c('dplyr', 'data.table')) %dopar% {
  #i = 14
  inter_all_c = inter_all_c_movingwin_l[[i]]
  t_windows = i
  inter_all_c$ra_m_i = as.numeric(inter_all_c$ra_m_i)
  inter_all_c$ra_m_j = as.numeric(inter_all_c$ra_m_j)
  inter_all_c_l_2 = split(inter_all_c, inter_all_c$f_p)
  
  dat_suc_sp_l = lapply(inter_all_c_l_2, function(x){
    #x = inter_all_c_l_2[[398]]
    trans_plot = x
    inv_suc = trans_plot %>% filter(stage_i != 'native' & stage_j == 'native') 
    inv_sp = unique(inv_suc$species_i)
    if (nrow(inv_suc) != 0) {
      dat_suc_sp = data.frame()
      for (j in 1:length(inv_sp)) {
        inv_suc_sp = inv_suc %>% filter(species_i == inv_sp[j])
        mnd = mean(inv_suc_sp$nd)
        mlgfd = mean(inv_suc_sp$lgfd)
        mpd = mean(inv_suc_sp$Phylo_dis)
        mfunc_d = mean(inv_suc_sp$Multi_traits)
        mconti_func_d = mean(inv_suc_sp$Multi_conti_traits)
        mnd.a = sum(inv_suc_sp$nd*inv_suc_sp$ra_m_fit_t_i*inv_suc_sp$ra_m_fit_t_j)/sum(inv_suc_sp$ra_m_fit_t_i*inv_suc_sp$ra_m_fit_t_j)
        mlgfd.a = sum(inv_suc_sp$lgfd*inv_suc_sp$ra_m_fit_t_i*inv_suc_sp$ra_m_fit_t_j)/sum(inv_suc_sp$ra_m_fit_t_i*inv_suc_sp$ra_m_fit_t_j)
        mpd.a = sum(inv_suc_sp$Phylo_dis*inv_suc_sp$ra_m_fit_t_i*inv_suc_sp$ra_m_fit_t_j)/sum(inv_suc_sp$ra_m_fit_t_i*inv_suc_sp$ra_m_fit_t_j)
        mfunc_d.a = sum(inv_suc_sp$Multi_traits*inv_suc_sp$ra_m_fit_t_i*inv_suc_sp$ra_m_fit_t_j)/sum(inv_suc_sp$ra_m_fit_t_i*inv_suc_sp$ra_m_fit_t_j)
        mconti_func_d.a = sum(inv_suc_sp$Multi_conti_traits*inv_suc_sp$ra_m_fit_t_i*inv_suc_sp$ra_m_fit_t_j)/sum(inv_suc_sp$ra_m_fit_t_i*inv_suc_sp$ra_m_fit_t_j)
        mnnd = min(inv_suc_sp$nd)
        mnlgfd = min(inv_suc_sp$lgfd)
        mntd = min(inv_suc_sp$Phylo_dis)
        mnfunc_d = min(inv_suc_sp$Multi_traits)
        mnconti_func_d = min(inv_suc_sp$Multi_conti_traits)
        
        dat_suc_sp_1 = data.frame(f_p = unique(trans_plot$f_p),
                                  plot = unique(trans_plot$f_p),
                                  field = unique(trans_plot$field),
                                  t_win = t_windows,
                                  species = inv_sp[j], 
                                  stage = unique(inv_suc_sp$stage_i),
                                  estab = unique(inv_suc_sp$estab_i),
                                  domin = unique(inv_suc_sp$domin_i),
                                  mnd = mnd, mlgfd = mlgfd, mpd = mpd, mfunc_d = mfunc_d,
                                  mconti_func_d = mconti_func_d,
                                  mnd.a = mnd.a, mlgfd.a = mlgfd.a, mpd.a = mpd.a,
                                  mfunc_d.a = mfunc_d.a, mconti_func_d.a = mconti_func_d.a,
                                  mnnd = mnnd, mnlgfd = mnlgfd, mntd = mntd,
                                  mnfunc_d = mnfunc_d, mnconti_func_d = mnconti_func_d)
        dat_suc_sp = rbind(dat_suc_sp, dat_suc_sp_1)
      } 
      return(dat_suc_sp) 
    } else {
      dat_suc_sp = data.frame()
      return(dat_suc_sp) 
    }
  })
  dat_suc_sp = rbindlist(dat_suc_sp_l)
  return(dat_suc_sp)
  
}
stopCluster(cl)

#### Fast start: load transformed data 
load("D:/R projects/BSS/results/fit_results/dat_suc_sp_movingwin.rdata")
dat_suc_sp_movingwin = rbindlist(dat_suc_sp_movingwin)

dat_suc_sp_movingwin$stage_level = NA
dat_suc_sp_movingwin[dat_suc_sp_movingwin$stage == 'introduce',]$stage_level = 1
dat_suc_sp_movingwin$stage_level = as.numeric(dat_suc_sp_movingwin$stage_level)
dat_suc_sp_movingwin[dat_suc_sp_movingwin$stage == 'establish',]$stage_level = 2
dat_suc_sp_movingwin[dat_suc_sp_movingwin$stage == 'dominant',]$stage_level = 3
dat_suc_sp_movingwin$stage_level = as.factor(dat_suc_sp_movingwin$stage_level)
str(dat_suc_sp_movingwin)

dat_suc_sp_movingwin$ff_p = as.factor(dat_suc_sp_movingwin$f_p)
dat_suc_sp_movingwin$fplot = as.factor(dat_suc_sp_movingwin$plot)
dat_suc_sp_movingwin$ffield = as.factor(dat_suc_sp_movingwin$field)
dat_suc_sp_movingwin$ft_win = as.factor(dat_suc_sp_movingwin$t_win)

library(ordinal)
library(lme4)
library(INLA)
library(brinla)

##### Analyse for mnd mfd abundance weighted mean
### Establishment ~ mnd,mfd,time_windows
pc_prior=list(prec=list("pc.prec",param=c(0.1,0.01)))
gammaprior=list(prec=list(prior="loggamma",param=c(0.01,0.01)))
dat_suc_sp_movingwin$t_win1 = dat_suc_sp_movingwin$t_win

mod_estab_sp_ab = inla(estab ~ mnd.a*t_win + mlgfd.a*t_win
                       + f(species, model="iid", hyper = pc_prior)
                       + f(field, model="iid", hyper = pc_prior)
                       + f(field:plot, model="iid", hyper = pc_prior)
                       + f(t_win1,model="ar1", hyper = pc_prior),
                       control.compute = list(dic=T,waic=T,cpo=T),
                     family="binomial", data=dat_suc_sp_movingwin)
mod_estab_sp = inla(estab ~ mnd*t_win + mlgfd*t_win
                         + f(species, model="iid", hyper = pc_prior)
                       + f(field, model="iid", hyper = pc_prior)
                       + f(field:plot, model="iid", hyper = pc_prior)
                       + f(t_win1,model="ar1", hyper = pc_prior),
                       control.compute = list(dic=T,waic=T,cpo=T),
                       family="binomial", data=dat_suc_sp_movingwin)
mod_estab_sp_mnd = inla(estab ~ mnnd*t_win + mnlgfd*t_win
                         + f(species, model="iid", hyper = pc_prior)
                       + f(field, model="iid", hyper = pc_prior)
                       + f(field:plot, model="iid", hyper = pc_prior)
                       + f(t_win1,model="ar1", hyper = pc_prior),
                       control.compute = list(dic=T,waic=T,cpo=T),
                       family="binomial", data=dat_suc_sp_movingwin)
c(mod_estab_sp$waic$waic, mod_estab_sp_ab$waic$waic, mod_estab_sp_mnd$waic$waic)
mod_estab_sp_1 = inla(estab ~ mnd*t_win + mlgfd*t_win
                         + f(field, model="iid", hyper = pc_prior)
                         + f(field:plot, model="iid", hyper = pc_prior)
                         + f(t_win1,model="ar1", hyper = pc_prior),
                         control.compute = list(config = T,dic=T,waic=T,cpo=T),
                       family="binomial", data=dat_suc_sp_movingwin)

c(mod_estab_sp$waic$waic, mod_estab_sp_1$waic$waic)
mod_estab_sp_1$summary.fixed
summary(mod_estab_sp_1)
mod_estab_sp$summary.fixed
summary(mod_estab_sp)

##### Plot interaction effects #####
tmp = inla.posterior.sample(100, mod_estab_sp_1)
tmp_r = inla.posterior.sample.eval(c('t_win1'), tmp)
tmp_1 = inla.posterior.sample.eval(c('mnd', 'mlgfd',
                                     't_win', "mnd:t_win",
                                     't_win:mlgfd'), tmp)
row.names(tmp_1) = c('mnd', 'mlgfd',
                     't_win', "mnd:t_win",
                     't_win:mlgfd')

# sequences of hypothetical trait values to predict off of
seq.t_win = seq(1, 14, length.out = 14) 

# Create a design matrix based on the equation: 1 = intercept, my.seq = trail values to predict off of
dm.t_win = cbind(1, seq.t_win)  

# Grab the appropriate slope parameters
mcmc.traits_mnd = cbind(tmp_1[1,], tmp_1[4,])
mcmc.traits_mlgfd = cbind(tmp_1[2,], tmp_1[5,])

# Some quick matrix math, multiplying the design matrix by the mcmc slope parameters
pred.t_win_mnd = mcmc.traits_mnd %*% t(dm.t_win)   
pred.t_win.random_mnd = mcmc.traits_mnd %*% t(dm.t_win) + t(tmp_r)

pred.t_win_mlgfd = mcmc.traits_mlgfd %*% t(dm.t_win)   
pred.t_win.random_mlgfd = mcmc.traits_mlgfd %*% t(dm.t_win) + t(tmp_r)

pm1_mnd = apply(pred.t_win.random_mnd, c(2), function(x) mean(x, na.rm=TRUE))   
cri1_mnd = apply(pred.t_win.random_mnd, c(2), function(x) quantile(x, prob = c(0.025, 0.975, 0.05, 0.95)))

data_random_mnd_estab = data.frame(
  "t_win" = dm.t_win[,2],
  "random_m" = pm1_mnd,
  "random_lower95" = c(t(cri1_mnd[3,])),
  "random_upper95" = c(t(cri1_mnd[4,]))
)

pm1_mlgfd = apply(pred.t_win.random_mlgfd, c(2), function(x) mean(x, na.rm=TRUE))   
cri1_mlgfd = apply(pred.t_win.random_mlgfd, c(2), function(x) quantile(x, prob = c(0.025, 0.975, 0.05, 0.95)))

data_random_mlgfd_estab = data.frame(
  "t_win" = dm.t_win[,2],
  "random_m" = pm1_mlgfd,
  "random_lower95" = c(t(cri1_mlgfd[3,])),
  "random_upper95" = c(t(cri1_mlgfd[4,]))
)

# Summarize the predictions. We really only need the median and CRIs
pred.t_win_mnd = apply(pred.t_win_mnd, 2, quantile, probs = c(0.025, 0.5, 0.975)) 

(dat.pred.estab_mnd = data.frame(
  "t_win" = seq.t_win,
  t(pred.t_win_mnd)
))
colnames(dat.pred.estab_mnd)[2:4] = c("lower95", "med", "upper95")

pred.t_win_mlgfd = apply(pred.t_win_mlgfd, 2, quantile, probs = c(0.025, 0.5, 0.975)) 

(dat.pred.estab_mlgfd = data.frame(
  "t_win" = seq.t_win,
  t(pred.t_win_mlgfd)
))
colnames(dat.pred.estab_mlgfd)[2:4] = c("lower95", "med", "upper95")


plot.mnd_effect_estab = ggplot() +
  theme_classic()+ 
  geom_ribbon(data = dat.pred.estab_mnd, aes(x = t_win, y = med, ymin = lower95, ymax = upper95), 
              fill = "grey30", alpha = 0.2)+
  geom_smooth(data = dat.pred.estab_mnd, aes(x = t_win, y = med), 
              se = T, color = "grey30") +
  geom_pointrange(data = data_random_mnd_estab, aes(x = t_win, y = random_m,
                                          ymin = random_lower95, ymax = random_upper95),
                  color = "grey30") +
  labs(x = "Year", y = "Mnd-establishment relationship") +
  scale_x_continuous(breaks = seq.t_win) +
  theme(
    #axis.text.x = element_blank(),
    #axis.text.y = element_blank(),
    #axis.ticks = element_blank(),
    axis.text.x = element_text(face = "bold", size = 14), 
    axis.text.y = element_text(face = "bold", size = 14), 
    axis.title.x = element_text(face = "bold", size = 18), 
    axis.title.y = element_text(face = "bold", size = 18))
plot.mnd_effect_estab

plot.mlgfd_effect_estab = ggplot() +
  theme_classic()+ 
  geom_ribbon(data = dat.pred.estab_mlgfd, aes(x = t_win, y = med, ymin = lower95, ymax = upper95), 
              fill = "grey30", alpha = 0.2)+
  geom_smooth(data = dat.pred.estab_mlgfd, aes(x = t_win, y = med), 
              se = T, color = "grey30") +
  geom_pointrange(data = data_random_mlgfd_estab, aes(x = t_win, y = random_m,
                                          ymin = random_lower95, ymax = random_upper95),
                  color = "grey30") +
  labs(x = "Year", y = "Mlgfd-establishment relationship") +
  scale_x_continuous(breaks = seq.t_win) +
  theme(
    #axis.text.x = element_blank(),
    #axis.text.y = element_blank(),
    #axis.ticks = element_blank(),
    axis.text.x = element_text(face = "bold", size = 14), 
    axis.text.y = element_text(face = "bold", size = 14), 
    axis.title.x = element_text(face = "bold", size = 18), 
    axis.title.y = element_text(face = "bold", size = 18))
plot.mlgfd_effect_estab

############ Dominant ~ mnd,mfd,time_windows ############
## tranform data
dat_dom_sp_movingwin = dat_suc_sp_movingwin %>% filter(stage %in% c('establish',
                                                                    'dominant'))
pc_prior=list(prec=list("pc.prec",param=c(0.1,0.01)))
gammaprior=list(prec=list(prior="loggamma",param=c(0.01,0.01)))
dat_dom_sp_movingwin$t_win1 = dat_dom_sp_movingwin$t_win
mod_domin_sp_ab = inla(domin ~ mnd.a*t_win + mlgfd.a*t_win
                       + f(species, model="iid", hyper = pc_prior)
                       + f(field, model="iid", hyper = pc_prior)
                       + f(field:plot, model="iid", hyper = pc_prior)
                       + f(t_win1,model="ar1", hyper = pc_prior),
                       control.compute = list(dic=T,waic=T,cpo=T),
                       family="binomial", data=dat_dom_sp_movingwin)
mod_domin_sp = inla(domin ~ mnd*t_win + mlgfd*t_win
                    + f(species, model="iid", hyper = pc_prior)
                    + f(field, model="iid", hyper = pc_prior)
                    + f(field:plot, model="iid", hyper = pc_prior)
                    + f(t_win1,model="ar1", hyper = pc_prior),
                    control.compute = list(dic=T,waic=T,cpo=T),
                    family="binomial", data=dat_dom_sp_movingwin)
mod_domin_sp_mnd = inla(domin ~ mnnd*t_win + mnlgfd*t_win
                        + f(species, model="iid", hyper = pc_prior)
                        + f(field, model="iid", hyper = pc_prior)
                        + f(field:plot, model="iid", hyper = pc_prior)
                        + f(t_win1,model="ar1", hyper = pc_prior),
                        control.compute = list(dic=T,waic=T,cpo=T),
                        family="binomial", data=dat_dom_sp_movingwin)
c(mod_domin_sp$waic$waic, mod_domin_sp_ab$waic$waic, mod_domin_sp_mnd$waic$waic)

mod_domin_sp_ab_1 = inla(domin ~ mnd.a*t_win + mlgfd.a*t_win 
                         + f(field, model="iid", hyper = pc_prior)
                         + f(field:plot, model="iid", hyper = pc_prior) 
                         + f(t_win1,model="ar1", hyper = pc_prior),
                         control.compute = list(config = T,dic=T,waic=T,cpo=T),
                         family="binomial", data=dat_dom_sp_movingwin)
summary(mod_domin_sp_ab_1)
c(mod_domin_sp_ab$waic$waic, mod_domin_sp_ab_1$waic$waic)

mod_domin_sp_ab$summary.fixed
summary(mod_domin_sp_ab)

##### Plot interaction effects #####
tmp = inla.posterior.sample(100, mod_domin_sp_ab_1)
tmp_r = inla.posterior.sample.eval(c('t_win1'), tmp)
tmp_1 = inla.posterior.sample.eval(c('mnd.a', 'mlgfd.a',
                                     't_win', "mnd.a:t_win",
                                     't_win:mlgfd.a'), tmp)
row.names(tmp_1) = c('mnd.a', 'mlgfd.a',
                     't_win', "mnd.a:t_win",
                     't_win:mlgfd.a')

# sequences of hypothetical trait values to predict off of
seq.t_win = seq(1, 14, length.out = 14) 

# Create a design matrix based on the equation: 1 = intercept, my.seq = trail values to predict off of
dm.t_win = cbind(1, seq.t_win)  

# Grab the appropriate slope parameters
mcmc.traits_mnd = cbind(tmp_1[1,], tmp_1[4,])
mcmc.traits_mlgfd = cbind(tmp_1[2,], tmp_1[5,])

# Some quick matrix math, multiplying the design matrix by the mcmc slope parameters
pred.t_win_mnd = mcmc.traits_mnd %*% t(dm.t_win)   
pred.t_win.random_mnd = mcmc.traits_mnd %*% t(dm.t_win) + t(tmp_r)

pred.t_win_mlgfd = mcmc.traits_mlgfd %*% t(dm.t_win)   
pred.t_win.random_mlgfd = mcmc.traits_mlgfd %*% t(dm.t_win) + t(tmp_r)

pm1_mnd = apply(pred.t_win.random_mnd, c(2), function(x) mean(x, na.rm=TRUE))   
cri1_mnd = apply(pred.t_win.random_mnd, c(2), function(x) quantile(x, prob = c(0.025, 0.975, 0.05, 0.95)))

data_random_mnd_domin = data.frame(
  "t_win" = dm.t_win[,2],
  "random_m" = pm1_mnd,
  "random_lower95" = c(t(cri1_mnd[3,])),
  "random_upper95" = c(t(cri1_mnd[4,]))
)

pm1_mlgfd = apply(pred.t_win.random_mlgfd, c(2), function(x) mean(x, na.rm=TRUE))   
cri1_mlgfd = apply(pred.t_win.random_mlgfd, c(2), function(x) quantile(x, prob = c(0.025, 0.975, 0.05, 0.95)))

data_random_mlgfd_domin = data.frame(
  "t_win" = dm.t_win[,2],
  "random_m" = pm1_mlgfd,
  "random_lower95" = c(t(cri1_mlgfd[3,])),
  "random_upper95" = c(t(cri1_mlgfd[4,]))
)

# Summarize the predictions. We really only need the median and CRIs
pred.t_win_mnd = apply(pred.t_win_mnd, 2, quantile, probs = c(0.025, 0.5, 0.975)) 

(dat.pred.domin_mnd = data.frame(
  "t_win" = seq.t_win,
  t(pred.t_win_mnd)
))
colnames(dat.pred.domin_mnd)[2:4] = c("lower95", "med", "upper95")

pred.t_win_mlgfd = apply(pred.t_win_mlgfd, 2, quantile, probs = c(0.025, 0.5, 0.975)) 

(dat.pred.domin_mlgfd = data.frame(
  "t_win" = seq.t_win,
  t(pred.t_win_mlgfd)
))
colnames(dat.pred.domin_mlgfd)[2:4] = c("lower95", "med", "upper95")


plot.mnd_effect_domin = ggplot() +
  theme_classic()+ 
  geom_ribbon(data = dat.pred.domin_mnd, aes(x = t_win, y = med, ymin = lower95, ymax = upper95), 
              fill = "grey30", alpha = 0.2)+
  geom_smooth(data = dat.pred.domin_mnd, aes(x = t_win, y = med), 
              se = F, color = "grey30") +
  geom_pointrange(data = data_random_mnd_domin, aes(x = t_win, y = random_m,
                                                    ymin = random_lower95, ymax = random_upper95),
                  color = "grey30") +
  labs(x = "Year", y = "Mnd.a-dominance relationship") +
  scale_x_continuous(breaks = seq.t_win) +
  theme(
    #axis.text.x = element_blank(),
    #axis.text.y = element_blank(),
    #axis.ticks = element_blank(),
    axis.text.x = element_text(face = "bold", size = 14), 
    axis.text.y = element_text(face = "bold", size = 14), 
    axis.title.x = element_text(face = "bold", size = 18), 
    axis.title.y = element_text(face = "bold", size = 18))
plot.mnd_effect_domin

plot.mlgfd_effect_domin = ggplot() +
  theme_classic()+ 
  geom_ribbon(data = dat.pred.domin_mlgfd, aes(x = t_win, y = med, ymin = lower95, ymax = upper95), 
              fill = "grey30", alpha = 0.2)+
  geom_smooth(data = dat.pred.domin_mlgfd, aes(x = t_win, y = med), 
              se = F, color = "grey30") +
  geom_pointrange(data = data_random_mlgfd_domin, aes(x = t_win, y = random_m,
                                                      ymin = random_lower95, ymax = random_upper95),
                  color = "grey30") +
  labs(x = "Year", y = "Mlgfd.a-dominlishment relationship") +
  scale_x_continuous(breaks = seq.t_win) +
  theme(
    #axis.text.x = element_blank(),
    #axis.text.y = element_blank(),
    #axis.ticks = element_blank(),
    axis.text.x = element_text(face = "bold", size = 14), 
    axis.text.y = element_text(face = "bold", size = 14), 
    axis.title.x = element_text(face = "bold", size = 18), 
    axis.title.y = element_text(face = "bold", size = 18))
plot.mlgfd_effect_domin

####### Merge plot for effect ~ mnd mfd at different time windows
gap = ggplot(NULL)+theme_void()
library(ggpubr)
Fig.1_split_stage = ggarrange(gap,Fig.1_compare_estab,gap,
                              Fig.1_compare_domin, nrow = 1, ncol = 4,
                              labels = c('','a)','','b)'), hjust = 1,
                              vjust = 1.5, widths = c(0.1, 1, 0.1, 1),
                              font.label = list(size = 150))

ggsave(plot = Fig.1_split_stage,
       'results/figures_alltime/Fig.1_split_stage.svg',
       width = 320,height = 100, dpi = 300, units = 'cm',
       limitsize = F)

library(ggpubr)
success.mnd.mlgfd.effects = ggarrange(gap, plot.mnd_effect_estab,
                                      gap, plot.mlgfd_effect_estab,
                                    gap, plot.mnd_effect_domin,
                                    gap, plot.mlgfd_effect_domin,
                                    nrow = 2, ncol = 4,
                                    hjust = 0.5,
                                    vjust = 1.5,
                                    labels = list('','a)', '', 'b)',
                                                  '', 'c)', '', 'd)'),
                                    widths = c(0.1, 1, 0.1, 1,
                                               0.1, 1, 0.1, 1),
                                    font.label = list(size = 80,
                                                      face = "bold",
                                                      color ="black"))

ggsave('results/figures_movingwin/success.mnd.mlgfd.effects.svg',
       width = 110, height = 100, units = c('cm'),
       dpi = 300, limitsize = F, plot = success.mnd.mlgfd.effects)
ggsave('results/figures_movingwin/success.mnd.mlgfd.effects_2.svg',
       width = 55, height = 55, units = c('cm'),
       dpi = 300, limitsize = F, plot = success.mnd.mlgfd.effects)
ggsave('results/figures_movingwin/success.mnd.mlgfd.effects_3.svg',
       width = 55, height = 55, units = c('cm'),
       dpi = 300, limitsize = F, plot = success.mnd.mlgfd.effects)
ggsave('results/figures_movingwin/success.mnd.mlgfd.effects.png',
       width = 110, height = 110, units = c('cm'),
       dpi = 150, limitsize = F, plot = success.mnd.mlgfd.effects)
ggsave('results/figures_movingwin/success.mnd.mlgfd.effects.pdf',
       width = 110, height = 110, units = c('cm'),
       dpi = 150, limitsize = F, plot = success.mnd.mlgfd.effects,
       device = 'pdf')



################ Invasion success~nd*t_windows+lgfd*t_windows+pd*t_windows+func_d*t_windows #########

#### load transformed data
load("D:/R projects/BSS/results/fit_results/dat_suc_sp_movingwin.rdata")
dat_suc_sp_movingwin = data.frame(rbindlist(dat_suc_sp_movingwin))
dat_suc_sp_movingwin$initial_age = dat_suc_sp_movingwin$t_win+(dat_suc_sp_movingwin$t_win-1)
dat_suc_sp_movingwin$stage_level = NA
dat_suc_sp_movingwin[dat_suc_sp_movingwin$stage == 'introduce',]$stage_level = 1
dat_suc_sp_movingwin$stage_level = as.numeric(dat_suc_sp_movingwin$stage_level)
dat_suc_sp_movingwin[dat_suc_sp_movingwin$stage == 'establish',]$stage_level = 2
dat_suc_sp_movingwin[dat_suc_sp_movingwin$stage == 'dominant',]$stage_level = 3
dat_suc_sp_movingwin$stage_level = as.factor(dat_suc_sp_movingwin$stage_level)
str(dat_suc_sp_movingwin)

dat_suc_sp_movingwin$ff_p = as.factor(dat_suc_sp_movingwin$f_p)
dat_suc_sp_movingwin$fplot = as.factor(dat_suc_sp_movingwin$plot)
dat_suc_sp_movingwin$ffield = as.factor(dat_suc_sp_movingwin$field)
dat_suc_sp_movingwin$finitial_age = as.factor(dat_suc_sp_movingwin$initial_age)

library(ordinal)
library(lme4)
library(INLA)
library(brinla)
library(ape)

##### Analyse for mnd mfd abundance weighted mean
### Establishment ~ nd*initial_agedows+lgfd*initial_agedows+pd*initial_agedows+func_d*initial_agedows ####
pc_prior=list(prec=list("pc.prec",param=c(0.1,0.01)))
gammaprior=list(prec=list(prior="loggamma",param=c(0.01,0.01)))
dat_suc_sp_movingwin$initial_age1 = dat_suc_sp_movingwin$initial_age

## Rescale data
numcols = grep("^m",names(dat_suc_sp_movingwin))
dat_suc_sp_movingwins = dat_suc_sp_movingwin
dat_suc_sp_movingwins[,numcols] = scale(dat_suc_sp_movingwins[,numcols])
dat_suc_sp_movingwins$initial_age = scale(dat_suc_sp_movingwins$initial_age)

dat_suc_sp_movingwins$species_1 = as.factor(dat_suc_sp_movingwins$species)
estab_sp_names = unique(dat_suc_sp_movingwins$species)

tree = read.tree('data/original data/phylo_tree332.txt')
estab_tree_fit = keep.tip(tree, estab_sp_names)
estab_vcv_tree = ape::vcv(estab_tree_fit, model = "Brownian", corr = FALSE)
estab_vcv_tree_sparse = inla.as.sparse(solve(estab_vcv_tree))

mod_estab_sp_ab = inla(estab ~ mnd.a*initial_age+mlgfd.a*initial_age+mpd.a*initial_age+mconti_func_d.a*initial_age+ 
                        # f(species, model="iid", hyper = pc_prior) + 
                       f(field, model="iid", hyper = pc_prior)
                       + f(field:plot, model="iid", hyper = pc_prior)
                       + f(initial_age1,model="ar1", hyper = pc_prior)+
                         f(species_1, model="generic0",
                           Cmatrix= estab_vcv_tree_sparse,
                           values = estab_sp_names, hyper=pc_prior),
                       control.compute = list(dic=T,waic=T,cpo=T),
                       family="binomial", data=dat_suc_sp_movingwins)

mod_estab_sp = inla(estab ~ mnd*initial_age+mlgfd*initial_age+mpd*initial_age+mconti_func_d*initial_age+ 
                       #f(species, model="iid", hyper = pc_prior)+ 
                       f(field, model="iid", hyper = pc_prior)
                       + f(field:plot, model="iid", hyper = pc_prior)+
                         f(initial_age1, model="ar1", hyper = pc_prior) +
                         f(species_1, model="generic0",
                        Cmatrix= estab_vcv_tree_sparse,
                        values = estab_sp_names, hyper=pc_prior),
                       control.compute = list(config = T,dic=T,waic=T,cpo=T),
                       family="binomial", data=dat_suc_sp_movingwins)

mod_estab_sp_mnd = inla(estab ~ mnnd*initial_age+mnlgfd*initial_age+mntd*initial_age+mnconti_func_d*initial_age+
                       # f(species, model="iid", hyper = pc_prior)+
                        f(field, model="iid", hyper = pc_prior)
                       + f(field:plot, model="iid", hyper = pc_prior)
                       + f(initial_age1,model="ar1", hyper = pc_prior)+
                         f(species_1, model="generic0",
                           Cmatrix= estab_vcv_tree_sparse,
                           values = estab_sp_names, hyper=pc_prior),
                       control.compute = list(dic=T,waic=T,cpo=T),
                       family="binomial", data=dat_suc_sp_movingwins)
c(mod_estab_sp$waic$waic, mod_estab_sp_ab$waic$waic, mod_estab_sp_mnd$waic$waic) # Chose mpd, the lowest waic
summary(mod_estab_sp)
summary(mod_estab_sp_ab)
summary(mod_estab_sp_mnd)

##### Calculate relative importance of ND FD PD FUNC_D along the time gradient #####



######### Plot interaction effects #########
tmp = inla.posterior.sample(100, mod_estab_sp)
tmp_r = inla.posterior.sample.eval(c('initial_age1'), tmp)
tmp_1 = inla.posterior.sample.eval(c('mnd', 'mlgfd',
                                     'initial_age', "mnd:initial_age",
                                     'initial_age:mlgfd'), tmp)
tmp_all = inla.posterior.sample.eval(c('mnd', 'mlgfd', 'mpd', 'mconti_func_d',
                                     'initial_age', "mnd:initial_age",
                                     'initial_age:mlgfd', "initial_age:mpd", "initial_age:mconti_func_d"), tmp)
row.names(tmp_all) = c('mnd', 'mlgfd', 'mpd', 'mconti_func_d',
                       'initial_age', "mnd:initial_age",
                       'initial_age:mlgfd', "initial_age:mpd", "initial_age:mconti_func_d")
row.names(tmp_1) = c('mnd', 'mlgfd',
                     'initial_age', "mnd:initial_age",
                     'initial_age:mlgfd')

# sequences of hypothetical trait values to predict off of
plot.initial_age = seq(1, 27, length.out = 14) 
seq.initial_age = seq(min(dat_suc_sp_movingwins$initial_age), max(dat_suc_sp_movingwins$initial_age),
                length.out = 14)

# Create a design matrix based on the equation: 1 = intercept, my.seq = trail values to predict off of
dm.initial_age = cbind(1, seq.initial_age)  

# Grab the appropriate slope parameters
mcmc.traits_mnd = cbind(tmp_1[1,], tmp_1[4,])
mcmc.traits_mlgfd = cbind(tmp_1[2,], tmp_1[5,])
mcmc.traits_mpd = cbind(tmp_all[3,], tmp_all[8,])
mcmc.traits_mconti_func_d = cbind(tmp_all[4,], tmp_all[9,])


# Some quick matrix math, multiplying the design matrix by the mcmc slope parameters
pred.initial_age_mnd = mcmc.traits_mnd %*% t(dm.initial_age)   
pred.initial_age.random_mnd = mcmc.traits_mnd %*% t(dm.initial_age) + t(tmp_r)

pred.initial_age_mlgfd = mcmc.traits_mlgfd %*% t(dm.initial_age)   
pred.initial_age.random_mlgfd = mcmc.traits_mlgfd %*% t(dm.initial_age) + t(tmp_r)

pred.initial_age_mpd = mcmc.traits_mpd %*% t(dm.initial_age)   
pred.initial_age.random_mpd = mcmc.traits_mpd %*% t(dm.initial_age) + t(tmp_r)

pred.initial_age_mconti_func_d = mcmc.traits_mconti_func_d %*% t(dm.initial_age)   
pred.initial_age.random_mconti_func_d = mcmc.traits_mconti_func_d %*% t(dm.initial_age) + t(tmp_r)

pm1_mnd = apply(pred.initial_age.random_mnd, c(2), function(x) mean(x, na.rm=TRUE))   
cri1_mnd = apply(pred.initial_age.random_mnd, c(2), function(x) quantile(x, prob = c(0.025, 0.975, 0.05, 0.95)))

data_random_mnd_estab = data.frame(
  "initial_age" = plot.initial_age,
  "random_m" = pm1_mnd,
  "random_lower95" = c(t(cri1_mnd[3,])),
  "random_upper95" = c(t(cri1_mnd[4,]))
)

pm1_mlgfd = apply(pred.initial_age.random_mlgfd, c(2), function(x) mean(x, na.rm=TRUE))   
cri1_mlgfd = apply(pred.initial_age.random_mlgfd, c(2), function(x) quantile(x, prob = c(0.025, 0.975, 0.05, 0.95)))

data_random_mlgfd_estab = data.frame(
  "initial_age" = plot.initial_age,
  "random_m" = pm1_mlgfd,
  "random_lower95" = c(t(cri1_mlgfd[3,])),
  "random_upper95" = c(t(cri1_mlgfd[4,]))
)

# Summarize the predictions. We really only need the median and CRIs
(dat.pred.estab_mnd = data.frame(
  "initial_age" = plot.initial_age,
  t(apply(pred.initial_age_mnd, 2, quantile, probs = c(0.025, 0.5, 0.975))) ,
  mean = apply(pred.initial_age_mnd, 2, mean) 
))
colnames(dat.pred.estab_mnd)[2:5] = c("lower95", "med", "upper95", "mean")

(dat.pred.estab_mlgfd = data.frame(
  "initial_age" = plot.initial_age,
  t(apply(pred.initial_age_mlgfd, 2, quantile, probs = c(0.025, 0.5, 0.975))),
  mean = apply(pred.initial_age_mlgfd, 2, mean) 
))
colnames(dat.pred.estab_mlgfd)[2:5] = c("lower95", "med", "upper95", "mean")

(dat.pred.estab_mpd = data.frame(
  "initial_age" = plot.initial_age,
  t(apply(pred.initial_age_mpd, 2, quantile, probs = c(0.025, 0.5, 0.975))),
  mean = apply(pred.initial_age_mpd, 2, mean)
))
colnames(dat.pred.estab_mpd)[2:5] = c("lower95", "med", "upper95", "mean")

pred.initial_age_mconti_func_d = apply(pred.initial_age_mconti_func_d, 2, quantile, probs = c(0.025, 0.5, 0.975)) 

(dat.pred.estab_mconti_func_d = data.frame(
  "initial_age" = plot.initial_age,
  t(apply(pred.initial_age_mconti_func_d, 2, quantile, probs = c(0.025, 0.5, 0.975))),
  mean = apply(pred.initial_age_mconti_func_d, 2, mean)
))
colnames(dat.pred.estab_mconti_func_d)[2:4] = c("lower95", "med", "upper95")

dat_pred.estab_all = as.data.frame(cbind(initial_age = dat.pred.estab_mnd$initial_age,
                           mnd = dat.pred.estab_mnd$mean,
                           mlgfd = dat.pred.estab_mlgfd$mean,
                           mpd = dat.pred.estab_mpd$mean,
                           mconti_func_d = dat.pred.estab_mconti_func_d$mean))
dat_pred.estab_all$sum = NA
dat_pred.estab_all$sum = rowSums(abs(dat_pred.estab_all[, c(2:5)]))
R2_per = as.data.frame(t(apply(dat_pred.estab_all[, c(2:6)], 1, function(x){abs(x)/x[5]})))
R2_per$sum = plot.initial_age
colnames(R2_per)[5] = 'initial_age'

library(tidyr)
R2_per_1 = R2_per %>% pivot_longer(cols = mnd:mconti_func_d, values_to = "R2")
(plot_R2_per_estab = ggplot()+ 
              theme_classic()+ 
              geom_col(data = R2_per_1, aes(x = initial_age, y = R2*100, fill = name))+
              scale_fill_viridis_d()+
              labs(x = "Initial age", y = expression(paste('Establishment-',"R" ^ "2", "%"))) +
              scale_x_continuous(breaks = plot.initial_age) +
              guides(colour=guide_legend(title=NULL,
                               override.aes = list(size = 8)))+
              theme(
                legend.title = element_blank(),
                legend.text = element_text(size = 30),
                axis.ticks = element_line(linewidth = 0.5),
                axis.ticks.length = unit(-0.4,'lines'),
                axis.text.x = element_text(face = "bold", size = 28), 
                axis.text.y = element_text(face = "bold", size = 28), 
                axis.title.x = element_text(face = "bold", size = 36), 
                axis.title.y = element_text(face = "bold", size = 36)
              ))

plot.mnd_effect_estab = ggplot() +
  theme_classic()+
  geom_ribbon(data = dat.pred.estab_mnd, aes(x = initial_age, y = med, ymin = lower95, ymax = upper95), 
              fill = "grey30", alpha = 0.2, size = 2)+
  geom_smooth(data = dat.pred.estab_mnd, aes(x = initial_age, y = med), 
              se = T, color = "grey30", linewidth = 4) +
  geom_pointrange(data = data_random_mnd_estab, aes(x = initial_age, y = random_m,
                                                    ymin = random_lower95, ymax = random_upper95),
                  color = "grey30", size = 2) +
  labs(x = "Initial age", y = "Mnd-establishment relationship") +
  scale_x_continuous(breaks = plot.initial_age) +
  theme(
    #axis.text.x = element_blank(),
    #axis.text.y = element_blank(),
    #axis.ticks = element_blank(),
    axis.ticks = element_line(linewidth = 0.5),
    axis.ticks.length = unit(-0.4,'lines'),
    axis.text.x = element_text(face = "bold", size = 28), 
    axis.text.y = element_text(face = "bold", size = 28), 
    axis.title.x = element_text(face = "bold", size = 36), 
    axis.title.y = element_text(face = "bold", size = 36)
  )
plot.mnd_effect_estab

plot.mlgfd_effect_estab = ggplot() +
  theme_classic()+
  geom_ribbon(data = dat.pred.estab_mlgfd, aes(x = initial_age, y = med, ymin = lower95, ymax = upper95), 
              fill = "grey30", alpha = 0.2, size = 2)+
  geom_smooth(data = dat.pred.estab_mlgfd, aes(x = initial_age, y = med), 
              se = T, color = "grey30", linewidth = 4) +
  geom_pointrange(data = data_random_mlgfd_estab, aes(x = initial_age, y = random_m,
                                                      ymin = random_lower95, ymax = random_upper95),
                  color = "grey30", size = 2) +
  labs(x = "Initial age", y = "Mlgfd-establishment relationship") +
  scale_x_continuous(breaks = plot.initial_age) +
  theme(
    #axis.text.x = element_blank(),
    #axis.text.y = element_blank(),
    #axis.ticks = element_blank(),
    axis.ticks = element_line(linewidth = 0.5),
    axis.ticks.length = unit(-0.4,'lines'),
    axis.text.x = element_text(face = "bold", size = 28), 
    axis.text.y = element_text(face = "bold", size = 28), 
    axis.title.x = element_text(face = "bold", size = 36), 
    axis.title.y = element_text(face = "bold", size = 36)
  )
plot.mlgfd_effect_estab


############ Dominant ~ nd*initial_agedows+lgfd*initial_agedows+pd*initial_agedows+func_d*initial_agedows ############
## tranform data
dat_dom_sp_movingwin = dat_suc_sp_movingwin %>% filter(stage %in% c('establish',
                                                                    'dominant'))
pc_prior=list(prec=list("pc.prec",param=c(0.1,0.01)))
gammaprior=list(prec=list(prior="loggamma",param=c(0.01,0.01)))
dat_dom_sp_movingwin$initial_age1 = dat_dom_sp_movingwin$initial_age

## Rescale data
numcols = grep("^m",names(dat_dom_sp_movingwin))
dat_dom_sp_movingwins = dat_dom_sp_movingwin
dat_dom_sp_movingwins[,numcols] = scale(dat_dom_sp_movingwins[,numcols])
dat_dom_sp_movingwins$initial_age = scale(dat_dom_sp_movingwins$initial_age)

dat_dom_sp_movingwins$species_1 = as.factor(dat_dom_sp_movingwins$species)
domin_sp_names = unique(dat_dom_sp_movingwins$species)

tree = read.tree('data/original data/phylo_tree332.txt')
domin_tree_fit = keep.tip(tree, domin_sp_names)
domin_vcv_tree = ape::vcv(domin_tree_fit, model = "Brownian", corr = FALSE)
domin_vcv_tree_sparse = inla.as.sparse(solve(domin_vcv_tree))

mod_domin_sp_ab = inla(domin ~ mnd.a*initial_age+mlgfd.a*initial_age+mpd.a*initial_age+mconti_func_d.a*initial_age+ 
                       #+ f(species, model="iid", hyper = pc_prior)+ 
                       f(field, model="iid", hyper = pc_prior)
                       + f(field:plot, model="iid", hyper = pc_prior)
                       + f(initial_age1,model="ar1", hyper = pc_prior)+
                         f(species_1, model="generic0",
                           Cmatrix= domin_vcv_tree_sparse,
                           values = domin_sp_names, hyper=pc_prior),
                       control.compute = list(config = T,dic=T,waic=T,cpo=T),
                       family="binomial", data=dat_dom_sp_movingwins)
mod_domin_sp = inla(domin ~ mnd*initial_age+mlgfd*initial_age+mpd.a*initial_age+mconti_func_d*initial_age+ 
                    #+ f(species, model="iid", hyper = pc_prior)+ 
                    f(field, model="iid", hyper = pc_prior)
                    + f(field:plot, model="iid", hyper = pc_prior)
                    + f(initial_age1,model="ar1", hyper = pc_prior)+
                      f(species_1, model="generic0",
                        Cmatrix= domin_vcv_tree_sparse,
                        values = domin_sp_names, hyper=pc_prior),
                    control.compute = list(dic=T,waic=T,cpo=T),
                    family="binomial", data=dat_dom_sp_movingwins)
mod_domin_sp_mnd = inla(domin ~ mnnd*initial_age+mnlgfd*initial_age+mntd*initial_age+mnconti_func_d*initial_age+ 
                       # + f(species, model="iid", hyper = pc_prior) + 
                        f(field, model="iid", hyper = pc_prior)
                        + f(field:plot, model="iid", hyper = pc_prior)
                        + f(initial_age1,model="ar1", hyper = pc_prior)+
                         f(species_1, model="generic0",
                           Cmatrix= domin_vcv_tree_sparse,
                           values = domin_sp_names, hyper=pc_prior),
                        control.compute = list(dic=T,waic=T,cpo=T),
                        family="binomial", data=dat_dom_sp_movingwins)
c(mod_domin_sp$waic$waic, mod_domin_sp_ab$waic$waic, mod_domin_sp_mnd$waic$waic)
## Chose mpd.a: the lowest waic
summary(mod_domin_sp_ab)

######### Plot interaction effects #########
tmp = inla.posterior.sample(100, mod_domin_sp_ab)
tmp_r = inla.posterior.sample.eval(c('initial_age1'), tmp)
tmp_1 = inla.posterior.sample.eval(c('mnd.a', 'mlgfd.a',
                                     'initial_age', "mnd.a:initial_age",
                                     'initial_age:mlgfd.a'), tmp)
tmp_all = inla.posterior.sample.eval(c('mnd.a', 'mlgfd.a', 'mpd.a', 'mconti_func_d.a',
                                       'initial_age', "mnd.a:initial_age",
                                       'initial_age:mlgfd.a', "initial_age:mpd.a", "initial_age:mconti_func_d.a"), tmp)
row.names(tmp_all) = c('mnd.a', 'mlgfd.a', 'mpd.a', 'mconti_func_d.a',
                       'initial_age', "mnd.a:initial_age",
                       'initial_age:mlgfd.a', "initial_age:mpd.a", "initial_age:mconti_func_d.a")
row.names(tmp_1) = c('mnd.a', 'mlgfd.a',
                     'initial_age', "mnd.a:initial_age",
                     'initial_age:mlgfd.a')

# sequences of hypothetical trait values to predict off of
plot.initial_age = seq(1, 27, length.out = 14) 
seq.initial_age = seq(min(dat_dom_sp_movingwins$initial_age), max(dat_dom_sp_movingwins$initial_age),
                length.out = 14)

# Create a design matrix based on the equation: 1 = intercept, my.seq = trail values to predict off of
dm.initial_age = cbind(1, seq.initial_age)  

# Grab the appropriate slope parameters
mcmc.traits_mnd = cbind(tmp_1[1,], tmp_1[4,])
mcmc.traits_mlgfd = cbind(tmp_1[2,], tmp_1[5,])
mcmc.traits_mpd = cbind(tmp_all[3,], tmp_all[8,])
mcmc.traits_mconti_func_d = cbind(tmp_all[4,], tmp_all[9,])


# Some quick matrix math, multiplying the design matrix by the mcmc slope parameters
pred.initial_age_mnd = mcmc.traits_mnd %*% t(dm.initial_age)   
pred.initial_age.random_mnd = mcmc.traits_mnd %*% t(dm.initial_age) + t(tmp_r)

pred.initial_age_mlgfd = mcmc.traits_mlgfd %*% t(dm.initial_age)   
pred.initial_age.random_mlgfd = mcmc.traits_mlgfd %*% t(dm.initial_age) + t(tmp_r)

pred.initial_age_mpd = mcmc.traits_mpd %*% t(dm.initial_age)   
pred.initial_age.random_mpd = mcmc.traits_mpd %*% t(dm.initial_age) + t(tmp_r)

pred.initial_age_mconti_func_d = mcmc.traits_mconti_func_d %*% t(dm.initial_age)   
pred.initial_age.random_mconti_func_d = mcmc.traits_mconti_func_d %*% t(dm.initial_age) + t(tmp_r)

pm1_mnd = apply(pred.initial_age.random_mnd, c(2), function(x) mean(x, na.rm=TRUE))   
cri1_mnd = apply(pred.initial_age.random_mnd, c(2), function(x) quantile(x, prob = c(0.025, 0.975, 0.05, 0.95)))

data_random_mnd_domin = data.frame(
  "initial_age" = plot.initial_age,
  "random_m" = pm1_mnd,
  "random_lower95" = c(t(cri1_mnd[3,])),
  "random_upper95" = c(t(cri1_mnd[4,]))
)

pm1_mlgfd = apply(pred.initial_age.random_mlgfd, c(2), function(x) mean(x, na.rm=TRUE))   
cri1_mlgfd = apply(pred.initial_age.random_mlgfd, c(2), function(x) quantile(x, prob = c(0.025, 0.975, 0.05, 0.95)))

data_random_mlgfd_domin = data.frame(
  "initial_age" = plot.initial_age,
  "random_m" = pm1_mlgfd,
  "random_lower95" = c(t(cri1_mlgfd[3,])),
  "random_upper95" = c(t(cri1_mlgfd[4,]))
)

# Summarize the predictions. We really only need the median and CRIs
(dat.pred.domin_mnd = data.frame(
  "initial_age" = plot.initial_age,
  t(apply(pred.initial_age_mnd, 2, quantile, probs = c(0.025, 0.5, 0.975))) ,
  mean = apply(pred.initial_age_mnd, 2, mean) 
))
colnames(dat.pred.domin_mnd)[2:5] = c("lower95", "med", "upper95", "mean")

(dat.pred.domin_mlgfd = data.frame(
  "initial_age" = plot.initial_age,
  t(apply(pred.initial_age_mlgfd, 2, quantile, probs = c(0.025, 0.5, 0.975))),
  mean = apply(pred.initial_age_mlgfd, 2, mean) 
))
colnames(dat.pred.domin_mlgfd)[2:5] = c("lower95", "med", "upper95", "mean")

(dat.pred.domin_mpd = data.frame(
  "initial_age" = plot.initial_age,
  t(apply(pred.initial_age_mpd, 2, quantile, probs = c(0.025, 0.5, 0.975))),
  mean = apply(pred.initial_age_mpd, 2, mean)
))
colnames(dat.pred.domin_mpd)[2:5] = c("lower95", "med", "upper95", "mean")

pred.initial_age_mconti_func_d = apply(pred.initial_age_mconti_func_d, 2, quantile, probs = c(0.025, 0.5, 0.975)) 

(dat.pred.domin_mconti_func_d = data.frame(
  "initial_age" = plot.initial_age,
  t(apply(pred.initial_age_mconti_func_d, 2, quantile, probs = c(0.025, 0.5, 0.975))),
  mean = apply(pred.initial_age_mconti_func_d, 2, mean)
))
colnames(dat.pred.domin_mconti_func_d)[2:4] = c("lower95", "med", "upper95")

dat_pred.domin_all = as.data.frame(cbind(initial_age = dat.pred.domin_mnd$initial_age,
                                         mnd.a = dat.pred.domin_mnd$mean,
                                         mlgfd.a = dat.pred.domin_mlgfd$mean,
                                         mpd.a = dat.pred.domin_mpd$mean,
                                         mconti_func_d.a = dat.pred.domin_mconti_func_d$mean))
dat_pred.domin_all$sum = NA
dat_pred.domin_all$sum = rowSums(abs(dat_pred.domin_all[, c(2:5)]))
R2_per = as.data.frame(t(apply(dat_pred.domin_all[, c(2:6)], 1, function(x){abs(x)/x[5]})))
R2_per$sum = plot.initial_age
colnames(R2_per)[5] = 'initial_age'

library(tidyr)
R2_per_1 = R2_per %>% pivot_longer(cols = mnd.a:mconti_func_d.a, values_to = "R2")
(plot_R2_per_domin = ggplot()+ 
    theme_classic()+ 
    geom_col(data = R2_per_1, aes(x = initial_age, y = R2*100, fill = name))+
    scale_fill_viridis_d()+
    labs(x = "Initial age", y = expression(paste('Dominance-',"R" ^ "2", "%"))) +
    scale_x_continuous(breaks = plot.initial_age) +
    guides(colour=guide_legend(title=NULL,
                               override.aes = list(size = 8)))+
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 30),
      axis.ticks = element_line(linewidth = 0.5),
      axis.ticks.length = unit(-0.4,'lines'),
      axis.text.x = element_text(face = "bold", size = 28), 
      axis.text.y = element_text(face = "bold", size = 28), 
      axis.title.x = element_text(face = "bold", size = 36), 
      axis.title.y = element_text(face = "bold", size = 36)
    ))


plot.mnd_effect_domin = ggplot() +
  theme_classic()+
  geom_ribbon(data = dat.pred.domin_mnd, aes(x = initial_age, y = med, ymin = lower95, ymax = upper95), 
              fill = "grey30", alpha = 0.2, size = 2)+
  geom_smooth(data = dat.pred.domin_mnd, aes(x = initial_age, y = med), 
              se = F, color = "grey30", linewidth = 4) +
  geom_pointrange(data = data_random_mnd_domin, aes(x = initial_age, y = random_m,
                                                    ymin = random_lower95, ymax = random_upper95),
                  color = "grey30", size = 2) +
  labs(x = "Initial age", y = "Mnd.a-dominance relationship") +
  scale_x_continuous(breaks = plot.initial_age) +
  theme(
    #axis.text.x = element_blank(),
    #axis.text.y = element_blank(),
    #axis.ticks = element_blank(),
    axis.ticks = element_line(linewidth = 0.5),
    axis.ticks.length = unit(-0.4,'lines'),
    axis.text.x = element_text(face = "bold", size = 28), 
    axis.text.y = element_text(face = "bold", size = 28), 
    axis.title.x = element_text(face = "bold", size = 36), 
    axis.title.y = element_text(face = "bold", size = 36)
  )
plot.mnd_effect_domin

plot.mlgfd_effect_domin = ggplot() +
  theme_classic()+
  geom_ribbon(data = dat.pred.domin_mlgfd, aes(x = initial_age, y = med, ymin = lower95, ymax = upper95), 
              fill = "grey30", alpha = 0.2, size = 2)+
  geom_smooth(data = dat.pred.domin_mlgfd, aes(x = initial_age, y = med), 
              se = F, color = "grey30", linewidth = 4) +
  geom_pointrange(data = data_random_mlgfd_domin, aes(x = initial_age, y = random_m,
                                                      ymin = random_lower95, ymax = random_upper95),
                  color = "grey30", size = 2) +
  labs(x = "Initial age", y = "Mlgfd.a-dominance relationship") +
  scale_x_continuous(breaks = plot.initial_age) +
  theme(
    #axis.text.x = element_blank(),
    #axis.text.y = element_blank(),
    #axis.ticks = element_blank(),
    axis.ticks = element_line(linewidth = 0.5),
    axis.ticks.length = unit(-0.4,'lines'),
    axis.text.x = element_text(face = "bold", size = 28), 
    axis.text.y = element_text(face = "bold", size = 28), 
    axis.title.x = element_text(face = "bold", size = 36), 
    axis.title.y = element_text(face = "bold", size = 36)
  )

plot.mlgfd_effect_domin


####### Merge plot for success ~ mnd mfd at different time windows #######
gap = ggplot(NULL)+theme_void()
library(ggpubr)
success.mnd.mlgfd.effects = ggarrange(gap, plot.mnd_effect_estab,
                                      gap, plot.mlgfd_effect_estab,
                                      gap, plot_R2_per_estab,
                                      gap, plot.mnd_effect_domin,
                                      gap, plot.mlgfd_effect_domin,
                                      gap, plot_R2_per_domin,
                                      nrow = 2, ncol = 6,
                                      hjust = 1,
                                      vjust = 1.5,
                                      labels = list('','a)', '', 'b)',
                                                    '', 'c)', '', 'd)',
                                                    '', 'e)', '', 'f)'),
                                      widths = c(0.1, 1, 0.1, 1, 0.1, 1.4,
                                                 0.1, 1, 0.1, 1, 0.1, 1.4),
                                      font.label = list(size = 50,
                                                        face = "bold",
                                                        color ="black"))

ggsave('results/figures_movingwin/success.mnd.mlgfd.effects_full_model.svg',
       width = 36, height = 20, units = c('in'),
       dpi = 300, limitsize = F, plot = success.mnd.mlgfd.effects)




############# Invasion impact ~ mnd mfd abundance weighted mean ############
load("D:/R projects/BSS/results/fit_results/dat_imp_sp_moving.rdata")
dat_imp_sp_moving = rbindlist(dat_imp_sp_moving)
dat_imp_sp_moving$gain_loss = as.numeric(dat_imp_sp_moving$gain_loss)

dat_imp_sp_moving$intro = 0
dat_imp_sp_moving = dat_imp_sp_moving %>% relocate(intro, .before = estab)
dat_imp_sp_moving[dat_imp_sp_moving$estab == 0 & dat_imp_sp_moving$domin == 0,]$intro = 1
dat_imp_sp_moving_c = dat_imp_sp_moving %>% filter(!is.na(gain_loss))
dat_imp_sp_moving_c$gain_loss_2 = dat_imp_sp_moving_c$gain_loss
dat_imp_sp_moving_c$mean.abs.change_2 = dat_imp_sp_moving_c$mean.abs.change*-1
dat_imp_sp_moving_c$mean.relative.change_2 = dat_imp_sp_moving_c$mean.relative.change*-1
dat_imp_sp_moving_c[dat_imp_sp_moving_c$gain_loss == 0,]$gain_loss_2 = 1
dat_imp_sp_moving_c[dat_imp_sp_moving_c$gain_loss == 1,]$gain_loss_2 = 0 

dat_imp_sp_moving_intro = dat_imp_sp_moving_c %>% filter(intro == 1)
dat_imp_sp_moving_estab = dat_imp_sp_moving_c %>% filter(estab == 1)
dat_imp_sp_moving_domin = dat_imp_sp_moving_c %>% filter(domin == 1)

#### Set up model for intro vs. native impact using INLA 
library(INLA)
library(brinla)

mod_intro_sp_imp_ab = inla(gain_loss_2 ~ mnd.a + mlgfd.a + t_win + 
                         + t_win:mnd.a + t_win:mlgfd.a +
                         + f(species, model="iid") + f(field, model="iid") +
                         f(field:plot, model="iid"),
                       control.compute = list(dic=T,waic=T,cpo=T),
                       family="binomial", data=dat_imp_sp_moving_intro)
mod_intro_sp_imp_ab_1 = inla(gain_loss_2 ~ mnd.a + mlgfd.a + t_win + 
                           + t_win:mnd.a + t_win:mlgfd.a 
                         + f(field, model="iid") +
                           f(field:plot, model="iid"),
                         control.compute = list(dic=T,waic=T,cpo=T),
                         family="binomial", data=dat_imp_sp_moving_intro)
c(mod_intro_sp_imp_ab$waic$waic, mod_intro_sp_imp_ab_1$waic$waic)
summary(mod_intro_sp_imp_ab)
summary(mod_intro_sp_imp_ab_1)
mod_intro_sp_imp_ab$summary.fixed
summary(mod_intro_sp_imp_ab)

mod_intro_sp_imp_ab$summary.random$field$fID = as.character(mod_intro_sp_imp_ab$summary.random$field$ID)
group_1 = length(mod_intro_sp_imp_ab$summary.random$species$mean)
group_2_1 = length(mod_intro_sp_imp_ab$summary.random$field$mean)
group_2_2 = length(mod_intro_sp_imp_ab$summary.random$`field:plot`$mean)
mod_intro_sp_imp_ab$summary.random$field = arrange(mod_intro_sp_imp_ab$summary.random$field, 
                                               mod_intro_sp_imp_ab$summary.random$field$ID)
mod_intro_sp_imp_ab$summary.random$`field:plot`$ID_1 = as.numeric(sapply(str_split(mod_intro_sp_imp_ab$summary.random$`field:plot`$ID,
                                                                               '_'), function(x){x[1]}))
mod_intro_sp_imp_ab$summary.random$`field:plot` = arrange(mod_intro_sp_imp_ab$summary.random$`field:plot`,
                                                      mod_intro_sp_imp_ab$summary.random$`field:plot`$ID_1,
                                                      mod_intro_sp_imp_ab$summary.random$`field:plot`$ID)

### mean nd effects on intro invasion impact
mnd.a.intro.imp.effect = mod_intro_sp_imp_ab$summary.fixed[2,1] +
  rep(mod_intro_sp_imp_ab$summary.random$species$mean,
      group_2_2)+
  rep(rep(mod_intro_sp_imp_ab$summary.random$field$mean,
          c(40, 48, 48, 48, 48, 40, 48, 48, 48, 48)), 
      each = group_1)+
  rep(mod_intro_sp_imp_ab$summary.random$`field:plot`$mean, each = group_1)+
  mod_intro_sp_imp_ab$summary.fixed[5,1]*(rep(summarise(group_by(dat_imp_sp_moving_intro,species),
                                                    t_win=mean(t_win))$t_win,group_2_2)  + 
                                        rep(summarise(group_by(dat_imp_sp_moving_intro,field,plot),
                                                      t_win=mean(t_win))$t_win, each = group_1))

mnd.a.intro.imp.effect.sd = rep(mod_intro_sp_imp_ab$summary.random$species$sd,
                            group_2_2)+
  rep(rep(mod_intro_sp_imp_ab$summary.random$field$sd,
          c(40, 48, 48, 48, 48, 40, 48, 48, 48, 48)), 
      each = group_1)+
  rep(mod_intro_sp_imp_ab$summary.random$`field:plot`$sd, each = group_1)

mnd.a.intro.imp.effect.df = data.frame(field = rep(rep(mod_intro_sp_imp_ab$summary.random$field$ID,
                                                   c(40, 48, 48, 48, 48, 40, 48, 48, 48, 48)), 
                                               each = group_1),
                                   plot = rep(mod_intro_sp_imp_ab$summary.random$`field:plot`$ID_1,
                                              each = group_1),
                                   species = rep(mod_intro_sp_imp_ab$summary.random$species$ID,
                                                 group_2_2),
                                   t_win = (rep(summarise(group_by(dat_imp_sp_moving_intro,species),
                                                          t_win=mean(t_win))$t_win,group_2_2)  + 
                                              rep(summarise(group_by(dat_imp_sp_moving_intro,field,plot),
                                                            t_win=mean(t_win))$t_win, each = group_1)),
                                   mnd.a.intro.imp.effect,
                                   sd=mnd.a.intro.imp.effect.sd)

plot.mnd.a.intro.imp.effect = ggplot(data=mnd.a.intro.imp.effect.df,
                                     aes(x=t_win,y=mnd.a.intro.imp.effect))+
  geom_point(size=3.5, alpha = 0.03)+
  geom_errorbar(aes(ymin=mnd.a.intro.imp.effect-sd,ymax=mnd.a.intro.imp.effect+sd),
                alpha = 0.03)+
  geom_abline(intercept=mod_intro_sp_imp_ab$summary.fixed[2,1],
              slope=mod_intro_sp_imp_ab$summary.fixed[5,1],
              linewidth = 2)+
  theme_custom()

### mean lgfd effects on intro invasion impact
mlgfd.a.intro.imp.effect = mod_intro_sp_imp_ab$summary.fixed[3,1] +
  rep(mod_intro_sp_imp_ab$summary.random$species$mean,
      group_2_2)+
  rep(rep(mod_intro_sp_imp_ab$summary.random$field$mean,
          c(40, 48, 48, 48, 48, 40, 48, 48, 48, 48)), 
      each = group_1)+
  rep(mod_intro_sp_imp_ab$summary.random$`field:plot`$mean, each = group_1)+
  mod_intro_sp_imp_ab$summary.fixed[6,1]*(rep(summarise(group_by(dat_imp_sp_moving_intro,species),
                                                        t_win=mean(t_win))$t_win,group_2_2)  + 
                                          rep(summarise(group_by(dat_imp_sp_moving_intro,field,plot),
                                                          t_win=mean(t_win))$t_win, each = group_1))

mlgfd.a.intro.imp.effect.sd = rep(mod_intro_sp_imp_ab$summary.random$species$sd,
                                group_2_2)+
  rep(rep(mod_intro_sp_imp_ab$summary.random$field$sd,
          c(40, 48, 48, 48, 48, 40, 48, 48, 48, 48)), 
      each = group_1)+
  rep(mod_intro_sp_imp_ab$summary.random$`field:plot`$sd, each = group_1)

mlgfd.a.intro.imp.effect.df = data.frame(field = rep(rep(mod_intro_sp_imp_ab$summary.random$field$ID,
                                                       c(40, 48, 48, 48, 48, 40, 48, 48, 48, 48)), 
                                                   each = group_1),
                                       plot = rep(mod_intro_sp_imp_ab$summary.random$`field:plot`$ID_1,
                                                  each = group_1),
                                       species = rep(mod_intro_sp_imp_ab$summary.random$species$ID,
                                                     group_2_2),
                                       t_win = (rep(summarise(group_by(dat_imp_sp_moving_intro,species),
                                                              t_win=mean(t_win))$t_win,group_2_2)  + 
                                                  rep(summarise(group_by(dat_imp_sp_moving_intro,field,plot),
                                                                t_win=mean(t_win))$t_win, each = group_1)),
                                       mlgfd.a.intro.imp.effect,
                                       sd=mlgfd.a.intro.imp.effect.sd)

plot.mlgfd.a.intro.imp.effect = ggplot(data=mlgfd.a.intro.imp.effect.df,
                                     aes(x=t_win,y=mlgfd.a.intro.imp.effect))+
  geom_point(size=3.5, alpha = 0.03)+
  geom_errorbar(aes(ymin=mlgfd.a.intro.imp.effect-sd,ymax=mlgfd.a.intro.imp.effect+sd),
                alpha = 0.03)+
  geom_abline(intercept=mod_intro_sp_imp_ab$summary.fixed[3,1],
              slope=mod_intro_sp_imp_ab$summary.fixed[6,1],
              linewidth = 2)+
  theme_custom()

#### Set up model for estab vs. native impact using INLA 
library(INLA)
library(brinla)

mod_estab_sp_imp_ab = inla(gain_loss_2 ~ mnd.a + mlgfd.a + t_win + 
                             + t_win:mnd.a + t_win:mlgfd.a +
                             + f(species, model="iid") + f(field, model="iid") +
                             f(field:plot, model="iid"),
                           control.compute = list(dic=T,waic=T,cpo=T),
                           family="binomial", data=dat_imp_sp_moving_estab)
mod_estab_sp_imp_ab_1 = inla(gain_loss_2 ~ mnd.a + mlgfd.a + t_win + 
                               + t_win:mnd.a + t_win:mlgfd.a 
                             + f(field, model="iid") +
                               f(field:plot, model="iid"),
                             control.compute = list(dic=T,waic=T,cpo=T),
                             family="binomial", data=dat_imp_sp_moving_estab)
c(mod_estab_sp_imp_ab$waic$waic, mod_estab_sp_imp_ab_1$waic$waic)
summary(mod_estab_sp_imp_ab)
summary(mod_estab_sp_imp_ab_1)
mod_estab_sp_imp_ab$summary.fixed
summary(mod_estab_sp_imp_ab)

group_1 = length(mod_estab_sp_imp_ab$summary.random$species$mean)
group_2_1 = length(mod_estab_sp_imp_ab$summary.random$field$mean)
group_2_2 = length(mod_estab_sp_imp_ab$summary.random$`field:plot`$mean)
mod_estab_sp_imp_ab$summary.random$field = arrange(mod_estab_sp_imp_ab$summary.random$field, 
                                                   mod_estab_sp_imp_ab$summary.random$field$ID)
mod_estab_sp_imp_ab$summary.random$`field:plot`$ID_1 = as.numeric(sapply(str_split(mod_estab_sp_imp_ab$summary.random$`field:plot`$ID,
                                                                                   '_'), function(x){x[1]}))
mod_estab_sp_imp_ab$summary.random$`field:plot` = arrange(mod_estab_sp_imp_ab$summary.random$`field:plot`,
                                                          mod_estab_sp_imp_ab$summary.random$`field:plot`$ID_1,
                                                          mod_estab_sp_imp_ab$summary.random$`field:plot`$ID)

### mean nd effects on estab invasion impact
mnd.a.estab.imp.effect = mod_estab_sp_imp_ab$summary.fixed[2,1] +
  rep(mod_estab_sp_imp_ab$summary.random$species$mean,
      group_2_2)+
  rep(rep(mod_estab_sp_imp_ab$summary.random$field$mean,
          c(40, 48, 48, 47, 48, 38, 48, 48, 48, 48)), 
      each = group_1)+
  rep(mod_estab_sp_imp_ab$summary.random$`field:plot`$mean, each = group_1)+
  mod_estab_sp_imp_ab$summary.fixed[5,1]*(rep(summarise(group_by(dat_imp_sp_moving_estab,species),
                                                        t_win=mean(t_win))$t_win,group_2_2)  + 
                                            rep(summarise(group_by(dat_imp_sp_moving_estab,field,plot),
                                                          t_win=mean(t_win))$t_win, each = group_1))

mnd.a.estab.imp.effect.sd = rep(mod_estab_sp_imp_ab$summary.random$species$sd,
                                group_2_2)+
  rep(rep(mod_estab_sp_imp_ab$summary.random$field$sd,
          c(40, 48, 48, 47, 48, 38, 48, 48, 48, 48)), 
      each = group_1)+
  rep(mod_estab_sp_imp_ab$summary.random$`field:plot`$sd, each = group_1)

mnd.a.estab.imp.effect.df = data.frame(field = rep(rep(mod_estab_sp_imp_ab$summary.random$field$ID,
                                                       c(40, 48, 48, 47, 48, 38, 48, 48, 48, 48)), 
                                                   each = group_1),
                                       plot = rep(mod_estab_sp_imp_ab$summary.random$`field:plot`$ID_1,
                                                  each = group_1),
                                       species = rep(mod_estab_sp_imp_ab$summary.random$species$ID,
                                                     group_2_2),
                                       t_win = (rep(summarise(group_by(dat_imp_sp_moving_estab,species),
                                                              t_win=mean(t_win))$t_win,group_2_2)  + 
                                                  rep(summarise(group_by(dat_imp_sp_moving_estab,field,plot),
                                                                t_win=mean(t_win))$t_win, each = group_1)),
                                       mnd.a.estab.imp.effect,
                                       sd=mnd.a.estab.imp.effect.sd)

plot.mnd.a.estab.imp.effect = ggplot(data=mnd.a.estab.imp.effect.df,
                                     aes(x=t_win,y=mnd.a.estab.imp.effect))+
  geom_point(size=3.5, alpha = 0.03)+
  geom_errorbar(aes(ymin=mnd.a.estab.imp.effect-sd,ymax=mnd.a.estab.imp.effect+sd),
                alpha = 0.03)+
  geom_abline(intercept=mod_estab_sp_imp_ab$summary.fixed[2,1],
              slope=mod_estab_sp_imp_ab$summary.fixed[5,1],
              linewidth = 2)+
  theme_custom()

### mean lgfd effects on estab invasion impact
mlgfd.a.estab.imp.effect = mod_estab_sp_imp_ab$summary.fixed[3,1] +
  rep(mod_estab_sp_imp_ab$summary.random$species$mean,
      group_2_2)+
  rep(rep(mod_estab_sp_imp_ab$summary.random$field$mean,
          c(40, 48, 48, 47, 48, 38, 48, 48, 48, 48)), 
      each = group_1)+
  rep(mod_estab_sp_imp_ab$summary.random$`field:plot`$mean, each = group_1)+
  mod_estab_sp_imp_ab$summary.fixed[6,1]*(rep(summarise(group_by(dat_imp_sp_moving_estab,species),
                                                        t_win=mean(t_win))$t_win,group_2_2)  + 
                                            rep(summarise(group_by(dat_imp_sp_moving_estab,field,plot),
                                                          t_win=mean(t_win))$t_win, each = group_1))

mlgfd.a.estab.imp.effect.sd = rep(mod_estab_sp_imp_ab$summary.random$species$sd,
                                  group_2_2)+
  rep(rep(mod_estab_sp_imp_ab$summary.random$field$sd,
          c(40, 48, 48, 47, 48, 38, 48, 48, 48, 48)), 
      each = group_1)+
  rep(mod_estab_sp_imp_ab$summary.random$`field:plot`$sd, each = group_1)

mlgfd.a.estab.imp.effect.df = data.frame(field = rep(rep(mod_estab_sp_imp_ab$summary.random$field$ID,
                                                         c(40, 48, 48, 47, 48, 38, 48, 48, 48, 48)), 
                                                     each = group_1),
                                         plot = rep(mod_estab_sp_imp_ab$summary.random$`field:plot`$ID_1,
                                                    each = group_1),
                                         species = rep(mod_estab_sp_imp_ab$summary.random$species$ID,
                                                       group_2_2),
                                         t_win = (rep(summarise(group_by(dat_imp_sp_moving_estab,species),
                                                                t_win=mean(t_win))$t_win,group_2_2)  + 
                                                    rep(summarise(group_by(dat_imp_sp_moving_estab,field,plot),
                                                                  t_win=mean(t_win))$t_win, each = group_1)),
                                         mlgfd.a.estab.imp.effect,
                                         sd=mlgfd.a.estab.imp.effect.sd)

plot.mlgfd.a.estab.imp.effect = ggplot(data=mlgfd.a.estab.imp.effect.df,
                                       aes(x=t_win,y=mlgfd.a.estab.imp.effect))+
  geom_point(size=3.5, alpha = 0.03)+
  geom_errorbar(aes(ymin=mlgfd.a.estab.imp.effect-sd,ymax=mlgfd.a.estab.imp.effect+sd),
                alpha = 0.03)+
  geom_abline(intercept=mod_estab_sp_imp_ab$summary.fixed[3,1],
              slope=mod_estab_sp_imp_ab$summary.fixed[6,1],
              linewidth = 2)+
  theme_custom()

#### Set up model for domin vs. native impact using INLA 
library(INLA)
library(brinla)

mod_domin_sp_imp_ab = inla(gain_loss_2 ~ mnd.a + mlgfd.a + t_win + 
                             + t_win:mnd.a + t_win:mlgfd.a +
                             + f(species, model="iid") + f(field, model="iid") +
                             f(field:plot, model="iid"),
                           control.compute = list(dic=T,waic=T,cpo=T),
                           family="binomial", data=dat_imp_sp_moving_domin)
mod_domin_sp_imp_ab_1 = inla(gain_loss_2 ~ mnd.a + mlgfd.a + t_win + 
                               + t_win:mnd.a + t_win:mlgfd.a 
                             + f(field, model="iid") +
                               f(field:plot, model="iid"),
                             control.compute = list(dic=T,waic=T,cpo=T),
                             family="binomial", data=dat_imp_sp_moving_domin)
c(mod_domin_sp_imp_ab$waic$waic, mod_domin_sp_imp_ab_1$waic$waic)
summary(mod_domin_sp_imp_ab)
summary(mod_domin_sp_imp_ab_1)
mod_domin_sp_imp_ab$summary.fixed
summary(mod_domin_sp_imp_ab)

group_1 = length(mod_domin_sp_imp_ab$summary.random$species$mean)
group_2_1 = length(mod_domin_sp_imp_ab$summary.random$field$mean)
group_2_2 = length(mod_domin_sp_imp_ab$summary.random$`field:plot`$mean)
mod_domin_sp_imp_ab$summary.random$field = arrange(mod_domin_sp_imp_ab$summary.random$field, 
                                                   mod_domin_sp_imp_ab$summary.random$field$ID)
mod_domin_sp_imp_ab$summary.random$`field:plot`$ID_1 = as.numeric(sapply(str_split(mod_domin_sp_imp_ab$summary.random$`field:plot`$ID,
                                                                                   '_'), function(x){x[1]}))
mod_domin_sp_imp_ab$summary.random$`field:plot` = arrange(mod_domin_sp_imp_ab$summary.random$`field:plot`,
                                                          mod_domin_sp_imp_ab$summary.random$`field:plot`$ID_1,
                                                          mod_domin_sp_imp_ab$summary.random$`field:plot`$ID)
table(mod_domin_sp_imp_ab$summary.random$`field:plot`$ID_1)

### mean nd effects on domin invasion impact
mnd.a.domin.imp.effect = mod_domin_sp_imp_ab$summary.fixed[2,1] +
  rep(mod_domin_sp_imp_ab$summary.random$species$mean,
      group_2_2)+
  rep(rep(mod_domin_sp_imp_ab$summary.random$field$mean,
          c(26, 42, 44, 36, 42, 26, 34, 39, 37, 39)), 
      each = group_1)+
  rep(mod_domin_sp_imp_ab$summary.random$`field:plot`$mean, each = group_1)+
  mod_domin_sp_imp_ab$summary.fixed[5,1]*(rep(summarise(group_by(dat_imp_sp_moving_domin,species),
                                                        t_win=mean(t_win))$t_win,group_2_2)  + 
                                            rep(summarise(group_by(dat_imp_sp_moving_domin,field,plot),
                                                          t_win=mean(t_win))$t_win, each = group_1))

mnd.a.domin.imp.effect.sd = rep(mod_domin_sp_imp_ab$summary.random$species$sd,
                                group_2_2)+
  rep(rep(mod_domin_sp_imp_ab$summary.random$field$sd,
          c(26, 42, 44, 36, 42, 26, 34, 39, 37, 39)), 
      each = group_1)+
  rep(mod_domin_sp_imp_ab$summary.random$`field:plot`$sd, each = group_1)

mnd.a.domin.imp.effect.df = data.frame(field = rep(rep(mod_domin_sp_imp_ab$summary.random$field$ID,
                                                       c(26, 42, 44, 36, 42, 26, 34, 39, 37, 39)), 
                                                   each = group_1),
                                       plot = rep(mod_domin_sp_imp_ab$summary.random$`field:plot`$ID_1,
                                                  each = group_1),
                                       species = rep(mod_domin_sp_imp_ab$summary.random$species$ID,
                                                     group_2_2),
                                       t_win = (rep(summarise(group_by(dat_imp_sp_moving_domin,species),
                                                              t_win=mean(t_win))$t_win,group_2_2)  + 
                                                  rep(summarise(group_by(dat_imp_sp_moving_domin,field,plot),
                                                                t_win=mean(t_win))$t_win, each = group_1)),
                                       mnd.a.domin.imp.effect,
                                       sd=mnd.a.domin.imp.effect.sd)

plot.mnd.a.domin.imp.effect = ggplot(data=mnd.a.domin.imp.effect.df,
                                     aes(x=t_win,y=mnd.a.domin.imp.effect))+
  geom_point(size=3.5, alpha = 0.03)+
  geom_errorbar(aes(ymin=mnd.a.domin.imp.effect-sd,ymax=mnd.a.domin.imp.effect+sd),
                alpha = 0.03)+
  geom_abline(intercept=mod_domin_sp_imp_ab$summary.fixed[2,1],
              slope=mod_domin_sp_imp_ab$summary.fixed[5,1],
              linewidth = 2)+
  theme_custom()

### mean lgfd effects on domin invasion impact
mlgfd.a.domin.imp.effect = mod_domin_sp_imp_ab$summary.fixed[3,1] +
  rep(mod_domin_sp_imp_ab$summary.random$species$mean,
      group_2_2)+
  rep(rep(mod_domin_sp_imp_ab$summary.random$field$mean,
          c(26, 42, 44, 36, 42, 26, 34, 39, 37, 39)), 
      each = group_1)+
  rep(mod_domin_sp_imp_ab$summary.random$`field:plot`$mean, each = group_1)+
  mod_domin_sp_imp_ab$summary.fixed[6,1]*(rep(summarise(group_by(dat_imp_sp_moving_domin,species),
                                                        t_win=mean(t_win))$t_win,group_2_2)  + 
                                            rep(summarise(group_by(dat_imp_sp_moving_domin,field,plot),
                                                          t_win=mean(t_win))$t_win, each = group_1))

mlgfd.a.domin.imp.effect.sd = rep(mod_domin_sp_imp_ab$summary.random$species$sd,
                                  group_2_2)+
  rep(rep(mod_domin_sp_imp_ab$summary.random$field$sd,
          c(26, 42, 44, 36, 42, 26, 34, 39, 37, 39)), 
      each = group_1)+
  rep(mod_domin_sp_imp_ab$summary.random$`field:plot`$sd, each = group_1)

mlgfd.a.domin.imp.effect.df = data.frame(field = rep(rep(mod_domin_sp_imp_ab$summary.random$field$ID,
                                                         c(26, 42, 44, 36, 42, 26, 34, 39, 37, 39)), 
                                                     each = group_1),
                                         plot = rep(mod_domin_sp_imp_ab$summary.random$`field:plot`$ID_1,
                                                    each = group_1),
                                         species = rep(mod_domin_sp_imp_ab$summary.random$species$ID,
                                                       group_2_2),
                                         t_win = (rep(summarise(group_by(dat_imp_sp_moving_domin,species),
                                                                t_win=mean(t_win))$t_win,group_2_2)  + 
                                                    rep(summarise(group_by(dat_imp_sp_moving_domin,field,plot),
                                                                  t_win=mean(t_win))$t_win, each = group_1)),
                                         mlgfd.a.domin.imp.effect,
                                         sd=mlgfd.a.domin.imp.effect.sd)

plot.mlgfd.a.domin.imp.effect = ggplot(data=mlgfd.a.domin.imp.effect.df,
                                       aes(x=t_win,y=mlgfd.a.domin.imp.effect))+
  geom_point(size=3.5, alpha = 0.03)+
  geom_errorbar(aes(ymin=mlgfd.a.domin.imp.effect-sd,ymax=mlgfd.a.domin.imp.effect+sd),
                alpha = 0.03)+
  geom_abline(intercept=mod_domin_sp_imp_ab$summary.fixed[3,1],
              slope=mod_domin_sp_imp_ab$summary.fixed[6,1],
              linewidth = 2)+
  theme_custom()

####### Merge plot for impact ~ mnd mfd at different time windows

library(ggpubr)
impact.mnd.mlgfd.effects = ggarrange(plot.mnd.a.intro.imp.effect,
                                    plot.mlgfd.a.intro.imp.effect,
                                    plot.mnd.a.estab.imp.effect,
                                    plot.mlgfd.a.estab.imp.effect,
                                    plot.mnd.a.domin.imp.effect,
                                    plot.mlgfd.a.domin.imp.effect,
                                    nrow = 3, ncol = 2, labels = list('a)',
                                    'b)', 'c)', 'd)', 'e)', 'f)'),
                                    font.label = list(size = 45,
                                                      face = "bold",
                                                      color ="black"))

ggsave('results/figures_movingwin/impact.mnd.mlgfd.effects.svg',
       width = 110, height = 160, units = c('cm'),
       dpi = 2, limitsize = F, plot = impact.mnd.mlgfd.effects)
ggsave('results/figures_movingwin/impact.mnd.mlgfd.effects_2.svg',
       width = 55, height = 80, units = c('cm'),
       dpi = 2, limitsize = F, plot = impact.mnd.mlgfd.effects)
ggsave('results/figures_movingwin/impact.mnd.mlgfd.effects.png',
       width = 110, height = 160, units = c('cm'),
       dpi = 150, limitsize = F, plot = impact.mnd.mlgfd.effects)
ggsave('results/figures_movingwin/impact.mnd.mlgfd.effects.pdf',
       width = 110, height = 160, units = c('cm'),
       dpi = 150, limitsize = F, plot = impact.mnd.mlgfd.effects)

