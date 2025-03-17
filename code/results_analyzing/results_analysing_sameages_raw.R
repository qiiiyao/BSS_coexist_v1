# Load packages
rm(list = ls())


############### Fast Start ####################
load("D:/R projects/BSS/results/fit_results/BSS_exculde_trees_raw/plot_sameages_filter_25/inter_all_c_alltime.rdata")
library(dplyr)
library(betareg)
library(stringr)
library(data.table)
library(ape)
library(reticulate)
inter_all_c = inter_all_c_alltime
# get the abs value of lgfd
inter_all_c$ablgfd = inter_all_c$lgfd
inter_all_c$ablgfd[inter_all_c$ablgfd < 0] = inter_all_c$ablgfd[inter_all_c$ablgfd < 0]*-1

####### N-F mapping ######
# loads the relevant python code
source_python("code/function/numerical_NFD.py")
source_python("code/function/scheme_plot.py")
source_python("code/function/scheme_plot_my.py")
# Source the basic functions for constructing the invasion graphs (original source: 10.5281/zenodo.7111753). 
source("code/function/invasion_graph_main_functions.R")
# Source the auxiliary invasion graph functions
source("code/function/invasion_graph_aux_functions.R")

plo = unique(unique(inter_all_c$f_p))
intra_l = list.files('results/fit_results/BSS_exculde_trees_raw/parameters/',
                     pattern = 'intra', full.names = T)

pars = list()
for (i in 1:length(plo)) {
  #i = 1
  inter_all_c_fp = inter_all_c %>% filter(f_p == plo[i])
  nd_m = dcast(inter_all_c, species_i ~ species_j, value.var = 'nd',
               fun.aggregate = sum, fill = 0)
  nd_m = nd_m[,-1]
  nd_m = as.matrix(nd_m)
  row.names(nd_m) = colnames(nd_m)
  ablgfd_m = dcast(inter_all_c, species_i ~ species_j, value.var = 'ablgfd',
                   fun.aggregate = sum, fill = 0)
  ablgfd_m = ablgfd_m[,-1]
  ablgfd_m = as.matrix(ablgfd_m)
  row.names(ablgfd_m) = colnames(ablgfd_m)
  ### Calculate the original functional distance
  tree = read.tree('data/original data/phylo_tree332.txt')
  sp_fit = unique(inter_all_c$species_i)
  tree_fit = keep.tip(tree, sp_fit)
  pd_m = cophenetic(tree_fit)
  mantel(pd_m, nd_m)
  mantel(pd_m, ablgfd_m)
  # cannot allocate too big vector
  # with(inter_all_c %>% filter(field == '1'), adonis2(nd ~ Phylo_dis,
  # strata = f_p, data = inter_all_c %>% filter(field == '1')))
  data(dune)
  data(dune.env)
  
  load(intra_l[grep(paste0(plo[i],'\\.'), intra_l)])
  native = inter_all_c_fp %>% filter(stage_ij == 'native_native')
  native_sp = unique(native$species_i)
  alien = inter_all_c_fp %>% filter(stage_j != 'native')
  alien_sp = unique(alien$species_j)
  pars[[i]] = list()
  for (z in 1:length(alien_sp)) {
    #z = 2
    native_onealien = inter_all_c_fp %>% filter(species_i %in% c(native_sp,
                                                                 alien_sp[z]) &
                                                  species_j %in% c(native_sp,
                                                                   alien_sp[z]))
    native_onealien = arrange(native_onealien, native_onealien$species_i)
    spe = unique(native_onealien$species_i)
    aii = unique(native_onealien$aii)
    a = dcast(native_onealien, species_i ~ species_j, value.var = 'aij',
              fun.aggregate = sum, fill = 0)
    
    
    a = a[,-1]
    a = as.matrix(a)
    row.names(a) = colnames(a)
    focal_sp = which(colnames(a) == alien_sp[z])
    for (j in 1:length(aii)) {
      a[colnames(a)[j], colnames(a)[j]] = aii[j]
    }
    result = tryCatch(
      {
        my_result = compute_NFD_LV(a)
        my_result
      },
      error = function(err) {
        
        cat("Error:", conditionMessage(err), "\n")
        NA  
      }
    )
    
    pars[[i]][[z]] = result
  }
  names(pars[[i]]) = alien_sp
}

names(pars) = plo[i]

#### Check the correlation between nd and fd
cor.test(inter_all_c$nd, inter_all_c$ablgfd) 
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
              inter_all_c$estab_i == 1&
              inter_all_c$domin_i == 0,]$stage_i_domin = 'estab.nodomin'
inter_all_c[inter_all_c$intro_i == 1&
              inter_all_c$estab_i == 1&
              inter_all_c$domin_i == 1,]$stage_i_domin = 'estab.domin'
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
library(ggpubr)
library(ggeffects)
library(stringr)
source('code/function/plot_func.R')

#### Just three stages #########
unique(inter_all_c$stage_ij)
inter_all_forcoexplot = inter_all_c %>%
  filter(stage_ij %in% c("intro_native",
                         "estab_native",
                         "domin_native"))
mod_nd_compare_1 = lmer(nd ~ stage_ij+(1|field/plot),
                        data = inter_all_forcoexplot, 
                        REML = F)
mod_nd_compare_lm = lm(nd ~ stage_ij,
                       data = inter_all_forcoexplot)
anova(mod_nd_compare_1, mod_nd_compare_lm)
mod_nd_compare = lmer(nd ~ stage_ij+(1|field/plot),
                      data = inter_all_forcoexplot, 
                      REML = T)

mod_fd_compare = lmer(lgfd ~ stage_ij+(1|field/plot),
                      data = inter_all_forcoexplot,
                      REML = T)
summary(mod_nd_compare)
summary(mod_fd_compare)


pre_mod_nd_compare = as.data.frame(ggpredict(mod_nd_compare,
                                             terms = 'stage_ij'))

pre_mod_fd_compare = as.data.frame(ggpredict(mod_fd_compare,
                                             terms = 'stage_ij'))
pre_ndfd = data.frame(nd = pre_mod_nd_compare$predicted,
                      lgfd = pre_mod_fd_compare$predicted,
                      stage = pre_mod_nd_compare$x,
                      nd_sd = pre_mod_nd_compare$std.error,
                      nd_low = pre_mod_nd_compare$conf.low,
                      nd_high = pre_mod_nd_compare$conf.high,
                      fd_sd = pre_mod_fd_compare$std.error,
                      fd_low = pre_mod_fd_compare$conf.low,
                      fd_high = pre_mod_fd_compare$conf.high)

spss.f = function(x) log10(1-x)
x = seq(-1, 1, 0.001)
spss.ff = function(x) -log10(1-x)

Fig.1 = ggplot(NULL) +
  geom_pointrange(data = pre_ndfd,
                  mapping = aes(x = nd,y = lgfd, ymin = fd_low,
                                ymax = fd_high,
                                color = stage),
                  fatten = 20,linewidth = 4)+
  geom_pointrange(data = pre_ndfd,
                  mapping = aes(y = lgfd,x = nd, xmin = nd_low,
                                xmax = nd_high,
                                color = stage),
                  fatten = 20,linewidth = 4)+
  scale_color_brewer(palette="Set2")+
  scale_x_continuous(limits=c(0, 0.6)) +
  scale_y_continuous(limits=c(-0.1, 0.1)) +
  geom_hline(yintercept=0, linetype="dashed", size=0.5) +
  geom_vline(xintercept=0, linetype="dashed", size=0.5) +
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
  #annotate(geom="text", x=-0.7, y=0.1, label="Priority effect", size=20) +
  guides(colour=guide_legend(title="Invasion stage")) +
  theme_custom_legend()

ggsave('results/figures_sameages_raw/Fig.1.svg',
       plot=Fig.1, 
       width=105, height=80, dpi=600, units='cm',
       limitsize=F)

########## Stage split two parts: estab noestab domin nodomin ##########
#### just mean and sd of raw data ####
inter_all_forexotics = inter_all_c %>%
  filter(stage_i != 'native', stage_j == 'native')
length(unique(inter_all_forexotics$species_i))
length(unique(inter_all_forexotics$species_j))

inter_all_forestab = inter_all_c %>%
  filter(stage_ij_estab %in% c("intro.estab_native",
                               "intro.noestab_native"))

length(unique(inter_all_forestab$sp_pair))
length(unique(inter_all_forestab$species_i))
length(unique(inter_all_forestab$species_j))

inter_all_forestab_msd = inter_all_forestab %>% 
  group_by(stage_ij_estab) %>% summarise_at(vars(nd, lgfd),
                                            c(mean, sd), na.rm = T)

colnames(inter_all_forestab_msd) = c('stage', 'nd_mean', 'lgfd_mean', 'nd_sd',
                                     'lgfd_sd')
spss.f = function(x) log10(1-x)
x = seq(-1, 1, 0.001)
spss.ff = function(x) -log10(1-x)


Fig.1_compare_estab_original = ggplot(NULL) +
  geom_point(data = inter_all_forestab, aes(x = nd, y = lgfd,
                                            color = stage_ij_estab),
             alpha = 0.1, size = 30)+
  geom_point(data = inter_all_forestab_msd,
             mapping = aes(x = nd_mean,y = lgfd_mean),
             color = 'black',
             size = 50)+
  geom_pointrange(data = inter_all_forestab_msd,
                  mapping = aes(x = nd_mean,y = lgfd_mean,
                                ymin = lgfd_mean-lgfd_sd,
                                ymax = lgfd_mean+lgfd_sd,
                                color = stage),
                  fatten = 70,linewidth = 6)+
  geom_pointrange(data = inter_all_forestab_msd,
                  mapping = aes(x = nd_mean,y = lgfd_mean,
                                xmin = nd_mean-nd_sd,
                                xmax = nd_mean+nd_sd,
                                color = stage),
                  fatten = 70,linewidth = 6)+
  scale_color_discrete(type = c('#efcc01',
                                '#025302'
  ))+
  #scale_color_brewer(palette="Set2")+
  scale_x_continuous(limits=c(-1, 1.5)) +
  scale_y_continuous(limits=c(-0.8, 0.8)) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=3) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=3) +
  xlab("Niche difference") +
  ylab(expression(paste(log[10],
                        "(Fitness difference)",
                        sep=' ')))+
  stat_function(fun=spss.f, colour="black", size = 5)+
  stat_function(fun=spss.ff, colour="black", size = 5)+
  #scale_x_continuous(breaks=c(-1, 1)) +
  annotate(geom="text", x=0.2, y=0.7, label=str_wrap("Invasion and exclude residents",
                                                     width = 20), size = 50) +
  annotate(geom="text", x=0.25, y=-0.7, label=str_wrap("Residents repel invasion",
                                                       width = 20), size = 50) +
  annotate(geom="text", x=1.2, y=0.7, label=str_wrap("Invasion and Coexistence",
                                                     width = 20), size = 50) +
  annotate(geom="text", x=1.2, y=-0.7, label=str_wrap("Invasion and Coexistence",
                                                      width = 20), size = 50) +
  #annotate(geom="text", x=-0.7, y=0.1, label="Priority effect", size=20) +
  guides(colour=guide_legend(title=NULL,
                             override.aes = list(size = 8))) +
  theme(panel.background = element_rect(fill = 'white',color = 'black',linewidth = 0.5),
        panel.grid = element_blank(),
        legend.position = c(0.18, 0.9),
        legend.text = element_text(size=120),
        #legend.title = element_text(size=80),
        plot.margin = margin(10,10,10,10),
        plot.background = element_blank(),
        text = element_text(size = 150),
        axis.ticks = element_line(linewidth = 3),
        axis.ticks.length = unit(-5,'lines'),
        axis.title.y = element_text(margin = margin(0,3,0,0),color = '#000000'),
        axis.title.x = element_text(margin = margin(5,0,0,0),color = '#000000'),
        axis.text.y = element_text(margin = margin(0,3,0,0),color = '#000000'),
        axis.text.x = element_text(margin = margin(5,0,0,0),color = '#000000'))


inter_all_fordomin = inter_all_c %>%
  filter(stage_ij_domin %in% c("estab.domin_native",
                               "estab.nodomin_native"))

inter_all_fordomin_msd = inter_all_fordomin %>% 
  group_by(stage_ij_domin) %>% summarise_at(vars(nd, lgfd),
                                            c(mean, sd), na.rm = T)

colnames(inter_all_fordomin_msd) = c('stage', 'nd_mean', 'lgfd_mean', 'nd_sd',
                                     'lgfd_sd')

spss.f = function(x) log10(1-x)
x = seq(-1, 1, 0.001)
spss.ff = function(x) -log10(1-x)

Fig.1_compare_domin_original = ggplot(NULL) +
  geom_point(data = inter_all_fordomin, aes(x = nd, y = lgfd,
                                            color = stage_ij_domin
  ), alpha = 0.1,
  size = 30)+
  geom_point(data = inter_all_fordomin_msd,
             mapping = aes(x = nd_mean,y = lgfd_mean),
             color = 'black',
             size = 50)+
  geom_pointrange(data = inter_all_fordomin_msd,
                  mapping = aes(x = nd_mean,y = lgfd_mean,
                                ymin = lgfd_mean-lgfd_sd,
                                ymax = lgfd_mean+lgfd_sd,
                                color = stage),
                  fatten = 70,linewidth = 6)+
  geom_pointrange(data = inter_all_fordomin_msd,
                  mapping = aes(x = nd_mean,y = lgfd_mean,
                                xmin = nd_mean-nd_sd,
                                xmax = nd_mean+nd_sd,
                                color = stage),
                  fatten = 70,linewidth = 6)+
  scale_color_discrete(type = c('#f40407',
                                '#0266be'
                                
  ))+
  #scale_color_brewer(palette="Set2")+
  scale_x_continuous(limits=c(-1, 1.5)) +
  scale_y_continuous(limits=c(-0.8, 0.8)) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=3) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=3) +
  labs(x = "Niche difference", 
       y = NULL) +
  stat_function(fun=spss.f, colour="black", size = 5)+
  stat_function(fun=spss.ff, colour="black", size = 5)+
  #scale_x_continuous(breaks=c(-1, 1)) +
  annotate(geom="text", x=0.3, y=0.7, label=str_wrap("Invasion and exclude residents",
                                                     width = 20), size = 50) +
  annotate(geom="text", x=0.3, y=-0.7, label=str_wrap("Residents repel invasion",
                                                      width = 20), size = 50) +
  annotate(geom="text", x=1.2, y=0.7, label=str_wrap("Invasion and Coexistence",
                                                     width = 20), size = 50) +
  annotate(geom="text", x=1.2, y=-0.7, label=str_wrap("Invasion and Coexistence",
                                                      width = 20), size = 50) +
  #annotate(geom="text", x=-0.7, y=0.1, label="Priority effect", size=20) +
  guides(colour=guide_legend(title=NULL,
                             override.aes = list(size = 8))) +
  theme(panel.background = element_rect(fill = 'white',color = 'black',linewidth = 0.5),
        panel.grid = element_blank(),
        legend.position = c(0.18, 0.9),
        legend.text = element_text(size = 120),
        plot.margin = margin(10,10,10,10),
        plot.background = element_blank(),
        text = element_text(size = 150),
        axis.ticks = element_line(linewidth = 3),
        axis.ticks.length = unit(-5,'lines'),
        axis.title.y = element_text(margin = margin(0,3,0,0),color = '#000000'),
        axis.title.x = element_text(margin = margin(5,0,0,0),color = '#000000'),
        axis.text.y = element_text(margin = margin(0,3,0,0),color = '#000000'),
        axis.text.x = element_text(margin = margin(5,0,0,0),color = '#000000'))
Fig.1_compare_domin_original

gap = ggplot(NULL)+theme_void()
library(ggpubr)
Fig.1_split_stage_original = ggarrange(gap,Fig.1_compare_estab_original, gap,
                                       Fig.1_compare_domin_original, nrow = 4, ncol = 1,
                                       labels = c('','a)','','b)'), hjust = -0.5,
                                       vjust = 0, heights = c(0.1, 1, 0.1, 1),
                                       font.label = list(size = 150))

ggsave(plot = Fig.1_split_stage_original,
       'results/figures_sameages_raw/Fig.1_split_stage_original.svg',
       width = 165,height = 340, dpi = 300, units = 'cm',
       limitsize = F)

library(ggpubr)
Fig.1_split_stage_original_forppt = ggarrange(gap,Fig.1_compare_estab_original, gap,
                                              Fig.1_compare_domin_original,
                                              nrow = 1, ncol = 4,
                                              labels = c('','a)','','b)'), hjust = 0,
                                              vjust = 1, widths = c(0.05, 1, 0.05, 1),
                                              font.label = list(size = 150))

ggsave(plot = Fig.1_split_stage_original_forppt,
       'results/figures_sameages_raw/Fig.1_split_stage_original_forppt.svg',
       width = 350,height = 165, dpi = 300, units = 'cm',
       limitsize = F)

#### estab noestab ####
### Normal methods ###
unique(inter_all_c$stage_ij_estab)
inter_all_forestab = inter_all_c %>%
  filter(stage_ij_estab %in% c("intro.estab_native",
                               "intro.noestab_native"))


mod_nd_compare_1 = lmer(nd ~ stage_ij_estab+(1|field/plot),
                        data = inter_all_forestab, 
                        REML = F)
mod_nd_compare_lm = lm(nd ~ stage_ij_estab,
                       data = inter_all_forestab)
anova(mod_nd_compare_1, mod_nd_compare_lm)
mod_nd_compare_estab_1 = lmer(nd ~ stage_ij_estab+(1|field/plot),
                              data = inter_all_forestab, 
                              REML = T)
mod_nd_compare_estab_2 = lmer(nd ~ stage_ij_estab+(1|field/plot)+(1|sp_pair),
                              data = inter_all_forestab, 
                              REML = T)
mod_nd_compare_estab = lmer(nd ~ stage_ij_estab+(1|field/plot)+
                              (1|species_i) +(1|species_j) +
                              (1|sp_pair),
                            data = inter_all_forestab, 
                            REML = T)
anova(mod_nd_compare_estab_1, mod_nd_compare_estab) ### select mod_nd_compare_estab
anova(mod_nd_compare_estab_2, mod_nd_compare_estab) ### select mod_nd_compare_estab

mod_fd_compare_estab_1 = lmer(lgfd ~ stage_ij_estab+(1|field/plot),
                              data = inter_all_forestab,
                              REML = T)
mod_fd_compare_estab = lmer(lgfd ~ stage_ij_estab+(1|field/plot)+
                              (1|species_i) + (1|species_j) +
                              (1|sp_pair),
                            data = inter_all_forestab,
                            REML = T)
mod_fd_compare_estab_2 = lmer(lgfd ~ stage_ij_estab+(1|field/plot)+
                                (1|sp_pair),
                              data = inter_all_forestab,
                              REML = T)

anova(mod_fd_compare_estab_1, mod_fd_compare_estab) ### select mod_fd_compare_estab
anova(mod_fd_compare_estab_2, mod_fd_compare_estab) ### select mod_fd_compare_estab

summary(mod_nd_compare_estab)
summary(mod_fd_compare_estab)

pre_mod_nd_compare_estab = as.data.frame(ggpredict(mod_nd_compare_estab,
                                                   terms = 'stage_ij_estab',
                                                   type = 're'))

pre_mod_fd_compare_estab = as.data.frame(ggpredict(mod_fd_compare_estab,
                                                   terms = 'stage_ij_estab',
                                                   type = 're'))

pre_ndfd_compare_estab = data.frame(nd = pre_mod_nd_compare_estab$predicted,
                                    lgfd = pre_mod_fd_compare_estab$predicted,
                                    stage = pre_mod_nd_compare_estab$x,
                                    nd_sd = pre_mod_nd_compare_estab$std.error,
                                    nd_low = pre_mod_nd_compare_estab$conf.low,
                                    nd_high = pre_mod_nd_compare_estab$conf.high,
                                    fd_sd = pre_mod_fd_compare_estab$std.error,
                                    fd_low = pre_mod_fd_compare_estab$conf.low,
                                    fd_high = pre_mod_fd_compare_estab$conf.high)

spss.f = function(x) log10(1-x)
x = seq(-1, 1, 0.001)
spss.ff = function(x) -log10(1-x)
# Define a collection of palettes to alter the default based on number of levels to encode
# Template function for creating densities grouped by a variable

Fig.1_compare_estab = ggplot(NULL) +
  geom_point(data = inter_all_forestab, aes(x = nd, y = lgfd,
                                            color = stage_ij_estab),
             alpha = 0.1, size = 30)+
  geom_point(data = pre_ndfd_compare_estab,
             mapping = aes(x = nd,y = lgfd),
             color = 'black',
             size = 50)+
  geom_pointrange(data = pre_ndfd_compare_estab,
                  mapping = aes(x = nd,y = lgfd, ymin = fd_low,
                                ymax = fd_high,
                                color = stage),
                  fatten = 70,linewidth = 3)+
  geom_pointrange(data = pre_ndfd_compare_estab,
                  mapping = aes(y = lgfd,x = nd, xmin = nd_low,
                                xmax = nd_high,
                                color = stage),
                  fatten = 70,linewidth = 3)+
  scale_color_discrete(type = c('#efcc01',
                                '#025302'
  ))+
  #scale_color_brewer(palette="Set2")+
  scale_x_continuous(limits=c(-1, 1.5)) +
  scale_y_continuous(limits=c(-0.8, 0.8)) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=3) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=3) +
  xlab("Niche difference") +
  ylab(expression(paste(log[10],
                        "(Fitness difference)",
                        sep=' ')))+
  stat_function(fun=spss.f, colour="black", size = 5)+
  stat_function(fun=spss.ff, colour="black", size = 5)+
  #scale_x_continuous(breaks=c(-1, 1)) +
  annotate(geom="text", x=0.2, y=0.7, label=str_wrap("Invasion and exclude residents",
                                                     width = 20), size = 50) +
  annotate(geom="text", x=0.25, y=-0.7, label=str_wrap("Residents repel invasion",
                                                       width = 20), size = 50) +
  annotate(geom="text", x=1.2, y=0.7, label=str_wrap("Invasion and Coexistence",
                                                     width = 20), size = 50) +
  annotate(geom="text", x=1.2, y=-0.7, label=str_wrap("Invasion and Coexistence",
                                                      width = 20), size = 50) +
  #annotate(geom="text", x=-0.7, y=0.1, label="Priority effect", size=20) +
  guides(colour=guide_legend(title=NULL,
                             override.aes = list(size = 8))) +
  theme(panel.background = element_rect(fill = 'white',color = 'black',linewidth = 0.5),
        panel.grid = element_blank(),
        legend.position = c(0.18, 0.9),
        legend.text = element_text(size=120),
        #legend.title = element_text(size=80),
        plot.margin = margin(10,10,10,10),
        plot.background = element_blank(),
        text = element_text(size = 150),
        axis.ticks = element_line(linewidth = 3),
        axis.ticks.length = unit(-5,'lines'),
        axis.title.y = element_text(margin = margin(0,3,0,0),color = '#000000',face = 'bold'),
        axis.title.x = element_text(margin = margin(5,0,0,0),color = '#000000',face = 'bold'),
        axis.text.y = element_text(margin = margin(0,3,0,0),color = '#000000'),
        axis.text.x = element_text(margin = margin(5,0,0,0),color = '#000000'))

Fig.1_compare_estab

#### domin nodomin ####
unique(inter_all_c$stage_ij_domin)
inter_all_fordomin = inter_all_c %>%
  filter(stage_ij_domin %in% c("estab.domin_native",
                               "estab.nodomin_native"))

mod_nd_compare_domin_1 = lmer(nd ~ stage_ij_domin+(1|field/plot),
                              data = inter_all_fordomin, 
                              REML = F)
mod_nd_compare_domin_2 = lmer(nd ~ stage_ij_domin+(1|field/plot)+(1|sp_pair),
                              data = inter_all_fordomin, 
                              REML = T)
mod_nd_compare_domin = lmer(nd ~ stage_ij_domin+(1|field/plot)+
                              (1|species_i) + (1|species_j) +
                              (1|sp_pair),
                            data = inter_all_fordomin, 
                            REML = T)
anova(mod_nd_compare_domin, mod_nd_compare_domin_1) ## select mod_nd_compare_domin
anova(mod_nd_compare_domin_2, mod_nd_compare_domin) ## select mod_nd_compare_domin

mod_fd_compare_domin_1 = lmer(lgfd ~ stage_ij_domin+(1|field/plot),
                              data = inter_all_fordomin,
                              REML = T)
mod_fd_compare_domin_2 = lmer(lgfd ~ stage_ij_domin+(1|field/plot)+(1|sp_pair),
                              data = inter_all_fordomin,
                              REML = T)
mod_fd_compare_domin = lmer(lgfd ~ stage_ij_domin+(1|field/plot)+
                              (1|species_i) + (1|species_j) +
                              (1|sp_pair),
                            data = inter_all_fordomin,
                            REML = T)
anova(mod_fd_compare_domin,mod_fd_compare_domin_1) ## select mod_fd_compare_domin
summary(mod_nd_compare_domin)
summary(mod_fd_compare_domin)


pre_mod_nd_compare_domin = as.data.frame(ggpredict(mod_nd_compare_domin,
                                                   terms = 'stage_ij_domin',
                                                   type = 're'))

pre_mod_fd_compare_domin = as.data.frame(ggpredict(mod_fd_compare_domin,
                                                   terms = 'stage_ij_domin',
                                                   type = 're'))

pre_ndfd_compare_domin = data.frame(nd = pre_mod_nd_compare_domin$predicted,
                                    lgfd = pre_mod_fd_compare_domin$predicted,
                                    stage = pre_mod_nd_compare_domin$x,
                                    nd_sd = pre_mod_nd_compare_domin$std.error,
                                    nd_low = pre_mod_nd_compare_domin$conf.low,
                                    nd_high = pre_mod_nd_compare_domin$conf.high,
                                    fd_sd = pre_mod_fd_compare_domin$std.error,
                                    fd_low = pre_mod_fd_compare_domin$conf.low,
                                    fd_high = pre_mod_fd_compare_domin$conf.high)

spss.f = function(x) log10(1-x)
x = seq(-1, 1, 0.001)
spss.ff = function(x) -log10(1-x)

Fig.1_compare_domin = ggplot(NULL) +
  geom_point(data = inter_all_fordomin, aes(x = nd, y = lgfd,
                                            color = stage_ij_domin
  ), alpha = 0.1,
  size = 30)+
  geom_point(data = pre_ndfd_compare_domin,
             mapping = aes(x = nd,y = lgfd),
             color = 'black',
             size = 50)+
  geom_pointrange(data = pre_ndfd_compare_domin,
                  mapping = aes(x = nd,y = lgfd, ymin = fd_low,
                                ymax = fd_high,
                                color = stage),
                  fatten = 70,linewidth = 3)+
  geom_pointrange(data = pre_ndfd_compare_domin,
                  mapping = aes(y = lgfd,x = nd, xmin = nd_low,
                                xmax = nd_high,
                                color = stage),
                  fatten = 70,linewidth = 3)+
  scale_color_discrete(type = c('#f40407',
                                '#0266be'
                                
  ))+
  #scale_color_brewer(palette="Set2")+
  scale_x_continuous(limits=c(-1, 1.5)) +
  scale_y_continuous(limits=c(-0.8, 0.8)) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=3) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=3) +
  xlab("Niche difference") +
  ylab(expression(paste(log[10],
                        "(Fitness difference)",
                        sep=' ')))+
  stat_function(fun=spss.f, colour="black", size = 5)+
  stat_function(fun=spss.ff, colour="black", size = 5)+
  #scale_x_continuous(breaks=c(-1, 1)) +
  annotate(geom="text", x=0.3, y=0.7, label=str_wrap("Invasion and exclude residents",
                                                     width = 20), size = 50) +
  annotate(geom="text", x=0.3, y=-0.7, label=str_wrap("Residents repel invasion",
                                                      width = 20), size = 50) +
  annotate(geom="text", x=1.2, y=0.7, label=str_wrap("Invasion and Coexistence",
                                                     width = 20), size = 50) +
  annotate(geom="text", x=1.2, y=-0.7, label=str_wrap("Invasion and Coexistence",
                                                      width = 20), size = 50) +
  #annotate(geom="text", x=-0.7, y=0.1, label="Priority effect", size=20) +
  guides(colour=guide_legend(title=NULL,
                             override.aes = list(size = 8))) +
  theme(panel.background = element_rect(fill = 'white',color = 'black',linewidth = 0.5),
        panel.grid = element_blank(),
        legend.position = c(0.18, 0.9),
        legend.text = element_text(size = 120),
        plot.margin = margin(10,10,10,10),
        plot.background = element_blank(),
        text = element_text(size = 150),
        axis.ticks = element_line(linewidth = 3),
        axis.ticks.length = unit(-5,'lines'),
        axis.title.y = element_text(margin = margin(0,3,0,0),color = '#000000',face = 'bold'),
        axis.title.x = element_text(margin = margin(5,0,0,0),color = '#000000',face = 'bold'),
        axis.text.y = element_text(margin = margin(0,3,0,0),color = '#000000'),
        axis.text.x = element_text(margin = margin(5,0,0,0),color = '#000000'))
Fig.1_compare_domin

#ggsave('results/figures_sameages_raw/Fig.1_compare_domin.svg',
#       plot=Fig.1_compare_domin, 
#       width=105, height=80, dpi=600, units='cm',
#       limitsize=F)

gap = ggplot(NULL)+theme_void()
library(ggpubr)
Fig.1_split_stage = ggarrange(gap,Fig.1_compare_estab,gap,
                              Fig.1_compare_domin, nrow = 1, ncol = 4,
                              labels = c('','a)','','b)'), hjust = 1,
                              vjust = 1.5, widths = c(0.1, 1, 0.1, 1),
                              font.label = list(size = 150))

ggsave(plot = Fig.1_split_stage,
       'results/figures_sameages_raw/Fig.1_split_stage.svg',
       width = 320,height = 100, dpi = 300, units = 'cm',
       limitsize = F)

Fig.1_split_stage_2 = ggarrange(gap,Fig.1_compare_estab,gap,
                                Fig.1_compare_domin, nrow = 4, ncol = 1,
                                labels = c('','a)','','b)'), hjust = -0.5,
                                vjust = 0, heights = c(0.1, 1, 0.1, 1),
                                font.label = list(size = 150))

ggsave(plot = Fig.1_split_stage_2,
       'results/figures_sameages_raw/Fig.1_split_stage_2.svg',
       width = 165,height = 340, dpi = 300, units = 'cm',
       limitsize = F)


######## Analyze relationship between ND & FD and Functional traits & phylogeny ########
library(sjPlot)
library(INLA)
library(inlabru)

inter_all_c_nd_posi = inter_all_c %>% filter(nd > 0)

inter_all_c_n_a_1 = inter_all_c %>% filter(stage_i != 'native' & stage_j == 'native')
inter_all_c_n_a_2 = inter_all_c %>% filter(stage_i == 'native' & stage_j != 'native')
inter_all_c_n_a = rbind(inter_all_c_n_a_1, inter_all_c_n_a_2)
sp_fit_n_a = unique(inter_all_c_n_a$species_i)
tree_fit_n_a = keep.tip(tree, sp_fit_n_a)
vcv_tree_n_a = ape::vcv(tree_fit_n_a, model = "Brownian", corr = FALSE)
vcv_tree_sparse_n_a = inla.as.sparse(solve(vcv_tree_n_a))
save(vcv_tree_n_a,
     file = 'code/results_analyzing/analysing_sameages_raw_data/vcv_tree_n_a.rdata')

### Calculate the original functional distance
tree = read.tree('data/original data/phylo_tree332.txt')
sp_fit = unique(inter_all_c$species_i)
tree_fit = keep.tip(tree, sp_fit)
vcv_tree = ape::vcv(tree_fit, model = "Brownian", corr = FALSE)
vcv_tree_sparse = inla.as.sparse(solve(vcv_tree))
save(vcv_tree_sparse, file = 'code/results_analyzing/analysing_sameages_raw_data/vcv_tree_sparse.rdata')


##### Nonlinear mixed effects model ######
###### nd/ablgfd ~ pd ######
library(nlme)
library(lme4)
library(minpack.lm)

max_e = function(a, b, dist){ 
  f = a*(1-exp(-(b*dist)))
  return(f)}
max_e_2 = function(b, dist){ 
  f = (1-exp(-(b*dist)))
  return(f)}
power = function(a, b, dist) {
  f = a*(dist^b)
  return(f)
}
hyperbola = function(a, b, dist) {
  f = dist/(a+b*dist)
  return(f)
}
logistic = function(a, b, c, dist) {
  f = a*dist/(1+exp(1)*((b-dist)/c))
}

dist = inter_all_c_nd_posi$Phylo_dis
#plot(dist, max_e(1, 0.01, dist))
#plot(dist, max_e_2(0.002740836, dist))
#plot(dist, power(1, 0.01, dist))
#plot(dist, hyperbola(1, 0.01, dist))
#plot(dist, logistic(a = 2, b = 1, c = 5,
#                    dist))

model_e1 = nlsLM(nd ~ max_e(a, b, dist = Phylo_dis),
                start=list(a = 1, b = 0.01),
                data = inter_all_c_nd_posi)
model_e2 = nlsLM(nd ~ max_e_2(b, dist = Phylo_dis),
                start=list(b = 0.01),
                data = inter_all_c_nd_posi)
model_power = nlsLM(nd ~ power(a, b, dist = Phylo_dis),
                 start=list(a = 1, b = 0.01),
                 data = inter_all_c_nd_posi)
model_hyperbola = nlsLM(nd ~ hyperbola(a, b, dist = Phylo_dis),
                 start=list(a = 1, b = 0.01),
                 data = inter_all_c_nd_posi)
model_logistic = nlsLM(nd ~ logistic(a, b, c,
                                     dist = Phylo_dis),
                 start=list(a = 2, b = 1, c = 5),
                 data = inter_all_c_nd_posi)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

plot(dist, max_e(nls_coff_e1[1,1],nls_coff_e1[2,1],dist))
plot(dist, max_e_2(nls_coff_e2[1,1],dist))
#options(error=recover)
model_final = nlme(nd~max_e(a, b, dist = Phylo_dis),
                   start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1]),
                   fixed = a+b~1,
                   groups = ~field/f_p,
                   random = a+b~1|field/f_p,
                   data=inter_all_c_nd_posi, 
                   na.action = na.omit, method = 'REML',
                   control = nlmeControl(maxIter = 100000, minScale = 0.0000001, pnlsTol = 1e-10,
                                         tolerance = 0.000001))
max_eform = ~a*(1-exp(-(b*Phylo_dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
              function.arg=c("Phylo_dist","a","b"))
(nm1 = nlmer(nd ~ max_e_fun(Phylo_dis, a, b) ~ a|field/f_p,
              inter_all_c_nd_posi,
              start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1]),
             nlmerControl(optimizer = "Nelder_Mead", tolPwrss = 1e-100,
                          optCtrl = list())))
summary(model_final)

model_final_2 = nlme(nd~max_e_2(b, dist = Phylo_dis),
                     start = c(b=nls_coff_e2[1,1]),
                     fixed = b~1,
                     groups = ~field/f_p,
                     random = b~1|field/f_p,
                     data=inter_all_c_nd_posi, 
                     na.action = na.omit, method = 'REML',
                     control = nlmeControl(maxIter = 10000))
summary(model_final_2)
anova(model_final, model_final_2)
nlm_coff = model_final$coefficients$fixed
plot(dist, max_e(nlm_coff[1],nlm_coff[2],dist))
nlm_coff_2 = model_final_2$coefficients$fixed
plot(dist, max_e_2(nlm_coff_2[1],dist))

new_dist_1 = data.frame(Phylo_dis = dist, field = factor(0), f_p = factor(0))
prdc = predict_nlme(model_final_2, newdata = new_dist_1,
                    interval = "confidence", level = 0.95)

library(ggplot2)
ggplot() + 
  geom_point(data = inter_all_c_nd_posi, aes(x = Phylo_dis, y = nd)) + 
  geom_line(aes(x = inter_all_c_nd_posi$Phylo_dis, y = prdc[,1])) + 
  geom_ribbon(aes(x = inter_all_c_nd_posi$Phylo_dis,
                  ymin = prdc[,3],
                  ymax = prdc[,4]), 
              fill = "purple", alpha = 0.2)

###### Inla methods #####
pc_prior = list(prec = list("pc.prec",param=c(0.1,0.01)))
inter_all_c_nd_posi$species_i_1 = factor(inter_all_c_nd_posi$species_i)
inter_all_c_nd_posi$species_j_1 = factor(inter_all_c_nd_posi$species_j)
sp_name = unique(inter_all_c_nd_posi$species_i)

##### nd ~ phylo_dist ####
lik_e1 = inlabru::like(family = "gaussian", data = inter_all_c_nd_posi,
  formula = nd ~ a*(1-exp(-(b*Phylo_dis)))+ field + f_p 
  + species_i + species_j + sp_pair + species_i_1 + species_j_1
             )

lik_e2 = inlabru::like(family = "gaussian", data = inter_all_c_nd_posi,
                       formula = nd ~ (1-exp(-(b*Phylo_dis)))+ field + f_p  
                       + species_i + species_j + sp_pair + species_i_1 + species_j_1
)

lik_power = inlabru::like(family = "gaussian", data = inter_all_c_nd_posi,
                       formula = nd ~ a*(Phylo_dis^b)+ field + f_p 
                       + species_i + species_j + sp_pair + species_i_1 + species_j_1
)
lik_hyperbola = inlabru::like(family = "gaussian", data = inter_all_c_nd_posi,
                       formula = nd ~ Phylo_dis/(a+b*Phylo_dis)+ field + f_p 
                       + species_i + species_j + sp_pair + species_i_1 + species_j_1
)

cmp = ~ a(1) + b(1) + 
        field(1, model="iid", hyper = pc_prior)+
        f_p(1, model="iid", hyper = pc_prior) +
        species_i(1, model="iid", hyper = pc_prior)+
        species_j(1, model="iid", hyper = pc_prior)+ 
        sp_pair(1, model="iid", hyper = pc_prior)+
        species_i_1(1, model="generic0", Cmatrix=vcv_tree_sparse, values = sp_name, hyper = pc_prior) + 
        species_j_1(1, model="generic0", Cmatrix=vcv_tree_sparse, values = sp_name, hyper = pc_prior)
cmp_e2 = ~ b(1) + 
  field(1, model="iid", hyper = pc_prior)+
  f_p(1, model="iid", hyper = pc_prior) +
  species_i(1, model="iid", hyper = pc_prior)+
  species_j(1, model="iid", hyper = pc_prior)+ 
  sp_pair(1, model="iid", hyper = pc_prior)+
  species_i_1(1, model="generic0", Cmatrix=vcv_tree_sparse, values = sp_name, hyper = pc_prior) + 
  species_j_1(1, model="generic0", Cmatrix=vcv_tree_sparse, values = sp_name, hyper = pc_prior)

fit_e1 = bru(components = cmp, lik_e1, 
             options = list(bru_initial = list(a = nls_coff_e1[1,1],
                                               b = nls_coff_e1[2,1]),
                       control.compute = list(dic=T, waic=T, cpo=T),
                       quantiles=c(0.025,0.5,0.975)))
summary(fit_e1)

fit_e2 = bru(components = cmp, lik_e2, 
             options = list(bru_initial = list(b = nls_coff_e2[1,1]),
                            control.compute = list(dic=T, waic=T, cpo=T),
                            quantiles=c(0.025,0.5,0.975)))
summary(fit_e2)
fit_e2$summary.fixed

fit_power = bru(components = cmp, lik_power, 
             options = list(bru_initial = list(a = nls_coff_power[1,1], 
                                              b = nls_coff_power[2,1]),
                            control.compute = list(dic=T, waic=T, cpo=T),
                            quantiles=c(0.025,0.5,0.975)))
summary(fit_power)

fit_hyperbola = bru(components = cmp, lik_hyperbola, 
             options = list(bru_initial = list(a = nls_coff_hyperbola[1,1], 
                                               b = nls_coff_hyperbola[2,1]),
                            control.compute = list(dic=T, waic=T, cpo=T),
                            quantiles=c(0.025,0.5,0.975)))
summary(fit_hyperbola)
c(fit_e1$waic$waic, fit_power$waic$waic, fit_hyperbola$waic$waic)


mod_nd_pd_fixed = fit_hyperbola$summary.fixed
plot(dist, power(mod_nd_pd_fixed[1,1],mod_nd_pd_fixed[2,1],dist))
plot(dist, hyperbola(mod_nd_pd_fixed[1,1],mod_nd_pd_fixed[2,1],dist))
plot(dist, max_e(mod_nd_pd_fixed[1,1],mod_nd_pd_fixed[2,1],dist))
plot(dist, max_e(0.156,0.765,dist))
plot(dist, max_e_2(mod_nd_pd_fixed[1,1],dist))
plot(dist, power(nls_coff_power[1,1],nls_coff_power[2,1],dist))
plot(dist, inter_all_c_nd_posi$nd)

pred_hyperbola = predict(fit_hyperbola, dist, ~ dist/(a+b*dist), n.samples = 1000)
pred_hyperbola
hnline = data.frame(distance = dist, p = pred_hyperbola$mean,
                     lower = pred_hyperbola$q0.025, upper = pred_hyperbola$q0.975)

ggplot() +
  geom_point(aes(x = inter_all_c_nd_posi$Phylo_dis, y = inter_all_c_nd_posi$nd), size = 3) +
  geom_line(data = hnline, aes(x = distance, y = p)) +
  geom_ribbon(
    data = hnline, aes(x = distance, ymin = lower, ymax = upper),
    alpha = 0.2, lty = 2
  )

ggplot() +
  gg(pred_hyperbola) +
  geom_point(data = cd, aes(x = x, y = count / exposure), cex = 2) +
  geom_point(data = true.lambda, aes(x, y), pch = "_", cex = 9, col = "blue") +
  coord_cartesian(xlim = c(0, 55), ylim = c(0, 6)) +
  xlab("x") +
  ylab("Intensity")

##### ablgfd ~ phylo_dist ####
lik_ablgfd = inlabru::like(family = "gaussian", data = inter_all_c_nd_posi,
                    formula = ablgfd ~ a*(1-exp(-(b*Phylo_dis)))+species_i + species_j +
                      sp_pair + field + f_p + species_i_1 + species_j_1
)
fit_fd = bru(components = cmp, lik_ablgfd, 
          options = list(bru_initial = list(a = 0.05, b = 21)),
          )
summary(fit_fd)


### Below model run in linux due to RAM limited!
#### nd ~ pd + fd ####
pc_prior=list(prec=list("pc.prec",param=c(0.1,0.01)))
gammaprior=list(prec=list(prior="loggamma",param=c(0.01,0.01)))
wishart1d = list(theta = list(prior = "wishart1d", 
                              param =c(0.01,0.01, 1)))
inter_all_c$species_i_1 = factor(inter_all_c$species_i)
inter_all_c$species_j_1 = factor(inter_all_c$species_j)
sp_name = unique(inter_all_c$species_i)
formula = nd ~ Phylo_dis + Multi_traits + 
  f(species_i, model="iid", hyper = pc_prior)+
  f(species_j, model="iid", hyper = pc_prior)+ 
  f(sp_pair, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) +
  f(species_i_1, model="generic0",
    Cmatrix=vcv_tree_sparse, values = sp_name, 
    hyper = pc_prior) + 
  f(species_j_1, model="generic0",
    Cmatrix=vcv_tree_sparse, values = sp_name,
    hyper = pc_prior)

mod_nd_pd_fd = inla(formula,
                    control.compute = list(dic=T, waic=T, cpo=T),
                    quantiles=c(0.025,0.5,0.975), data=inter_all_c)
mod_nd_pd_fd = mod_nd_pd
summary(mod_nd_pd_fd)
mod_nd_pd_fd$summary.fixed

##### lgfd ~ pd +fd ####
# Just retain one pair of speices whose lgfd > 0 at the least phygenetic distance
inter_all_c_sp_pair = split(inter_all_c, inter_all_c$sp_pair)
inter_all_c_forfd = lapply(inter_all_c_sp_pair,
                           function(x){y = x %>% arrange(Phylo_dis)
                           if (y[1,]$lgfd > 0) {return(y)} else (return(NULL))
                           })
inter_all_c_forfd = rbindlist(inter_all_c_forfd)
length(unique(inter_all_c_forfd$sp_pair))
length(unique(inter_all_c$sp_pair))

mod_lgfd_pd_fd = inla(lgfd ~ Phylo_dis + Multi_traits + f(species_i, model="iid", hyper = gammaprior)+
                        f(species_j, model="iid", hyper = gammaprior)+ 
                        f(sp_pair, model="iid", hyper = gammaprior)+
                        f(field, model="iid", hyper = gammaprior)+
                        f(f_p, model="iid", hyper = gammaprior) +
                        f(species_i_1, model="generic0",
                          Cmatrix=vcv_tree_sparse, values = sp_name, 
                          hyper = gammaprior) + 
                        f(species_j_1, model="generic0",
                          Cmatrix=vcv_tree_sparse, values = sp_name,
                          hyper = gammaprior),
                      control.compute = list(dic=T,waic=T,cpo=T),
                      quantiles=c(0.025,0.5,0.975), data=inter_all_c_forfd)
summary(mod_lgfd_pd_fd)
mod_lgfd_pd_fd$summary.fixed

#### Plot nd lgfd ~ pd+fd ####
load("D:/R projects/BSS/results/fit_results/BSS_exculde_trees_raw/mod_nd_pd_fd_n_a_2_phy_combs.rdata")
load("D:/R projects/BSS/results/fit_results/BSS_exculde_trees_raw/mod_nd_pd_fd_n_a_2_trait_combs.rdata")
load("D:/R projects/BSS/results/fit_results/BSS_exculde_trees_raw/mod_ablgfd_pd_fd_n_a_2_phy_combs.rdata")
load("D:/R projects/BSS/results/fit_results/BSS_exculde_trees_raw/mod_ablgfd_pd_fd_n_a_2_trait_combs.rdata")
mod_nd_pd_fd_fixed = mod_nd_pd_fd$summary.fixed
mod_ablgfd_pd_fd_fixed = mod_ablgfd_pd_fd$summary.fixed
#### Draw predicted line and real data #######
mod_nd_pd_fd_n_a_2_trait_combs$summary.fixed
newdata_for_phy = data.frame(Phylo_dis=seq(min(inter_all_c_n_a$Phylo_dis),
                                           max(inter_all_c_n_a$Phylo_dis),length=100),
                             Multi_traits=mean(inter_all_c_n_a$Multi_traits))
newdata_for_trait = data.frame(Multi_traits=seq(min(inter_all_c_n_a$Multi_traits),
                                                max(inter_all_c_n_a$Multi_traits),length=100),
                               Phylo_dis=mean(inter_all_c_n_a$Phylo_dis))

nd_predicted.trait = mod_nd_pd_fd_n_a_2_trait_combs$summary.lincomb.derived
nd_newdata.lmm_trait = data.frame(newdata_for_trait,
                                  mean=nd_predicted.trait$mean,
                                  lower=nd_predicted.trait$'0.025quant',
                                  upper=nd_predicted.trait$'0.975quant')
nd_predicted.phy = mod_nd_pd_fd_n_a_2_phy_combs$summary.lincomb.derived
nd_newdata.lmm_phy = data.frame(newdata_for_phy,
                                mean=nd_predicted.phy$mean,
                                lower=nd_predicted.phy$'0.025quant',
                                upper=nd_predicted.phy$'0.975quant')

mod_ablgfd_pd_fd_n_a_2_trait_combs$summary.fixed
ablgfd_predicted.trait = mod_ablgfd_pd_fd_n_a_2_trait_combs$summary.lincomb.derived
ablgfd_newdata.lmm_trait = data.frame(newdata_for_trait,
                                      mean=ablgfd_predicted.trait$mean,
                                      lower=ablgfd_predicted.trait$'0.025quant',
                                      upper=ablgfd_predicted.trait$'0.975quant')
ablgfd_predicted.phy = mod_ablgfd_pd_fd_n_a_2_phy_combs$summary.lincomb.derived
ablgfd_newdata.lmm_phy = data.frame(newdata_for_phy,
                                    mean=ablgfd_predicted.phy$mean,
                                    lower=ablgfd_predicted.phy$'0.025quant',
                                    upper=ablgfd_predicted.phy$'0.975quant')

(plo_nd_phy = ggplot()+
    geom_point(position = position_jitter(width = 2, height = 0),
               shape = 1, alpha = 0.7, color = 'grey', size = 25,
               data=inter_all_c_n_a, aes(x=Phylo_dis, y = nd))+
    geom_line(data=nd_newdata.lmm_phy,aes(x=Phylo_dis,y=mean), linetype = "dashed",
              size = 12, color = 'blue')+
    geom_ribbon(data=nd_newdata.lmm_phy,
                aes(x=Phylo_dis,ymin=lower,ymax=upper),fill="purple",alpha=0.3)+
    labs(x = NULL, y = 'Miche difference')+
    theme_custom())

(plo_nd_trait = ggplot()+
    geom_point(position = position_jitter(width = 0.001, height = 0),
               shape = 1, alpha = 0.7, color = 'grey', size = 25,
               data=inter_all_c_n_a,aes(x=Multi_traits, y = nd))+
    geom_line(data=nd_newdata.lmm_trait,aes(x=Multi_traits,y=mean), linetype = "dashed",
              size = 12, color = 'blue')+
    geom_ribbon(data=nd_newdata.lmm_trait,
                aes(x=Multi_traits,ymin=lower,ymax=upper),fill="purple",alpha=0.3)+
    labs(x = NULL, y = NULL)+
    theme_custom())

(plo_ablgfd_phy = ggplot()+
    geom_point(position = position_jitter(width = 2, height = 0),
               shape = 1, alpha = 0.7, color = 'grey', size = 25,
               data=inter_all_c_n_a, aes(x=Phylo_dis, y = ablgfd))+
    geom_line(data=ablgfd_newdata.lmm_phy,aes(x=Phylo_dis,y=mean), linetype = "dashed",
              size = 12, color = 'blue')+
    geom_ribbon(data=ablgfd_newdata.lmm_phy,
                aes(x=Phylo_dis,ymin=lower,ymax=upper),fill="purple",alpha=0.3)+
    labs(x = 'Phylogenetic distance', y = 'abs(log10(CD))')+
    theme_custom())

(plo_ablgfd_trait = ggplot()+
    geom_point(position = position_jitter(width = 0.001, height = 0),
               shape = 1, alpha = 0.7, color = 'grey', size = 25,
               data=inter_all_c_n_a,aes(x=Multi_traits, y = ablgfd))+
    geom_line(data=ablgfd_newdata.lmm_trait,aes(x=Multi_traits,y=mean),
              size = 12, color = 'blue')+
    geom_ribbon(data=ablgfd_newdata.lmm_trait,
                aes(x=Multi_traits,ymin=lower,ymax=upper),fill="purple",alpha=0.3)+
    labs(x = 'Functional distance', y = NULL)+
    theme_custom())

gap = ggplot(NULL)+theme_void()
library(ggpubr)
plo_ndfd_pd_func_d = ggarrange(plo_nd_phy, gap, plo_nd_trait,
                               plo_ablgfd_phy, gap, plo_ablgfd_trait,
                               nrow = 2, ncol = 3,
                               widths = c(1, 0.1, 1),
                               font.label = list(size = 100))

ggsave(plot = plo_ndfd_pd_func_d,
       'results/figures_sameages_raw/plo_ndfd_pd_func_d.jpg',
       width = 200,height = 170,dpi = 300, units = 'cm',
       limitsize = F)

coeffs_nd_pdgower = data.frame(stage = c('PD', 'FD'),
                               nd = c(mod_nd_pd_fd_fixed[1,1],
                                      mod_nd_pd_fd_fixed[2,1]),
                               lt = c(mod_nd_pd_fd_fixed[1,5],
                                      mod_nd_pd_fd_fixed[2,5]),
                               ut = c(mod_nd_pd_fd_fixed[1,3],
                                      mod_nd_pd_fd_fixed[2,3]))

coeffs_ablgfd_pdgower = data.frame(stage = c('PD', 'FD'),
                                   lgfd = c(mod_ablgfd_pd_fd_fixed[1,1],
                                            mod_ablgfd_pd_fd_fixed[2,1]),
                                   lt = c(mod_ablgfd_pd_fd_fixed[1,5],
                                          mod_ablgfd_pd_fd_fixed[2,5]),
                                   ut = c(mod_ablgfd_pd_fd_fixed[1,3],
                                          mod_ablgfd_pd_fd_fixed[2,3]))
library(ggplot2)
library(ggplot2)
p_nd_pdgower = ggplot(coeffs_nd_pdgower, aes(x=nd, y=stage)) + 
  geom_pointrange(aes(xmin=lt, xmax=ut), linewidth = 2,
                  fatten = 40) + 
  geom_vline(aes(xintercept=0), colour="#BB0000", linetype="dashed", size = 2) +
  ggtitle('ND')+
  labs(x = 'Coefficients estimated', y = '') +
  scale_color_manual(labels=c(0, 1), values = c(2, 4))+
  scale_x_continuous(limits = c(-0.0015, 0.003), breaks = seq(-0.002, 0.003, by = 0.001))+
  theme(plot.title = element_text(hjust = 0.5,
                                  vjust=0.5, size=150, face='bold'), 
        plot.title.position = "plot",
        panel.background = element_rect(fill = 'white', color = 'black'),
        legend.key = element_rect(fill = "white"),
        panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        text = element_text(size = 150),
        axis.ticks.y = element_line(linetype=1,color="black",size=1),
        axis.ticks.x = element_line(linetype=1,color="black",size=1),
        axis.ticks.length.y = unit(1,'cm'),
        axis.ticks.length.x = unit(1,'cm'),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.7, vjust = 0.7,
                                   color = 'black'),
        axis.text.y = element_text(angle = 0, vjust = 0.5,color = 'black'),
        axis.title.x = element_text(margin=margin(t=20),
                                    color = 'black'))

p_ablgfd_pdgower = ggplot(coeffs_ablgfd_pdgower, aes(x=lgfd, y=stage)) + 
  geom_pointrange(aes(xmin=lt, xmax=ut), linewidth = 2,
                  fatten = 40) + 
  geom_vline(aes(xintercept=0), colour="#BB0000", linetype="dashed", size = 2)+
  ggtitle('abs(log10(FD))')+
  labs(x = 'Coefficients estimated', y = '') +
  scale_color_manual(labels=c(0, 1), values = c(2, 4))+
  #ylim(1, 3) +
  theme(plot.title = element_text(hjust = 0.5,
                                  vjust=0.5, size=150, face='bold'), 
        plot.title.position = "plot",
        panel.background = element_rect(fill = 'white', color = 'black'),
        legend.key = element_rect(fill = "white"),
        panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        text = element_text(size = 150),
        axis.ticks.y = element_line(linetype=1,color="black",size=1),
        axis.ticks.x = element_line(linetype=1,color="black",size=1),
        axis.ticks.length.y = unit(1,'cm'),
        axis.ticks.length.x = unit(1,'cm'),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.7, vjust = 0.7,color="black"),
        axis.text.y = element_text(angle = 0, vjust = 0.5,color="black"),
        axis.title.x = element_text(margin=margin(t=20),
                                    color = 'black'))
library(ggpubr)
p_ndfd_pdgower = ggarrange(p_nd_pdgower, p_ablgfd_pdgower,
                           nrow = 1, ncol = 2,
                           labels = c('a)', 'b)'),
                           font.label = list(size = 100))

ggsave(plot = p_ndfd_pdgower,
       'results/figures_sameages_raw/p_ndfd_pdgower.svg',
       width = 200,height = 80,dpi = 300, units = 'cm',
       limitsize = F)
ggsave(plot = p_ndfd_pdgower,
       'results/figures_sameages_raw/p_ndfd_pdgower_forppt.svg',
       width = 180,height = 80,dpi = 300, units = 'cm',
       limitsize = F)

######## Analyse invasion success ########
inter_all_c$ra_m_real_t_i = as.numeric(inter_all_c$ra_m_real_t_i)
inter_all_c$ra_m_real_t_j = as.numeric(inter_all_c$ra_m_real_t_j)
inter_all_c_l_2 = split(inter_all_c, inter_all_c$f_p)

##### Species level #####
dat_suc_sp = data.frame()
for (i in 1:length(inter_all_c_l_2)) {
  #i = 1
  trans_plot = inter_all_c_l_2[[i]]
  inv_suc = trans_plot %>% filter(stage_i != 'native' & stage_j == 'native') 
  inv_sp = unique(inv_suc$species_i)
  for (j in 1:length(inv_sp)) {
    inv_suc_sp = inv_suc %>% filter(species_i == inv_sp[j])
    mnd = mean(inv_suc_sp$nd)
    mlgfd = mean(inv_suc_sp$lgfd)
    mablgfd = mean(inv_suc_sp$ablgfd)
    mpd = mean(inv_suc_sp$Phylo_dis)
    mfunc_d = mean(inv_suc_sp$Multi_traits)
    mconti_func_d = mean(inv_suc_sp$Multi_conti_traits)
    mgrowth = mean(inv_suc_sp$growth)
    mspan = mean(inv_suc_sp$span)
    mpollination = mean(inv_suc_sp$pollination)
    mdispersal = mean(inv_suc_sp$dispersal)
    mclonality = mean(inv_suc_sp$clonality)
    mheight = mean(inv_suc_sp$height)
    mldmc = mean(inv_suc_sp$ldmc)
    msla = mean(inv_suc_sp$sla)
    mseedmass = mean(inv_suc_sp$seedmass)
    mnd.a = sum(inv_suc_sp$nd*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mlgfd.a = sum(inv_suc_sp$lgfd*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mablgfd.a = sum(inv_suc_sp$ablgfd*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mpd.a = sum(inv_suc_sp$Phylo_dis*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mfunc_d.a = sum(inv_suc_sp$Multi_traits*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mconti_func_d.a = sum(inv_suc_sp$Multi_conti_traits*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mgrowth.a = sum(inv_suc_sp$growth*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mspan.a = sum(inv_suc_sp$span*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mpollination.a = sum(inv_suc_sp$pollination*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mdispersal.a = sum(inv_suc_sp$dispersal*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mclonality.a = sum(inv_suc_sp$clonality*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mheight.a = sum(inv_suc_sp$height*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mldmc.a = sum(inv_suc_sp$ldmc*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    msla.a = sum(inv_suc_sp$sla*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mseedmass.a = sum(inv_suc_sp$seedmass*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mnnd = min(inv_suc_sp$nd)
    mnlgfd = min(inv_suc_sp$lgfd)
    mnablgfd = min(inv_suc_sp$ablgfd)
    mntd = min(inv_suc_sp$Phylo_dis)
    mnfunc_d = min(inv_suc_sp$Multi_traits)
    mnconti_func_d = min(inv_suc_sp$Multi_conti_traits)
    mngrowth = min(inv_suc_sp$growth)
    mnspan = min(inv_suc_sp$span)
    mnpollination = min(inv_suc_sp$pollination)
    mndispersal = min(inv_suc_sp$dispersal)
    mnclonality = min(inv_suc_sp$clonality)
    mnheight = min(inv_suc_sp$height)
    mnldmc = min(inv_suc_sp$ldmc)
    mnsla = min(inv_suc_sp$sla)
    mnseedmass = min(inv_suc_sp$seedmass)
    
    dat_suc_sp_1 = data.frame(f_p = unique(trans_plot$f_p), plot = unique(trans_plot$f_p), field = unique(trans_plot$field),
                              species = inv_sp[j], 
                              stage = unique(inv_suc_sp$stage_i),
                              estab = unique(inv_suc_sp$estab_i),
                              domin = unique(inv_suc_sp$domin_i),
                              mnd = mnd, mlgfd = mlgfd, mpd = mpd, mfunc_d = mfunc_d,
                              mconti_func_d = mconti_func_d, mgrowth = mgrowth, mspan = mspan, 
                              mpollination = mpollination, mdispersal = mdispersal, 
                              mclonality = mclonality, mheight = mheight, mldmc = mldmc,
                              msla = msla, mseedmass = mseedmass, 
                              mnd.a = mnd.a, mlgfd.a = mlgfd.a, mablgfd.a, mpd.a = mpd.a,
                              mfunc_d.a = mfunc_d.a, mconti_func_d.a = mconti_func_d.a,
                              mgrowth.a = mgrowth.a, mspan.a = mspan.a,mpollination.a = mpollination.a,
                              mdispersal.a = mdispersal.a, mclonality.a = mclonality.a,
                              mheight.a = mheight.a, mldmc.a = mldmc.a, msla.a = msla.a,
                              mseedmass.a = mseedmass.a, 
                              mnnd = mnnd, mnlgfd = mnlgfd, mntd = mntd,
                              mnfunc_d = mnfunc_d, mnconti_func_d = mnconti_func_d, 
                              mngrowth = mngrowth, mnspan = mnspan, 
                              mnpollination = mnpollination, mndispersal = mndispersal, 
                              mnclonality = mnclonality, mnheight = mnheight, mnldmc = mnldmc,
                              mnsla = mnsla, mnseedmass = mnseedmass)
    dat_suc_sp = rbind(dat_suc_sp, dat_suc_sp_1)
  }
}

dat_suc_sp$stage_level = NA
dat_suc_sp[dat_suc_sp$stage == 'introduce',]$stage_level = 1
dat_suc_sp[dat_suc_sp$stage == 'establish',]$stage_level = 2
dat_suc_sp[dat_suc_sp$stage == 'dominant',]$stage_level = 3
dat_suc_sp$stage_level = ordered(dat_suc_sp$stage_level)
str(dat_suc_sp)
save(dat_suc_sp,
     file = 'code/results_analyzing/analysing_sameages_raw_data/dat_suc_sp.rdata')

############ Fast start for analyzing invasion success probability ~ nd+fd #########
load('code/results_analyzing/analysing_sameages_raw_data/dat_suc_sp.rdata')
library(ordinal)
library(lme4)
library(INLA)
library(brinla)
library(brms)
library(ape)

numcols = grep("^m",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])
dat_suc_sps$species_1 = as.factor(dat_suc_sps$species)
estab_sp_names = unique(dat_suc_sps$species)

tree = read.tree('data/original data/phylo_tree332.txt')
estab_tree_fit = keep.tip(tree, estab_sp_names)
estab_vcv_tree = ape::vcv(estab_tree_fit, model = "Brownian", corr = FALSE)
estab_vcv_tree_sparse = inla.as.sparse(solve(estab_vcv_tree))
dat_suc_sps$species_1 = dat_suc_sps$species

## Check the co-linearity
car::vif(
  glmer(estab ~ mnd.a + mlgfd.a + (1|species) + (1|f_p),
        family=binomial,data=dat_suc_sps)
)
## VIF < 3
with(dat_suc_sps,cor(mnd.a,mlgfd.a))

## VIF < 3
with(dat_suc_sps,cor(mnd.a,mlgfd.a))

# cor_coefficients = -0.56<0.7
### Analyse for mnd mfd abundance weighted mean #############
#### Establishment
pc_prior = list(prec=list("pc.prec", param=c(0.1,0.01)))

mod_suc_sp_ab_md.a = inla(estab ~ mnd.a+mlgfd.a+mnd.a:mlgfd.a+
                            f(species, model = "iid", hyper = pc_prior)+
                            f(field, model="iid", hyper = pc_prior) +
                            f(f_p, model="iid", hyper = pc_prior)+
                            f(species_1, model="generic0",
                              Cmatrix= estab_vcv_tree_sparse,
                              values = estab_sp_names, hyper=pc_prior),
                          control.compute = list(dic=T,waic=T,cpo=T,
                                                 config = TRUE),
                          quantiles=c(0.025,0.5,0.975),
                          family="binomial", data=dat_suc_sps)
summary(mod_suc_sp_ab_md.a)

mod_suc_sp_ab_md = inla(estab ~ mnd + mlgfd + mnd:mlgfd+
                          f(species, model = "iid", hyper = pc_prior)+
                          f(field, model="iid", hyper = pc_prior) +
                          f(f_p, model="iid", hyper = pc_prior) + 
                          f(species_1, model="generic0",
                            Cmatrix= estab_vcv_tree_sparse,
                            values = estab_sp_names, hyper=pc_prior),
                        control.compute = list(dic=T,waic=T,cpo=T,
                                               config = TRUE),
                        quantiles=c(0.025,0.5,0.975),
                        family="binomial", data=dat_suc_sps)
summary(mod_suc_sp_ab_md)

mod_suc_sp_ab_mnd = inla(estab ~ mnnd + mnlgfd +
                           f(species, model = "iid", hyper = pc_prior)+
                           f(field, model="iid", hyper = pc_prior) +
                           f(f_p, model="iid", hyper = pc_prior)+
                           f(species_1, model="generic0",
                             Cmatrix= estab_vcv_tree_sparse,
                             values = estab_sp_names, hyper=pc_prior),
                         control.compute = list(dic=T,waic=T,cpo=T,
                                                config = TRUE),
                         quantiles=c(0.025,0.5,0.975),
                         family="binomial", data=dat_suc_sps)
c(mod_suc_sp_ab_md$waic$waic, mod_suc_sp_ab_md.a$waic$waic,
  mod_suc_sp_ab_mnd$waic$waic) 
# choose mnd.a because there is just 13 differences and the calculation way is more reasonable
summary(mod_suc_sp_ab_md)
summary(mod_suc_sp_ab_md.a)
summary(mod_suc_sp_ab_mnd)
mod_suc_sp_ab_md_fixed = mod_suc_sp_ab_md$summary.fixed
mod_suc_sp_ab_md.a_fixed = mod_suc_sp_ab_md.a$summary.fixed
mod_suc_sp_ab_mnd_fixed = mod_suc_sp_ab_mnd$summary.fixed

mod_suc_sp_ab_mpd_mfunc_d = inla(estab ~ mpd + mconti_func_d +
                           f(species, model = "iid", hyper = pc_prior)+
                           f(field, model="iid", hyper = pc_prior) +
                           f(f_p, model="iid", hyper = pc_prior)+
                           f(species_1, model="generic0",
                             Cmatrix= estab_vcv_tree_sparse,
                             values = estab_sp_names, hyper=pc_prior),
                         control.compute = list(dic=T,waic=T,cpo=T,
                                                config = TRUE),
                         quantiles=c(0.025,0.5,0.975),
                         family="binomial", data=dat_suc_sps)
mod_suc_sp_ab_mpd.a_mfunc_d.a = inla(estab ~ mpd.a + mconti_func_d.a +
                                   f(species, model = "iid", hyper = pc_prior)+
                                   f(field, model="iid", hyper = pc_prior) +
                                   f(f_p, model="iid", hyper = pc_prior)+
                                   f(species_1, model="generic0",
                                     Cmatrix= estab_vcv_tree_sparse,
                                     values = estab_sp_names, hyper=pc_prior),
                                 control.compute = list(dic=T,waic=T,cpo=T,
                                                        config = TRUE),
                                 quantiles=c(0.025,0.5,0.975),
                                 family="binomial", data=dat_suc_sps)

mod_suc_sp_ab_mntd_mfunc_d = inla(estab ~ mntd + mnconti_func_d +
                                   f(species, model = "iid", hyper = pc_prior)+
                                   f(field, model="iid", hyper = pc_prior) +
                                   f(f_p, model="iid", hyper = pc_prior)+
                                   f(species_1, model="generic0",
                                     Cmatrix= estab_vcv_tree_sparse,
                                     values = estab_sp_names, hyper=pc_prior),
                                 control.compute = list(dic=T,waic=T,cpo=T,
                                                        config = TRUE),
                                 quantiles=c(0.025,0.5,0.975),
                                 family="binomial", data=dat_suc_sps)
c(mod_suc_sp_ab_mpd_mfunc_d$waic$waic, mod_suc_sp_ab_mpd.a_mfunc_d.a$waic$waic,
  mod_suc_sp_ab_mntd_mfunc_d$waic$waic)

c(mod_suc_sp_ab_mpd.a_mfunc_d.a$waic$waic, mod_suc_sp_ab_md.a$waic$waic)
## [1] 2677.738 2597.872

##### Dominant
dat_dom_sp = dat_suc_sps %>% filter(stage %in% c('establish', 'dominant'))

domin_sp_names = unique(dat_dom_sp$species)

tree = read.tree('data/original data/phylo_tree332.txt')
domin_tree_fit = keep.tip(tree, domin_sp_names)
domin_vcv_tree = ape::vcv(domin_tree_fit, model = "Brownian", corr = FALSE)
domin_vcv_tree_sparse = inla.as.sparse(solve(domin_vcv_tree))
dat_dom_sp$species_1 = dat_dom_sp$species

## Check the co-linearity
car::vif(
  glmer(domin ~ mnd.a + mlgfd.a + (1|species) + (1|f_p),
        family=binomial, data=dat_dom_sp)
)
## VIF < 3
with(dat_dom_sp,cor(mnd.a,mlgfd.a))
# cor_coefficients = -0.44<0.7

mod_dom_sp_ab_md.a = inla(domin ~ mnd.a+ mlgfd.a+
                            f(species, model = "iid", hyper = pc_prior)+
                            f(field, model="iid", hyper = pc_prior) +
                            f(f_p, model="iid", hyper = pc_prior)+
                            f(species_1, model="generic0",
                              Cmatrix= domin_vcv_tree_sparse,
                              values = domin_sp_names, hyper=pc_prior),
                          control.compute = list(dic=T,waic=T,cpo=T,
                                                 config = TRUE),
                          quantiles=c(0.025,0.5,0.975),
                          family="binomial", data=dat_dom_sp)
mod_dom_sp_ab_md = inla(domin ~ mnd+ mlgfd+
                          f(species, model = "iid", hyper = pc_prior)+
                          f(field, model="iid", hyper = pc_prior) +
                          f(f_p, model="iid", hyper = pc_prior)+
                          f(species_1, model="generic0",
                            Cmatrix= domin_vcv_tree_sparse,
                            values = domin_sp_names, hyper=pc_prior),
                        control.compute = list(dic=T,waic=T,cpo=T,
                                               config = TRUE),
                        quantiles=c(0.025,0.5,0.975),
                        family="binomial", data=dat_dom_sp)
mod_dom_sp_ab_mnd = inla(domin ~ mnnd+ mnlgfd+
                           f(species, model = "iid", hyper = pc_prior)+
                           f(field, model="iid", hyper = pc_prior) +
                           f(f_p, model="iid", hyper = pc_prior)+
                           f(species_1, model="generic0",
                             Cmatrix= domin_vcv_tree_sparse,
                             values = domin_sp_names, hyper=pc_prior),
                         control.compute = list(dic=T,waic=T,cpo=T,
                                                config = TRUE),
                         quantiles=c(0.025,0.5,0.975),
                         family="binomial", data=dat_dom_sp)
c(mod_dom_sp_ab_md$waic$waic, mod_dom_sp_ab_md.a$waic$waic,
  mod_dom_sp_ab_mnd$waic$waic) # choose mnd.a, same reason with establishment 
summary(mod_dom_sp_ab_md)
summary(mod_dom_sp_ab_md.a)
summary(mod_dom_sp_ab_mnd)
mod_dom_sp_ab_md_fixed = mod_dom_sp_ab_md$summary.fixed
mod_dom_sp_ab_md.a_fixed = mod_dom_sp_ab_md.a$summary.fixed
mod_dom_sp_ab_mnd_fixed = mod_dom_sp_ab_mnd$summary.fixed

### Phy and Fun analyse
mod_dom_sp_ab_mpd_mfunc_d = inla(domin ~ mpd + mconti_func_d +
                                   f(species, model = "iid", hyper = pc_prior)+
                                   f(field, model="iid", hyper = pc_prior) +
                                   f(f_p, model="iid", hyper = pc_prior)+
                                   f(species_1, model="generic0",
                                     Cmatrix= domin_vcv_tree_sparse,
                                     values = domin_sp_names, hyper=pc_prior),
                                 control.compute = list(dic=T,waic=T,cpo=T,
                                                        config = TRUE),
                                 quantiles=c(0.025,0.5,0.975),
                                 family="binomial", data=dat_dom_sp)
mod_dom_sp_ab_mpd.a_mfunc_d.a = inla(domin ~ mpd.a + mconti_func_d.a +
                                       f(species, model = "iid", hyper = pc_prior)+
                                       f(field, model="iid", hyper = pc_prior) +
                                       f(f_p, model="iid", hyper = pc_prior)+
                                       f(species_1, model="generic0",
                                         Cmatrix= domin_vcv_tree_sparse,
                                         values = domin_sp_names, hyper=pc_prior),
                                     control.compute = list(dic=T,waic=T,cpo=T,
                                                            config = TRUE),
                                     quantiles=c(0.025,0.5,0.975),
                                     family="binomial", data=dat_dom_sp)

mod_dom_sp_ab_mntd_mfunc_d = inla(domin ~ mntd + mnconti_func_d +
                                    f(species, model = "iid", hyper = pc_prior)+
                                    f(field, model="iid", hyper = pc_prior) +
                                    f(f_p, model="iid", hyper = pc_prior)+
                                    f(species_1, model="generic0",
                                      Cmatrix= domin_vcv_tree_sparse,
                                      values = domin_sp_names, hyper=pc_prior),
                                  control.compute = list(dic=T,waic=T,cpo=T,
                                                         config = TRUE),
                                  quantiles=c(0.025,0.5,0.975),
                                  family="binomial", data=dat_dom_sp)
c(mod_dom_sp_ab_mpd_mfunc_d$waic$waic, mod_dom_sp_ab_mpd.a_mfunc_d.a$waic$waic,
  mod_dom_sp_ab_mntd_mfunc_d$waic$waic)

## Compare model for ndcd and pdfd
c(mod_dom_sp_ab_mpd_mfunc_d$waic$waic, mod_dom_sp_ab_md.a$waic$waic)
# [1] 1400.040 1362.904

#### All stages
# Normal method
cor.test(dat_suc_sps$mpd, dat_suc_sps$mlgfd)
mod_sp_all_ab = clm(stage_level ~ mnd.a + mlgfd.a,
                    data = dat_suc_sps)
mod_sp_all_ab_1 = clmm(stage_level ~ mnd.a + mlgfd.a+ (1|species)+
                         (1|ffield/fplot),
                       data = dat_suc_sps)
summary(mod_sp_all_ab)
summary(mod_sp_all_ab_1)

# Bayesian method given the phylogenetic independence 
options(mc.cores = parallel::detectCores(logical = F)) 
mod_sp_all_md_brms = brm(stage_level ~ mnd + mlgfd + (1|species)+
                           (1|field/plot), data=dat_suc_sps,
                         cov_ranef = list(species = estab_vcv_tree), # phylo
                         family=cumulative("logit", link_disc = "log",
                                           threshold = "flexible"),
                         chains = 4,warmup = 1000, iter = 3000, refresh = 100,   
                         control = list(max_treedepth = 20,
                                        adapt_delta = 0.99), cores = 10)
save(mod_sp_all_md_brms,
     file = 'code/results_analyzing/analysing_sameages_raw_data/mod_sp_all_md_brms.rdata')

mod_sp_all_md.a_brms = brm(stage_level ~ mnd.a + mlgfd.a + (1|species)+
                             (1|field/plot), data=dat_suc_sps,
                           cov_ranef = list(species = estab_vcv_tree), # phylo
                           family=cumulative("logit", link_disc = "log",
                                             threshold = "flexible"),
                           chains = 4,warmup = 1000, iter = 3000, refresh = 100,   
                           control = list(max_treedepth = 20,
                                          adapt_delta = 0.99), cores = 10)
save(mod_sp_all_md.a_brms,
     file = 'code/results_analyzing/analysing_sameages_raw_data/mod_sp_all_md.a_brms.rdata')

mod_sp_all_mnd_brms = brm(stage_level ~ mnnd + mnlgfd + (1|species)+
                            (1|field/plot), data=dat_suc_sps,
                          cov_ranef = list(species = estab_vcv_tree), # phylo
                          family=cumulative("logit", link_disc = "log",
                                            threshold = "flexible"),
                          chains = 4,warmup = 1000, iter = 3000, refresh = 100,   
                          control = list(max_treedepth = 20,
                                         adapt_delta = 0.99), cores = 10)
save(mod_sp_all_mnd_brms,
     file = 'code/results_analyzing/analysing_sameages_raw_data/mod_sp_all_mnd_brms.rdata')

mod_sp_all_md.a_fixed = fixef(mod_sp_all_md.a_brms,
                              probs = c(0.025, 0.975),
                              summary = T,
                              robust = T)
mod_sp_all_md_fixed = fixef(mod_sp_all_md_brms,
                            probs = c(0.025, 0.975),
                            summary = T,
                            robust = T)
mod_sp_all_mnd_fixed = fixef(mod_sp_all_mnd_brms,
                            probs = c(0.025, 0.975),
                            summary = T,
                            robust = T)
####################### Plot ##########################
# plot invasion success mnd mnlgfd

coeffs_mnd = data.frame(stage = c('establishment', 'dominant', 'overall'),
                        nd = c(mod_suc_sp_ab_md_fixed[2,1],
                               mod_dom_sp_ab_md_fixed[2,1],
                               mod_sp_all_md_fixed[3,1]),
                        ut = c(mod_suc_sp_ab_md_fixed[2,5],
                               mod_dom_sp_ab_md_fixed[2,5],
                               mod_sp_all_md_fixed[3,4]),
                        lt = c(mod_suc_sp_ab_md_fixed[2,3],
                               mod_dom_sp_ab_md_fixed[2,3],
                               mod_sp_all_md_fixed[3,3]))
coeffs_mlgfd = data.frame(stage = c('establishment', 'dominant', 'overall'),
                          lgfd = c(mod_suc_sp_ab_md_fixed[3,1],
                                   mod_dom_sp_ab_md_fixed[3,1],
                                   mod_sp_all_md_fixed[4,1]),
                          ut = c(mod_suc_sp_ab_md_fixed[3,5],
                                 mod_dom_sp_ab_md_fixed[3,5],
                                 mod_sp_all_md_fixed[4,4]),
                          lt = c(mod_suc_sp_ab_md_fixed[3,3],
                                 mod_dom_sp_ab_md_fixed[3,3],
                                 mod_sp_all_md_fixed[4,3]))
library(ggplot2)
(p_mod_nd = ggplot(coeffs_mnd, aes(x=nd, y=stage)) + 
    geom_pointrange(aes(xmin=lt, xmax=ut), linewidth = 3,
                    fatten = 40) + 
    geom_vline(aes(xintercept=0), colour="#BB0000", linetype="dashed", size = 2) +
    ggtitle('MND')+
    labs(x = 'Coefficients estimated', y = '') +
    scale_color_manual(labels=c(0, 1), values = c(2, 4))+
    #ylim(1, 3) +
    scale_y_discrete(limits = c('overall', 'dominant', 'establishment'))+
    theme_custom()+
    theme(axis.ticks = element_line(linewidth = 3),
          axis.ticks.length = unit(-2,'cm')))

(p_mod_lgfd = ggplot(coeffs_mlgfd, aes(x=lgfd, y=stage)) + 
    geom_pointrange(aes(xmin=lt, xmax=ut), linewidth = 3,
                    fatten = 40) + 
    geom_vline(aes(xintercept=0), colour="#BB0000", linetype="dashed", size = 2)+
    ggtitle('MCD')+
    labs(x = 'Coefficients estimated', y = '') +
    scale_color_manual(labels=c(0, 1), values = c(2, 4))+
    scale_y_discrete(limits = c('overall', 'dominant', 'establishment'))+
    theme_custom()+
    theme(axis.ticks = element_line(linewidth = 3),
          axis.ticks.length = unit(-2,'cm')))

# plot invasion success mnd.a mnlgfd.a
coeffs_mnd_ab = data.frame(stage = c('establishment', 'dominant', 'overall'),
                           nd = c(mod_suc_sp_ab_md.a_fixed[2,1],
                                  mod_dom_sp_ab_md.a_fixed[2,1],
                                  mod_sp_all_md.a_fixed[3,1]),
                           ut = c(mod_suc_sp_ab_md.a_fixed[2,5],
                                  mod_dom_sp_ab_md.a_fixed[2,5],
                                  mod_sp_all_md.a_fixed[3,4]),
                           lt = c(mod_suc_sp_ab_md.a_fixed[2,3],
                                  mod_dom_sp_ab_md.a_fixed[2,3],
                                  mod_sp_all_md.a_fixed[3,3]))
coeffs_mlgfd_ab = data.frame(stage = c('establishment', 'dominant', 'overall'),
                             lgfd = c(mod_suc_sp_ab_md.a_fixed[3,1],
                                      mod_dom_sp_ab_md.a_fixed[3,1],
                                      mod_sp_all_md.a_fixed[4,1]),
                             ut = c(mod_suc_sp_ab_md.a_fixed[3,5],
                                    mod_dom_sp_ab_md.a_fixed[3,5],
                                    mod_sp_all_md.a_fixed[4,4]),
                             lt = c(mod_suc_sp_ab_md.a_fixed[3,3],
                                    mod_dom_sp_ab_md.a_fixed[3,3],
                                    mod_sp_all_md.a_fixed[4,3]))
library(ggplot2)
(p_mod_nd_ab = ggplot(coeffs_mnd_ab, aes(x=nd, y=stage)) + 
  geom_pointrange(aes(xmin=lt, xmax=ut), linewidth = 2,
                  fatten = 40) + 
  geom_vline(aes(xintercept=0), colour="#BB0000", linetype="dashed", size = 2) +
  ggtitle(expression(MND[ab]))+
  labs(x = 'Coefficients estimated', y = '') +
  scale_color_manual(labels=c(0, 1), values = c(2, 4))+
  #ylim(1, 3) +
  scale_y_discrete(limits = c('overall', 'dominant', 'establishment'))+
  theme_custom()+
  theme(axis.ticks = element_line(linewidth = 3),
          axis.ticks.length = unit(-2,'cm')))

(p_mod_lgfd_ab = ggplot(coeffs_mlgfd_ab, aes(x=lgfd, y=stage))+ 
  geom_pointrange(aes(xmin=lt, xmax=ut),
                  linewidth = 2,
                  fatten = 40) + 
  geom_vline(aes(xintercept=0), colour="#BB0000",
             linetype="dashed", size = 2)+
  ggtitle(expression(MCD[ab]))+
  labs(x = 'Coefficients estimated', y = '') +
  scale_color_manual(labels=c(0, 1), values = c(2, 4))+
  #ylim(1, 3) +
  scale_y_discrete(limits = c('overall', 'dominant', 'establishment'))+
  theme_custom()+
  theme(axis.ticks = element_line(linewidth = 3),
          axis.ticks.length = unit(-2,'cm')))


# plot invasion success mnnd mnlgfd
coeffs_mnd_mntd = data.frame(stage = c('establishment', 'dominant', 'overall'),
                             nd = c(mod_suc_sp_ab_mnd_fixed[2,1],
                                    mod_dom_sp_ab_mnd_fixed[2,1],
                                    mod_sp_all_mnd_fixed[3,1]),
                             ut = c(mod_suc_sp_ab_mnd_fixed[2,5],
                                    mod_dom_sp_ab_mnd_fixed[2,5],
                                    mod_sp_all_mnd_fixed[3,4]),
                             lt = c(mod_suc_sp_ab_mnd_fixed[2,3],
                                    mod_dom_sp_ab_mnd_fixed[2,3],
                                    mod_sp_all_mnd_fixed[3,3]))

coeffs_mlgfd_mntd = data.frame(stage = c('establishment', 'dominant', 'overall'),
                               lgfd = c(mod_suc_sp_ab_mnd_fixed[3,1],
                                        mod_dom_sp_ab_mnd_fixed[3,1],
                                        mod_sp_all_mnd_fixed[4,1]),
                               ut = c(mod_suc_sp_ab_mnd_fixed[3,5],
                                      mod_dom_sp_ab_mnd_fixed[3,5],
                                      mod_sp_all_mnd_fixed[4,4]),
                               lt = c(mod_suc_sp_ab_mnd_fixed[3,3],
                                      mod_dom_sp_ab_mnd_fixed[3,3],
                                      mod_sp_all_mnd_fixed[4,3]))
library(ggplot2)
(p_mod_nd_mntd = ggplot(coeffs_mnd_mntd, aes(x=nd, y=stage)) + 
  geom_pointrange(aes(xmin=lt, xmax=ut), linewidth = 2,
                  fatten = 40) + 
  geom_vline(aes(xintercept=0), colour="#BB0000", linetype="dashed", size = 2) +
  ggtitle('MNND')+
  labs(x = 'Coefficients estimated', y = '') +
  scale_color_manual(labels=c(0, 1), values = c(2, 4))+
  #ylim(1, 3) +
  scale_y_discrete(limits = c('overall', 'dominant', 'establishment'))+
  theme_custom()+
  theme(axis.ticks = element_line(linewidth = 3),
          axis.ticks.length = unit(-2,'cm')))

(p_mod_lgfd_mntd = ggplot(coeffs_mlgfd_mntd, aes(x=lgfd, y=stage))+ 
  geom_pointrange(aes(xmin=lt, xmax=ut),
                  linewidth = 2,
                  fatten = 40) + 
  geom_vline(aes(xintercept=0), colour="#BB0000",
             linetype="dashed", size = 2)+
  ggtitle('MNCD')+
  labs(x = 'Coefficients estimated', y = '') +
  scale_color_manual(labels=c(0, 1), values = c(2, 4))+
  #ylim(1, 3) +
  scale_y_discrete(limits = c('overall', 'dominant', 'establishment'))+
  theme_custom()+
  theme(axis.ticks = element_line(linewidth = 3),
        axis.ticks.length = unit(-2,'cm')))

library(ggpubr)
Fig.2_ab = ggarrange(p_mod_nd_ab, p_mod_lgfd_ab, nrow = 1, ncol = 2,
                     labels = c('a)', 'b)'), font.label = list(size = 200))
Fig.2_mnd = ggarrange(p_mod_nd, p_mod_lgfd, nrow = 1, ncol = 2,
                      labels = c('a)', 'b)'), hjust = -2, font.label = list(size = 200))
Fig.2_mntd = ggarrange(p_mod_nd_mntd, p_mod_lgfd_mntd, nrow = 1, ncol = 2,
                       labels = c('a)', 'b)'), font.label = list(size = 200))

ggsave(plot = Fig.2_ab,
       'results/figures_sameages_raw/Fig.2_ab.svg',
       width = 200,height = 80, dpi = 300, units = 'cm',
       limitsize = F)
ggsave(plot = Fig.2_mnd,
       'results/figures_sameages_raw/Fig.2_mnd.svg',
       width = 240,height = 80, dpi = 300, units = 'cm',
       limitsize = F)
ggsave(plot = Fig.2_mntd,
       'results/figures_sameages_raw/Fig.2_mnnd.svg',
       width = 200,height = 80, dpi = 300, units = 'cm',
       limitsize = F)


############### Fast start for analyzing mnd ~ mpd, mnd ~ mfunc_d, mfd ~ mpd, mfd ~ mfunc_d ###########
require(INLA)
require(inlabru)
require(tibble)
require(lme4)
require(nlme)
require(lmerTest)
require(minpack.lm)

load('code/results_analyzing/analysing_sameages_raw_data/dat_suc_sp.rdata')
plot(dat_suc_sp$mpd.a, dat_suc_sp$mnd.a)
plot(dat_suc_sp$mpd, dat_suc_sp$mnd)
plot(dat_suc_sp$mpd, dat_suc_sp$mnd)

max_e = function(a, b, dist){ 
  f = a*(1-exp(-(b*dist)))
  return(f)}
max_e_2 = function(b, dist){ 
  f = (1-exp(-(b*dist)))
  return(f)}
power = function(a, b, dist) {
  f = a*(dist^b)
  return(f)
}
hyperbola = function(a, b, dist) {
  f = dist/(a+b*dist)
  return(f)
}
logistic = function(a, b, c, dist) {
  f = a*dist/(1+exp(1)*((b-dist)/c))
}


####### mnd.a ~ mpd.a ######
mpd.a = dat_suc_sp$mpd.a
#plot(dist, max_e(1, 0.01, dist))
#plot(dist, max_e_2(0.002740836, dist))
#plot(dist, power(1, 0.01, dist))
#plot(dist, hyperbola(1, 0.01, dist))
#plot(dist, logistic(a = 2, b = 1, c = 5,
#                    dist))

model_e1 = nlsLM(mnd.a ~ max_e(a, b, dist = mpd.a),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mnd.a ~ max_e_2(b, dist = mpd.a),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mnd.a ~ power(a, b, dist = mpd.a),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mnd.a ~ hyperbola(a, b, dist = mpd.a),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mnd.a ~ logistic(a, b, c,
                                     dist = mpd.a),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

plot(mpd.a, max_e(nls_coff_e1[1,1],nls_coff_e1[2,1],mpd.a))
plot(mpd.a, max_e_2(nls_coff_e2[1,1],mpd.a))
#options(error=recover)

## lm
(mod_lmer = lmer(mnd.a ~ 0 + mpd.a + (1|field/f_p),
                 data = dat_suc_sp, REML = TRUE))
summary(mod_lmer)
ggpredict(mod_lmer, terms = 'mpd.a')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                   function.arg=c("dist","a","b"))
(mod_e = nlmer(mnd.a ~ max_e_fun(dist = mpd.a, a, b) ~ (a|field/f_p) + (a|species),
              dat_suc_sp,
              start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_e2 = nlmer(mnd.a ~ max_e2_fun(dist = mpd.a, b) ~ (b|field/f_p) + (b|species),
                dat_suc_sp,
                start = c(b = nls_coff_e2[1,1])))
summary(mod_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                   function.arg=c("dist", "a", "b"))
(mod_power = nlmer(mnd.a ~ power_fun(dist = mpd.a, a, b) ~ (a|field/f_p) + (a|species),
                dat_suc_sp,
                start = c(a = nls_coff_power[1,1],
                          b = nls_coff_power[2,1])))
summary(mod_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_hyperbola = nlmer(mnd.a ~ hyperbola_fun(dist = mpd.a, a, b) ~ a|field/f_p,
                dat_suc_sp,
                start = c(a = nls_coff_hyperbola[1,1],
                          b = nls_coff_hyperbola[2,1])))
summary(mod_hyperbola)

anova(mod_lmer, mod_e, mod_e2, mod_power, mod_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_lmer, mod_e, mod_e2, mod_power, mod_hyperbola)
require(ggeffects)
ggpredict(mod_lmer, terms = 'mpd.a')
ggpredict(mod_e, terms = 'mpd.a')
mod_e_coff = summary(mod_e)$coefficients

new_dist_1 = data.frame(mpd.a = mpd.a, field = factor(0), f_p = factor(0))
prdc = predict_nlme(mod_e, newdata = new_dist_1,
                    interval = "confidence", level = 0.95)

library(ggplot2)
plo_mnd.a_mpd.a = ggplot() + 
  geom_point(data = dat_suc_sp, aes(x = mpd.a, y = mnd.a),
             shape = 1, alpha = 1, color = 'grey', size = 16) + 
  geom_line(aes(x = seq(0, max(mpd.a), length.out = 100),
                y = max_e(mod_e_coff[1,1],mod_e_coff[2,1],seq(0, max(mpd.a),
                                                              length.out = 100))), 
            col = 'blue', linewidth = 4) + 
  theme_custom() #+ 
  #geom_ribbon(aes(x = inter_all_c_nd_posi$Phylo_dis,
              #    ymin = prdc[,3],
              #    ymax = prdc[,4]), 
           #   fill = "purple", alpha = 0.2)

### Bayesian inla methods to inculde phylogenetic independence
pc_prior = list(prec = list("pc.prec",param=c(0.1,0.01)))
dat_suc_sp$species_1 = factor(dat_suc_sp$species)
sp_name = unique(dat_suc_sp$species)

### Calculate the original functional distance
tree = read.tree('data/original data/phylo_tree332.txt')
sp_md = unique(dat_suc_sp$species)
tree_md = keep.tip(tree, sp_md)
vcv_tree_md = ape::vcv(tree_md, model = "Brownian", corr = FALSE)
vcv_tree_md_sparse = inla.as.sparse(solve(vcv_tree_md))

### lmer ###
formula = mnd.a ~ 0 + mpd.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
   hyper = pc_prior)

mod_mnd.a_mpd.a = inla(formula,
                    control.compute = list(dic=T, waic=T, cpo=T),
                    quantiles=c(0.025,0.5,0.975), data=dat_suc_sp)

summary(mod_mnd.a_mpd.a)
mod_mnd.a_mpd.a$summary.fixed

### nlmer ###
lik_e1 = inlabru::like(family = "gaussian", data = dat_suc_sp,
                       formula = mnd.a ~ a*(1-exp(-(b*mpd.a)))+ field + f_p
                       + species + species_1)
lik_e2 = inlabru::like(family = "gaussian", data = dat_suc_sp,
                       formula = mnd.a ~ (1-exp(-(b*mpd.a)))+ field + f_p 
                       + species + species_1)
lik_power = inlabru::like(family = "gaussian", data = dat_suc_sp,
                          formula = mnd.a ~ a*(mpd.a^b)+ field + f_p 
                          + species + species_1)
lik_hyperbola = inlabru::like(family = "gaussian", data = dat_suc_sp,
                              formula = mnd.a ~ mpd.a/(a+b*mpd.a)+ field + f_p 
                              + species + species_1)
cmp = ~ a(1) + b(1) + 
  field(1, model="iid", hyper = pc_prior)+
  f_p(1, model="iid", hyper = pc_prior) +
  species(1, model="iid", hyper = pc_prior)+
  species_1(1, model="generic0", Cmatrix=vcv_tree_md_sparse, values = sp_md, hyper = pc_prior) 

cmp_e2 = ~ b(1) + 
  field(1, model="iid", hyper = pc_prior)+
  f_p(1, model="iid", hyper = pc_prior) +
  species(1, model="iid", hyper = pc_prior)+
  species_1(1, model="generic0", Cmatrix=vcv_tree_md_sparse, values = sp_md, hyper = pc_prior) 

fit_e1 = bru(components = cmp, lik_e1, 
             options = list(bru_initial = list(a = nls_coff_e1[1,1],
                                               b = nls_coff_e1[2,1]),
                            control.compute = list(dic=T, waic=T, cpo=T),
                            quantiles=c(0.025,0.5,0.975)))
summary(fit_e1)
fit_e1$summary.fixed
plot(mpd.a, max_e(fit_e1$summary.fixed[1,1],fit_e1$summary.fixed[2,1],mpd.a))

fit_e2 = bru(components = cmp, lik_e2, 
             options = list(bru_initial = list(b = nls_coff_e2[1,1]),
                            control.compute = list(dic=T, waic=T, cpo=T),
                            quantiles=c(0.025,0.5,0.975)))
summary(fit_e2)
fit_e2$summary.fixed
plot(mpd.a, max_e_2(fit_e2$summary.fixed[1,1],mpd.a))

fit_power = bru(components = cmp, lik_power, 
                options = list(bru_initial = list(a = nls_coff_power[1,1], 
                                                  b = nls_coff_power[2,1]),
                               control.compute = list(dic=T, waic=T, cpo=T),
                               quantiles=c(0.025,0.5,0.975)))
summary(fit_power)
fit_power$summary.fixed
plot(mpd.a, power(fit_power$summary.fixed[2,1], fit_power$summary.fixed[1,1], mpd.a))

fit_hyperbola = bru(components = cmp, lik_hyperbola, 
                    options = list(bru_initial = list(a = nls_coff_hyperbola[1,1], 
                                                      b = nls_coff_hyperbola[2,1]),
                                   control.compute = list(dic=T, waic=T, cpo=T),
                                   quantiles=c(0.025,0.5,0.975)))
summary(fit_hyperbola)
plot(mpd.a, hyperbola(fit_hyperbola$summary.fixed[2,1], fit_hyperbola$summary.fixed[1,1], mpd.a))

c(mod_mnd.a_mpd.a$waic$waic, fit_e1$waic$waic, fit_e2$waic$waic,
  fit_power$waic$waic, fit_hyperbola$waic$waic) ### lmer is the best after considering phy_independece



####### mlgfd.a ~ mpd.a ######
mpd.a = dat_suc_sp$mpd.a
#plot(dist, max_e(1, 0.01, dist))
#plot(dist, max_e_2(0.002740836, dist))
#plot(dist, power(1, 0.01, dist))
#plot(dist, hyperbola(1, 0.01, dist))
#plot(dist, logistic(a = 2, b = 1, c = 5,
#                    dist))

model_e1 = nlsLM(mlgfd.a ~ max_e(a, b, dist = mpd.a),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mlgfd.a ~ max_e_2(b, dist = mpd.a),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mlgfd.a ~ power(a, b, dist = mpd.a),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mlgfd.a ~ hyperbola(a, b, dist = mpd.a),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mlgfd.a ~ logistic(a, b, c,
                                        dist = mpd.a),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

plot(mpd.a, max_e(nls_coff_e1[1,1],nls_coff_e1[2,1],mpd.a))
plot(mpd.a, max_e(-5,6.8,mpd.a))
plot(mpd.a, max_e_2(nls_coff_e2[1,1],mpd.a))
#options(error=recover)

### Bayesian inla methods to inculde phylogenetic independence
pc_prior = list(prec = list("pc.prec",param=c(0.1,0.01)))
dat_suc_sp$species_1 = factor(dat_suc_sp$species)

### Calculate the original functional distance
tree = read.tree('data/original data/phylo_tree332.txt')
sp_md = unique(dat_suc_sp$species)
tree_md = keep.tip(tree, sp_md)
vcv_tree_md = ape::vcv(tree_md, model = "Brownian", corr = FALSE)
vcv_tree_md_sparse = inla.as.sparse(solve(vcv_tree_md))

### lmer ###
formula = mlgfd.a ~ 0 + mpd.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)

mod_mlgfd.a_mpd.a = inla(formula,
                       control.compute = list(dic=T, waic=T, cpo=T),
                       quantiles=c(0.025,0.5,0.975), data=dat_suc_sp)

summary(mod_mlgfd.a_mpd.a)
mod_mlgfd.a_mpd.a$summary.fixed

### nlmer ###
lik_e1 = inlabru::like(family = "gaussian", data = dat_suc_sp,
                       formula = mlgfd.a ~ a*(1-exp(-(b*mpd.a)))+ field + f_p
                       + species + species_1)
lik_e2 = inlabru::like(family = "gaussian", data = dat_suc_sp,
                       formula = mlgfd.a ~ (1-exp(-(b*mpd.a)))+ field + f_p 
                       + species + species_1)
lik_power = inlabru::like(family = "gaussian", data = dat_suc_sp,
                          formula = mlgfd.a ~ a*(mpd.a^b)+ field + f_p 
                          + species + species_1)
lik_hyperbola = inlabru::like(family = "gaussian", data = dat_suc_sp,
                              formula = mlgfd.a ~ mpd.a/(a+b*mpd.a)+ field + f_p 
                              + species + species_1)
cmp = ~ a(1) + b(1) + 
  field(1, model="iid", hyper = pc_prior)+
  f_p(1, model="iid", hyper = pc_prior) +
  species(1, model="iid", hyper = pc_prior)+
  species_1(1, model="generic0", Cmatrix=vcv_tree_md_sparse, values = sp_md, hyper = pc_prior) 

cmp_e2 = ~ b(1) + 
  field(1, model="iid", hyper = pc_prior)+
  f_p(1, model="iid", hyper = pc_prior) +
  species(1, model="iid", hyper = pc_prior)+
  species_1(1, model="generic0", Cmatrix=vcv_tree_md_sparse, values = sp_md, hyper = pc_prior) 

fit_e1 = bru(components = cmp, lik_e1, 
             options = list(bru_initial = list(a = nls_coff_e1[1,1],
                                               b = nls_coff_e1[2,1]),
                            control.compute = list(dic=T, waic=T, cpo=T),
                            quantiles=c(0.025,0.5,0.975)))
summary(fit_e1)
fit_e1$summary.fixed
plot(mpd.a, max_e(fit_e1$summary.fixed[1,1],fit_e1$summary.fixed[2,1],mpd.a))

fit_e2 = bru(components = cmp, lik_e2, 
             options = list(bru_initial = list(b = nls_coff_e2[1,1]),
                            control.compute = list(dic=T, waic=T, cpo=T),
                            quantiles=c(0.025,0.5,0.975)))
summary(fit_e2)
fit_e2$summary.fixed
plot(mpd.a, max_e_2(fit_e2$summary.fixed[1,1],mpd.a))

fit_power = bru(components = cmp, lik_power, 
                options = list(bru_initial = list(a = nls_coff_power[1,1], 
                                                  b = nls_coff_power[2,1]),
                               control.compute = list(dic=T, waic=T, cpo=T),
                               quantiles=c(0.025,0.5,0.975)))
summary(fit_power)
fit_power$summary.fixed
plot(mpd.a, power(fit_power$summary.fixed[2,1], fit_power$summary.fixed[1,1], mpd.a))

fit_hyperbola = bru(components = cmp, lik_hyperbola, 
                    options = list(bru_initial = list(a = nls_coff_hyperbola[1,1], 
                                                      b = nls_coff_hyperbola[2,1]),
                                   control.compute = list(dic=T, waic=T, cpo=T),
                                   quantiles=c(0.025,0.5,0.975)))
summary(fit_hyperbola)
plot(mpd.a, hyperbola(fit_hyperbola$summary.fixed[2,1], fit_hyperbola$summary.fixed[1,1], mpd.a))

c(fit_e1$waic$waic, fit_e2$waic$waic,
  fit_power$waic$waic, fit_hyperbola$waic$waic) ### lmer is the best after considering phy_independece


####### mnd.a ~ mconti_func_d.a ######
mconti_func_d.a = dat_suc_sp$mconti_func_d.a
plot(mconti_func_d.a, dat_suc_sp$mnd.a)
plot(dat_suc_sp$mgrowth.a, dat_suc_sp$mnd.a)
plot(dat_suc_sp$mspan.a, dat_suc_sp$mnd.a)
plot(dat_suc_sp$mpollination.a, dat_suc_sp$mnd.a)
plot(dat_suc_sp$mdispersal.a, dat_suc_sp$mnd.a)
plot(dat_suc_sp$mclonality.a, dat_suc_sp$mnd.a)
plot(dat_suc_sp$mheight.a, dat_suc_sp$mnd.a)
plot(dat_suc_sp$mldmc.a, dat_suc_sp$mnd.a)
plot(dat_suc_sp$msla.a, dat_suc_sp$mnd.a)

plot(mconti_func_d.a, dat_suc_sp$mablgfd.a)
plot(dat_suc_sp$mgrowth.a, dat_suc_sp$mablgfd.a)
plot(dat_suc_sp$mspan.a, dat_suc_sp$mablgfd.a)
plot(dat_suc_sp$mpollination.a, dat_suc_sp$mablgfd.a)
plot(dat_suc_sp$mdispersal.a, dat_suc_sp$mablgfd.a)
plot(dat_suc_sp$mclonality.a, dat_suc_sp$mablgfd.a)
plot(dat_suc_sp$mheight.a, dat_suc_sp$mablgfd.a)
plot(dat_suc_sp$mldmc.a, dat_suc_sp$mablgfd.a)
plot(dat_suc_sp$msla.a, dat_suc_sp$mablgfd.a)

#### normal lmer 
(mod_lmer = lmer(mnd.a ~ 0 + mconti_func_d.a + (1|field/f_p) + (1|species),
                 data = dat_suc_sp, REML = TRUE))
summary(mod_lmer)

### Bayesian inla methods to inculde phylogenetic independence
pc_prior = list(prec = list("pc.prec",param=c(0.1,0.01)))
dat_suc_sp$species_1 = factor(dat_suc_sp$species)
sp_name = unique(dat_suc_sp$species)

numcols = grep("^m",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])

### Calculate the original functional distance
tree = read.tree('data/original data/phylo_tree332.txt')
sp_md = unique(dat_suc_sp$species)
tree_md = keep.tip(tree, sp_md)
vcv_tree_md = ape::vcv(tree_md, model = "Brownian", corr = FALSE)
vcv_tree_md_sparse = inla.as.sparse(solve(vcv_tree_md))

newdata_for_trait = data.frame(mconti_func_d.a=seq(min(dat_suc_sp$mconti_func_d.a),
                                                   max(dat_suc_sp$mconti_func_d.a),length=100))
trait_matrix=model.matrix(~mconti_func_d.a-1,data=newdata_for_trait)
trait_matrix=as.data.frame(trait_matrix)
trait_matrix_lincombs=inla.make.lincombs(trait_matrix)

### lmer ###
formula = mnd.a ~ mconti_func_d.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mnd.a_mconti_func_d.a = inla(formula,
                                 control.compute = list(dic=T, waic=T, cpo=T),
                                 quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mnd.a_mconti_func_d.a)
mod_mnd.a_mconti_func_d.a$summary.fixed

mod_mnd.a_mconti_func_d.a_combs = inla(formula,
                                       control.compute = list(dic=T, waic=T, cpo=T),
                                       quantiles=c(0.025,0.5,0.975), data=dat_suc_sps,
                                       lincomb=phy_matrix_lincombs)

# mgrowth.a
formula = mnd.a ~ mgrowth.a +
    f(species, model="iid", hyper = pc_prior)+
    f(field, model="iid", hyper = pc_prior)+
    f(f_p, model="iid", hyper = pc_prior) + 
    f(species_1, model="generic0",
      Cmatrix=vcv_tree_md_sparse, values = sp_md, 
      hyper = pc_prior)
mod_mnd.a_mgrowth.a = inla(formula,
                             control.compute = list(dic=T, waic=T, cpo=T),
                             quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mnd.a_mgrowth.a)
mod_mnd.a_mgrowth.a$summary.fixed

# mspan.a
formula = mnd.a ~ mspan.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mnd.a_mspan.a = inla(formula,
                           control.compute = list(dic=T, waic=T, cpo=T),
                           quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mnd.a_mspan.a)
mod_mnd.a_mspan.a$summary.fixed

# mpollination.a
formula = mnd.a ~ mpollination.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mnd.a_mpollination.a = inla(formula,
                           control.compute = list(dic=T, waic=T, cpo=T),
                           quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mnd.a_mpollination.a)
mod_mnd.a_mpollination.a$summary.fixed

# mdispersal.a
formula = mnd.a ~ mdispersal.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mnd.a_mdispersal.a = inla(formula,
                           control.compute = list(dic=T, waic=T, cpo=T),
                           quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mnd.a_mdispersal.a)
mod_mnd.a_mdispersal.a$summary.fixed

# mclonality.a
formula = mnd.a ~ mclonality.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mnd.a_mclonality.a = inla(formula,
                              control.compute = list(dic=T, waic=T, cpo=T),
                              quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mnd.a_mclonality.a)
mod_mnd.a_mclonality.a$summary.fixed

# mheight.a
formula = mnd.a ~ mheight.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mnd.a_mheight.a = inla(formula,
                              control.compute = list(dic=T, waic=T, cpo=T),
                              quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mnd.a_mheight.a)
mod_mnd.a_mheight.a$summary.fixed

# mldmc.a
formula = mnd.a ~ mldmc.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mnd.a_mldmc.a = inla(formula,
                           control.compute = list(dic=T, waic=T, cpo=T),
                           quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mnd.a_mldmc.a)
mod_mnd.a_mldmc.a$summary.fixed

# msla.a
formula = mnd.a ~ msla.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mnd.a_msla.a = inla(formula,
                         control.compute = list(dic=T, waic=T, cpo=T),
                         quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mnd.a_msla.a)
mod_mnd.a_msla.a$summary.fixed

# mseedmass.a
formula = mnd.a ~ mseedmass.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mnd.a_mseedmass.a = inla(formula,
                        control.compute = list(dic=T, waic=T, cpo=T),
                        quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mnd.a_mseedmass.a)
mod_mnd.a_mseedmass.a$summary.fixed

dat_trait_mnd.a = data.frame(
  Trait = names(dat_suc_sps)[28:36],
  Value = c(mod_mnd.a_mgrowth.a$summary.fixed$mean[2],
            mod_mnd.a_mspan.a$summary.fixed$mean[2],
            mod_mnd.a_mpollination.a$summary.fixed$mean[2],
            mod_mnd.a_mdispersal.a$summary.fixed$mean[2],
            mod_mnd.a_mclonality.a$summary.fixed$mean[2],
            mod_mnd.a_mheight.a$summary.fixed$mean[2],
            mod_mnd.a_mldmc.a$summary.fixed$mean[2],
            mod_mnd.a_msla.a$summary.fixed$mean[2],
            mod_mnd.a_mseedmass.a$summary.fixed$mean[2]),
  lw = c(mod_mnd.a_mgrowth.a$summary.fixed$'0.025quant'[2],
         mod_mnd.a_mspan.a$summary.fixed$'0.025quant'[2],
         mod_mnd.a_mpollination.a$summary.fixed$'0.025quant'[2],
         mod_mnd.a_mdispersal.a$summary.fixed$'0.025quant'[2],
         mod_mnd.a_mclonality.a$summary.fixed$'0.025quant'[2],
         mod_mnd.a_mheight.a$summary.fixed$'0.025quant'[2],
         mod_mnd.a_mldmc.a$summary.fixed$'0.025quant'[2],
         mod_mnd.a_msla.a$summary.fixed$'0.025quant'[2],
         mod_mnd.a_mseedmass.a$summary.fixed$'0.025quant'[2]),
  uw = c(mod_mnd.a_mgrowth.a$summary.fixed$'0.975quant'[2],
         mod_mnd.a_mspan.a$summary.fixed$'0.975quant'[2],
         mod_mnd.a_mpollination.a$summary.fixed$'0.975quant'[2],
         mod_mnd.a_mdispersal.a$summary.fixed$'0.975quant'[2],
         mod_mnd.a_mclonality.a$summary.fixed$'0.975quant'[2],
         mod_mnd.a_mheight.a$summary.fixed$'0.975quant'[2],
         mod_mnd.a_mldmc.a$summary.fixed$'0.975quant'[2],
         mod_mnd.a_msla.a$summary.fixed$'0.975quant'[2],
         mod_mnd.a_mseedmass.a$summary.fixed$'0.975quant'[2])
)
dat_trait_mnd.a = arrange(dat_trait_mnd.a, dat_trait_mnd.a$Trait)

(p_mnd.a_traits =  ggplot(dat_trait_mnd.a, aes(x = Trait, y = Value, group = 1)) +
  geom_polygon(aes(y = uw), fill = "grey50", alpha = 0.5) +
  geom_polygon(aes(y = lw), fill = "grey99", alpha = 0.7) +
  geom_polygon(fill = NA, colour = "purple", size = 2) +
  coord_polar() +
  geom_hline(yintercept = c(-0.3, 0.3), colour = "black", size = 2) +
  geom_vline(xintercept = c(0:8), colour = "black", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 2)  +
  scale_x_discrete(labels=c("clonality","dispersal","growth","height",
                            "ldmc","pollination","seedmass","sla",
                              "span")) +
  scale_y_continuous(limits = c(-0.3, 0.3),
                    breaks = seq(-0.3, 0.3, by = 0.3),
                    labels = seq(-0.3, 0.3, by = 0.3))+
  theme_custom() +
  theme(panel.background = element_blank(),
    panel.border = element_blank(),
          panel.grid  = element_blank(),
          legend.key = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.length = unit(-1, "lines"),
          axis.ticks.margin = unit(1.3,"lines"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          axis.line=element_line(),
          axis.line.x=element_blank(),
          axis.line.y = element_blank(),
          plot.title = element_text(hjust = 0, size = 20))+
  labs(x = "", y = ""))



####### mlgfd.a ~ mconti_func_d.a ######
mconti_func_d.a = dat_suc_sp$mconti_func_d.a
plot(mconti_func_d.a, dat_suc_sp$mlgfd.a)
plot(dat_suc_sp$mgrowth.a, dat_suc_sp$mlgfd.a)
plot(dat_suc_sp$mspan.a, dat_suc_sp$mlgfd.a)
plot(dat_suc_sp$mpollination.a, dat_suc_sp$mlgfd.a)
plot(dat_suc_sp$mdispersal.a, dat_suc_sp$mlgfd.a)
plot(dat_suc_sp$mclonality.a, dat_suc_sp$mlgfd.a)
plot(dat_suc_sp$mheight.a, dat_suc_sp$mlgfd.a)
plot(dat_suc_sp$mldmc.a, dat_suc_sp$mlgfd.a)
plot(dat_suc_sp$msla.a, dat_suc_sp$mlgfd.a)

#### normal lmer 
(mod_lmer = lmer(mlgfd.a ~ 0 + mconti_func_d.a + (1|field/f_p) + (1|species),
                 data = dat_suc_sp, REML = TRUE))
summary(mod_lmer)

### Bayesian inla methods to inculde phylogenetic independence
pc_prior = list(prec = list("pc.prec",param=c(0.1,0.01)))
dat_suc_sp$species_1 = factor(dat_suc_sp$species)
sp_name = unique(dat_suc_sp$species)

numcols = grep("^m",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])

### Calculate the original functional distance
tree = read.tree('data/original data/phylo_tree332.txt')
sp_md = unique(dat_suc_sp$species)
tree_md = keep.tip(tree, sp_md)
vcv_tree_md = ape::vcv(tree_md, model = "Brownian", corr = FALSE)
vcv_tree_md_sparse = inla.as.sparse(solve(vcv_tree_md))

newdata_for_trait = data.frame(mconti_func_d.a=seq(min(dat_suc_sp$mconti_func_d.a),
                                                   max(dat_suc_sp$mconti_func_d.a),length=100))
trait_matrix=model.matrix(~mconti_func_d.a-1,data=newdata_for_trait)
trait_matrix=as.data.frame(trait_matrix)
trait_matrix_lincombs=inla.make.lincombs(trait_matrix)

### lmer ###
formula = mlgfd.a ~ mconti_func_d.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mlgfd.a_mconti_func_d.a = inla(formula,
                                 control.compute = list(dic=T, waic=T, cpo=T),
                                 quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mlgfd.a_mconti_func_d.a)
mod_mlgfd.a_mconti_func_d.a$summary.fixed

mod_mlgfd.a_mconti_func_d.a_combs = inla(formula,
                                       control.compute = list(dic=T, waic=T, cpo=T),
                                       quantiles=c(0.025,0.5,0.975), data=dat_suc_sps,
                                       lincomb=phy_matrix_lincombs)

# mgrowth.a
formula = mlgfd.a ~ mgrowth.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mlgfd.a_mgrowth.a = inla(formula,
                           control.compute = list(dic=T, waic=T, cpo=T),
                           quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mlgfd.a_mgrowth.a)
mod_mlgfd.a_mgrowth.a$summary.fixed

# mspan.a
formula = mlgfd.a ~ mspan.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mlgfd.a_mspan.a = inla(formula,
                         control.compute = list(dic=T, waic=T, cpo=T),
                         quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mlgfd.a_mspan.a)
mod_mlgfd.a_mspan.a$summary.fixed

# mpollination.a
formula = mlgfd.a ~ mpollination.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mlgfd.a_mpollination.a = inla(formula,
                                control.compute = list(dic=T, waic=T, cpo=T),
                                quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mlgfd.a_mpollination.a)
mod_mlgfd.a_mpollination.a$summary.fixed

# mdispersal.a
formula = mlgfd.a ~ mdispersal.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mlgfd.a_mdispersal.a = inla(formula,
                              control.compute = list(dic=T, waic=T, cpo=T),
                              quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mlgfd.a_mdispersal.a)
mod_mlgfd.a_mdispersal.a$summary.fixed

# mclonality.a
formula = mlgfd.a ~ mclonality.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mlgfd.a_mclonality.a = inla(formula,
                              control.compute = list(dic=T, waic=T, cpo=T),
                              quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mlgfd.a_mclonality.a)
mod_mlgfd.a_mclonality.a$summary.fixed

# mheight.a
formula = mlgfd.a ~ mheight.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mlgfd.a_mheight.a = inla(formula,
                           control.compute = list(dic=T, waic=T, cpo=T),
                           quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mlgfd.a_mheight.a)
mod_mlgfd.a_mheight.a$summary.fixed

# mldmc.a
formula = mlgfd.a ~ mldmc.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mlgfd.a_mldmc.a = inla(formula,
                         control.compute = list(dic=T, waic=T, cpo=T),
                         quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mlgfd.a_mldmc.a)
mod_mlgfd.a_mldmc.a$summary.fixed

# msla.a
formula = mlgfd.a ~ msla.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mlgfd.a_msla.a = inla(formula,
                        control.compute = list(dic=T, waic=T, cpo=T),
                        quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mlgfd.a_msla.a)
mod_mlgfd.a_msla.a$summary.fixed

# mseedmass.a
formula = mlgfd.a ~ mseedmass.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mlgfd.a_mseedmass.a = inla(formula,
                             control.compute = list(dic=T, waic=T, cpo=T),
                             quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mlgfd.a_mseedmass.a)
mod_mlgfd.a_mseedmass.a$summary.fixed

dat_trait_mlgfd.a = data.frame(
  Trait = names(dat_suc_sps)[28:36],
  Value = c(mod_mlgfd.a_mgrowth.a$summary.fixed$mean[2],
            mod_mlgfd.a_mspan.a$summary.fixed$mean[2],
            mod_mlgfd.a_mpollination.a$summary.fixed$mean[2],
            mod_mlgfd.a_mdispersal.a$summary.fixed$mean[2],
            mod_mlgfd.a_mclonality.a$summary.fixed$mean[2],
            mod_mlgfd.a_mheight.a$summary.fixed$mean[2],
            mod_mlgfd.a_mldmc.a$summary.fixed$mean[2],
            mod_mlgfd.a_msla.a$summary.fixed$mean[2],
            mod_mlgfd.a_mseedmass.a$summary.fixed$mean[2]),
  lw = c(mod_mlgfd.a_mgrowth.a$summary.fixed$'0.025quant'[2],
         mod_mlgfd.a_mspan.a$summary.fixed$'0.025quant'[2],
         mod_mlgfd.a_mpollination.a$summary.fixed$'0.025quant'[2],
         mod_mlgfd.a_mdispersal.a$summary.fixed$'0.025quant'[2],
         mod_mlgfd.a_mclonality.a$summary.fixed$'0.025quant'[2],
         mod_mlgfd.a_mheight.a$summary.fixed$'0.025quant'[2],
         mod_mlgfd.a_mldmc.a$summary.fixed$'0.025quant'[2],
         mod_mlgfd.a_msla.a$summary.fixed$'0.025quant'[2],
         mod_mlgfd.a_mseedmass.a$summary.fixed$'0.025quant'[2]),
  uw = c(mod_mlgfd.a_mgrowth.a$summary.fixed$'0.975quant'[2],
         mod_mlgfd.a_mspan.a$summary.fixed$'0.975quant'[2],
         mod_mlgfd.a_mpollination.a$summary.fixed$'0.975quant'[2],
         mod_mlgfd.a_mdispersal.a$summary.fixed$'0.975quant'[2],
         mod_mlgfd.a_mclonality.a$summary.fixed$'0.975quant'[2],
         mod_mlgfd.a_mheight.a$summary.fixed$'0.975quant'[2],
         mod_mlgfd.a_mldmc.a$summary.fixed$'0.975quant'[2],
         mod_mlgfd.a_msla.a$summary.fixed$'0.975quant'[2],
         mod_mlgfd.a_mseedmass.a$summary.fixed$'0.975quant'[2])
)
dat_trait_mlgfd.a = arrange(dat_trait_mlgfd.a, dat_trait_mlgfd.a$Trait)

(p_mlgfd.a_traits =  ggplot(dat_trait_mlgfd.a, aes(x = Trait, y = Value, group = 1)) +
  geom_polygon(aes(y = uw), fill = "grey50", alpha = 0.5) +
  geom_polygon(aes(y = lw), fill = "grey99", alpha = 0.7) +
  geom_polygon(fill = NA, colour = "purple", size = 2) +
  coord_polar() +
  geom_hline(yintercept = c(-0.3, 0.3), colour = "black", size = 2) +
  geom_vline(xintercept = c(0:8), colour = "black", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 2)  +
  scale_x_discrete(labels=c("clonality","dispersal","growth","height",
                            "ldmc","pollination","seedmass","sla",
                            "span")) +
  scale_y_continuous(limits = c(-0.3, 0.3),
                     breaks = seq(-0.3, 0.3, by = 0.3),
                     labels = seq(-0.3, 0.3, by = 0.3))+
  theme_custom() +
  theme(panel.background = element_blank(),
    panel.border = element_blank(),
        panel.grid  = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(-1, "lines"),
        axis.ticks.margin = unit(1.3,"lines"),
        axis.text.x = element_text(size = 23),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.line=element_line(),
        axis.line.x=element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(hjust = 0, size = 20))+
  labs(x = "", y = ""))


####### mablgfd.a ~ mconti_func_d.a ######
mconti_func_d.a = dat_suc_sp$mconti_func_d.a
plot(mconti_func_d.a, dat_suc_sp$mablgfd.a)
plot(dat_suc_sp$mgrowth.a, dat_suc_sp$mablgfd.a)
plot(dat_suc_sp$mspan.a, dat_suc_sp$mablgfd.a)
plot(dat_suc_sp$mpollination.a, dat_suc_sp$mablgfd.a)
plot(dat_suc_sp$mdispersal.a, dat_suc_sp$mablgfd.a)
plot(dat_suc_sp$mclonality.a, dat_suc_sp$mablgfd.a)
plot(dat_suc_sp$mheight.a, dat_suc_sp$mablgfd.a)
plot(dat_suc_sp$mldmc.a, dat_suc_sp$mablgfd.a)
plot(dat_suc_sp$msla.a, dat_suc_sp$mablgfd.a)

#### normal lmer 
(mod_lmer = lmer(mablgfd.a ~ 0 + mconti_func_d.a + (1|field/f_p) + (1|species),
                 data = dat_suc_sp, REML = TRUE))
summary(mod_lmer)

### Bayesian inla methods to inculde phylogenetic independence
pc_prior = list(prec = list("pc.prec",param=c(0.1,0.01)))
dat_suc_sp$species_1 = factor(dat_suc_sp$species)
sp_name = unique(dat_suc_sp$species)

numcols = grep("^m",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])

### Calculate the original functional distance
tree = read.tree('data/original data/phylo_tree332.txt')
sp_md = unique(dat_suc_sp$species)
tree_md = keep.tip(tree, sp_md)
vcv_tree_md = ape::vcv(tree_md, model = "Brownian", corr = FALSE)
vcv_tree_md_sparse = inla.as.sparse(solve(vcv_tree_md))

newdata_for_trait = data.frame(mconti_func_d.a=seq(min(dat_suc_sp$mconti_func_d.a),
                                                   max(dat_suc_sp$mconti_func_d.a),length=100))
trait_matrix=model.matrix(~mconti_func_d.a-1,data=newdata_for_trait)
trait_matrix=as.data.frame(trait_matrix)
trait_matrix_lincombs=inla.make.lincombs(trait_matrix)

### lmer ###
formula = mablgfd.a ~ mconti_func_d.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mablgfd.a_mconti_func_d.a = inla(formula,
                                   control.compute = list(dic=T, waic=T, cpo=T),
                                   quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mablgfd.a_mconti_func_d.a)
mod_mablgfd.a_mconti_func_d.a$summary.fixed

mod_mablgfd.a_mconti_func_d.a_combs = inla(formula,
                                         control.compute = list(dic=T, waic=T, cpo=T),
                                         quantiles=c(0.025,0.5,0.975), data=dat_suc_sps,
                                         lincomb=phy_matrix_lincombs)

# mgrowth.a
formula = mablgfd.a ~ mgrowth.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mablgfd.a_mgrowth.a = inla(formula,
                             control.compute = list(dic=T, waic=T, cpo=T),
                             quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mablgfd.a_mgrowth.a)
mod_mablgfd.a_mgrowth.a$summary.fixed

# mspan.a
formula = mablgfd.a ~ mspan.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mablgfd.a_mspan.a = inla(formula,
                           control.compute = list(dic=T, waic=T, cpo=T),
                           quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mablgfd.a_mspan.a)
mod_mablgfd.a_mspan.a$summary.fixed

# mpollination.a
formula = mablgfd.a ~ mpollination.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mablgfd.a_mpollination.a = inla(formula,
                                  control.compute = list(dic=T, waic=T, cpo=T),
                                  quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mablgfd.a_mpollination.a)
mod_mablgfd.a_mpollination.a$summary.fixed

# mdispersal.a
formula = mablgfd.a ~ mdispersal.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mablgfd.a_mdispersal.a = inla(formula,
                                control.compute = list(dic=T, waic=T, cpo=T),
                                quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mablgfd.a_mdispersal.a)
mod_mablgfd.a_mdispersal.a$summary.fixed

# mclonality.a
formula = mablgfd.a ~ mclonality.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mablgfd.a_mclonality.a = inla(formula,
                                control.compute = list(dic=T, waic=T, cpo=T),
                                quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mablgfd.a_mclonality.a)
mod_mablgfd.a_mclonality.a$summary.fixed

# mheight.a
formula = mablgfd.a ~ mheight.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mablgfd.a_mheight.a = inla(formula,
                             control.compute = list(dic=T, waic=T, cpo=T),
                             quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mablgfd.a_mheight.a)
mod_mablgfd.a_mheight.a$summary.fixed

# mldmc.a
formula = mablgfd.a ~ mldmc.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mablgfd.a_mldmc.a = inla(formula,
                           control.compute = list(dic=T, waic=T, cpo=T),
                           quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mablgfd.a_mldmc.a)
mod_mablgfd.a_mldmc.a$summary.fixed

# msla.a
formula = mablgfd.a ~ msla.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mablgfd.a_msla.a = inla(formula,
                          control.compute = list(dic=T, waic=T, cpo=T),
                          quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mablgfd.a_msla.a)
mod_mablgfd.a_msla.a$summary.fixed

# mseedmass.a
formula = mablgfd.a ~ mseedmass.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mablgfd.a_mseedmass.a = inla(formula,
                               control.compute = list(dic=T, waic=T, cpo=T),
                               quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mablgfd.a_mseedmass.a)
mod_mablgfd.a_mseedmass.a$summary.fixed

## all traits
formula = mablgfd.a ~ mgrowth.a + mspan.a + mpollination.a + mdispersal.a + 
  mclonality.a + mheight.a + mldmc.a + msla.a + mseedmass.a + 
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mablgfd.a_alltraits = inla(formula,
                                 control.compute = list(dic=T, waic=T, cpo=T),
                                 quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mablgfd.a_alltraits)
mod_mablgfd.a_alltraits$summary.fixed

dat_trait_mablgfd.a = data.frame(
  Trait = names(dat_suc_sps)[28:36],
  Value = c(mod_mablgfd.a_mgrowth.a$summary.fixed$mean[2],
            mod_mablgfd.a_mspan.a$summary.fixed$mean[2],
            mod_mablgfd.a_mpollination.a$summary.fixed$mean[2],
            mod_mablgfd.a_mdispersal.a$summary.fixed$mean[2],
            mod_mablgfd.a_mclonality.a$summary.fixed$mean[2],
            mod_mablgfd.a_mheight.a$summary.fixed$mean[2],
            mod_mablgfd.a_mldmc.a$summary.fixed$mean[2],
            mod_mablgfd.a_msla.a$summary.fixed$mean[2],
            mod_mablgfd.a_mseedmass.a$summary.fixed$mean[2]),
  lw = c(mod_mablgfd.a_mgrowth.a$summary.fixed$'0.025quant'[2],
         mod_mablgfd.a_mspan.a$summary.fixed$'0.025quant'[2],
         mod_mablgfd.a_mpollination.a$summary.fixed$'0.025quant'[2],
         mod_mablgfd.a_mdispersal.a$summary.fixed$'0.025quant'[2],
         mod_mablgfd.a_mclonality.a$summary.fixed$'0.025quant'[2],
         mod_mablgfd.a_mheight.a$summary.fixed$'0.025quant'[2],
         mod_mablgfd.a_mldmc.a$summary.fixed$'0.025quant'[2],
         mod_mablgfd.a_msla.a$summary.fixed$'0.025quant'[2],
         mod_mablgfd.a_mseedmass.a$summary.fixed$'0.025quant'[2]),
  uw = c(mod_mablgfd.a_mgrowth.a$summary.fixed$'0.975quant'[2],
         mod_mablgfd.a_mspan.a$summary.fixed$'0.975quant'[2],
         mod_mablgfd.a_mpollination.a$summary.fixed$'0.975quant'[2],
         mod_mablgfd.a_mdispersal.a$summary.fixed$'0.975quant'[2],
         mod_mablgfd.a_mclonality.a$summary.fixed$'0.975quant'[2],
         mod_mablgfd.a_mheight.a$summary.fixed$'0.975quant'[2],
         mod_mablgfd.a_mldmc.a$summary.fixed$'0.975quant'[2],
         mod_mablgfd.a_msla.a$summary.fixed$'0.975quant'[2],
         mod_mablgfd.a_mseedmass.a$summary.fixed$'0.975quant'[2])
)
dat_trait_mablgfd.a = arrange(dat_trait_mablgfd.a, dat_trait_mablgfd.a$Trait)
require(ggplot2)
(p_mablgfd.a_traits =  ggplot(dat_trait_mablgfd.a, aes(x = Trait, y = Value, group = 1)) +
    geom_polygon(aes(y = uw), fill = "grey50", alpha = 0.5) +
    geom_polygon(aes(y = lw), fill = "grey99", alpha = 0.7) +
    geom_polygon(fill = NA, colour = "purple", size = 2) +
    coord_polar() +
    geom_hline(yintercept = c(-0.3, 0.3), colour = "black", size = 2) +
    geom_vline(xintercept = c(0:8), colour = "black", size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 2)  +
    scale_x_discrete(labels=c("clonality","dispersal","growth","height",
                              "ldmc","pollination","seedmass","sla",
                              "span")) +
    scale_y_continuous(limits = c(-0.3, 0.3),
                       breaks = seq(-0.3, 0.3, by = 0.3),
                       labels = seq(-0.3, 0.3, by = 0.3))+
    theme_custom() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid  = element_blank(),
          legend.key = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.length = unit(-1, "lines"),
          axis.ticks.margin = unit(1.3,"lines"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          axis.line=element_line(),
          axis.line.x=element_blank(),
          axis.line.y = element_blank(),
          plot.title = element_text(hjust = 0, size = 20))+
    labs(x = "", y = ""))


## Merge plots mpd.a results for mnd.a and mablgfd.a ##
gap = ggplot(NULL)+theme_void()
require(ggpubr)
plo_mnd.a_mablgfd.a_traits = ggarrange(p_mnd.a_traits, p_mablgfd.a_traits,
                                      nrow = 1, ncol = 2,
                                      widths = c(1, 1), labels = c('a)', 'b)'),
                                      hjust = -3, #vjust = 0.1,
                                      font.label = list(size = 100))

ggsave(plot = plo_mnd.a_mablgfd.a_traits,
       'results/figures_sameages_raw/plo_mnd.a_mablgfd.a_traits.jpg',
       width = 120, height = 40,dpi = 300, units = 'cm',
       limitsize = F)


### nlmer ###
lik_e1 = inlabru::like(family = "gaussian", data = dat_suc_sp,
                       formula = mnd.a ~ a*(1-exp(-(b*mconti_func_d.a)))+ field + f_p
                       + species + species_1)
lik_e2 = inlabru::like(family = "gaussian", data = dat_suc_sp,
                       formula = mnd.a ~ (1-exp(-(b*mconti_func_d.a)))+ field + f_p 
                       + species + species_1)
lik_power = inlabru::like(family = "gaussian", data = dat_suc_sp,
                          formula = mnd.a ~ a*(mconti_func_d.a^b)+ field + f_p 
                          + species + species_1)
lik_hyperbola = inlabru::like(family = "gaussian", data = dat_suc_sp,
                              formula = mnd.a ~ mconti_func_d.a/(a+b*mconti_func_d.a)+ field + f_p 
                              + species + species_1)
cmp = ~ a(1) + b(1) + 
  field(1, model="iid", hyper = pc_prior)+
  f_p(1, model="iid", hyper = pc_prior) +
  species(1, model="iid", hyper = pc_prior)+
  species_1(1, model="generic0", Cmatrix=vcv_tree_md_sparse, values = sp_md, hyper = pc_prior) 

cmp_e2 = ~ b(1) + 
  field(1, model="iid", hyper = pc_prior)+
  f_p(1, model="iid", hyper = pc_prior) +
  species(1, model="iid", hyper = pc_prior)+
  species_1(1, model="generic0", Cmatrix=vcv_tree_md_sparse, values = sp_md, hyper = pc_prior) 

fit_e1 = bru(components = cmp, lik_e1, 
             options = list(bru_initial = list(a = nls_coff_e1[1,1],
                                               b = nls_coff_e1[2,1]),
                            control.compute = list(dic=T, waic=T, cpo=T),
                            quantiles=c(0.025,0.5,0.975)))
summary(fit_e1)
fit_e1$summary.fixed
plot(mconti_func_d.a, max_e(fit_e1$summary.fixed[1,1],fit_e1$summary.fixed[2,1],mconti_func_d.a))

fit_e2 = bru(components = cmp, lik_e2, 
             options = list(bru_initial = list(b = nls_coff_e2[1,1]),
                            control.compute = list(dic=T, waic=T, cpo=T),
                            quantiles=c(0.025,0.5,0.975)))
summary(fit_e2)
fit_e2$summary.fixed
plot(mconti_func_d.a, max_e_2(fit_e2$summary.fixed[1,1],mconti_func_d.a))

fit_power = bru(components = cmp, lik_power, 
                options = list(bru_initial = list(a = nls_coff_power[1,1], 
                                                  b = nls_coff_power[2,1]),
                               control.compute = list(dic=T, waic=T, cpo=T),
                               quantiles=c(0.025,0.5,0.975)))
summary(fit_power)
fit_power$summary.fixed
plot(mconti_func_d.a, power(fit_power$summary.fixed[2,1], fit_power$summary.fixed[1,1], mconti_func_d.a))

fit_hyperbola = bru(components = cmp, lik_hyperbola, 
                    options = list(bru_initial = list(a = nls_coff_hyperbola[1,1], 
                                                      b = nls_coff_hyperbola[2,1]),
                                   control.compute = list(dic=T, waic=T, cpo=T),
                                   quantiles=c(0.025,0.5,0.975)))
summary(fit_hyperbola)
plot(mconti_func_d.a, hyperbola(fit_hyperbola$summary.fixed[2,1],
                                fit_hyperbola$summary.fixed[1,1],
                                mconti_func_d.a))

c(mod_mnd.a_mconti_func_d.a$waic$waic, fit_e1$waic$waic, fit_e2$waic$waic,
  fit_power$waic$waic, fit_hyperbola$waic$waic) ### lmer is the best after considering phy_independece

mod_nd_pd_fd_fixed = mod_nd_pd_fd$summary.fixed
mod_ablgfd_pd_fd_fixed = mod_ablgfd_pd_fd$summary.fixed

require(ggplot2)
mnd.a_predicted.phy = mod_mnd.a_mconti_func_d.a_combs$summary.lincomb.derived
mnd.a_newdata.lmm_phy = data.frame(newdata_for_phy,
                                   mean=mnd.a_predicted.phy$mean,
                                   lower=mnd.a_predicted.phy$'0.025quant',
                                   upper=mnd.a_predicted.phy$'0.975quant')

plo_mnd.a_mconti_func_d.a = ggplot() + 
  geom_point(data = dat_suc_sp, aes(x = mconti_func_d.a, y = mnd.a),
             shape = 1, alpha = 1, color = 'grey', size = 16) + 
  
  geom_line(aes(x = seq(0, max(mconti_func_d.a), length.out = 100),
                y = max_e(mod_e_coff[1,1],mod_e_coff[2,1],seq(0, max(mconti_func_d.a),
                                                              length.out = 100))), 
            col = 'blue', linewidth = 4) + 
  geom_line(data=mnd.a_newdata.lmm_phy, aes(x=mconti_func_d.a,y=mean), linetype = "dashed",
            size = 12, color = 'blue')+
  geom_ribbon(data=mnd.a_newdata.lmm_phy,
              aes(x=mconti_func_d.a,ymin=lower,ymax=upper),fill="purple",alpha=0.3)+
  theme_custom() 

############ Fast start for analyzing invasion success probability ~ mnd+mfd+mpd+mfunc_d #########
library(INLA)
library(tibble)
library(lme4)
load('code/results_analyzing/analysing_sameages_raw_data/dat_suc_sp.rdata')
numcols = grep("^m",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])

## Check the co-linearity
car::vif(
  glmer(estab ~ mnd + mlgfd + mpd + mconti_func_d + (1|species) + (1|f_p),
        family = binomial, data = dat_suc_sps)
)

car::vif(
  glmer(estab ~ mnd.a + mlgfd.a + mpd.a + mconti_func_d.a + (1|species) + (1|f_p),
        family = binomial, data = dat_suc_sps)
)

car::vif(
  glmer(estab ~ mnnd + mnlgfd + mntd + mnconti_func_d + (1|species) + (1|f_p),
        family=binomial,data=dat_suc_sps)
)

### no problem
#### estab ####
pc_prior = list(prec=list("pc.prec", param=c(0.1,0.01)))
dat_suc_sps$species_1 = as.factor(dat_suc_sps$species)
estab_sp_names = unique(dat_suc_sps$species)

tree = read.tree('data/original data/phylo_tree332.txt')
estab_tree_fit = keep.tip(tree, estab_sp_names)
estab_vcv_tree = ape::vcv(estab_tree_fit, model = "Brownian", corr = FALSE)
estab_vcv_tree_sparse = inla.as.sparse(solve(estab_vcv_tree))

estab_inla.model_pd_fund1 = inla(estab~mpd+mconti_func_d+
                                   f(species, model = "iid", hyper = pc_prior)+
                                   f(field, model="iid", hyper = pc_prior) +
                                   f(f_p, model="iid", hyper = pc_prior)+
                                   f(species_1, model="generic0",
                                     Cmatrix= estab_vcv_tree_sparse,
                                     values = estab_sp_names, hyper=pc_prior),
                                 family="binomial",data=dat_suc_sps,
                                 control.compute=list(dic=T,waic=T,cpo=T),
                                 quantiles=c(0.025,0.15,0.5,0.85,0.975))
summary(estab_inla.model_pd_fund1)

estab_inla.model_all1 = inla(estab~mnd+mlgfd+mpd+mconti_func_d+
                               f(species, model = "iid", hyper = pc_prior)+
                               f(field, model="iid", hyper = pc_prior) +
                               f(f_p, model="iid", hyper = pc_prior)+
                               f(species_1, model="generic0",
                                 Cmatrix= estab_vcv_tree_sparse,
                                 values = estab_sp_names, hyper=pc_prior),
                             family="binomial",data=dat_suc_sps,
                             control.compute=list(dic=T,waic=T,cpo=T),
                             quantiles=c(0.025,0.15,0.5,0.85,0.975))
summary(estab_inla.model_all1)
c(estab_inla.model_pd_fund1$waic$waic, estab_inla.model_all1$waic$waic) ### Using the inla.model_all1

estab_inla.model_all1_md.a = inla(estab~mnd.a+mlgfd.a+mpd.a+mconti_func_d.a+
                                    f(species, model = "iid", hyper = pc_prior)+
                                    f(field, model="iid", hyper = pc_prior) +
                                    f(f_p, model="iid", hyper = pc_prior)+
                                    f(species_1, model="generic0",
                                      Cmatrix=estab_vcv_tree_sparse,
                                      values = estab_sp_names, hyper=pc_prior),
                                  family="binomial",data=dat_suc_sps,
                                  control.compute=list(dic=T,waic=T,cpo=T),
                                  quantiles=c(0.025,0.15,0.5,0.85,0.975))
summary(estab_inla.model_all1_md.a)

estab_inla.model_all1_mnd = inla(estab~mnnd+mnlgfd+mntd+mnconti_func_d+
                                   f(species, model = "iid", hyper = pc_prior)+
                                   f(field, model="iid", hyper = pc_prior) +
                                   f(f_p, model="iid", hyper = pc_prior)+
                                   f(species_1, model="generic0",
                                     Cmatrix=estab_vcv_tree_sparse,
                                     values = estab_sp_names, hyper=pc_prior),
                                 family="binomial",data=dat_suc_sps,
                                 control.compute=list(dic=T,waic=T,cpo=T),
                                 quantiles=c(0.025,0.15,0.5,0.85,0.975))
summary(estab_inla.model_all1_mnd)

c(estab_inla.model_all1$waic$waic, estab_inla.model_all1_md.a$waic$waic,
  estab_inla.model_all1_mnd$waic$waic) ### Using the inla.model_all1

# plot 
estab_data.inla.all.varied_intercept1 = estab_inla.model_all1$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower1="0.025quant",lower2="0.15quant",median="0.5quant",upper2="0.85quant",upper1="0.975quant")

estab_data.inla.all.varied_intercept = estab_data.inla.all.varied_intercept1%>%
  mutate(rowname=c("mnd","mlgfd","mpd","mconti_func_d"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

# point + effect size
estab_inla.varied.intercept.plot =
  ggplot(data=estab_data.inla.all.varied_intercept,aes(x=mean,y=rowname,color=rowname))+
  geom_point(size=80)+
  ggtitle('Establishment')+
  geom_linerange(aes(xmin=lower1,xmax=upper1),size=4.8)+
  geom_linerange(aes(xmin=lower2,xmax=upper2),size=12)+
  geom_vline(xintercept=0,linetype=2,color="grey40",linewidth=6)+
  scale_y_discrete(limits=arrange(estab_data.inla.all.varied_intercept,percent)$rowname)+
  scale_x_continuous(name="Standardized effects",limits=c(-1,1))+
  scale_color_viridis_d()+
  guides(color="none")+
  theme(panel.background = element_rect(fill = 'white',color = 'black',linewidth = 2),
        panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.18, 0.9),
        legend.text = element_text(size = 120),
        plot.margin = unit(c(0.4,0,0,0.6),units="lines"),
        plot.background = element_blank(),
        plot.title = element_text(color = '#000000',
                                  size = 220, face = 'bold'),
        text = element_text(size = 180),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size = 220,margin = margin(15,0,0,0),
                                  color = '#000000'),
        axis.ticks = element_line(linewidth = 3),
        axis.ticks.length = unit(-5,'lines'),
        axis.text.y = element_text(margin = margin(0,8,0,0),color = '#000000',
                                   size = 220),
        axis.text.x = element_text(margin = margin(8,0,0,0),color = '#000000',
                                   size = 220))


# R2%
estab_inla.varied.intercept.R2.plot = 
  ggplot(data=estab_data.inla.all.varied_intercept,aes(percent,rowname,fill=rowname))+
  geom_bar(stat="identity",width=0.8)+
  geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),
            hjust=-0.1,size=80)+
  theme_void()+
  scale_y_discrete(limits=arrange(estab_data.inla.all.varied_intercept,percent)$rowname)+
  scale_fill_viridis_d()+
  theme(plot.margin=unit(c(0.4,0,1.8,-0.3),units="lines"))+
  xlim(0,0.8)+
  guides(fill="none")

# Merge effect size + R2 
Established_all = ggarrange(estab_inla.varied.intercept.plot,
                            estab_inla.varied.intercept.R2.plot,
                            widths=c(2,1.2))

#### domin ####
dat_dom_sp = dat_suc_sps %>% filter(stage %in% c('establish', 'dominant'))


## Check the co-linearity
car::vif(
  glmer(domin ~ mnd + mlgfd + mpd + mconti_func_d + (1|species) + (1|f_p),
        family=binomial,data=dat_dom_sp)
)

car::vif(
  glmer(domin ~ mnd.a + mlgfd.a + mpd.a + mconti_func_d.a + (1|species) + (1|f_p),
        family=binomial,data=dat_dom_sp)
)

car::vif(
  glmer(domin ~ mnnd + mnlgfd + mntd + mnconti_func_d + (1|species) + (1|f_p),
        family=binomial,data=dat_dom_sp)
)
### no problem

pc_prior = list(prec=list("pc.prec", param=c(0.1,0.01)))
dat_dom_sp$species_1 = as.factor(dat_dom_sp$species)
domin_sp_names = unique(dat_dom_sp$species)

tree = read.tree('data/original data/phylo_tree332.txt')
domin_tree_fit = keep.tip(tree, domin_sp_names)
domin_vcv_tree = ape::vcv(domin_tree_fit, model = "Brownian", corr = FALSE)
domin_vcv_tree_sparse = inla.as.sparse(solve(domin_vcv_tree))

domin_inla.model_all1 = inla(domin~mnd+mlgfd+mpd+mconti_func_d+
                               f(species, model = "iid", hyper = pc_prior)+
                               f(field, model="iid", hyper = pc_prior) +
                               f(f_p, model="iid", hyper = pc_prior)+
                               f(species_1, model="generic0",
                                 Cmatrix=domin_vcv_tree_sparse,
                                 values = domin_sp_names, hyper=pc_prior),
                             family="binomial",data=dat_dom_sp,
                             control.compute=list(dic=T,waic=T,cpo=T),
                             quantiles=c(0.025,0.15,0.5,0.85,0.975))
summary(domin_inla.model_all1)
domin_inla.model_all_pd_funcd1 = inla(domin~mpd+mconti_func_d+
                                        f(species, model = "iid", hyper = pc_prior)+
                                        f(field, model="iid", hyper = pc_prior) +
                                        f(f_p, model="iid", hyper = pc_prior)+
                                        f(species_1, model="generic0",
                                          Cmatrix=domin_vcv_tree_sparse,
                                          values = domin_sp_names, hyper=pc_prior),
                                      family="binomial",data=dat_dom_sp,
                                      control.compute=list(dic=T,waic=T,cpo=T),
                                      quantiles=c(0.025,0.15,0.5,0.85,0.975))
c(domin_inla.model_all1$waic$waic, domin_inla.model_all_pd_funcd1$waic$waic)
#domin_inla.model_all1_md.a = inla(domin~mnd.a+mlgfd.a+mpd.a+mconti_func_d.a+
#                             f(species, model = "iid", hyper = pc_prior)+
#                             f(field, model="iid", hyper = pc_prior) +
#                             f(f_p, model="iid", hyper = pc_prior)+
#                           f(species_1, model="generic0",
#                              Cmatrix=domin_vcv_tree_sparse,
#                              values = domin_sp_names, hyper=pc_prior),
#                          family="binomial",data=dat_dom_sp,
#                          control.compute=list(dic=T,waic=T,cpo=T),
#                          quantiles=c(0.025,0.15,0.5,0.85,0.975))
summary(domin_inla.model_all1_md.a)

#domin_inla.model_all1_mnd = inla(domin~mnnd+mnlgfd+mntd+mnconti_func_d+
#                                   f(species, model = "iid", hyper = pc_prior)+
#                                   f(field, model="iid", hyper = pc_prior) +
#                                   f(f_p, model="iid", hyper = pc_prior)+
##                                   f(species_1, model="generic0",
#                                    Cmatrix=domin_vcv_tree_sparse,
#                                     values = domin_sp_names, hyper=pc_prior),
#                                family="binomial",data=dat_dom_sp,
#                                control.compute=list(dic=T,waic=T,cpo=T),
#                                quantiles=c(0.025,0.15,0.5,0.85,0.975))
#summary(domin_inla.model_all1_mnd)
c(domin_inla.model_all1$waic$waic, domin_inla.model_all1_md.a$waic$waic,
  domin_inla.model_all1_mnd$waic$waic) ### Using the inla.model_all1

# plot 
domin_data.inla.all.varied_intercept1 = domin_inla.model_all1$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower1="0.025quant",lower2="0.15quant",median="0.5quant",upper2="0.85quant",upper1="0.975quant")

domin_data.inla.all.varied_intercept = domin_data.inla.all.varied_intercept1%>%
  mutate(rowname=c("mnd","mlgfd","mpd","mconti_func_d"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

# point + effect size
domin_inla.varied.intercept.plot =
  ggplot(data=domin_data.inla.all.varied_intercept,aes(x=mean,y=rowname,color=rowname))+
  geom_point(size=80)+
  ggtitle('Dominance')+
  geom_linerange(aes(xmin=lower1,xmax=upper1),size=4.8)+
  geom_linerange(aes(xmin=lower2,xmax=upper2),size=12)+
  geom_vline(xintercept=0,linetype=2,color="grey40",linewidth=6)+
  scale_y_discrete(limits=arrange(domin_data.inla.all.varied_intercept,percent)$rowname)+
  scale_x_continuous(name="Standardized effects",limits=c(-1,1))+
  scale_color_viridis_d()+
  guides(color="none")+
  theme(panel.background = element_rect(fill = 'white',color = 'black',linewidth = 2),
        panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.18, 0.9),
        legend.text = element_text(size = 120),
        plot.margin = unit(c(0.4,0,0,0.6),units="lines"),
        plot.background = element_blank(),
        plot.title = element_text(face = 'bold',color = '#000000',
                                  size = 220),
        text = element_text(size = 180),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size = 220,margin = margin(15,0,0,0),
                                  color = '#000000'),
        axis.ticks = element_line(linewidth = 3),
        axis.ticks.length = unit(-5,'lines'),
        axis.text.y = element_text(margin = margin(0,8,0,0),color = '#000000',
                                   size = 220),
        axis.text.x = element_text(margin = margin(8,0,0,0),color = '#000000',
                                   size = 220))

# R2%
domin_inla.varied.intercept.R2.plot = 
  ggplot(data=domin_data.inla.all.varied_intercept,aes(percent,rowname,fill=rowname))+
  geom_bar(stat="identity",width=0.8)+
  geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),
            hjust=-0.1,size=80)+
  theme_void()+
  scale_y_discrete(limits=arrange(domin_data.inla.all.varied_intercept,percent)$rowname)+
  scale_fill_viridis_d()+
  theme(plot.margin=unit(c(0.4,0,1.8,-0.3),units="lines"),
        text = element_text(size = 150))+
  xlim(0,0.7)+
  guides(fill="none")


# Merge effect size + R2 
Dominant_all = ggarrange(domin_inla.varied.intercept.plot,
                         domin_inla.varied.intercept.R2.plot,
                         widths=c(2,1.2))

### Merge two plots
gap = ggplot(NULL)+theme_void()
library(ggpubr)
Estab_dominant_all = ggarrange(gap,Established_all,gap,
                               Dominant_all, nrow = 4, ncol = 1,
                               labels = c('','a)','','b)'),
                               heights = c(0.1, 1, 0.1, 1),
                               font.label = list(size = 220))
ggsave(plot = Estab_dominant_all,
       'results/figures_sameages_raw/Estab_dominant_all.svg',
       width = 300,height = 350, dpi = 300, units = 'cm',
       limitsize = F)

gap = ggplot(NULL)+theme_void()
library(ggpubr)
Estab_dominant_all_forppt = ggarrange(gap,Established_all,gap,
                                      Dominant_all, nrow = 1, ncol = 4,
                                      labels = c('','a)','','b)'),
                                      widths = c(0.02, 1, 0.02, 1),
                                      font.label = list(size = 220))
ggsave(plot = Estab_dominant_all_forppt,
       'results/figures_sameages_raw/Estab_dominant_all_forppt.svg',
       width = 800,height = 300, dpi = 300, units = 'cm',
       limitsize = F)


###### SEM: Success ~ mND mRFD ~ mfd mpd #####
### Original model ###
library(tibble)
library(lme4)
library(optimx)
library(piecewiseSEM)

load('code/results_analyzing/analysing_sameages_raw_data/dat_suc_sp.rdata')
numcols = grep("^m",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])
dat_suc_sps$species_1 = as.factor(dat_suc_sps$species)
estab_sp_names = unique(dat_suc_sps$species)
summary(lmer(mnd.a~mpd.a+mconti_func_d.a+
       (1|species)+
       (1|field)+
       (1|f_p), data=dat_suc_sps))
summary(lmer(mlgfd.a~mpd.a+mconti_func_d.a+
               (1|species)+
               (1|field)+
               (1|f_p), data=dat_suc_sps))

estab_sem1_md.a = psem(
  glmer(estab~mnd.a+mlgfd.a+mpd.a+mconti_func_d.a+
          (1|species)+
          (1|field)+
          (1|f_p), family=binomial, data=dat_suc_sps,
        control=glmerControl(optimizer ='optimx',
                             optCtrl=list(method='nlminb'))),
  lmer(mnd.a~mpd.a+mconti_func_d.a+
         (1|species)+
         (1|field)+
         (1|f_p), data=dat_suc_sps),
  lmer(mlgfd.a~mpd.a+mconti_func_d.a+
         (1|species)+
         (1|field)+
         (1|f_p), data=dat_suc_sps),
  #mlgfd.a %~~% mnd.a,
  #estab %~~% mpd,
  #estab %~~% mconti_func_d,
  data=dat_suc_sps)

summary(estab_sem1_mpd)
print(fisherC(estab_sem1_mpd), digits = 9)
AIC(estab_sem1_mpd)

estab_sem2_md = psem(
  glmer(estab~mnd+mlgfd+mpd+mconti_func_d+
          (1|species)+
          (1|field)+
          (1|f_p), family=binomial, data=dat_suc_sps,
        control=glmerControl(optimizer ='optimx',
                             optCtrl=list(method='nlminb'))),
  lmer(mnd~mpd+mconti_func_d+
         (1|species)+
         (1|field)+
         (1|f_p), data=dat_suc_sps),
  lmer(mlgfd~mpd+mconti_func_d+
         (1|species)+
         (1|field)+
         (1|f_p), data=dat_suc_sps),
  #mlgfd.a %~~% mnd.a,
  #estab %~~% mpd,
  #estab %~~% mconti_func_d,
  data=dat_suc_sps)

summary(estab_sem2_md)
print(fisherC(estab_sem2_md), digits = 9)
AIC(estab_sem2_md)

estab_sem2_mnd = psem(
  glmer(estab~mnnd+mnlgfd+mntd+mnconti_func_d+
          (1|species)+
          (1|field)+
          (1|f_p), family=binomial, data=dat_suc_sps,
        control=glmerControl(optimizer ='optimx',
                             optCtrl=list(method='nlminb'))),
  lmer(mnnd~mntd+mnconti_func_d+
         (1|species)+
         (1|field)+
         (1|f_p), data=dat_suc_sps),
  lmer(mnlgfd~mntd+mnconti_func_d+
         (1|species)+
         (1|field)+
         (1|f_p), data=dat_suc_sps),
  #mlgfd.a %~~% mnd.a,
  #estab %~~% mpd,
  #estab %~~% mconti_func_d,
  data=dat_suc_sps)

summary(estab_sem2_mnd)
print(fisherC(estab_sem2_mnd), digits = 9)
AIC(estab_sem2_mnd)



####### MND MFD change with time #######
load("D:/R projects/BSS/code/data preparation/transformed data/fit_fp_same_ages_0.25_raw.RData")
inter_all_c_l_2 = split(inter_all_c, inter_all_c$f_p)
plots = names(fit_fp_same_ages_0.25_raw)
dat_suc_sp_timeall = data.frame()

for (i in 1:480) {
  #i = 1
  plot = fit_fp_same_ages_0.25_raw[[which(names(fit_fp_same_ages_0.25_raw) == plots[i])]]
  plot_1 = plot$re_cover_ab_f
  inter_all_plot = inter_all_c_l_2[[which(names(inter_all_c_l_2) == plots[i])]]
  
  ages = sort(unique(plot_1$fake_age))
  dat_suc_sp_plo = data.frame()
  for (j in 1:length(ages)) {
    #j = 1
    age = ages[j]
    plot_age = plot_1 %>% filter(fake_age == age)
    plot_age_mat = plot_age[,10:ncol(plot_age)]
    sps = names(plot_age_mat[which(plot_age_mat != 1e-06)])
    inter_all_plot_age = inter_all_plot %>% filter(species_i %in% sps & species_j %in% sps)
    if ('native' %in% unique(inter_all_plot_age$stage_i) &
        length(unique(inter_all_plot_age$stage_i)) > 1){
    inv_suc = inter_all_plot_age %>% filter(stage_i != 'native' & stage_j == 'native') 
    inv_sp = unique(inv_suc$species_i)
    dat_suc_sp_t = data.frame()
    
    for (j in 1:length(inv_sp)) {
      inv_suc_sp = inv_suc %>% filter(species_i == inv_sp[j])
      mnd = mean(inv_suc_sp$nd)
      mlgfd = mean(inv_suc_sp$lgfd)
      mablgfd = mean(inv_suc_sp$ablgfd)
      mpd = mean(inv_suc_sp$Phylo_dis)
      mfunc_d = mean(inv_suc_sp$Multi_traits)
      mconti_func_d = mean(inv_suc_sp$Multi_conti_traits)
      mgrowth = mean(inv_suc_sp$growth)
      mspan = mean(inv_suc_sp$span)
      mpollination = mean(inv_suc_sp$pollination)
      mdispersal = mean(inv_suc_sp$dispersal)
      mclonality = mean(inv_suc_sp$clonality)
      mheight = mean(inv_suc_sp$height)
      mldmc = mean(inv_suc_sp$ldmc)
      msla = mean(inv_suc_sp$sla)
      mseedmass = mean(inv_suc_sp$seedmass)
      mnd.a = sum(inv_suc_sp$nd*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
      mlgfd.a = sum(inv_suc_sp$lgfd*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
      mablgfd.a = sum(inv_suc_sp$ablgfd*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
      mpd.a = sum(inv_suc_sp$Phylo_dis*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
      mfunc_d.a = sum(inv_suc_sp$Multi_traits*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
      mconti_func_d.a = sum(inv_suc_sp$Multi_conti_traits*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
      mgrowth.a = sum(inv_suc_sp$growth*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
      mspan.a = sum(inv_suc_sp$span*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
      mpollination.a = sum(inv_suc_sp$pollination*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
      mdispersal.a = sum(inv_suc_sp$dispersal*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
      mclonality.a = sum(inv_suc_sp$clonality*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
      mheight.a = sum(inv_suc_sp$height*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
      mldmc.a = sum(inv_suc_sp$ldmc*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
      msla.a = sum(inv_suc_sp$sla*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
      mseedmass.a = sum(inv_suc_sp$seedmass*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
      mnnd = min(inv_suc_sp$nd)
      mnlgfd = min(inv_suc_sp$lgfd)
      mnablgfd = min(inv_suc_sp$ablgfd)
      mntd = min(inv_suc_sp$Phylo_dis)
      mnfunc_d = min(inv_suc_sp$Multi_traits)
      mnconti_func_d = min(inv_suc_sp$Multi_conti_traits)
      mngrowth = min(inv_suc_sp$growth)
      mnspan = min(inv_suc_sp$span)
      mnpollination = min(inv_suc_sp$pollination)
      mndispersal = min(inv_suc_sp$dispersal)
      mnclonality = min(inv_suc_sp$clonality)
      mnheight = min(inv_suc_sp$height)
      mnldmc = min(inv_suc_sp$ldmc)
      mnsla = min(inv_suc_sp$sla)
      mnseedmass = min(inv_suc_sp$seedmass)
      
      dat_suc_sp_1 = data.frame(f_p = unique(inter_all_plot_age$f_p),
                                plot = unique(inter_all_plot_age$f_p),
                                field = unique(inter_all_plot_age$field),
                                year = unique(plot_age$Year),
                                age = unique(plot_age$Age),
                                fake_age = unique(plot_age$fake_age),
                                species = inv_sp[j], 
                                stage = unique(inv_suc_sp$stage_i),
                                estab = unique(inv_suc_sp$estab_i),
                                domin = unique(inv_suc_sp$domin_i),
                                mnd = mnd, mlgfd = mlgfd, mpd = mpd, mfunc_d = mfunc_d,
                                mconti_func_d = mconti_func_d, mgrowth = mgrowth, mspan = mspan, 
                                mpollination = mpollination, mdispersal = mdispersal, 
                                mclonality = mclonality, mheight = mheight, mldmc = mldmc,
                                msla = msla, mseedmass = mseedmass, 
                                mnd.a = mnd.a, mlgfd.a = mlgfd.a, mablgfd.a, mpd.a = mpd.a,
                                mfunc_d.a = mfunc_d.a, mconti_func_d.a = mconti_func_d.a,
                                mgrowth.a = mgrowth.a, mspan.a = mspan.a,mpollination.a = mpollination.a,
                                mdispersal.a = mdispersal.a, mclonality.a = mclonality.a,
                                mheight.a = mheight.a, mldmc.a = mldmc.a, msla.a = msla.a,
                                mseedmass.a = mseedmass.a, 
                                mnnd = mnnd, mnlgfd = mnlgfd, mntd = mntd,
                                mnfunc_d = mnfunc_d, mnconti_func_d = mnconti_func_d, 
                                mngrowth = mngrowth, mnspan = mnspan, 
                                mnpollination = mnpollination, mndispersal = mndispersal, 
                                mnclonality = mnclonality, mnheight = mnheight, mnldmc = mnldmc,
                                mnsla = mnsla, mnseedmass = mnseedmass)
      dat_suc_sp_t = rbind(dat_suc_sp_t, dat_suc_sp_1)
    }
      } else {
      print('cannot calculate exotic-native differences')
     }
    
    dat_suc_sp_plo = rbind(dat_suc_sp_plo, dat_suc_sp_t)
    
  } 
  dat_suc_sp_timeall = rbind(dat_suc_sp_timeall, dat_suc_sp_plo)
}

dat_suc_sp_timeall$stage_level = NA
dat_suc_sp_timeall[dat_suc_sp_timeall$stage == 'introduce',]$stage_level = 1
dat_suc_sp_timeall[dat_suc_sp_timeall$stage == 'establish',]$stage_level = 2
dat_suc_sp_timeall[dat_suc_sp_timeall$stage == 'dominant',]$stage_level = 3
dat_suc_sp_timeall$stage_level = ordered(dat_suc_sp_timeall$stage_level)
str(dat_suc_sp_timeall)
save(dat_suc_sp_timeall,
     file = 'code/results_analyzing/analysing_sameages_raw_data/dat_suc_sp_timeall.rdata')

require(INLA)
numcols = grep("^m",names(dat_suc_sp_timeall))
dat_suc_sp_timealls = dat_suc_sp_timeall
dat_suc_sp_timealls[,c(6, numcols)] = scale(dat_suc_sp_timealls[,c(6, numcols)])
dat_suc_sp_timealls$species_1 = as.factor(dat_suc_sp_timealls$species)
estab_sp_names = unique(dat_suc_sp_timealls$species)

tree = read.tree('data/original data/phylo_tree332.txt')
estab_tree_fit = keep.tip(tree, estab_sp_names)
estab_vcv_tree = ape::vcv(estab_tree_fit, model = "Brownian", corr = FALSE)
estab_vcv_tree_sparse = inla.as.sparse(solve(estab_vcv_tree))
dat_suc_sp_timealls$species_1 = dat_suc_sp_timealls$species

dat_suc_sp_timeall_1 = dat_suc_sp_timeall %>% filter(f_p == '1_4')
plot(dat_suc_sp_timeall_1$fake_age, dat_suc_sp_timeall_1$mnd)
plot(dat_suc_sp_timeall_1$year, dat_suc_sp_timeall_1$mlgfd)

mod_mnd_age = lmer(mnd~fake_age+(1|field)+(1|f_p), data = dat_suc_sp_timealls)
summary(mod_mnd_age)

mod_mlgfd_age = lmer(mlgfd~fake_age+(1|field)+(1|f_p), data = dat_suc_sp_timealls)
summary(mod_mlgfd_age)

mod_mnd.a_age = lmer(mnd.a~fake_age+(1|field)+(1|f_p), data = dat_suc_sp_timealls)
summary(mod_mnd.a_age)

mod_mlgfd.a_age = lmer(mlgfd.a~fake_age+(1|field)+(1|f_p), data = dat_suc_sp_timealls)
summary(mod_mlgfd.a_age)

## Check the co-linearity
car::vif(
  glmer(estab ~ mnd.a + mlgfd.a + (1|species) + (1|f_p),
        family=binomial,data=dat_suc_sp_timealls)
)
## VIF < 3
with(dat_suc_sp_timealls,cor(mnd.a,mlgfd.a))

## VIF < 3
with(dat_suc_sp_timealls,cor(mnd.a,mlgfd.a))

# cor_coefficients = -0.56<0.7
### Analyse for mnd mfd abundance weighted mean #############
#### Establishment
pc_prior = list(prec=list("pc.prec", param=c(0.1,0.01)))
dat_suc_sp_timealls$fake_age_1 = as.factor(dat_suc_sp_timealls$fake_age)
mod_suc_sp_ab_md.a_age = inla(estab ~ mnd.a*fake_age + mlgfd.a*fake_age +
                            f(species, model = "iid", hyper = pc_prior) +
                            f(field, model="iid", hyper = pc_prior) +
                            f(fake_age_1,model="ar1", hyper = pc_prior)+
                            f(f_p, model="iid", hyper = pc_prior)+
                            f(species_1, model="generic0",
                              Cmatrix= estab_vcv_tree_sparse,
                              values = estab_sp_names, hyper=pc_prior),
                          control.compute = list(dic=T,waic=T,cpo=T,
                                                 config = TRUE),
                          quantiles=c(0.025,0.5,0.975),
                          family="binomial", data=dat_suc_sp_timealls)
summary(mod_suc_sp_ab_md.a_age)

mod_suc_sp_ab_md_age = inla(estab ~ mnd*fake_age+mlgfd*fake_age+
                          f(species, model = "iid", hyper = pc_prior)+
                          f(field, model="iid", hyper = pc_prior) + 
                          f(f_p, model="iid", hyper = pc_prior) + 
                            f(fake_age_1,model="ar1", hyper = pc_prior)+
                          f(species_1, model="generic0",
                            Cmatrix= estab_vcv_tree_sparse,
                            values = estab_sp_names, hyper=pc_prior),
                        control.compute = list(dic=T,waic=T,cpo=T,
                                               config = TRUE),
                        quantiles=c(0.025,0.5,0.975),
                        family="binomial", data=dat_suc_sp_timealls)
summary(mod_suc_sp_ab_md_age)

mod_suc_sp_ab_mnd_age = inla(estab ~ mnnd*fake_age+mnlgfd*fake_age+
                              f(species, model = "iid", hyper = pc_prior)+
                              f(field, model="iid", hyper = pc_prior) + 
                              f(f_p, model="iid", hyper = pc_prior) +
                               f(fake_age_1,model="ar1", hyper = pc_prior)+
                              f(species_1, model="generic0",
                                Cmatrix= estab_vcv_tree_sparse,
                                values = estab_sp_names, hyper=pc_prior),
                            control.compute = list(dic=T,waic=T,cpo=T,
                                                   config = TRUE),
                            quantiles=c(0.025,0.5,0.975),
                            family="binomial", data=dat_suc_sp_timealls)
summary(mod_suc_sp_ab_mnd_age)

c(mod_suc_sp_ab_md_age$waic$waic,
  mod_suc_sp_ab_md.a_age$waic$waic,
  mod_suc_sp_ab_mnd_age$waic$waic)

######### Plot interaction effects #########
tmp = inla.posterior.sample(100, mod_suc_sp_ab_md.a_age)
tmp_r = inla.posterior.sample.eval(c('fake_age_1'), tmp)
tmp_1 = inla.posterior.sample.eval(c('mnd.a', 'mlgfd.a',
                                     'fake_age', "mnd.a:fake_age",
                                     'fake_age:mlgfd.a'), tmp)
tmp_all = inla.posterior.sample.eval(c('mnd.a', 'mlgfd.a', 
                                       'fake_age', "mnd.a:fake_age",
                                       'fake_age:mlgfd.a'), tmp)
row.names(tmp_all) = c('mnd.a', 'mlgfd.a', 
                       'fake_age', "mnd.a:fake_age",
                       'fake_age:mlgfd.a')
row.names(tmp_1) = c('mnd.a', 'mlgfd.a',
                     'fake_age', "mnd.a:fake_age",
                     'fake_age:mlgfd.a')

# sequences of hypothetical trait values to predict off of
plot.fake_age = seq(1,51,2)
seq.fake_age = sort(unique(dat_suc_sp_timealls$fake_age))

# Create a design matrix based on the equation: 1 = intercept, my.seq = trail values to predict off of
dm.fake_age = cbind(1, seq.fake_age)  

# Grab the appropriate slope parameters
mcmc.traits_mnd.a = cbind(tmp_1[1,], tmp_1[4,])
mcmc.traits_mlgfd.a = cbind(tmp_1[2,], tmp_1[5,])

# Some quick matrix math, multiplying the design matrix by the mcmc slope parameters
pred.fake_age_mnd.a = mcmc.traits_mnd.a %*% t(dm.fake_age)   
pred.fake_age.random_mnd.a = mcmc.traits_mnd.a %*% t(dm.fake_age) + t(tmp_r)

pred.fake_age_mlgfd.a = mcmc.traits_mlgfd.a %*% t(dm.fake_age)   
pred.fake_age.random_mlgfd.a = mcmc.traits_mlgfd.a %*% t(dm.fake_age) + t(tmp_r)

pm1_mnd.a = apply(pred.fake_age_mnd.a, c(2), function(x) mean(x, na.rm=TRUE))   
cri1_mnd.a = apply(pred.fake_age_mnd.a, c(2), function(x) quantile(x, prob = c(0.025, 0.975, 0.05, 0.95)))

data_random_mnd.a_suc = data.frame(
  "fake_age" = plot.fake_age,
  "random_m" = pm1_mnd.a,
  "random_lower95" = c(t(cri1_mnd.a[3,])),
  "random_upper95" = c(t(cri1_mnd.a[4,]))
)

pm1_mlgfd.a = apply(pred.fake_age.random_mlgfd.a, c(2), function(x) mean(x, na.rm=TRUE))   
cri1_mlgfd.a = apply(pred.fake_age.random_mlgfd.a, c(2), function(x) quantile(x, prob = c(0.025, 0.975, 0.05, 0.95)))

data_random_mlgfd.a_suc = data.frame(
  "fake_age" = plot.fake_age,
  "random_m" = pm1_mlgfd.a,
  "random_lower95" = c(t(cri1_mlgfd.a[3,])),
  "random_upper95" = c(t(cri1_mlgfd.a[4,]))
)

# Summarize the predictions. We really only need the median and CRIs
(dat.pred.suc_mnd.a = data.frame(
  "fake_age" = plot.fake_age,
  t(apply(pred.fake_age_mnd.a, 2, quantile, probs = c(0.025, 0.5, 0.975))) ,
  mean = apply(pred.fake_age_mnd.a, 2, mean) 
))
colnames(dat.pred.suc_mnd.a)[2:5] = c("lower95", "med", "upper95", "mean")

(dat.pred.suc_mlgfd.a = data.frame(
  "fake_age" = plot.fake_age,
  t(apply(pred.fake_age_mlgfd.a, 2, quantile, probs = c(0.025, 0.5, 0.975))),
  mean = apply(pred.fake_age_mlgfd.a, 2, mean) 
))
colnames(dat.pred.suc_mlgfd.a)[2:5] = c("lower95", "med", "upper95", "mean")

plot.mnd.a_effect_suc = ggplot() +
  theme_classic()+
  geom_ribbon(data = dat.pred.suc_mnd.a, aes(x = fake_age, y = med, ymin = lower95, ymax = upper95), 
              fill = "grey30", alpha = 0.2, size = 2)+
  geom_smooth(data = dat.pred.suc_mnd.a, aes(x = fake_age, y = med), 
              se = F, color = "grey30", linewidth = 4) +
  geom_pointrange(data = data_random_mnd.a_suc, aes(x = fake_age, y = random_m,
                                                    ymin = random_lower95, ymax = random_upper95),
                  color = "grey30", size = 2) +
  labs(x = "Age", y = "mnd.a-establishment relationship") +
  scale_x_continuous(breaks = plot.fake_age) +
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
plot.mnd.a_effect_suc


plot.mlgfd.a_effect_suc = ggplot() +
  theme_classic()+
  geom_ribbon(data = dat.pred.suc_mlgfd.a, aes(x = fake_age, y = med, ymin = lower95, ymax = upper95), 
              fill = "grey30", alpha = 0.2, size = 2)+
  geom_smooth(data = dat.pred.suc_mlgfd.a, aes(x = fake_age, y = med), 
              se = F, color = "grey30", linewidth = 4) +
  geom_pointrange(data = data_random_mlgfd.a_suc, aes(x = fake_age, y = random_m,
                                                      ymin = random_lower95, ymax = random_upper95),
                  color = "grey30", size = 2) +
  labs(x = "Age", y = "mlgfd.a-establishment relationship") +
  scale_x_continuous(breaks = plot.fake_age) +
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

plot.mlgfd.a_effect_suc

gap = ggplot(NULL)+theme_void()
library(ggpubr)
success.mnd.mlgfd.effects = ggarrange(gap, plot.mnd.a_effect_suc,
                                      gap, plot.mlgfd.a_effect_suc,
                                     # gap, plot_R2_per_estab,
                                      #gap, plot.mnd_effect_domin,
                                     # gap, plot.mlgfd_effect_domin,
                                     # gap, plot_R2_per_domin,
                                      nrow = 1, ncol = 4,
                                      hjust = 1,
                                      vjust = 1.5,
                                      labels = list('','a)', '', 'b)'
                                                    #'', 'c)', '', 'd)',
                                                    #'', 'e)', '', 'f)'
                                                    ),
                                      widths = c(0.1, 1, 0.1, 1),
                                      font.label = list(size = 50,
                                                        face = "bold",
                                                        color ="black"))

ggsave('results/figures_sameages_raw/success.mnd.mlgfd.effects_model.svg',
       width = 36, height = 10, units = c('in'),
       dpi = 300, limitsize = F, plot = success.mnd.mlgfd.effects)
#################### Analyse invasion impact ####################
inter_all_c_l_2 = split(inter_all_c, inter_all_c$f_p)

##### Species level #####
dat_imp_sp_sum = data.frame()

for (i in 1:length(inter_all_c_l_2)) {
  #i = 1
  trans_plot = inter_all_c_l_2[[i]]
  imp_intro = trans_plot %>% 
    filter(stage_i == 'introduce' & stage_j == 'native')
  imp_estab = trans_plot %>% filter(stage_i %in% c('establish') &
                                      stage_j == 'native')
  imp_domin = trans_plot %>% 
    filter(stage_i == 'dominant' & stage_j == 'native')
  imp_all = trans_plot %>% 
    filter(stage_i != 'native' & stage_j == 'native')
  
  nati_sp = unique(imp_all$species_j)
  for (j in 1:length(nati_sp)) {
    #j = 1
    imp_intro_sp = imp_intro %>% filter(species_j == nati_sp[j])
    imp_estab_sp = imp_estab %>% filter(species_j == nati_sp[j])
    imp_domin_sp = imp_domin %>% filter(species_j == nati_sp[j])
    imp_all_sp = imp_all %>% filter(species_j == nati_sp[j])
    if (nrow(imp_intro_sp) > 0){
      mnd = mean(imp_intro_sp$nd)
      mlgfd = mean(imp_intro_sp$lgfd)
      mpd = mean(imp_intro_sp$Phylo_dis)
      mfunc_d = mean(imp_intro_sp$Multi_traits)
      mconti_func_d = mean(imp_intro_sp$Multi_conti_traits)
      mnd.a = sum(imp_intro_sp$nd*imp_intro_sp$ra_m_fit_t_i*imp_intro_sp$ra_m_fit_t_j)/sum(imp_intro_sp$ra_m_fit_t_i*imp_intro_sp$ra_m_fit_t_j)
      mlgfd.a = sum(imp_intro_sp$lgfd*imp_intro_sp$ra_m_fit_t_i*imp_intro_sp$ra_m_fit_t_j)/sum(imp_intro_sp$ra_m_fit_t_i*imp_intro_sp$ra_m_fit_t_j)
      mpd.a = sum(imp_intro_sp$Phylo_dis*imp_intro_sp$ra_m_fit_t_i*imp_intro_sp$ra_m_fit_t_j)/sum(imp_intro_sp$ra_m_fit_t_i*imp_intro_sp$ra_m_fit_t_j)
      mfunc_d.a = sum(imp_intro_sp$Multi_traits*imp_intro_sp$ra_m_fit_t_i*imp_intro_sp$ra_m_fit_t_j)/sum(imp_intro_sp$ra_m_fit_t_i*imp_intro_sp$ra_m_fit_t_j)
      mconti_func_d.a = sum(imp_intro_sp$Multi_conti_traits*imp_intro_sp$ra_m_fit_t_i*imp_intro_sp$ra_m_fit_t_j)/sum(imp_intro_sp$ra_m_fit_t_i*imp_intro_sp$ra_m_fit_t_j)
      mnnd = min(imp_intro_sp$nd)
      mnlgfd = min(imp_intro_sp$lgfd)
      mntd = min(imp_intro_sp$Phylo_dis)
      mnfunc_d = min(imp_intro_sp$Multi_traits)
      mnconti_func_d = min(imp_intro_sp$Multi_conti_traits)
      dat_imp_sp_intro = data.frame(f_p = unique(trans_plot$f_p),
                                    plot = unique(trans_plot$f_p),
                                    field = unique(trans_plot$field),
                                    species = nati_sp[j],
                                    gain_loss = unique(imp_all_sp$gain_loss.fittedtime_j),
                                    mean.relative.change = unique(imp_all_sp$mean.relative.change.fittedtime_j),
                                    mean.abs.change = unique(imp_all_sp$mean.absolute.change.fittedtime_j),
                                    stage = unique(imp_intro_sp$stage_i),
                                    estab = unique(imp_intro_sp$estab_i),
                                    domin = unique(imp_intro_sp$domin_i),
                                    mnd = mnd, mlgfd = mlgfd, mpd = mpd, mfunc_d = mfunc_d,
                                    mconti_func_d = mconti_func_d,
                                    mnd.a = mnd.a, mlgfd.a = mlgfd.a, mpd.a = mpd.a,
                                    mfunc_d.a = mfunc_d.a, mconti_func_d.a = mconti_func_d.a,
                                    mnnd = mnnd, mnlgfd = mnlgfd, mntd = mntd,
                                    mnfunc_d = mnfunc_d, mnconti_func_d = mnconti_func_d)
    } else {dat_imp_sp_intro = NULL}
    if (nrow(imp_estab_sp) > 0){
      mnd = mean(imp_estab_sp$nd)
      mlgfd = mean(imp_estab_sp$lgfd)
      mpd = mean(imp_estab_sp$Phylo_dis)
      mfunc_d = mean(imp_estab_sp$Multi_traits)
      mconti_func_d = mean(imp_estab_sp$Multi_conti_traits)
      mnd.a = sum(imp_estab_sp$nd*imp_estab_sp$ra_m_fit_t_i*imp_estab_sp$ra_m_fit_t_j)/sum(imp_estab_sp$ra_m_fit_t_i*imp_estab_sp$ra_m_fit_t_j)
      mlgfd.a = sum(imp_estab_sp$lgfd*imp_estab_sp$ra_m_fit_t_i*imp_estab_sp$ra_m_fit_t_j)/sum(imp_estab_sp$ra_m_fit_t_i*imp_estab_sp$ra_m_fit_t_j)
      mpd.a = sum(imp_estab_sp$Phylo_dis*imp_estab_sp$ra_m_fit_t_i*imp_estab_sp$ra_m_fit_t_j)/sum(imp_estab_sp$ra_m_fit_t_i*imp_estab_sp$ra_m_fit_t_j)
      mfunc_d.a = sum(imp_estab_sp$Multi_traits*imp_estab_sp$ra_m_fit_t_i*imp_estab_sp$ra_m_fit_t_j)/sum(imp_estab_sp$ra_m_fit_t_i*imp_estab_sp$ra_m_fit_t_j)
      mconti_func_d.a = sum(imp_estab_sp$Multi_conti_traits*imp_estab_sp$ra_m_fit_t_i*imp_estab_sp$ra_m_fit_t_j)/sum(imp_estab_sp$ra_m_fit_t_i*imp_estab_sp$ra_m_fit_t_j)
      mnnd = min(imp_estab_sp$nd)
      mnlgfd = min(imp_estab_sp$lgfd)
      mntd = min(imp_estab_sp$Phylo_dis)
      mnfunc_d = min(imp_estab_sp$Multi_traits)
      mnconti_func_d = min(imp_estab_sp$Multi_conti_traits)
      dat_imp_sp_estab = data.frame(f_p = unique(trans_plot$f_p),
                                    plot = unique(trans_plot$f_p),
                                    field = unique(trans_plot$field),
                                    species = nati_sp[j],
                                    gain_loss = unique(imp_all_sp$gain_loss.fittedtime_j),
                                    mean.relative.change = unique(imp_all_sp$mean.relative.change.fittedtime_j),
                                    mean.abs.change = unique(imp_all_sp$mean.absolute.change.fittedtime_j),
                                    stage = unique(imp_estab_sp$stage_i),
                                    estab = unique(imp_estab_sp$estab_i),
                                    domin = unique(imp_estab_sp$domin_i),
                                    mnd = mnd, mlgfd = mlgfd, mpd = mpd, mfunc_d = mfunc_d,
                                    mconti_func_d = mconti_func_d,
                                    mnd.a = mnd.a, mlgfd.a = mlgfd.a, mpd.a = mpd.a,
                                    mfunc_d.a = mfunc_d.a, mconti_func_d.a = mconti_func_d.a,
                                    mnnd = mnnd, mnlgfd = mnlgfd, mntd = mntd,
                                    mnfunc_d = mnfunc_d, mnconti_func_d = mnconti_func_d)
    } else {dat_imp_sp_estab = NULL}
    if (nrow(imp_domin_sp) > 0){
      mnd = mean(imp_domin_sp$nd)
      mlgfd = mean(imp_domin_sp$lgfd)
      mpd = mean(imp_domin_sp$Phylo_dis)
      mfunc_d = mean(imp_domin_sp$Multi_traits)
      mconti_func_d = mean(imp_domin_sp$Multi_conti_traits)
      mnd.a = sum(imp_domin_sp$nd*imp_domin_sp$ra_m_fit_t_i*imp_domin_sp$ra_m_fit_t_j)/sum(imp_domin_sp$ra_m_fit_t_i*imp_domin_sp$ra_m_fit_t_j)
      mlgfd.a = sum(imp_domin_sp$lgfd*imp_domin_sp$ra_m_fit_t_i*imp_domin_sp$ra_m_fit_t_j)/sum(imp_domin_sp$ra_m_fit_t_i*imp_domin_sp$ra_m_fit_t_j)
      mpd.a = sum(imp_domin_sp$Phylo_dis*imp_domin_sp$ra_m_fit_t_i*imp_domin_sp$ra_m_fit_t_j)/sum(imp_domin_sp$ra_m_fit_t_i*imp_domin_sp$ra_m_fit_t_j)
      mfunc_d.a = sum(imp_domin_sp$Multi_traits*imp_domin_sp$ra_m_fit_t_i*imp_domin_sp$ra_m_fit_t_j)/sum(imp_domin_sp$ra_m_fit_t_i*imp_domin_sp$ra_m_fit_t_j)
      mconti_func_d.a = sum(imp_domin_sp$Multi_conti_traits*imp_domin_sp$ra_m_fit_t_i*imp_domin_sp$ra_m_fit_t_j)/sum(imp_domin_sp$ra_m_fit_t_i*imp_domin_sp$ra_m_fit_t_j)
      mnnd = min(imp_domin_sp$nd)
      mnlgfd = min(imp_domin_sp$lgfd)
      mntd = min(imp_domin_sp$Phylo_dis)
      mnfunc_d = min(imp_domin_sp$Multi_traits)
      mnconti_func_d = min(imp_domin_sp$Multi_conti_traits)
      dat_imp_sp_domin = data.frame(f_p = unique(trans_plot$f_p),
                                    plot = unique(trans_plot$f_p),
                                    field = unique(trans_plot$field),
                                    species = nati_sp[j],
                                    gain_loss = unique(imp_all_sp$gain_loss.fittedtime_j),
                                    mean.relative.change = unique(imp_all_sp$mean.relative.change.fittedtime_j),
                                    mean.abs.change = unique(imp_all_sp$mean.absolute.change.fittedtime_j),
                                    stage = unique(imp_domin_sp$stage_i),
                                    estab = unique(imp_domin_sp$estab_i),
                                    domin = unique(imp_domin_sp$domin_i),
                                    mnd = mnd, mlgfd = mlgfd, mpd = mpd, mfunc_d = mfunc_d,
                                    mconti_func_d = mconti_func_d,
                                    mnd.a = mnd.a, mlgfd.a = mlgfd.a, mpd.a = mpd.a,
                                    mfunc_d.a = mfunc_d.a, mconti_func_d.a = mconti_func_d.a,
                                    mnnd = mnnd, mnlgfd = mnlgfd, mntd = mntd,
                                    mnfunc_d = mnfunc_d, mnconti_func_d = mnconti_func_d)
    } else {dat_imp_sp_domin = NULL}
    if (nrow(imp_all_sp) > 0){
      mnd = mean(imp_all_sp$nd)
      mlgfd = mean(imp_all_sp$lgfd)
      mpd = mean(imp_all_sp$Phylo_dis)
      mfunc_d = mean(imp_all_sp$Multi_traits)
      mconti_func_d = mean(imp_all_sp$Multi_conti_traits)
      mnd.a = sum(imp_all_sp$nd*imp_all_sp$ra_m_fit_t_i*imp_all_sp$ra_m_fit_t_j)/sum(imp_all_sp$ra_m_fit_t_i*imp_all_sp$ra_m_fit_t_j)
      mlgfd.a = sum(imp_all_sp$lgfd*imp_all_sp$ra_m_fit_t_i*imp_all_sp$ra_m_fit_t_j)/sum(imp_all_sp$ra_m_fit_t_i*imp_all_sp$ra_m_fit_t_j)
      mpd.a = sum(imp_all_sp$Phylo_dis*imp_all_sp$ra_m_fit_t_i*imp_all_sp$ra_m_fit_t_j)/sum(imp_all_sp$ra_m_fit_t_i*imp_all_sp$ra_m_fit_t_j)
      mfunc_d.a = sum(imp_all_sp$Multi_traits*imp_all_sp$ra_m_fit_t_i*imp_all_sp$ra_m_fit_t_j)/sum(imp_all_sp$ra_m_fit_t_i*imp_all_sp$ra_m_fit_t_j)
      mconti_func_d.a = sum(imp_all_sp$Multi_conti_traits*imp_all_sp$ra_m_fit_t_i*imp_all_sp$ra_m_fit_t_j)/sum(imp_all_sp$ra_m_fit_t_i*imp_all_sp$ra_m_fit_t_j)
      mnnd = min(imp_all_sp$nd)
      mnlgfd = min(imp_all_sp$lgfd)
      mntd = min(imp_all_sp$Phylo_dis)
      mnfunc_d = min(imp_all_sp$Multi_traits)
      mnconti_func_d = min(imp_all_sp$Multi_conti_traits)
      dat_imp_sp_all = data.frame(f_p = unique(trans_plot$f_p),
                                  plot = unique(trans_plot$f_p),
                                  field = unique(trans_plot$field),
                                  species = nati_sp[j],
                                  gain_loss = unique(imp_all_sp$gain_loss.fittedtime_j),
                                  mean.relative.change = unique(imp_all_sp$mean.relative.change.fittedtime_j),
                                  mean.abs.change = unique(imp_all_sp$mean.absolute.change.fittedtime_j),
                                  stage = paste(unique(imp_all_sp$stage_i), collapse = '/'),
                                  estab = paste(unique(imp_all_sp$estab_i), collapse = '/'),
                                  domin = paste(unique(imp_all_sp$domin_i), collapse = '/'),
                                  mnd = mnd, mlgfd = mlgfd, mpd = mpd, mfunc_d = mfunc_d,
                                  mconti_func_d = mconti_func_d,
                                  mnd.a = mnd.a, mlgfd.a = mlgfd.a, mpd.a = mpd.a,
                                  mfunc_d.a = mfunc_d.a, mconti_func_d.a = mconti_func_d.a,
                                  mnnd = mnnd, mnlgfd = mnlgfd, mntd = mntd,
                                  mnfunc_d = mnfunc_d, mnconti_func_d = mnconti_func_d)
    } else {dat_imp_sp_all = NULL}
    
    dat_imp_sp_1 = rbind(dat_imp_sp_all, dat_imp_sp_domin,
                         dat_imp_sp_estab, dat_imp_sp_intro)
    dat_imp_sp_sum = rbind(dat_imp_sp_sum, dat_imp_sp_1)
  }
}

save(dat_imp_sp_sum,
     file = 'code/results_analyzing/analysing_alltime_data/dat_imp_sp_sum.rdata')

####### Fast load data #######
load("code/results_analyzing/analysing_alltime_data/dat_imp_sp_sum.rdata")
dat_imp_sp_sum$intro = 0
dat_imp_sp_sum = dat_imp_sp_sum %>% relocate(intro, .before = estab)
dat_imp_sp_sum[dat_imp_sp_sum$estab == 0 & dat_imp_sp_sum$domin == 0,]$intro = 1
dat_imp_sp_sum_c = dat_imp_sp_sum %>% filter(!is.na(gain_loss))
dat_imp_sp_sum_c$gain_loss_2 = dat_imp_sp_sum_c$gain_loss
dat_imp_sp_sum_c$mean.abs.change_2 = dat_imp_sp_sum_c$mean.abs.change*-1
dat_imp_sp_sum_c$mean.relative.change_2 = dat_imp_sp_sum_c$mean.relative.change*-1
dat_imp_sp_sum_c[dat_imp_sp_sum_c$gain_loss == 0,]$gain_loss_2 = 1
dat_imp_sp_sum_c[dat_imp_sp_sum_c$gain_loss == 1,]$gain_loss_2 = 0

########### Gain_loss 01 impact #############
library(ordinal)
library(lme4)

### Analyse for mnd mfd abundance weighted mean
dat_imp_sp_intro = dat_imp_sp_sum_c %>% filter(intro == 1)
mod_imp_sp_intro_ab = glmer(gain_loss_2 ~ mnd.a + mlgfd.a +
                              (1|field/plot) + (1|species),
                            data = dat_imp_sp_intro, family = binomial)
#mod_imp_sp_intro_ab_2 = glmer(gain_loss_2 ~ mnd.a + mlgfd.a +
#                             (1|field/plot),
#                          data = dat_imp_sp_intro, family = binomial)
#anova(mod_imp_sp_intro_ab, mod_imp_sp_intro_ab_2)
#mod_imp_sp_intro_ab_1 = glm(gain_loss_2 ~ mnd.a + mlgfd.a,
#                            data = dat_imp_sp_intro, family = binomial)
#anova(mod_imp_sp_intro_ab, mod_imp_sp_intro_ab_1)
summary(mod_imp_sp_intro_ab)

mod_imp_sp_intro = glmer(gain_loss_2 ~ mnd + mlgfd + (1|field/plot) + 
                           (1|species),
                         data = dat_imp_sp_intro, family = binomial)
#mod_imp_sp_intro_2 = glmer(gain_loss_2 ~ mnd + mlgfd + (1|field/plot),
#                         data = dat_imp_sp_intro, family = binomial)
#mod_imp_sp_intro_1 = glm(gain_loss_2 ~ mnd + mlgfd,
#                         data = dat_imp_sp_intro, family = binomial)
#anova(mod_imp_sp_intro, mod_imp_sp_intro_1)
#anova(mod_imp_sp_intro_2, mod_imp_sp_intro)
summary(mod_imp_sp_intro)

dat_imp_sp_estab = dat_imp_sp_sum_c %>%
  filter(estab == 1 & domin == 0 & !is.na(gain_loss_2))
mod_imp_sp_estab_ab = glmer(gain_loss_2 ~ mnd.a + mlgfd.a + (1|field/plot) + 
                              (1|species),
                            data = dat_imp_sp_estab, family = binomial)
#mod_imp_sp_estab_ab_2 = glmer(gain_loss_2 ~ mnd.a + mlgfd.a + (1|field/plot),
#                            data = dat_imp_sp_estab, family = binomial)
#mod_imp_sp_estab_ab_1 = glm(gain_loss_2 ~ mnd.a + mlgfd.a,
#                            data = dat_imp_sp_estab, family = binomial)
#anova(mod_imp_sp_estab_ab, mod_imp_sp_estab_ab_1)
#anova(mod_imp_sp_estab_ab_2, mod_imp_sp_estab_ab)
summary(mod_imp_sp_estab_ab)

mod_imp_sp_estab = glmer(gain_loss_2 ~ mnd + mlgfd + (1|field/plot) +
                           + (1|species),
                         data = dat_imp_sp_estab, family = binomial)
#mod_imp_sp_estab_1 = glm(gain_loss_2 ~ mnd + mlgfd,
#                         data = dat_imp_sp_estab, family = binomial)
#mod_imp_sp_estab_2 = glmer(gain_loss_2 ~ mnd + mlgfd + (1|field/plot),
#                           data = dat_imp_sp_estab, family = binomial)
#anova(mod_imp_sp_estab, mod_imp_sp_estab_1)
#anova(mod_imp_sp_estab, mod_imp_sp_estab_2)
summary(mod_imp_sp_estab)

dat_imp_sp_domin = dat_imp_sp_sum_c %>% filter(domin == 1) 
mod_imp_sp_domin_ab = glmer(gain_loss_2 ~ mnd.a + mlgfd.a + (1|field/plot) + 
                              (1|species),
                            data = dat_imp_sp_domin, family = binomial)
#mod_imp_sp_domin_ab_1 = glm(gain_loss_2 ~ mnd.a + mlgfd.a,
#                            data = dat_imp_sp_domin, family = binomial)
#mod_imp_sp_domin_ab_2 = glmer(gain_loss_2 ~ mnd.a + mlgfd.a + (1|field/plot),
#                            data = dat_imp_sp_domin, family = binomial)
#anova(mod_imp_sp_domin_ab, mod_imp_sp_domin_ab_1)
#anova(mod_imp_sp_domin_ab, mod_imp_sp_domin_ab_2)
summary(mod_imp_sp_domin_ab)
summary(mod_imp_sp_domin_ab_2)

mod_imp_sp_domin = glmer(gain_loss_2 ~ mnd + mlgfd + (1|field/plot) + 
                           (1|species),
                         data = dat_imp_sp_domin, family = binomial)
#mod_imp_sp_domin_2 = glmer(gain_loss_2 ~ mnd + mlgfd + (1|field/plot),
#                           data = dat_imp_sp_domin, family = binomial)
#mod_imp_sp_domin_1 = glm(gain_loss_2 ~ mnd + mlgfd,
#                         data = dat_imp_sp_domin, family = binomial)
anova(mod_imp_sp_domin, mod_imp_sp_domin_1)
anova(mod_imp_sp_domin, mod_imp_sp_domin_2)
summary(mod_imp_sp_domin)
summary(mod_imp_sp_domin_1)
summary(mod_imp_sp_domin_2)

dat_imp_sp_all = dat_imp_sp_sum_c %>% filter(stage != 'introduce' &
                                               stage != 'establish' &
                                               stage != 'dominant') 
mod_imp_sp_all_ab = glmer(gain_loss_2 ~ mnd.a + mlgfd.a + (1|field/plot) +  
                            (1|species),
                          data = dat_imp_sp_all, family = binomial) ## singular fit
library(brms)
mod_imp_sp_all_ab_brm = brm(gain_loss_2 ~ mnd.a + mlgfd.a + (1|field/plot) +  
                              (1|species), data = dat_imp_sp_all,
                            family = bernoulli,
                            chains = 2, cores = 2)
summary(mod_imp_sp_all_ab_brm)

mod_imp_sp_all_ab_1 = glm(gain_loss_2 ~ mnd.a + mlgfd.a,
                          data = dat_imp_sp_all, family = binomial)
mod_imp_sp_all_ab_2 = glmer(gain_loss_2 ~ mnd.a + mlgfd.a + (1|field/plot),
                            data = dat_imp_sp_all, family = binomial)## singular fit
anova(mod_imp_sp_all_ab, mod_imp_sp_all_ab_1)
anova(mod_imp_sp_all_ab, mod_imp_sp_all_ab_2)
summary(mod_imp_sp_all_ab)
summary(mod_imp_sp_all_ab_1)

mod_imp_sp_all = glmer(gain_loss_2 ~ mnd + mlgfd + (1|field/plot) +  
                         (1|species),
                       data = dat_imp_sp_all, family = binomial)
mod_imp_sp_all_2 = glmer(gain_loss_2 ~ mnd + mlgfd + (1|field/plot),
                         data = dat_imp_sp_all, family = binomial)
mod_imp_sp_all_1 = glm(gain_loss_2 ~ mnd + mlgfd,
                       data = dat_imp_sp_all, family = binomial)
anova(mod_imp_sp_all, mod_imp_sp_all_1)
anova(mod_imp_sp_all, mod_imp_sp_all_2)
summary(mod_imp_sp_all)
summary(mod_imp_sp_all_2)

# plot invasion impact
coeffs_mnd_imp = data.frame(stage = c('introduce', 
                                      'establishment',
                                      'dominant',
                                      'overall'),
                            nd = c(summary(mod_imp_sp_intro)$coefficients[2,1],
                                   summary(mod_imp_sp_estab)$coefficients[2,1],
                                   summary(mod_imp_sp_domin)$coefficients[2,1],
                                   summary(mod_imp_sp_all)$coefficients[2,1]),
                            sd = c(summary(mod_imp_sp_intro)$coefficients[2,2]*1.96,
                                   summary(mod_imp_sp_estab)$coefficients[2,2]*1.96,
                                   summary(mod_imp_sp_domin)$coefficients[2,2]*1.96,
                                   summary(mod_imp_sp_all)$coefficients[2,2]*1.96))

coeffs_mlgfd_imp = data.frame(stage = c('introduce', 
                                        'establishment',
                                        'dominant',
                                        'overall'),
                              lgfd = c(summary(mod_imp_sp_intro)$coefficients[3,1],
                                       summary(mod_imp_sp_estab)$coefficients[3,1],
                                       summary(mod_imp_sp_domin)$coefficients[3,1],
                                       summary(mod_imp_sp_all)$coefficients[3,1]),
                              sd = c(summary(mod_imp_sp_intro)$coefficients[3,2]*1.96,
                                     summary(mod_imp_sp_estab)$coefficients[3,2]*1.96,
                                     summary(mod_imp_sp_domin)$coefficients[3,2]*1.96,
                                     summary(mod_imp_sp_all)$coefficients[3,2]*1.96))


library(ggplot2)
p_mod_imp_nd = ggplot(coeffs_mnd_imp, aes(x=nd, y=stage)) + 
  geom_pointrange(aes(xmin=nd-sd, xmax=nd+sd), linewidth = 2,
                  fatten = 40) + 
  geom_vline(aes(xintercept=0), colour="#BB0000", linetype="dashed", size = 2) +
  ggtitle('Species level: Mean ND')+
  labs(x = 'Coefficients estimated', y = '') +
  scale_color_manual(labels=c(0, 1), values = c(2, 4))+
  #ylim(1, 3) +
  theme(plot.title = element_text(hjust = 0.5,
                                  vjust=0.5, size=150, face='bold'), 
        panel.background = element_rect(fill = 'white', color = 'black'),
        legend.key = element_rect(fill = "white"),
        panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        text = element_text(size = 200),
        axis.ticks.y = element_line(linetype=1,color="black",size=0.5),
        axis.ticks.x = element_line(linetype=1,color="black",size=0.5),
        axis.ticks.length.y = unit(0.3,'cm'),
        axis.ticks.length.x = unit(0.3,'cm'),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.7, vjust = 0.7),
        axis.text.y = element_text(angle = 0, vjust = 0.5),
        axis.title.x = element_text(margin=margin(t=20)))

p_mod_imp_lgfd = ggplot(coeffs_mlgfd_imp, aes(x=lgfd, y=stage)) + 
  geom_pointrange(aes(xmin=lgfd-sd, xmax=lgfd+sd), linewidth = 2,
                  fatten = 40) + 
  geom_vline(aes(xintercept=0), colour="#BB0000", linetype="dashed", size = 2)+
  ggtitle('Species level: Mean log10(FD)')+
  labs(x = 'Coefficients estimated', y = '') +
  scale_color_manual(labels=c(0, 1), values = c(2, 4))+
  #ylim(1, 3) +
  theme(plot.title = element_text(hjust = 0.5,
                                  vjust=0.5, size=150, face='bold'), 
        panel.background = element_rect(fill = 'white', color = 'black'),
        legend.key = element_rect(fill = "white"),
        panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        text = element_text(size = 200),
        axis.ticks.y = element_line(linetype=1,color="black",size=0.5),
        axis.ticks.x = element_line(linetype=1,color="black",size=0.5),
        axis.ticks.length.y = unit(0.3,'cm'),
        axis.ticks.length.x = unit(0.3,'cm'),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.7, vjust = 0.7),
        axis.text.y = element_text(angle = 0, vjust = 0.5),
        axis.title.x = element_text(margin=margin(t=20)))

# plot invasion success ab 
coeffs_mnd_imp_ab = data.frame(stage = c('introduce', 
                                         'establishment',
                                         'dominant',
                                         'overall'),
                               nd = c(summary(mod_imp_sp_intro_ab)$coefficients[2,1],
                                      summary(mod_imp_sp_estab_ab)$coefficients[2,1],
                                      summary(mod_imp_sp_domin_ab)$coefficients[2,1],
                                      summary(mod_imp_sp_all_ab)$coefficients[2,1]),
                               sd = c(summary(mod_imp_sp_intro_ab)$coefficients[2,2]*1.96,
                                      summary(mod_imp_sp_estab_ab)$coefficients[2,2]*1.96,
                                      summary(mod_imp_sp_domin_ab)$coefficients[2,2]*1.96,
                                      summary(mod_imp_sp_all_ab)$coefficients[2,2]*1.96))

coeffs_mlgfd_imp_ab = data.frame(stage = c('introduce', 
                                           'establishment',
                                           'dominant',
                                           'overall'),
                                 lgfd = c(summary(mod_imp_sp_intro_ab)$coefficients[3,1],
                                          summary(mod_imp_sp_estab_ab)$coefficients[3,1],
                                          summary(mod_imp_sp_domin_ab)$coefficients[3,1],
                                          summary(mod_imp_sp_all_ab)$coefficients[3,1]),
                                 sd = c(summary(mod_imp_sp_intro_ab)$coefficients[3,2]*1.96,
                                        summary(mod_imp_sp_estab_ab)$coefficients[3,2]*1.96,
                                        summary(mod_imp_sp_domin_ab)$coefficients[3,2]*1.96,
                                        summary(mod_imp_sp_all_ab)$coefficients[3,2]*1.96))
library(ggplot2)
p_mod_imp_nd_ab = ggplot(coeffs_mnd_imp_ab, aes(x=nd, y=stage)) + 
  geom_pointrange(aes(xmin=nd-sd, xmax=nd+sd), linewidth = 2,
                  fatten = 40) + 
  geom_vline(aes(xintercept=0), colour="#BB0000", linetype="dashed", size = 2) +
  ggtitle(paste('Species level: Mean ',
                expression(ND[ab])))+
  labs(x = 'Coefficients estimated', y = '') +
  scale_color_manual(labels=c(0, 1), values = c(2, 4))+
  #ylim(1, 3) +
  theme(plot.title = element_text(hjust = 0.5,
                                  vjust=0.5, size=150, face='bold'), 
        panel.background = element_rect(fill = 'white', color = 'black'),
        legend.key = element_rect(fill = "white"),
        panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        text = element_text(size = 200),
        axis.ticks.y = element_line(linetype=1,color="black",size=0.5),
        axis.ticks.x = element_line(linetype=1,color="black",size=0.5),
        axis.ticks.length.y = unit(0.3,'cm'),
        axis.ticks.length.x = unit(0.3,'cm'),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.7, vjust = 0.7),
        axis.text.y = element_text(angle = 0, vjust = 0.5),
        axis.title.x = element_text(margin=margin(t=20)))

p_mod_imp_lgfd_ab = ggplot(coeffs_mlgfd_imp_ab, aes(x=lgfd, y=stage)) + 
  geom_pointrange(aes(xmin=lgfd-sd, xmax=lgfd+sd), linewidth = 2,
                  fatten = 40) + 
  geom_vline(aes(xintercept=0), colour="#BB0000", linetype="dashed", size = 2)+
  ggtitle(paste('Species level: Mean ',
                expression(log10(FD)[ab])))+
  labs(x = 'Coefficients estimated', y = '') +
  scale_color_manual(labels=c(0, 1), values = c(2, 4))+
  #ylim(1, 3) +
  theme(plot.title = element_text(hjust = 0.5,
                                  vjust=0.5, size=150, face='bold'), 
        panel.background = element_rect(fill = 'white', color = 'black'),
        legend.key = element_rect(fill = "white"),
        panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        text = element_text(size = 200),
        axis.ticks.y = element_line(linetype=1,color="black",size=0.5),
        axis.ticks.x = element_line(linetype=1,color="black",size=0.5),
        axis.ticks.length.y = unit(0.3,'cm'),
        axis.ticks.length.x = unit(0.3,'cm'),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.7, vjust = 0.7),
        axis.text.y = element_text(angle = 0, vjust = 0.5),
        axis.title.x = element_text(margin=margin(t=20)))

## Merge plot and save them ##
library(ggpubr)

Fig.3 = ggarrange(p_mod_imp_nd, p_mod_imp_lgfd, nrow = 2, ncol = 1,
                  labels = c('a)', 'b)'), font.label = list(size = 200))
Fig.3.ab = ggarrange(p_mod_imp_nd_ab, p_mod_imp_lgfd_ab, nrow = 2, ncol = 1,
                     labels = c('a)', 'b)'), font.label = list(size = 200))

ggsave(plot = Fig.3,
       'results/figures_sameages_raw/Fig.3.svg',
       width = 160,height = 240,dpi = 300, units = 'cm',
       limitsize = F)

ggsave(plot = Fig.3.ab,
       'results/figures_sameages_raw/Fig.3.ab.svg',
       width = 160,height = 240,dpi = 300, units = 'cm',
       limitsize = F)

########### Abundance change: continuous impact #############
library(ordinal)
library(lme4)
library(INLA)

## Rescale data
numcols = grep("^m.+d",names(dat_imp_sp_sum_c))
dat_imp_sp_sum_cs = dat_imp_sp_sum_c
dat_imp_sp_sum_cs[,numcols] = scale(dat_imp_sp_sum_cs[,numcols])

#### estab ####
dat_imp_sp_sum_cs$species_1 = as.factor(dat_imp_sp_sum_cs$species)

tree = read.tree('data/original data/phylo_tree332.txt')

pc_prior=list(prec=list("pc.prec",param=c(0.1,0.01)))
gammaprior=list(prec=list(prior="loggamma",param=c(0.01,0.01)))

## Intro ##
### Analyse for mnd mfd abundance weighted mean
dat_imp_sp_intro = dat_imp_sp_sum_cs %>% filter(intro == 1)
intro_sp_names = unique(dat_imp_sp_intro$species)
dat_imp_sp_intro$species_1 = as.factor(dat_imp_sp_intro$species)

intro_tree_fit = keep.tip(tree, intro_sp_names)
intro_vcv_tree = ape::vcv(intro_tree_fit, model = "Brownian", corr = FALSE)
intro_vcv_tree_sparse = inla.as.sparse(solve(intro_vcv_tree))

### Singular fit using lmer, so choose bayesian methods ###
intro_inla.model_all1 = inla(mean.relative.change_2~mnd+mlgfd+mpd+mconti_func_d+
                               f(species, model = "iid", hyper = pc_prior)+
                               f(field, model="iid", hyper = pc_prior) +
                               f(f_p, model="iid", hyper = pc_prior)+
                               f(species_1, model="generic0",
                                 Cmatrix= intro_vcv_tree_sparse,
                                 values = intro_sp_names, hyper=pc_prior),
                             family="gaussian",data=dat_imp_sp_intro,
                             control.compute=list(dic=T,waic=T,cpo=T),
                             quantiles=c(0.025,0.15,0.5,0.85,0.975))
summary(intro_inla.model_all1)

intro_inla.model_all1_md.a = inla(mean.relative.change_2~mnd.a+mlgfd.a+mpd.a+mconti_func_d.a+
                                    f(species, model = "iid", hyper = pc_prior)+
                                    f(field, model="iid", hyper = pc_prior) +
                                    f(f_p, model="iid", hyper = pc_prior)+
                                    f(species_1, model="generic0",
                                      Cmatrix = intro_vcv_tree_sparse,
                                      values = intro_sp_names, hyper=pc_prior),
                                  family="gaussian",data=dat_imp_sp_intro,
                                  control.compute=list(dic=T,waic=T,cpo=T),
                                  quantiles=c(0.025,0.15,0.5,0.85,0.975))
summary(intro_inla.model_all1_md.a)

intro_inla.model_all1_mnd = inla(mean.relative.change_2~mnnd+mnlgfd+mntd+mnconti_func_d+
                                   f(species, model = "iid", hyper = pc_prior)+
                                   f(field, model="iid", hyper = pc_prior) +
                                   f(f_p, model="iid", hyper = pc_prior)+
                                   f(species_1, model="generic0",
                                     Cmatrix = intro_vcv_tree_sparse,
                                     values = intro_sp_names, hyper=pc_prior),
                                 family="gaussian",data=dat_imp_sp_intro,
                                 control.compute=list(dic=T,waic=T,cpo=T),
                                 quantiles=c(0.025,0.15,0.5,0.85,0.975))
summary(intro_inla.model_all1_mnd)
c(intro_inla.model_all1$waic$waic, intro_inla.model_all1_md.a$waic$waic,
  intro_inla.model_all1_mnd$waic$waic) # no significant difference
## thus use the lowest: intro_inla.model_all1_mnd


# plot 
intro_imp_data.inla.all.varied_intercept1 = intro_inla.model_all1_mnd$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower1="0.025quant",lower2="0.15quant",median="0.5quant",upper2="0.85quant",upper1="0.975quant")

intro_imp_data.inla.all.varied_intercept = intro_imp_data.inla.all.varied_intercept1%>%
  mutate(rowname=c("mnd","mlgfd","mpd","mconti_func_d"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

# point + effect size
intro_imp_inla.varied.intercept.plot =
  ggplot(data=intro_imp_data.inla.all.varied_intercept,aes(x=mean,y=rowname,color=rowname))+
  geom_point(size=80)+
  ggtitle('intro_implishment')+
  geom_linerange(aes(xmin=lower1,xmax=upper1),size=4.8)+
  geom_linerange(aes(xmin=lower2,xmax=upper2),size=12)+
  geom_vline(xintercept=0,linetype=2,color="grey40",linewidth=6)+
  scale_y_discrete(limits=arrange(intro_imp_data.inla.all.varied_intercept,percent)$rowname)+#f	g'%R2gfe:fc