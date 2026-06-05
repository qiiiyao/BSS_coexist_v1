#### SEM: Success ~ mND + mRFD ~ mfd + mpd for all species ####
### Original model ###
rm(list = ls())
setwd("D:/R projects/BSS_coexist")

require(glmmTMB)
require(piecewiseSEM)
require(optimx)
require(dplyr)
require(svglite)
require(ggpubr)
require(cowplot)
source('code/function/plot_func.R')


colors_4d_dim = c(
  "#CDB867",  # muted mustard
  "#DD9368",  # muted orange
  "#6F8FCD",  # muted blue
  "#74B383"   # muted green
)

### SEM for establishment

load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/allRealTraits/dat_suc_sp.rdata')

numcols = grep("^m.",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])
dat_suc_sps$species_1 = as.factor(dat_suc_sps$species)
nrow(dat_suc_sps)

estab_sem1_md.a_all = psem(
  glmmTMB(estab~mnd.a+mlgfd.a
          +mpd.a_all
          +mconti_func_d.a_all
          #+(1|species)
          +(1|field)
          #+(1|f_p)
          , family = binomial, data = dat_suc_sps,
          control = glmmTMBControl(optimizer=optim,
                                 optArgs=list(method="BFGS"))),
  glmmTMB(mnd.a~mpd.a_all+
            mconti_func_d.a_all+
            #(1|species)+
            #(1|field)+
            (1|f_p)
          , family=gaussian, data = dat_suc_sps),
  glmmTMB(mlgfd.a~#mpd.a_all+
            mconti_func_d.a_all+
            (1|species)+
            #(1|field)#+
            (1|f_p)
          , family=gaussian, data=dat_suc_sps),
  mpd.a_all %~~% mconti_func_d.a_all,
  mlgfd.a %~~% mnd.a,
  #estab %~~% mpd,
  #estab %~~% mconti_func_d,
  data=dat_suc_sps)
summary(estab_sem1_md.a_all)

### SEM for Dominance ###
dat_dom_sp = dat_suc_sp %>% filter(stage %in% c("establish", "dominant"))
numcols = grep("^m.",names(dat_dom_sp))
dat_dom_sps = dat_dom_sp
dat_dom_sps[,numcols] = scale(dat_dom_sps[,numcols])
dat_dom_sps$species_1 = as.factor(dat_dom_sps$species)
nrow(dat_dom_sps)

domin_sem1_md.a_all = psem(
  glmmTMB(domin~mnd.a+
            mlgfd.a+
            mpd.a_all+
            mconti_func_d.a_all+
            #(1|species)
            (1|field)
          +(1|f_p)
          , family=binomial, data=dat_dom_sps,
          control=glmmTMBControl(optimizer=optim,
                                 optArgs=list(method="BFGS"))),
  glmmTMB(mnd.a~mpd.a_all+
            mconti_func_d.a_all+
            #(1|species)+
            #(1|field)
            +(1|f_p)
          , family=gaussian, data = dat_dom_sps),
  glmmTMB(mlgfd.a~#mpd.a_all+
            mconti_func_d.a_all+  
            #(1|species)+
            # (1|field)
            +(1|f_p)
          , family=gaussian, data=dat_dom_sps),
  mpd.a_all %~~% mconti_func_d.a_all,
  mlgfd.a %~~% mnd.a,
  #domin %~~% mpd,
  #domin %~~% mconti_func_d,
  data=dat_dom_sps)
summary(domin_sem1_md.a_all)


##### Draw plot for md.a ##### 
require(ggh4x)
### The followed code is adapted from Xu_et_al_2022_GCB
### Establishment corrected std coefficients in the glmer, following nakagawa 2013 MEE
sd.y_estab_sem1_md.a_all = sqrt(predict(estab_sem1_md.a_all[[1]],
                                        re.form=NA)%>%var+sum(unlist(VarCorr(estab_sem1_md.a_all[[1]])))+pi^2/3)
sd.y_estab_sem1_md.a_all
std.estab_sem1_md.a_all=unlist(fixef(estab_sem1_md.a_all[[1]]))[-1]*
  c(sd(dat_suc_sps$mnd.a),
    sd(dat_suc_sps$mlgfd.a),
    sd(dat_suc_sps$mpd.a_all),
    sd(dat_suc_sps$mconti_func_d.a_all))/sd.y_estab_sem1_md.a_all
std.estab_sem1_md.a_all

# merge the new corrected cofficients into the original cofficients data frame
coefs_estab_sem1_md.a_all=coefs(estab_sem1_md.a_all)# original coefficients of piecewiseSEM
coefs_estab_sem1_md.a_all[1:4,8]=std.estab_sem1_md.a_all

### Dominance corrected std coefficients in the glmer, following nakagawa 2013 MEE
sd.y_domin_sem1_md.a_all=sqrt(predict(domin_sem1_md.a_all[[1]],re.form=NA)%>%var+sum(unlist(VarCorr(domin_sem1_md.a_all[[1]])))+pi^2/3)
sd.y_domin_sem1_md.a_all
std.domin_sem1_md.a_all=unlist(fixef(domin_sem1_md.a_all[[1]]))[-1]*
  c(sd(dat_dom_sps$mnd.a),
    sd(dat_dom_sps$mlgfd.a),
    sd(dat_dom_sps$mpd.a_all),
    sd(dat_dom_sps$mconti_func_d.a_all))/sd.y_domin_sem1_md.a_all
std.domin_sem1_md.a_all

# merge the new corrected cofficients into the original cofficients data frame
coefs_domin_sem1_md.a_all=coefs(domin_sem1_md.a_all)# original coefficients of piecewiseSEM
coefs_domin_sem1_md.a_all[1:4,8]=std.domin_sem1_md.a_all
coefs_domin_sem1_md.a_all

### Calculate direct and indirect effects of SEM following Xu Meng (2022) GCB
### Estab
coefs_estab_sem1_md.a_all = coefs_estab_sem1_md.a_all[-c(nrow(coefs_estab_sem1_md.a_all),
                                                         (nrow(coefs_estab_sem1_md.a_all)-1)),]

sem_estab_effects_md.a_all = data.frame(predictor = rep(c('mnd.a', 'mlgfd.a',
                                                          'mpd.a_all', 'mconti_func_d.a_all'), 3),
                                        type = rep(c('direct', 'indirect', 'total'),
                                                   each = 4),
                                        value = 0,
                                        significant = 'yes')

coefs_estab_sem1_md.a_all_c = coefs_estab_sem1_md.a_all[coefs_estab_sem1_md.a_all$P.Value < 0.05,]
for (i in 1:length(unique(coefs_estab_sem1_md.a_all_c$Predictor))) {
  #i = 4
  predictor_1 = unique(coefs_estab_sem1_md.a_all_c$Predictor)[i]
  dat = coefs_estab_sem1_md.a_all_c[coefs_estab_sem1_md.a_all_c$Predictor == predictor_1,]
  dat_2 = coefs_estab_sem1_md.a_all_c[coefs_estab_sem1_md.a_all_c$Predictor %in% dat$Response,]
  dat_3 = rbind(dat, dat_2)
  if (nrow(dat_3) == 1) {
    sem_estab_effects_md.a_all[sem_estab_effects_md.a_all$predictor == predictor_1&
                                 sem_estab_effects_md.a_all$type == 'direct',]$value = dat_3$Std.Estimate
    sem_estab_effects_md.a_all[sem_estab_effects_md.a_all$predictor == predictor_1&
                                 sem_estab_effects_md.a_all$type == 'total',]$value = dat_3$Std.Estimate
  } else if(nrow(dat_3) != 1) {
    if (nrow(dat_3[dat_3$Response == 'estab' &
                   dat_3$Predictor == predictor_1,]) > 0) {
      sem_estab_effects_md.a_all[sem_estab_effects_md.a_all$predictor == predictor_1&
                                   sem_estab_effects_md.a_all$type == 'direct',]$value = dat_3[dat_3$Response == 'estab' &
                                                                                                 dat_3$Predictor == predictor_1,]$Std.Estimate
    }
    
    sem_estab_effects_md.a_all[sem_estab_effects_md.a_all$predictor == predictor_1&
                                 sem_estab_effects_md.a_all$type == 'indirect',]$value = sum(dat_3[dat_3$Response != 'estab',]$Std.Estimate*dat_3[dat_3$Response == 'estab' &
                                                                                                                                                    dat_3$Predictor != predictor_1,]$Std.Estimate)
    sem_estab_effects_md.a_all[sem_estab_effects_md.a_all$predictor == predictor_1&
                                 sem_estab_effects_md.a_all$type == 'total',]$value  = sem_estab_effects_md.a_all[sem_estab_effects_md.a_all$predictor == predictor_1&
                                                                                                                    sem_estab_effects_md.a_all$type == 'direct',]$value + sem_estab_effects_md.a_all[sem_estab_effects_md.a_all$predictor == predictor_1&
                                                                                                                                                                                                       sem_estab_effects_md.a_all$type == 'indirect',]$value
  }
}

### Domin
coefs_domin_sem1_md.a_all = coefs_domin_sem1_md.a_all[-c(nrow(coefs_domin_sem1_md.a_all),
                                                         (nrow(coefs_domin_sem1_md.a_all)-1)),]
sem_domin_effects_md.a_all = data.frame(predictor = rep(c('mnd.a', 'mlgfd.a',
                                                          'mpd.a_all', 'mconti_func_d.a_all'), 3),
                                        type = rep(c('direct', 'indirect', 'total'),
                                                   each = 4),
                                        value = 0,
                                        significant = 'yes')

coefs_domin_sem1_md.a_all_c = coefs_domin_sem1_md.a_all[coefs_domin_sem1_md.a_all$P.Value < 0.05,]
for (i in 1:length(unique(coefs_domin_sem1_md.a_all_c$Predictor))) {
  #i = 4
  predictor_1 = unique(coefs_domin_sem1_md.a_all_c$Predictor)[i]
  dat = coefs_domin_sem1_md.a_all_c[coefs_domin_sem1_md.a_all_c$Predictor == predictor_1,]
  dat_2 = coefs_domin_sem1_md.a_all_c[coefs_domin_sem1_md.a_all_c$Predictor %in% dat$Response,]
  dat_3 = rbind(dat, dat_2)
  if (nrow(dat_3) == 1) {
    sem_domin_effects_md.a_all[sem_domin_effects_md.a_all$predictor == predictor_1&
                                 sem_domin_effects_md.a_all$type == 'direct',]$value = dat_3$Std.Estimate
    sem_domin_effects_md.a_all[sem_domin_effects_md.a_all$predictor == predictor_1&
                                 sem_domin_effects_md.a_all$type == 'total',]$value = dat_3$Std.Estimate
  } else if(nrow(dat_3) != 1) {
    sem_domin_effects_md.a_all[sem_domin_effects_md.a_all$predictor == predictor_1&
                                 sem_domin_effects_md.a_all$type == 'direct',]$value = dat_3[dat_3$Response == 'domin' &
                                                                                               dat_3$Predictor == predictor_1,]$Std.Estimate
    sem_domin_effects_md.a_all[sem_domin_effects_md.a_all$predictor == predictor_1&
                                 sem_domin_effects_md.a_all$type == 'indirect',]$value = sum(dat_3[dat_3$Response != 'domin',]$Std.Estimate*dat_3[dat_3$Response == 'domin' &
                                                                                                                                                    dat_3$Predictor != predictor_1,]$Std.Estimate)
    sem_domin_effects_md.a_all[sem_domin_effects_md.a_all$predictor == predictor_1&
                                 sem_domin_effects_md.a_all$type == 'total',]$value  = sem_domin_effects_md.a_all[sem_domin_effects_md.a_all$predictor == predictor_1&
                                                                                                                    sem_domin_effects_md.a_all$type == 'direct',]$value + sem_domin_effects_md.a_all[sem_domin_effects_md.a_all$predictor == predictor_1&
                                                                                                                                                                                                       sem_domin_effects_md.a_all$type == 'indirect',]$value
  }
}



#-----------Draw plots for SEMs' relative total effects, direct and indirect effects 
### Estab
require(ggthemes)
require(ggplot2)
require(devEMF)
sem_estab_effects_md.a_all_total = sem_estab_effects_md.a_all %>% filter(type == 'total')
sem_estab_effects_md.a_all.percent = sem_estab_effects_md.a_all_total%>%
  mutate(percent=abs(value)/sum(abs(value)),value=NULL)
sem_estab_effects_md.a_all.longer = sem_estab_effects_md.a_all %>% filter(type != 'total' 
                                                                          #& value != 0
)

# Relative total effect 
#ggthemr::ggthemr(palette = "fresh", layout = "clean")
sem_estab_effects_md.a_all.percent.plot = 
  ggplot(data=sem_estab_effects_md.a_all.percent,aes(percent,predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_vline(xintercept = 0, linetype = 3, color="grey40", linewidth = 0.8)+
  geom_vline(xintercept = -0.05, linetype = 2, color="#93785B", linewidth = 0.6)+
  geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),hjust=-0.1)+
  ggtitle('Relative total effect')+
  scale_x_continuous(limits = c(-0.05, 0.6), breaks = seq(0, 0.5, 0.25),
                     labels = function(x) sprintf("%.2f", x))+
  scale_y_discrete(limits=c("mconti_func_d.a_all", "mpd.a_all", "mlgfd.a", "mnd.a"),
                   labels = '' ,
                   position="right")+
  scale_fill_manual(values = colors_4d_dim)+# color palatte in ggthemes
  xlab("")+
  ylab("")+
  #theme_test()+
  theme_regular_2()+
  annotate(
    geom = 'segment',
    y = Inf,
    yend = Inf,
    x = -Inf,
    xend = Inf,
    , color="#93785B", linewidth = 1)+
  theme(plot.title = element_text(family = 'Arial', face = 'plain', hjust = 0.5,
                                  colour = "black", size = 14),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color="#93785B", linewidth = 0.6),
        plot.margin = margin(t = 0.05, r = 0.1, b = -0.5, l = -0.1, unit = "cm"))

sem_estab_effects_md.a_all.percent.plot

# Direct and indirect effects
#ggthemr::ggthemr(palette = "fresh", layout = "clean")
sem_estab_effects_md.a_all.dir.mpd_all.plot = 
  ggplot((data=sem_estab_effects_md.a_all.longer%>% 
            filter(type == 'direct')),aes(x=value,y=predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_vline(xintercept = 0, linetype = 3, color="grey40", linewidth = 0.8)+
  scale_x_continuous(limits = c(-0.44, 0.5), breaks = seq(-0.3, 0.3, 0.3),
                     labels = function(x) sprintf("%.2f", x))+
  scale_y_discrete(limits=c("mconti_func_d.a_all", "mpd.a_all", "mlgfd.a", "mnd.a"),
                   labels=#c(bquote('I-N MFD'[ab]),bquote('I-N MPD'[ab]),bquote('I-N MRFD'[ab]),bquote('I-N MND'[ab]))
                     c(expression(MFD[ab]),
                       expression(MPD[ab]),
                       expression(RFD[ab]),
                       expression(ND[ab])) ,
                   position="left")+
  scale_fill_manual(values = colors_4d_dim)+# color palatte in ggthemes
  ggtitle("Direct effects")+
  xlab("")+
  ylab("")+
  #theme_test()+
  theme_regular_2() +
  annotate(
    geom = 'segment',
    y = Inf,
    yend = Inf,
    x = -Inf,
    xend = Inf,
    , color="#93785B", linewidth = 1
  ) +
  theme(plot.title = element_text(family = 'Arial', face = 'plain', hjust = 0.5,
                                  colour = "black", size = 14),
        panel.border = element_blank(),
        axis.line = element_line(color="#93785B", linewidth = 0.6),
        plot.margin = margin(t = 0.05, r = -0.05, b = -0.5, l = 0.15, unit = "cm"))
sem_estab_effects_md.a_all.dir.mpd_all.plot

sem_estab_effects_md.a_all.indir.mpd_all.plot = 
  ggplot((data=sem_estab_effects_md.a_all.longer %>% 
            filter(type == 'indirect')),aes(x=value,y=predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_vline(xintercept = 0, linetype = 3, color="grey40", linewidth = 0.8)+
  geom_vline(xintercept = -0.031, linetype = 2, color="#93785B", linewidth = 0.6)+
  scale_x_continuous(limits = c(-0.031, 0.031), breaks = seq(-0.02, 0.02, 0.02),
                     labels = function(x) sprintf("%.2f", x))+
  scale_y_discrete(limits=c("mconti_func_d.a_all", "mpd.a_all", "mlgfd.a", "mnd.a"),
                   position="right")+
  scale_fill_manual(values = colors_4d_dim)+# color palatte in ggthemes
  ggtitle("Indirect effects")+
  xlab("")+
  ylab("")+
  #theme_test()+
  theme_regular_2() +
  theme(plot.title = element_text(family = 'Arial', face = 'plain', hjust = 0.5,
                                  colour = "black", size = 14),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.length.y = unit(0, "pt"),
        axis.ticks.length.y.left = unit(0, "pt"),
        axis.ticks.length.y.right = unit(0, "pt"))  + 
  theme(plot.margin = margin(t = 0.05, r = 0, b = -0.5, l = -0.2, unit = "cm")) + 
  annotate(
    geom = 'segment',
    y = Inf,
    yend = Inf,
    x = -Inf,
    xend = Inf,
    , color="#93785B", linewidth = 1
  ) + 
  theme(panel.border = element_blank(),
        axis.line = element_line(color="#93785B", linewidth = 0.6))

### Space for SEM
sem_estab_md.a_all_space1 = ggplot(NULL)+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(0,1))+
  geom_hline(yintercept = 0.98, linetype = 1, color="#93785B", linewidth = 0.6)+
  theme_regular_2()+
  theme(plot.margin = margin(t = 0.05, r = 0.2, b = 0.05, l = 0.2, unit = "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

### domin
sem_domin_effects_md.a_all_total = sem_domin_effects_md.a_all %>% filter(type == 'total')
sem_domin_effects_md.a_all.percent = sem_domin_effects_md.a_all_total%>%
  mutate(percent=abs(value)/sum(abs(value)),value=NULL)
sem_domin_effects_md.a_all.longer = sem_domin_effects_md.a_all %>% filter(type != 'total' 
                                                                          #&value != 0
)

# Relative total effect
#ggthemr::ggthemr(palette = "fresh", layout = "clean")
sem_domin_effects_md.a_all.percent.plot = 
  ggplot(data=sem_domin_effects_md.a_all.percent,aes(percent,predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_vline(xintercept = 0, linetype = 3, color="grey40", linewidth = 0.8)+
  geom_vline(xintercept = -0.05, linetype = 2, color="#93785B", linewidth = 0.6)+
  geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),hjust=-0.1)+
  ggtitle("Relative total effect")+
  scale_x_continuous(limits = c(-0.05, 0.72), breaks = seq(0, 0.5, 0.25),
                     labels = function(x) sprintf("%.2f", x))+
  scale_y_discrete(limits=c("mconti_func_d.a_all", "mpd.a_all", "mlgfd.a", "mnd.a"),
                   labels=#c(bquote('E-N MFD'[ab]),bquote('E-N MPD'[ab]),bquote('E-N MRFD'[ab]),bquote('E-N MND'[ab]))
                     c(expression(MFD[ab]),
                       expression(MPD[ab]),
                       expression(RFD[ab]),
                       expression(ND[ab])),
                   position="right")+
  scale_fill_manual(values = colors_4d_dim)+# color palatte in ggthemes
  xlab("")+
  ylab("")+
  #theme_test()+
  theme_regular_2() +
  annotate(
    geom = 'segment',
    y = Inf,
    yend = Inf,
    x = -Inf,
    xend = Inf,
    , color="#93785B", linewidth = 1
  ) + 
  theme(plot.title = element_text(family = 'Arial', face = 'plain', hjust = 0.5,
                                  colour = "black", size = 14),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color="#93785B", linewidth = 0.6),
        plot.margin = margin(t = 0.05, r = 0.1, b = -0.5, l = -0.1, unit = "cm"))
sem_domin_effects_md.a_all.percent.plot

# Direct and indirect effects
sem_domin_effects_md.a_all.dir.mpd_all.plot = 
  ggplot((data=sem_domin_effects_md.a_all.longer%>% 
            filter(type == 'direct')),aes(x=value,y=predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_vline(xintercept = 0, linetype = 3, color="grey40", linewidth = 0.8)+
  scale_x_continuous(limits = c(-0.14, 0.14), breaks = seq(-0.1, 0.1, 0.1),
                     labels = function(x) sprintf("%.2f", x))+
  scale_y_discrete(limits = c("mconti_func_d.a_all", "mpd.a_all", "mlgfd.a", "mnd.a"),
                   labels = NULL,
                   position = "left")+
  scale_fill_manual(values = colors_4d_dim)+# color palatte in ggthemes
  ggtitle("Direct effects")+
  xlab("")+
  ylab("")+
  #theme_test()+
  theme_regular_2() +
  annotate(
    geom = 'segment',
    y = Inf,
    yend = Inf,
    x = -Inf,
    xend = Inf,
    , color="#93785B", linewidth = 1
  ) +
  theme(plot.title = element_text(family = 'Arial', face = 'plain', hjust = 0.5,
                                  colour = "black", size = 14),
        panel.border = element_blank(),
        axis.line = element_line(color="#93785B", linewidth = 0.6),
        plot.margin = margin(t = 0.05, r = 0.85, b = -0.5, l = -0.75, unit = "cm"))
sem_domin_effects_md.a_all.dir.mpd_all.plot

sem_domin_effects_md.a_all.indir.mpd_all.plot = 
  ggplot(data=(sem_domin_effects_md.a_all.longer %>% 
                 filter(type == 'indirect')),aes(x=value,y=predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_vline(xintercept = 0, linetype = 3, color="grey40", linewidth = 0.8)+
  geom_vline(xintercept = -0.14, linetype = 2, color="#93785B", linewidth = 0.6)+
  scale_x_continuous(limits = c(-0.14, 0.14), breaks = c(0),
                     labels = function(x) sprintf("%.2f", x))+
  scale_y_discrete(limits=c("mconti_func_d.a_all", "mpd.a_all", "mlgfd.a", "mnd.a"),
                   position="right")+
  scale_fill_manual(values = colors_4d_dim)+# color palatte in ggthemes
  guides(fill="none")+
  ggtitle("Indirect effects")+
  #theme_test()+
  xlab("")+
  ylab("")+
  #theme_test()+
  theme_regular_2() +
  theme(plot.title = element_text(family = 'Arial', face = 'plain', hjust = 0.5,
                                  colour = "black", size = 14),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.length.y = unit(0, "pt"),
        axis.ticks.length.y.left = unit(0, "pt"),
        axis.ticks.length.y.right = unit(0, "pt"))  + 
  theme(plot.margin = margin(t = 0.05, r = 0, b = -0.5, l = -0.2, unit = "cm"))+ 
  annotate(
    geom = 'segment',
    y = Inf,
    yend = Inf,
    x = -Inf,
    xend = Inf,
    , color="#93785B", linewidth = 1
  ) + 
  theme(panel.border = element_blank(),
        axis.line = element_line(color="#93785B", linewidth = 0.6))
sem_domin_effects_md.a_all.indir.mpd_all.plot

### Space for SEM
sem_domin_md.a_all_space1 = ggplot(NULL)+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(0,1))+
  geom_hline(yintercept = 0.98, linetype = 1, color="#93785B", linewidth = 0.6)+
  theme_regular_2()+
  theme(plot.margin = margin(t = 0.05, r = 0.2, b = 0.05, l = 0, unit = "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Merge and export eight plots
sem_md.a_all_plot_all = ggdraw() +
  
  # Upper SEM plots
  draw_plot(sem_estab_md.a_all_space1,
            x = 0.045, y = 0.195,
            width = 0.4, height = 0.525) +
  
  draw_plot(sem_domin_md.a_all_space1,
            x = 0.438, y = 0.195,
            width = 0.4, height = 0.525) +
  
  
  # =========================
# Lower effect plots
# =========================
draw_plot(sem_estab_effects_md.a_all.dir.mpd_all.plot,
          x = -0.006, y = 0.03,
          width = 0.18, height = 0.16) +
  
  draw_plot(sem_estab_effects_md.a_all.indir.mpd_all.plot,
            x = 0.170, y = 0.03,
            width = 0.12, height = 0.16) +
  
  draw_plot(sem_estab_effects_md.a_all.percent.plot,
            x = 0.280, y = 0.03,
            width = 0.17, height = 0.16) +
  
  draw_plot(sem_domin_effects_md.a_all.dir.mpd_all.plot,
            x = 0.439, y = 0.03,
            width = 0.15, height = 0.16) +
  
  draw_plot(sem_domin_effects_md.a_all.indir.mpd_all.plot,
            x = 0.564, y = 0.03,
            width = 0.12, height = 0.16) +
  
  draw_plot(sem_domin_effects_md.a_all.percent.plot,
            x = 0.674, y = 0.03,
            width = 0.1605, height = 0.16) +
  
  
  # =========================
# Labels
# =========================
draw_label("A   Establishment",
           x = 0.055, y = 0.708,
           size = 15,
           fontface = "plain",
           hjust = 0, vjust = 1) +
  
  draw_label("B   Dominance",
             x = 0.442, y = 0.7075,
             size = 15,
             fontface = "plain",
             hjust = 0, vjust = 1)

# Print the figure
sem_md.a_all_plot_all

#emf('results/figures/sem_md.a_all_plot_all.emf',
# width=16 * 1.2, height=24 * 1.2, coordDPI = 1200, 
# emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
# ensure text looks correct on the viewing system
# units = 'cm')
#sem_md.a_all_plot_all
#dev.off() #turn off device and finalize file

ggsave('results/figures/sem_md.a_all_plot_allRealTraits.svg',
       width=20 * 1.6, height=22 * 1.6, 
       units = 'cm')

