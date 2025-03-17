#### SEM: Success ~ mND + mRFD ~ mfd + mpd for all species ####
### Original model ###
rm(list = ls())
require(glmmTMB)
library(piecewiseSEM)
require(optimx)
require(dplyr)

### SEM for establishment
setwd("D:/R projects/BSS")
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_data/dat_suc_sp.rdata')

numcols = grep("^m.",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])
dat_suc_sps$species_1 = as.factor(dat_suc_sps$species)


estab_sem1_md.a_all = psem(
  glmmTMB(estab~mnd.a+mlgfd.a
          +mpd.a_all
          +mconti_func_d.a_all
          #+(1|species)
          +(1|field)
          #+(1|f_p)
          , family=binomial, data=dat_suc_sps,
          control=glmmTMBControl(optimizer=optim,
                                 optArgs=list(method="BFGS"))),
  glmmTMB(mnd.a~mpd.a_all+
            #mconti_func_d.a_all+
            #(1|species)+
            #+(1|field)
            (1|f_p)
          , family=gaussian, data = dat_suc_sps),
  glmmTMB(mlgfd.a~mpd.a_all+
            mconti_func_d.a_all+
            #(1|species)+
            (1|field)
          +(1|f_p)
          , family=gaussian, data=dat_suc_sps),
  mlgfd.a %~~% mnd.a,
  #estab %~~% mpd,
  #estab %~~% mconti_func_d,
  data=dat_suc_sps)
summary(estab_sem1_md.a_all)

### SEM for Dominance ###
dat_dom_sps = dat_suc_sps %>% filter(stage %in% c("establish", "dominant"))

domin_sem1_md.a_all = psem(
  glmmTMB(domin~mnd.a+
            mlgfd.a+
            mpd.a_all+
            mconti_func_d.a_all+
            #(1|species)+
            (1|field)
          +(1|f_p)
          , family=binomial, data=dat_dom_sps,
          control=glmmTMBControl(optimizer=optim,
                                 optArgs=list(method="BFGS"))),
  glmmTMB(mnd.a~mpd.a_all+
            mconti_func_d.a_all+
            #(1|species)+
            (1|field)
          #+(1|f_p)
          , family=gaussian, data = dat_dom_sps),
  glmmTMB(mlgfd.a~mpd.a_all+
            #mconti_func_d.a_all+
            (1|species)+
            (1|field)
          #+(1|f_p)
          , family=gaussian, data=dat_dom_sps),
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
  c(sd(dat_suc_sps$mpd.a_all),
    sd(dat_suc_sps$mconti_func_d.a_all),
    sd(dat_suc_sps$mnd),
    sd(dat_suc_sps$mlgfd))/sd.y_estab_sem1_md.a_all
std.estab_sem1_md.a_all

# merge the new corrected cofficients into the original cofficients data frame
coefs_estab_sem1_md.a_all=coefs(estab_sem1_md.a_all)# original coefficients of piecewiseSEM
#coefs_estab_sem1_md.a_all[1:4,8]=std.estab_sem1_md.a_all

### Dominance corrected std coefficients in the glmer, following nakagawa 2013 MEE
sd.y_domin_sem1_md.a_all=sqrt(predict(domin_sem1_md.a_all[[1]],re.form=NA)%>%var+sum(unlist(VarCorr(domin_sem1_md.a_all[[1]])))+pi^2/3)
sd.y_domin_sem1_md.a_all
std.domin_sem1_md.a_all=unlist(fixef(domin_sem1_md.a_all[[1]]))[-1]*
  c(sd(dat_dom_sps$mpd.a_all),
    sd(dat_dom_sps$mconti_func_d.a_all),
    sd(dat_dom_sps$mnd),
    sd(dat_dom_sps$mlgfd))/sd.y_domin_sem1_md.a_all
std.domin_sem1_md.a_all

# merge the new corrected cofficients into the original cofficients data frame
coefs_domin_sem1_md.a_all=coefs(domin_sem1_md.a_all)# original coefficients of piecewiseSEM
#coefs_domin_sem1_md.a_all[1:4,8]=std.domin_sem1_md.a_all
coefs_domin_sem1_md.a_all

### Calculate direct and indirect effects of SEM following Xu Meng (2022) GCB
### Estab
coefs_estab_sem1_md.a_all = coefs_estab_sem1_md.a_all[-nrow(coefs_estab_sem1_md.a_all),]
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
coefs_domin_sem1_md.a_all = coefs_domin_sem1_md.a_all[-nrow(coefs_domin_sem1_md.a_all),]
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
sem_estab_effects_md.a_all.percent.plot = 
  ggplot(data=sem_estab_effects_md.a_all.percent,aes(percent,predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),hjust=-0.1,size=3)+
  facet_wrap(~"Relative total effect")+
  scale_y_discrete(limits=c("mconti_func_d.a_all", "mpd.a_all", "mlgfd.a", "mnd.a"),
                   labels=#c(bquote('I-N MFD'[ab]),bquote('I-N MPD'[ab]),bquote('I-N MRFD'[ab]),bquote('I-N MND'[ab]))
                     c(expression(MFD[ab]),
                       expression(MPD[ab]),
                       expression(RFD[ab]),
                       expression(ND[ab])) ,
                   position="right")+
  scale_fill_manual(values = c("#ccb974ff", "#dd8452ff", "#4c72b0ff", "#55a868ff"))+# color palatte in ggthemes
  xlim(0,0.5)+
  xlab("")+
  #theme_test()+
  theme(axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        strip.text=element_text(face="bold"),
        plot.margin=unit(c(0.4,0.4,0.2,-0.2),units="lines"))+
  guides(fill="none")
sem_estab_effects_md.a_all.percent.plot

# Direct and indirect effects
sem_estab_effects_md.a_all.dir_indir.mpd_all.plot = 
  ggplot(data=sem_estab_effects_md.a_all.longer,aes(x=value,y=predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_vline(xintercept = 0, linetype = 2)+
  facet_wrap(~type,labeller=as_labeller(c("direct"="Direct effects","indirect"="Indirect effects")),
             scales = 'free')+
  ggh4x::facetted_pos_scales(
    x = list(
      `direct` = scale_x_continuous(limits = c(-0.44, 0.44), breaks = seq(-0.4, 0.4, 0.2)),
      `indirect` = scale_x_continuous(limits = c(-0.14, 0.14), breaks = seq(-0.1, 0.1, 0.1))  
    )
  )+
  scale_y_discrete(limits=c("mconti_func_d.a_all", "mpd.a_all", "mlgfd.a", "mnd.a"),
                   position="right")+
  scale_fill_manual(values = c("#ccb974ff", "#dd8452ff", "#4c72b0ff", "#55a868ff"))+# color palatte in ggthemes
  guides(fill="none")+
  #theme_test()+
  theme(axis.title=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        strip.text=element_text(face="bold"))
sem_estab_effects_md.a_all.dir_indir.mpd_all.plot

### Space for SEM
sem_estab_md.a_all_space1 = ggplot(NULL)+theme_void()
sem_estab_md.a_all_space2 = ggplot(NULL)+theme_void()

### domin
sem_domin_effects_md.a_all_total = sem_domin_effects_md.a_all %>% filter(type == 'total')
sem_domin_effects_md.a_all.percent = sem_domin_effects_md.a_all_total%>%
  mutate(percent=abs(value)/sum(abs(value)),value=NULL)
sem_domin_effects_md.a_all.longer = sem_domin_effects_md.a_all %>% filter(type != 'total' 
                                                                          #&value != 0
)

# Relative total effect
sem_domin_effects_md.a_all.percent.plot = 
  ggplot(data=sem_domin_effects_md.a_all.percent,aes(percent,predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),hjust=-0.1,size=3)+
  facet_wrap(~"Relative total effect")+
  scale_y_discrete(limits=c("mconti_func_d.a_all", "mpd.a_all", "mlgfd.a", "mnd.a"),
                   labels=#c(bquote('E-N MFD'[ab]),bquote('E-N MPD'[ab]),bquote('E-N MRFD'[ab]),bquote('E-N MND'[ab]))
                     c(expression(MFD[ab]),
                       expression(MPD[ab]),
                       expression(RFD[ab]),
                       expression(ND[ab])),
                   position="right")+
  scale_fill_manual(values = c("#ccb974ff", "#dd8452ff", "#4c72b0ff", "#55a868ff"))+# color palatte in ggthemes
  xlim(0,0.72)+
  xlab("")+
  #theme_test()+
  theme(axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        strip.text=element_text(face="bold"),
        plot.margin=unit(c(0.4,0.4,0.2,-0.2),units="lines"))+
  guides(fill="none")

# Direct and indirect effects
sem_domin_effects_md.a_all.dir_indir.mpd_all.plot = 
  ggplot(data=sem_domin_effects_md.a_all.longer,aes(x=value,y=predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_vline(xintercept = 0, linetype = 2)+
  facet_wrap(~type,labeller=as_labeller(c("direct"="Direct effects",
                                          "indirect"="Indirect effects")),
             scales = 'free')+
  ggh4x::facetted_pos_scales(
    x = list(
      `direct` = scale_x_continuous(limits = c(-0.44, 0.44), breaks = seq(-0.4, 0.4, 0.2)),
      `indirect` = scale_x_continuous(limits = c(-0.24, 0.24), breaks = seq(-0.2, 0.2, 0.1))  
    )
  )+
  scale_y_discrete(limits=c("mconti_func_d.a_all", "mpd.a_all", "mlgfd.a", "mnd.a"),
                   position="right")+
  scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2),
                     limits = c(-0.4, 0.4))+
  scale_fill_manual(values = c("#ccb974ff", "#dd8452ff", "#4c72b0ff", "#55a868ff"))+# color palatte in ggthemes
  guides(fill="none")+
  #theme_test()+
  theme(axis.title=element_blank(),axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text=element_text(face="bold"))

### Space for SEM
sem_domin_md.a_all_space1 = ggplot(NULL)+theme_void()
sem_domin_md.a_all_space2 = ggplot(NULL)+theme_void()

# Merge and export eight plots
sem_md.a_all_plot_all = ggarrange(sem_estab_md.a_all_space1, sem_estab_md.a_all_space2,
                                  sem_estab_effects_md.a_all.dir_indir.mpd_all.plot,
                                  sem_estab_effects_md.a_all.percent.plot,
                                  sem_domin_md.a_all_space1,sem_domin_md.a_all_space2,
                                  sem_domin_effects_md.a_all.dir_indir.mpd_all.plot,
                                  sem_domin_effects_md.a_all.percent.plot,
                                  nrow = 4, ncol = 2,
                                  heights = c(0.35, 0.15, 0.35, 0.15),
                                  labels = c('(a)', '', '','',
                                             '(b)', '', '',''),
                                  vjust = 1.8)

emf('results/figures_ages1_35_top50_equal_interval_bh/sem_md.a_all_plot_all.emf',
    width=16, height=24, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
sem_md.a_all_plot_all
dev.off() #turn off device and finalize file
