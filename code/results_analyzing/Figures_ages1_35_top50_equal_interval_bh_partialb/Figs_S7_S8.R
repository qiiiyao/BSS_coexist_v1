#### Fast start for analyzing mnd ~ mpd, mnd ~ mconti_func_d, mfd ~ mpd, mfd ~ mconti_func_d ####
rm(list = ls())

require(INLA)
require(phyr)
require(inlabru)
require(tibble)
require(lme4)
require(nlme)
require(lmerTest)
require(minpack.lm)

## Functions for plot 
source('code/function/plot_func.R')

load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/dat_suc_sp.rdata')

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




#### Relationships among abundance weighted mean differences ####
##### predictive curve for mnd.a ~ mpd.a_all #####
lincombs.data.mnd.a.mpd.a_all = data.frame(mpd.a_all=seq(min(dat_suc_sp$mpd.a_all),
                                                         max(dat_suc_sp$mpd.a_all),
                                                         length=100),
                                           mconti_func_d.a_all=mean(dat_suc_sp$mconti_func_d.a_all))

lincombs.matrix.mnd.a.mpd.a_all=model.matrix(~mpd.a_all,
                                             data=lincombs.data.mnd.a.mpd.a_all)
lincombs.matrix.mnd.a.mpd.a_all=as.data.frame(lincombs.matrix.mnd.a.mpd.a_all)
lincombs.mnd.a.mpd.a_all=inla.make.lincombs(lincombs.matrix.mnd.a.mpd.a_all)

inla.model_lincombs.mnd.a.mpd.a_all = pglmm(mnd.a ~ mpd.a_all +#(1|species) + 
                                              (1|f_p) + (1|field), data = dat_suc_sp,
                                            family = "gaussian", cov_ranef = list(species = tree),
                                            bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                        config = TRUE),
                                                                 quantiles=c(0.025,0.5,0.975),
                                                                 lincomb=lincombs.mnd.a.mpd.a_all,
                                                                 control.predictor=list(compute=T)),
                                            bayes = T)

inla.model_lincombs.mnd.a.mpd.a_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mnd.a.mpd.a_all$predicted.value=inla.model_lincombs.mnd.a.mpd.a_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mnd.a.mpd.a_all$lower=inla.model_lincombs.mnd.a.mpd.a_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mnd.a.mpd.a_all$upper=inla.model_lincombs.mnd.a.mpd.a_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mnd.a.mpd.a_all

save(lincombs.data.mnd.a.mpd.a_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnd.a.mpd.a_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnd.a.mpd.a_all.rdata')
(Fig.S7_mnd.a_mpd.a_all = ggplot(data=lincombs.data.mnd.a.mpd.a_all, aes(x=mpd.a_all, y=predicted.value))+
    geom_point(data=dat_suc_sp, aes(x=mpd.a_all, y=mnd.a),
               color = 'grey',
               alpha = 0.5,
               position=position_jitter(height=0.001))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='darkgray',alpha=0.2)+
    geom_line(color='grey1', linetype = 1, linewidth = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #scale_y_continuous(breaks = seq(0, 0.8, length.out = 5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    labs(x = expression('MPD'[ab]),
         y = expression('ND'[ab]))+
    annotate(geom="text",
             x=c((min(dat_suc_sp$mpd.a_all)+max(dat_suc_sp$mpd.a_all))/2,
                 (min(dat_suc_sp$mpd.a_all)+max(dat_suc_sp$mpd.a_all))/2),
             y=c(-0.68,-0.8),
             label=c("italic(β)['MPD'[ab]] == '-0.00002'",
                     "'95%CI' == '[-0.00004, -0.00001]'"),
             parse=T,size=4)+
    theme_regular()
)
# mpd.a_all   -0.00002   -0.00004   -0.00001

##### predictive curve for mnd.a ~ mconti_func_d.a_all #####
lincombs.data.mnd.a.mconti_func_d.a_all = data.frame(mconti_func_d.a_all=seq(min(dat_suc_sp$mconti_func_d.a_all),
                                                                 max(dat_suc_sp$mconti_func_d.a_all),
                                                                 length=100),
                                               mpd.a_all=mean(dat_suc_sp$mpd.a_all))

lincombs.matrix.mnd.a.mconti_func_d.a_all=model.matrix(~mconti_func_d.a_all,
                                                 data=lincombs.data.mnd.a.mconti_func_d.a_all)
lincombs.matrix.mnd.a.mconti_func_d.a_all=as.data.frame(lincombs.matrix.mnd.a.mconti_func_d.a_all)
lincombs.mnd.a.mconti_func_d.a_all=inla.make.lincombs(lincombs.matrix.mnd.a.mconti_func_d.a_all)

inla.model_lincombs.mnd.a.mconti_func_d.a_all = pglmm(mnd.a ~ mconti_func_d.a_all+#(1|species) + 
                                                  (1|f_p) + (1|field), data = dat_suc_sp,
                                                family = "gaussian", cov_ranef = list(species = tree),
                                                bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                            config = TRUE),
                                                                     quantiles=c(0.025,0.5,0.975),
                                                                     lincomb=lincombs.mnd.a.mconti_func_d.a_all,
                                                                     control.predictor=list(compute=T)),
                                                bayes = T)

inla.model_lincombs.mnd.a.mconti_func_d.a_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mnd.a.mconti_func_d.a_all$predicted.value=inla.model_lincombs.mnd.a.mconti_func_d.a_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mnd.a.mconti_func_d.a_all$lower=inla.model_lincombs.mnd.a.mconti_func_d.a_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mnd.a.mconti_func_d.a_all$upper=inla.model_lincombs.mnd.a.mconti_func_d.a_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mnd.a.mconti_func_d.a_all

save(lincombs.data.mnd.a.mconti_func_d.a_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnd.a.mconti_func_d.a_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnd.a.mconti_func_d.a_all.rdata')
(Fig.S7_mnd.a_mconti_func_d.a_all = ggplot(data=lincombs.data.mnd.a.mconti_func_d.a_all, 
                                     aes(x=mconti_func_d.a_all, y=predicted.value))+
    geom_point(data=dat_suc_sp, aes(x=mconti_func_d.a_all, y=mnd.a),
               #shape=1,
               alpha = 0.5,
               color = 'grey',
               position=position_jitter(height=0.001))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='darkgray',alpha=0.2)+
    geom_line(color='grey1', linetype = 1, linewidth = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #scale_y_continuous(breaks = seq(0, 0.8, length.out = 5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    
    labs(x= expression('MFD'[ab]),
         y= ' ')+
    annotate(geom="text",
             x=c((min(dat_suc_sp$mconti_func_d.a_all)+max(dat_suc_sp$mconti_func_d.a_all))/2,
                 (min(dat_suc_sp$mconti_func_d.a_all)+max(dat_suc_sp$mconti_func_d.a_all))/2),
             y=c(-0.68,-0.8),
             label=c("italic(β)['MFD'[ab]] == -0.00223",
                     "'95%CI' == '[-0.00429, -0.00017]'"),
             parse=T,size=4)+
    theme_regular()
)
#mconti_func_d.a_all -0.00223   -0.00429   -0.00017

##### predictive curve for mlgfd.a ~ mpd.a_all #####
lincombs.data.mlgfd.a.mpd.a_all = data.frame(mpd.a_all=seq(min(dat_suc_sp$mpd.a_all),
                                                           max(dat_suc_sp$mpd.a_all),
                                                           length=100),
                                             mconti_func_d.a_all=mean(dat_suc_sp$mconti_func_d.a_all))

lincombs.matrix.mlgfd.a.mpd.a_all=model.matrix(~mpd.a_all,
                                               data=lincombs.data.mlgfd.a.mpd.a_all)
lincombs.matrix.mlgfd.a.mpd.a_all=as.data.frame(lincombs.matrix.mlgfd.a.mpd.a_all)
lincombs.mlgfd.a.mpd.a_all=inla.make.lincombs(lincombs.matrix.mlgfd.a.mpd.a_all)

inla.model_lincombs.mlgfd.a.mpd.a_all = pglmm(mlgfd.a ~ mpd.a_all+#(1|species) + 
                                                (1|f_p) + (1|field), data = dat_suc_sp,
                                              family = "gaussian", cov_ranef = list(species = tree),
                                              bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                          config = TRUE),
                                                                   quantiles=c(0.025,0.5,0.975),
                                                                   lincomb=lincombs.mlgfd.a.mpd.a_all,
                                                                   control.predictor=list(compute=T)),
                                              bayes = T)


inla.model_lincombs.mlgfd.a.mpd.a_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mlgfd.a.mpd.a_all$predicted.value=inla.model_lincombs.mlgfd.a.mpd.a_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mlgfd.a.mpd.a_all$lower=inla.model_lincombs.mlgfd.a.mpd.a_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mlgfd.a.mpd.a_all$upper=inla.model_lincombs.mlgfd.a.mpd.a_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mlgfd.a.mpd.a_all

save(lincombs.data.mlgfd.a.mpd.a_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mlgfd.a.mpd.a_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mlgfd.a.mpd.a_all.rdata')
(Fig.S7_mlgfd.a_mpd.a_all = ggplot(data=lincombs.data.mlgfd.a.mpd.a_all, aes(x=mpd.a_all, y=predicted.value))+
    geom_point(data=dat_suc_sp, aes(x=mpd.a_all, y=mlgfd.a),
               alpha = 0.5,
               color = 'grey',
               position=position_jitter(height=0.001))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='darkgray',alpha=0.2)+
    geom_line(color='grey1', linetype = 1, linewidth = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    scale_y_continuous(breaks = seq(-2.5, 2, 0.5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    labs(x= ' ',
         y=expression('RFD'[ab]))+
    annotate(geom="text",
             x=c((min(dat_suc_sp$mpd.a_all)+max(dat_suc_sp$mpd.a_all))/2,
                 (min(dat_suc_sp$mpd.a_all)+max(dat_suc_sp$mpd.a_all))/2),
             y=c(-1.7,-2.0),
             label=c("italic(β)['MPD'[ab]] == 0.00024",
                     "'95%CI' == '[-0.00031, -0.00017]'"),
             parse=T,size=4)+
    theme_regular()    
)
# mpd.a_all   -0.00024   -0.00031   -0.00017

##### predictive curve for mlgfd.a ~ mconti_func_d.a_all #####
lincombs.data.mlgfd.a.mconti_func_d.a_all = data.frame(mconti_func_d.a_all=seq(
                                                                   min(dat_suc_sp$mconti_func_d.a_all),
                                                                   max(dat_suc_sp$mconti_func_d.a_all),
                                                                   length=100),
                                                 mpd.a_all=mean(dat_suc_sp$mpd.a_all))

lincombs.matrix.mlgfd.a.mconti_func_d.a_all=model.matrix(~mconti_func_d.a_all,
                                                   data=lincombs.data.mlgfd.a.mconti_func_d.a_all)
lincombs.matrix.mlgfd.a.mconti_func_d.a_all=as.data.frame(lincombs.matrix.mlgfd.a.mconti_func_d.a_all)
lincombs.mlgfd.a.mconti_func_d.a_all=inla.make.lincombs(lincombs.matrix.mlgfd.a.mconti_func_d.a_all)

inla.model_lincombs.mlgfd.a.mconti_func_d.a_all = pglmm(mlgfd.a ~ mconti_func_d.a_all+#(1|species) + 
                                                    (1|f_p) + (1|field), data = dat_suc_sp,
                                                  family = "gaussian", cov_ranef = list(species = tree),
                                                  bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                              config = TRUE),
                                                                       quantiles=c(0.025,0.5,0.975),
                                                                       lincomb=lincombs.mlgfd.a.mconti_func_d.a_all,
                                                                       control.predictor=list(compute=T)),
                                                  bayes = T)

inla.model_lincombs.mlgfd.a.mconti_func_d.a_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mlgfd.a.mconti_func_d.a_all$predicted.value=inla.model_lincombs.mlgfd.a.mconti_func_d.a_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mlgfd.a.mconti_func_d.a_all$lower=inla.model_lincombs.mlgfd.a.mconti_func_d.a_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mlgfd.a.mconti_func_d.a_all$upper=inla.model_lincombs.mlgfd.a.mconti_func_d.a_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mlgfd.a.mconti_func_d.a_all

save(lincombs.data.mlgfd.a.mconti_func_d.a_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mlgfd.a.mconti_func_d.a_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mlgfd.a.mconti_func_d.a_all.rdata')
(Fig.S7_mlgfd.a_mconti_func_d.a_all = ggplot(data=lincombs.data.mlgfd.a.mconti_func_d.a_all, 
                                       aes(x=mconti_func_d.a_all, y=predicted.value))+
    geom_point(data=dat_suc_sp, aes(x=mconti_func_d.a_all, y=mlgfd.a),
               alpha = 0.5,
               color = 'grey',
               position=position_jitter(height=0.001))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='darkgray',alpha=0.2)+
    geom_line(color='grey1', linetype = 1, linewidth = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    scale_y_continuous(breaks = seq(-2.5, 2, 0.5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    labs(x= ' ',
         y= ' ')+
    annotate(geom="text",
             x=c((min(dat_suc_sp$mconti_func_d.a_all)+max(dat_suc_sp$mconti_func_d.a_all))/2,
                 (min(dat_suc_sp$mconti_func_d.a_all)+max(dat_suc_sp$mconti_func_d.a_all))/2),
             y=c(-1.7,-2.0),
             label=c("italic(β)['MFD'[ab]] == -0.02736",
                     "'95%CI' == '[-0.03487, -0.01985]'"),
             parse=T,size=4)+
    theme_regular()    
)
#mconti_func_d.a_all -0.02736   -0.03487   -0.01985

###### merge plots for Fig. S7 ######
library(cowplot)

Fig.S7 = ggdraw() +
  draw_plot(Fig.S7_mlgfd.a_mpd.a_all, 0.02, 0.5, .48, .46) +
  draw_plot(Fig.S7_mlgfd.a_mconti_func_d.a_all, 0.5, 0.5, .48, .46) +
  draw_plot(Fig.S7_mnd.a_mpd.a_all, 0.02, 0.08, .48, .46) +
  draw_plot(Fig.S7_mnd.a_mconti_func_d.a_all, 0.5, 0.08, .48, .46) +
  draw_plot_label(c("a", "b", "c", "d"),
                  c(0.015+0.01, 0.495+0.01,
                    0.015+0.01, 0.495+0.01),
                  c(0.97, 0.97, 0.55, 0.55)
                  ,size = 18)

emf('results/Figures_ages1_35_top50_equal_interval_bh_partialb/Fig.S7.emf',
    width = 20, height = 20, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.S7
dev.off() #turn off device and finalize file


#### Relationships among mean differences ####
##### predictive curve for mnd ~ mpd_all #####
lincombs.data.mnd.mpd_all = data.frame(mpd_all=seq(min(dat_suc_sp$mpd_all),
                                                         max(dat_suc_sp$mpd_all),
                                                         length=100),
                                           mconti_func_d_all=mean(dat_suc_sp$mconti_func_d_all))

lincombs.matrix.mnd.mpd_all=model.matrix(~mpd_all,
                                             data=lincombs.data.mnd.mpd_all)
lincombs.matrix.mnd.mpd_all=as.data.frame(lincombs.matrix.mnd.mpd_all)
lincombs.mnd.mpd_all=inla.make.lincombs(lincombs.matrix.mnd.mpd_all)

inla.model_lincombs.mnd.mpd_all = pglmm(mnd ~ mpd_all +#(1|species) + 
                                              (1|f_p) + (1|field), data = dat_suc_sp,
                                            family = "gaussian", cov_ranef = list(species = tree),
                                            bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                        config = TRUE),
                                                                 quantiles=c(0.025,0.5,0.975),
                                                                 lincomb=lincombs.mnd.mpd_all,
                                                                 control.predictor=list(compute=T)),
                                            bayes = T)

inla.model_lincombs.mnd.mpd_all$inla.model$summary.fixed[c(1,3,5)]%>%round(7)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mnd.mpd_all$predicted.value=inla.model_lincombs.mnd.mpd_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mnd.mpd_all$lower=inla.model_lincombs.mnd.mpd_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mnd.mpd_all$upper=inla.model_lincombs.mnd.mpd_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mnd.mpd_all

save(lincombs.data.mnd.mpd_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnd.mpd_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnd.mpd_all.rdata')
(Fig.S8_mnd_mpd_all = ggplot(data=lincombs.data.mnd.mpd_all, aes(x=mpd_all, y=predicted.value))+
    geom_point(data=dat_suc_sp, aes(x=mpd_all, y=mnd),
               color = 'grey',
               alpha = 0.5,
               position=position_jitter(height=0.001))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='darkgray',alpha=0.2)+
    geom_line(color='grey1', linetype = 1, linewidth = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    labs(x = expression('MPD'),
         y = expression('ND'))+
    annotate(geom="text",
             x=c((min(dat_suc_sp$mpd_all)+max(dat_suc_sp$mpd_all))/2,
                 (min(dat_suc_sp$mpd_all)+max(dat_suc_sp$mpd_all))/2),
             y=c(-0.43,-0.5),
             label=c("italic(β)['MPD'] == '-0.0000166'",
                     "'95%CI' == '[-0.0000330, -0.0000002]'"),
             parse=T,size=4)+
    theme_regular()
)
#mpd_all     -0.00002   -0.00003    0.00000

##### predictive curve for mnd ~ mconti_func_d_all #####
lincombs.data.mnd.mconti_func_d_all = data.frame(mconti_func_d_all=seq(min(dat_suc_sp$mconti_func_d_all),
                                                                             max(dat_suc_sp$mconti_func_d_all),
                                                                             length=100),
                                                     mpd_all=mean(dat_suc_sp$mpd_all))

lincombs.matrix.mnd.mconti_func_d_all=model.matrix(~mconti_func_d_all,
                                                       data=lincombs.data.mnd.mconti_func_d_all)
lincombs.matrix.mnd.mconti_func_d_all=as.data.frame(lincombs.matrix.mnd.mconti_func_d_all)
lincombs.mnd.mconti_func_d_all=inla.make.lincombs(lincombs.matrix.mnd.mconti_func_d_all)

inla.model_lincombs.mnd.mconti_func_d_all = pglmm(mnd ~ mconti_func_d_all+#(1|species) + 
                                                        (1|f_p) + (1|field), data = dat_suc_sp,
                                                      family = "gaussian", cov_ranef = list(species = tree),
                                                      bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                                  config = TRUE),
                                                                           quantiles=c(0.025,0.5,0.975),
                                                                           lincomb=lincombs.mnd.mconti_func_d_all,
                                                                           control.predictor=list(compute=T)),
                                                      bayes = T)

inla.model_lincombs.mnd.mconti_func_d_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mnd.mconti_func_d_all$predicted.value=inla.model_lincombs.mnd.mconti_func_d_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mnd.mconti_func_d_all$lower=inla.model_lincombs.mnd.mconti_func_d_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mnd.mconti_func_d_all$upper=inla.model_lincombs.mnd.mconti_func_d_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mnd.mconti_func_d_all

save(lincombs.data.mnd.mconti_func_d_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnd.mconti_func_d_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnd.mconti_func_d_all.rdata')
(Fig.S8_mnd_mconti_func_d_all = ggplot(data=lincombs.data.mnd.mconti_func_d_all, 
                                           aes(x=mconti_func_d_all, y=predicted.value))+
    geom_point(data=dat_suc_sp, aes(x=mconti_func_d_all, y=mnd),
               #shape=1,
               alpha = 0.5,
               color = 'grey',
               position=position_jitter(height=0.001))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='darkgray',alpha=0.2)+
    geom_line(color='grey1', linetype = 1, linewidth = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    labs(x= expression('MFD'),
         y= ' ')+
    annotate(geom="text",
             x=c((min(dat_suc_sp$mconti_func_d_all)+max(dat_suc_sp$mconti_func_d_all))/2,
                 (min(dat_suc_sp$mconti_func_d_all)+max(dat_suc_sp$mconti_func_d_all))/2),
             y=c(-0.43,-0.5),
             label=c("italic(β)['MFD'] == '0.00270'",
                     "'95%CI' == '[0.00108, 0.00432]'"),
             parse=T,size=4)+
    theme_regular()
)
#mconti_func_d_all  0.00270    0.00108    0.00432

##### predictive curve for mlgfd ~ mpd_all #####
lincombs.data.mlgfd.mpd_all = data.frame(mpd_all=seq(min(dat_suc_sp$mpd_all),
                                                           max(dat_suc_sp$mpd_all),
                                                           length=100),
                                             mconti_func_d_all=mean(dat_suc_sp$mconti_func_d_all))

lincombs.matrix.mlgfd.mpd_all=model.matrix(~mpd_all,
                                               data=lincombs.data.mlgfd.mpd_all)
lincombs.matrix.mlgfd.mpd_all=as.data.frame(lincombs.matrix.mlgfd.mpd_all)
lincombs.mlgfd.mpd_all=inla.make.lincombs(lincombs.matrix.mlgfd.mpd_all)

inla.model_lincombs.mlgfd.mpd_all = pglmm(mlgfd ~ mpd_all+#(1|species) + 
                                                (1|f_p) + (1|field), data = dat_suc_sp,
                                              family = "gaussian", cov_ranef = list(species = tree),
                                              bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                          config = TRUE),
                                                                   quantiles=c(0.025,0.5,0.975),
                                                                   lincomb=lincombs.mlgfd.mpd_all,
                                                                   control.predictor=list(compute=T)),
                                              bayes = T)


inla.model_lincombs.mlgfd.mpd_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mlgfd.mpd_all$predicted.value=inla.model_lincombs.mlgfd.mpd_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mlgfd.mpd_all$lower=inla.model_lincombs.mlgfd.mpd_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mlgfd.mpd_all$upper=inla.model_lincombs.mlgfd.mpd_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mlgfd.mpd_all

save(lincombs.data.mlgfd.mpd_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mlgfd.mpd_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mlgfd.mpd_all.rdata')
(Fig.S8_mlgfd_mpd_all = ggplot(data=lincombs.data.mlgfd.mpd_all, aes(x=mpd_all, y=predicted.value))+
    geom_point(data=dat_suc_sp, aes(x=mpd_all, y=mlgfd),
               alpha = 0.5,
               color = 'grey',
               position=position_jitter(height=0.001))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='darkgray',alpha=0.2)+
    geom_line(color='grey1', linetype = 1, linewidth = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    scale_y_continuous(breaks = seq(-0.5, 1.5, 0.5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    labs(x= ' ',
         y=expression('RFD'))+
    annotate(geom="text",
             x=c((min(dat_suc_sp$mpd_all)+max(dat_suc_sp$mpd_all))/2,
                 (min(dat_suc_sp$mpd_all)+max(dat_suc_sp$mpd_all))/2),
             y=c(-0.83,-1),
             label=c("italic(β)['MPD'] == -'0.00026'",
                     "'95%CI' == '[-0.00036, -0.00017]'"),
             parse=T,size=4)+
    theme_regular()    
)
# mpd_all  -0.00026   -0.00036   -0.00017

##### predictive curve for mlgfd ~ mconti_func_d_all #####
lincombs.data.mlgfd.mconti_func_d_all = data.frame(mconti_func_d_all=seq(
  min(dat_suc_sp$mconti_func_d_all),
  max(dat_suc_sp$mconti_func_d_all),
  length=100),
  mpd_all=mean(dat_suc_sp$mpd_all))

lincombs.matrix.mlgfd.mconti_func_d_all=model.matrix(~mconti_func_d_all,
                                                         data=lincombs.data.mlgfd.mconti_func_d_all)
lincombs.matrix.mlgfd.mconti_func_d_all=as.data.frame(lincombs.matrix.mlgfd.mconti_func_d_all)
lincombs.mlgfd.mconti_func_d_all=inla.make.lincombs(lincombs.matrix.mlgfd.mconti_func_d_all)

inla.model_lincombs.mlgfd.mconti_func_d_all = pglmm(mlgfd ~ mconti_func_d_all+#(1|species) + 
                                                          (1|f_p) + (1|field), data = dat_suc_sp,
                                                        family = "gaussian", cov_ranef = list(species = tree),
                                                        bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                                    config = TRUE),
                                                                             quantiles=c(0.025,0.5,0.975),
                                                                             lincomb=lincombs.mlgfd.mconti_func_d_all,
                                                                             control.predictor=list(compute=T)),
                                                        bayes = T)

inla.model_lincombs.mlgfd.mconti_func_d_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mlgfd.mconti_func_d_all$predicted.value=inla.model_lincombs.mlgfd.mconti_func_d_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mlgfd.mconti_func_d_all$lower=inla.model_lincombs.mlgfd.mconti_func_d_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mlgfd.mconti_func_d_all$upper=inla.model_lincombs.mlgfd.mconti_func_d_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mlgfd.mconti_func_d_all

save(lincombs.data.mlgfd.mconti_func_d_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mlgfd.mconti_func_d_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mlgfd.mconti_func_d_all.rdata')
(Fig.S8_mlgfd_mconti_func_d_all = ggplot(data=lincombs.data.mlgfd.mconti_func_d_all, 
                                             aes(x=mconti_func_d_all, y=predicted.value))+
    geom_point(data=dat_suc_sp, aes(x=mconti_func_d_all, y=mlgfd),
               alpha = 0.5,
               color = 'grey',
               position=position_jitter(height=0.001))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='darkgray',alpha=0.2)+
    geom_line(color='grey1', linetype = 2, linewidth = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    scale_y_continuous(breaks = seq(-0.5, 1.5, 0.5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    labs(x= ' ',
         y= ' ')+
    annotate(geom="text",
             x=c((min(dat_suc_sp$mconti_func_d_all)+max(dat_suc_sp$mconti_func_d_all))/2,
                 (min(dat_suc_sp$mconti_func_d_all)+max(dat_suc_sp$mconti_func_d_all))/2),
             y=c(-0.83,-1),
             label=c("italic(β)['MFD'] == '0.00050'",
                     "'95%CI' == '[-0.00833, 0.00933]'"),
             parse=T,size=4)+
    theme_regular()    
)
# mconti_func_d_all  0.00050   -0.00833    0.00933

###### merge plots for Fig. S8 ######
library(cowplot)

Fig.S8 = ggdraw() +
  draw_plot(Fig.S8_mlgfd_mpd_all, 0.02, 0.5, .48, .46) +
  draw_plot(Fig.S8_mlgfd_mconti_func_d_all, 0.5, 0.5, .48, .46) +
  draw_plot(Fig.S8_mnd_mpd_all, 0.02, 0.08, .48, .46) +
  draw_plot(Fig.S8_mnd_mconti_func_d_all, 0.5, 0.08, .48, .46) +
  draw_plot_label(c("a", "b", "c", "d"),
                  c(0.015+0.01, 0.495+0.01,
                    0.015+0.01, 0.495+0.01),
                  c(0.97, 0.97, 0.55, 0.55)
                  ,size = 18
  )
emf('results/Figures_ages1_35_top50_equal_interval_bh_partialb/Fig.S8.emf',
    width = 20, height = 20, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.S8
dev.off() #turn off device and finalize file


#### Relationships among mean nearest differences ####
##### predictive curve for mnnd ~ mntd_all #####
lincombs.data.mnnd.mntd_all = data.frame(mntd_all=seq(min(dat_suc_sp$mntd_all),
                                                         max(dat_suc_sp$mntd_all),
                                                         length=100),
                                           mnconti_func_d_all=mean(dat_suc_sp$mnconti_func_d_all))

lincombs.matrix.mnnd.mntd_all=model.matrix(~mntd_all,
                                             data=lincombs.data.mnnd.mntd_all)
lincombs.matrix.mnnd.mntd_all=as.data.frame(lincombs.matrix.mnnd.mntd_all)
lincombs.mnnd.mntd_all=inla.make.lincombs(lincombs.matrix.mnnd.mntd_all)

inla.model_lincombs.mnnd.mntd_all = pglmm(mnnd ~ mntd_all +#(1|species) + 
                                              (1|f_p) + (1|field), data = dat_suc_sp,
                                            family = "gaussian", cov_ranef = list(species = tree),
                                            bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                        config = TRUE),
                                                                 quantiles=c(0.025,0.5,0.975),
                                                                 lincomb=lincombs.mnnd.mntd_all,
                                                                 control.predictor=list(compute=T)),
                                            bayes = T)

inla.model_lincombs.mnnd.mntd_all$inla.model$summary.fixed[c(1,3,5)]%>%round(6)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mnnd.mntd_all$predicted.value=inla.model_lincombs.mnnd.mntd_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mnnd.mntd_all$lower=inla.model_lincombs.mnnd.mntd_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mnnd.mntd_all$upper=inla.model_lincombs.mnnd.mntd_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mnnd.mntd_all

save(lincombs.data.mnnd.mntd_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnnd.mntd_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnnd.mntd_all.rdata')
#mntd_all   0.00043    0.00034    0.00052
(Fig.S9_mnnd_mntd_all = ggplot(data=lincombs.data.mnnd.mntd_all, aes(x=mntd_all, y=predicted.value))+
    geom_point(data=dat_suc_sp, aes(x=mntd_all, y=mnnd),
               color = 'grey',
               alpha = 0.5,
               position=position_jitter(height=0.001))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='darkgray',alpha=0.2)+
    geom_line(color='grey1', linetype = 1, linewidth = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #scale_y_continuous(breaks = seq(0, 0.8, length.out = 5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    labs(x = expression('NPD'),
         y = expression('ND'[nearest]))+
    annotate(geom="text",
             x=c((min(dat_suc_sp$mntd_all)+max(dat_suc_sp$mntd_all))/2,
                 (min(dat_suc_sp$mntd_all)+max(dat_suc_sp$mntd_all))/2),
             y=c(-1.1,-1.2),
             label=c("italic(β)['NPD'] == '-0.000013'",
                     "'95%CI' == '[-0.000021, -0.000004]'"),
             parse=T,size=4)+
    theme_regular()
)
# mntd_all -0.000013  -0.000021  -0.000004

##### predictive curve for mnnd ~ mnconti_func_d_all #####
lincombs.data.mnnd.mnconti_func_d_all = data.frame(mnconti_func_d_all=seq(min(dat_suc_sp$mnconti_func_d_all),
                                                                             max(dat_suc_sp$mnconti_func_d_all),
                                                                             length=100),
                                                     mntd_all=mean(dat_suc_sp$mntd_all))

lincombs.matrix.mnnd.mnconti_func_d_all=model.matrix(~mnconti_func_d_all,
                                                       data=lincombs.data.mnnd.mnconti_func_d_all)
lincombs.matrix.mnnd.mnconti_func_d_all=as.data.frame(lincombs.matrix.mnnd.mnconti_func_d_all)
lincombs.mnnd.mnconti_func_d_all=inla.make.lincombs(lincombs.matrix.mnnd.mnconti_func_d_all)

inla.model_lincombs.mnnd.mnconti_func_d_all = pglmm(mnnd ~ mnconti_func_d_all+#(1|species) + 
                                                        (1|f_p) + (1|field), data = dat_suc_sp,
                                                      family = "gaussian", cov_ranef = list(species = tree),
                                                      bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                                  config = TRUE),
                                                                           quantiles=c(0.025,0.5,0.975),
                                                                           lincomb=lincombs.mnnd.mnconti_func_d_all,
                                                                           control.predictor=list(compute=T)),
                                                      bayes = T)

inla.model_lincombs.mnnd.mnconti_func_d_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mnnd.mnconti_func_d_all$predicted.value=inla.model_lincombs.mnnd.mnconti_func_d_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mnnd.mnconti_func_d_all$lower=inla.model_lincombs.mnnd.mnconti_func_d_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mnnd.mnconti_func_d_all$upper=inla.model_lincombs.mnnd.mnconti_func_d_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mnnd.mnconti_func_d_all

save(lincombs.data.mnnd.mnconti_func_d_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnnd.mnconti_func_d_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnnd.mnconti_func_d_all.rdata')
(Fig.S9_mnnd_mnconti_func_d_all = ggplot(data=lincombs.data.mnnd.mnconti_func_d_all, 
                                           aes(x=mnconti_func_d_all, y=predicted.value))+
    geom_point(data=dat_suc_sp, aes(x=mnconti_func_d_all, y=mnnd),
               #shape=1,
               alpha = 0.5,
               color = 'grey',
               position=position_jitter(height=0.001))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='darkgray',alpha=0.2)+
    geom_line(color='grey1', linetype = 2, linewidth = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #scale_y_continuous(breaks = seq(0, 0.8, length.out = 5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    
    labs(x= expression('NFD'),
         y= ' ')+
    annotate(geom="text",
             x=c((min(dat_suc_sp$mnconti_func_d_all)+max(dat_suc_sp$mnconti_func_d_all))/2,
                 (min(dat_suc_sp$mnconti_func_d_all)+max(dat_suc_sp$mnconti_func_d_all))/2),
             y=c(-1.1,-1.2),
             label=c("italic(β)['NFD'] == '0.00182'",
                     "'95%CI' == '[-0.00039, 0.00403]'"),
             parse=T,size=4)+
    theme_regular()
)
#mnconti_func_d_all  0.00182   -0.00039    0.00403

##### predictive curve for mnlgfd ~ mntd_all #####
lincombs.data.mnlgfd.mntd_all = data.frame(mntd_all=seq(min(dat_suc_sp$mntd_all),
                                                           max(dat_suc_sp$mntd_all),
                                                           length=100),
                                           mnconti_func_d_all=mean(dat_suc_sp$mnconti_func_d_all))

lincombs.matrix.mnlgfd.mntd_all=model.matrix(~mntd_all,
                                               data=lincombs.data.mnlgfd.mntd_all)
lincombs.matrix.mnlgfd.mntd_all=as.data.frame(lincombs.matrix.mnlgfd.mntd_all)
lincombs.mnlgfd.mntd_all=inla.make.lincombs(lincombs.matrix.mnlgfd.mntd_all)

inla.model_lincombs.mnlgfd.mntd_all = pglmm(mnlgfd ~ mntd_all+#(1|species) + 
                                                (1|f_p) + (1|field), data = dat_suc_sp,
                                              family = "gaussian", cov_ranef = list(species = tree),
                                              bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                          config = TRUE),
                                                                   quantiles=c(0.025,0.5,0.975),
                                                                   lincomb=lincombs.mnlgfd.mntd_all,
                                                                   control.predictor=list(compute=T)),
                                              bayes = T)


inla.model_lincombs.mnlgfd.mntd_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mnlgfd.mntd_all$predicted.value=inla.model_lincombs.mnlgfd.mntd_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mnlgfd.mntd_all$lower=inla.model_lincombs.mnlgfd.mntd_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mnlgfd.mntd_all$upper=inla.model_lincombs.mnlgfd.mntd_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mnlgfd.mntd_all

save(lincombs.data.mnlgfd.mntd_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnlgfd.mntd_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnlgfd.mntd_all.rdata')
(Fig.S9_mnlgfd_mntd_all = ggplot(data=lincombs.data.mnlgfd.mntd_all, aes(x=mntd_all, y=predicted.value))+
    geom_point(data=dat_suc_sp, aes(x=mntd_all, y=mnlgfd),
               alpha = 0.5,
               color = 'grey',
               position=position_jitter(height=0.001))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='darkgray',alpha=0.2)+
    geom_line(color='grey1', linetype = 1, linewidth = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    scale_y_continuous(breaks = seq(-2.5, 2, 0.5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    labs(x= ' ',
         y=expression('RFD'[nearest]))+
    annotate(geom="text",
             x=c((min(dat_suc_sp$mntd_all)+max(dat_suc_sp$mntd_all))/2,
                 (min(dat_suc_sp$mntd_all)+max(dat_suc_sp$mntd_all))/2),
             y=c(-2.6,-2.9),
             label=c("italic(β)['NPD'] == '-0.00008'",
                     "'95%CI' == '[-0.00014, -0.00001]'"),
             parse=T,size=4)+
    theme_regular()    
)
# mntd_all -0.00003   -0.00008    0.00002

##### predictive curve for mnlgfd ~ mnconti_func_d_all #####
lincombs.data.mnlgfd.mnconti_func_d_all = data.frame(mnconti_func_d_all=seq(
  min(dat_suc_sp$mnconti_func_d_all),
  max(dat_suc_sp$mnconti_func_d_all),
  length=100),
  mntd_all=mean(dat_suc_sp$mntd_all))

lincombs.matrix.mnlgfd.mnconti_func_d_all=model.matrix(~mnconti_func_d_all,
                                                         data=lincombs.data.mnlgfd.mnconti_func_d_all)
lincombs.matrix.mnlgfd.mnconti_func_d_all=as.data.frame(lincombs.matrix.mnlgfd.mnconti_func_d_all)
lincombs.mnlgfd.mnconti_func_d_all=inla.make.lincombs(lincombs.matrix.mnlgfd.mnconti_func_d_all)

inla.model_lincombs.mnlgfd.mnconti_func_d_all = pglmm(mnlgfd ~ mnconti_func_d_all+#(1|species) + 
                                                          (1|f_p) + (1|field), data = dat_suc_sp,
                                                        family = "gaussian", cov_ranef = list(species = tree),
                                                        bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                                    config = TRUE),
                                                                             quantiles=c(0.025,0.5,0.975),
                                                                             lincomb=lincombs.mnlgfd.mnconti_func_d_all,
                                                                             control.predictor=list(compute=T)),
                                                        bayes = T)

inla.model_lincombs.mnlgfd.mnconti_func_d_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mnlgfd.mnconti_func_d_all$predicted.value=inla.model_lincombs.mnlgfd.mnconti_func_d_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mnlgfd.mnconti_func_d_all$lower=inla.model_lincombs.mnlgfd.mnconti_func_d_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mnlgfd.mnconti_func_d_all$upper=inla.model_lincombs.mnlgfd.mnconti_func_d_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mnlgfd.mnconti_func_d_all

save(lincombs.data.mnlgfd.mnconti_func_d_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnlgfd.mnconti_func_d_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnlgfd.mnconti_func_d_all.rdata')
(Fig.S9_mnlgfd_mnconti_func_d_all = ggplot(data=lincombs.data.mnlgfd.mnconti_func_d_all, 
                                             aes(x=mnconti_func_d_all, y=predicted.value))+
    geom_point(data=dat_suc_sp, aes(x=mnconti_func_d_all, y=mnlgfd),
               alpha = 0.5,
               color = 'grey',
               position=position_jitter(height=0.001))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='darkgray',alpha=0.2)+
    geom_line(color='grey1', linetype = 1, linewidth = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    scale_y_continuous(breaks = seq(-2.5, 2, 0.5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    labs(x= ' ',
         y= ' ')+
    annotate(geom="text",
             x=c((min(dat_suc_sp$mnconti_func_d_all)+max(dat_suc_sp$mnconti_func_d_all))/2,
                 (min(dat_suc_sp$mnconti_func_d_all)+max(dat_suc_sp$mnconti_func_d_all))/2),
             y=c(-2.6,-2.9),
             label=c("italic(β)['NFD'] == 0.04343",
                     "'95%CI' == '[0.03121, 0.05565]'"),
             parse=T,size=4)+
    theme_regular()    
)
# mnconti_func_d_all  0.04343    0.03121    0.05565


###### merge plots for Fig. S9 ######
library(cowplot)

Fig.S9 = ggdraw() +
  draw_plot(Fig.S9_mnlgfd_mntd_all, 0.02, 0.5, .48, .46) +
  draw_plot(Fig.S9_mnlgfd_mnconti_func_d_all, 0.5, 0.5, .48, .46) +
  draw_plot(Fig.S9_mnnd_mntd_all, 0.02, 0.08, .48, .46) +
  draw_plot(Fig.S9_mnnd_mnconti_func_d_all, 0.5, 0.08, .48, .46) +
  draw_plot_label(c("a", "b", "c", "d"),
                  c(0.015+0.01, 0.495+0.01,
                    0.015+0.01, 0.495+0.01),
                  c(0.97, 0.97, 0.55, 0.55)
                  ,size = 18
  )

emf('results/Figures_ages1_35_top50_equal_interval_bh_partialb/Fig.S9.emf',
    width = 20, height = 20, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.S9
dev.off() #turn off device and finalize file


