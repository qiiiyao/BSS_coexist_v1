#### Fast start for analyzing mnd ~ mpd, mnd ~ mfunc_d, mfd ~ mpd, mfd ~ mfunc_d ####
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
theme_regular = function(x){
  ggplot2::theme_test() + 
    theme(text = element_text(family = 'Arial'),
          axis.title = element_text(face="bold"),
          axis.title.x = element_text(margin = margin(t = 0.5, r = 0,
                                                      b = 0, l = 0,
                                                      unit = "cm")),
          axis.title.y = element_text(margin = margin(t = 0, r = 0.5,
                                                      b = 0, l = 0,
                                                      unit = "cm")),
          plot.margin = margin(t = 0.2, r = 0.2, b = 0.5, l = 0.2, 
                               unit = "cm"))
}


load('code/results_analyzing/analysing_sameages_top50_data/dat_suc_sp.rdata')

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
lincombs.data.mnd.a.mpd.a_all = data.frame(mpd.a_all=seq(0.0001,max(dat_suc_sp$mpd.a_all),length=100),
                                     mfunc_d.a_all=mean(dat_suc_sp$mfunc_d.a_all))

lincombs.matrix.mnd.a.mpd.a_all=model.matrix(~mpd.a_all+mfunc_d.a_all,
                                       data=lincombs.data.mnd.a.mpd.a_all)
lincombs.matrix.mnd.a.mpd.a_all=as.data.frame(lincombs.matrix.mnd.a.mpd.a_all)
lincombs.mnd.a.mpd.a_all=inla.make.lincombs(lincombs.matrix.mnd.a.mpd.a_all)

inla.model_lincombs.mnd.a.mpd.a_all = pglmm(mnd.a ~ mpd.a_all+mfunc_d.a_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mnd.a.mpd.a_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mnd.a.mpd.a_all.rdata')
(Fig.S7_mnd.a_mpd.a_all = ggplot(data=lincombs.data.mnd.a.mpd.a_all, aes(x=mpd.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='grey1',alpha=0.2)+
    geom_line(color='grey1', linetype = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #scale_y_continuous(breaks = seq(0, 0.8, length.out = 5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mpd.a_all, y=mnd.a),
               shape=1,
               alpha = 0.1,
               position=position_jitter(height=0.001))+
    labs(x = expression(bold('MPD'[ab])),
         y = expression(bold('MND'[ab])))+
    annotate(geom="text",
             x=c((0.0001+max(dat_suc_sp$mpd.a_all))/2,
                 (0.0001+max(dat_suc_sp$mpd.a_all))/2),
             y=c(-0.1,-0.2),
             label=c("italic(β)['MPD'[ab]] == 0.00031", "'95%CI' == '[0.00012, 0.00051]'"),
             parse=T,size=3.5)+
    theme_regular()
  )

##### predictive curve for mnd.a ~ mfunc_d.a_all #####
lincombs.data.mnd.a.mfunc_d.a_all = data.frame(mfunc_d.a_all=seq(0.0001,
                                                                 max(dat_suc_sp$mfunc_d.a_all),
                                                                 length=100),
                                               mpd.a_all=mean(dat_suc_sp$mpd.a_all))

lincombs.matrix.mnd.a.mfunc_d.a_all=model.matrix(~mpd.a_all+mfunc_d.a_all,
                                             data=lincombs.data.mnd.a.mfunc_d.a_all)
lincombs.matrix.mnd.a.mfunc_d.a_all=as.data.frame(lincombs.matrix.mnd.a.mfunc_d.a_all)
lincombs.mnd.a.mfunc_d.a_all=inla.make.lincombs(lincombs.matrix.mnd.a.mfunc_d.a_all)

inla.model_lincombs.mnd.a.mfunc_d.a_all = pglmm(mnd.a ~ mpd.a_all+mfunc_d.a_all+(1|species) + 
                                              (1|f_p) + (1|field), data = dat_suc_sp,
                                            family = "gaussian", cov_ranef = list(species = tree),
                                            bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                        config = TRUE),
                                                                 quantiles=c(0.025,0.5,0.975),
                                                                 lincomb=lincombs.mnd.a.mfunc_d.a_all,
                                                                 control.predictor=list(compute=T)),
                                            bayes = T)

inla.model_lincombs.mnd.a.mfunc_d.a_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mnd.a.mfunc_d.a_all$predicted.value=inla.model_lincombs.mnd.a.mfunc_d.a_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mnd.a.mfunc_d.a_all$lower=inla.model_lincombs.mnd.a.mfunc_d.a_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mnd.a.mfunc_d.a_all$upper=inla.model_lincombs.mnd.a.mfunc_d.a_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mnd.a.mfunc_d.a_all

save(lincombs.data.mnd.a.mfunc_d.a_all, 
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mnd.a.mfunc_d.a_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mnd.a.mfunc_d.a_all.rdata')
(Fig.S7_mnd.a_mfunc_d.a_all = ggplot(data=lincombs.data.mnd.a.mfunc_d.a_all, 
                                     aes(x=mfunc_d.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='grey1',alpha=0.2)+
    geom_line(color='grey1', linetype = 1)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #scale_y_continuous(breaks = seq(0, 0.8, length.out = 5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mfunc_d.a_all, y=mnd.a),
               shape=1,
               alpha = 0.1,
               position=position_jitter(height=0.001))+
    labs(x= expression(bold('MFD'[ab])),
         y= ' ')+
    annotate(geom="text",
             x=c((0.0001+max(dat_suc_sp$mfunc_d.a_all))/2,
                 (0.0001+max(dat_suc_sp$mfunc_d.a_all))/2),
             y=c(-0.1,-0.2),
             label=c("italic(β)['MFD'[ab]] == -0.08748", "'95%CI' == '[-0.15936, -0.01571]'"),
             parse=T,size=3.5)+
    theme_regular()
  )


##### predictive curve for mlgfd.a ~ mpd.a_all #####
lincombs.data.mlgfd.a.mpd.a_all = data.frame(mpd.a_all=seq(0.0001,max(dat_suc_sp$mpd.a_all),length=100),
                                           mfunc_d.a_all=mean(dat_suc_sp$mfunc_d.a_all))

lincombs.matrix.mlgfd.a.mpd.a_all=model.matrix(~mpd.a_all+mfunc_d.a_all,
                                             data=lincombs.data.mlgfd.a.mpd.a_all)
lincombs.matrix.mlgfd.a.mpd.a_all=as.data.frame(lincombs.matrix.mlgfd.a.mpd.a_all)
lincombs.mlgfd.a.mpd.a_all=inla.make.lincombs(lincombs.matrix.mlgfd.a.mpd.a_all)

inla.model_lincombs.mlgfd.a.mpd.a_all = pglmm(mlgfd.a ~ mpd.a_all+mfunc_d.a_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mlgfd.a.mpd.a_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mlgfd.a.mpd.a_all.rdata')
(Fig.S7_mlgfd.a_mpd.a_all = ggplot(data=lincombs.data.mlgfd.a.mpd.a_all, aes(x=mpd.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='grey1',alpha=0.2)+
    geom_line(color='grey1', linetype = 1)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #scale_y_continuous(breaks = seq(0, 0.8, length.out = 5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mpd.a_all, y=mlgfd.a),
               shape=1,
               alpha = 0.1,
               position=position_jitter(height=0.001))+
    labs(x= ' ',
         y=expression(bold('MRFD'[ab])))+
    annotate(geom="text",
             x=c((0.0001+max(dat_suc_sp$mpd.a_all))/2,
                 (0.0001+max(dat_suc_sp$mpd.a_all))/2),
             y=c(0.5,0.4),
             label=c("italic(β)['MPD'[ab]] == -0.00027", "'95%CI' == '[-0.00045, -0.00008]'"),
             parse=T,size=3.5)+
    theme_regular()    
  )

##### predictive curve for mlgfd.a ~ mfunc_d.a_all #####
lincombs.data.mlgfd.a.mfunc_d.a_all = data.frame(mfunc_d.a_all=seq(0.0001,
                                                                 max(dat_suc_sp$mfunc_d.a_all),
                                                                 length=100),
                                               mpd.a_all=mean(dat_suc_sp$mpd.a_all))

lincombs.matrix.mlgfd.a.mfunc_d.a_all=model.matrix(~mpd.a_all+mfunc_d.a_all,
                                                 data=lincombs.data.mlgfd.a.mfunc_d.a_all)
lincombs.matrix.mlgfd.a.mfunc_d.a_all=as.data.frame(lincombs.matrix.mlgfd.a.mfunc_d.a_all)
lincombs.mlgfd.a.mfunc_d.a_all=inla.make.lincombs(lincombs.matrix.mlgfd.a.mfunc_d.a_all)

inla.model_lincombs.mlgfd.a.mfunc_d.a_all = pglmm(mlgfd.a ~ mpd.a_all+mfunc_d.a_all+(1|species) + 
                                                  (1|f_p) + (1|field), data = dat_suc_sp,
                                                family = "gaussian", cov_ranef = list(species = tree),
                                                bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                            config = TRUE),
                                                                     quantiles=c(0.025,0.5,0.975),
                                                                     lincomb=lincombs.mlgfd.a.mfunc_d.a_all,
                                                                     control.predictor=list(compute=T)),
                                                bayes = T)

inla.model_lincombs.mlgfd.a.mfunc_d.a_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mlgfd.a.mfunc_d.a_all$predicted.value=inla.model_lincombs.mlgfd.a.mfunc_d.a_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mlgfd.a.mfunc_d.a_all$lower=inla.model_lincombs.mlgfd.a.mfunc_d.a_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mlgfd.a.mfunc_d.a_all$upper=inla.model_lincombs.mlgfd.a.mfunc_d.a_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mlgfd.a.mfunc_d.a_all

save(lincombs.data.mlgfd.a.mfunc_d.a_all, 
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mlgfd.a.mfunc_d.a_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mlgfd.a.mfunc_d.a_all.rdata')
(Fig.S7_mlgfd.a_mfunc_d.a_all = ggplot(data=lincombs.data.mlgfd.a.mfunc_d.a_all, 
                                     aes(x=mfunc_d.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='grey1',alpha=0.2)+
    geom_line(color='grey1', linetype = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #scale_y_continuous(breaks = seq(0, 0.8, length.out = 5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mfunc_d.a_all, y=mlgfd.a),
               shape=1,
               alpha = 0.1,
               position=position_jitter(height=0.001))+
    labs(x= ' ',
         y= ' ')+
    annotate(geom="text",
             x=c((0.0001+max(dat_suc_sp$mfunc_d.a_all))/2,
                 (0.0001+max(dat_suc_sp$mfunc_d.a_all))/2),
             y=c(0.5,0.4),
             label=c("italic(β)['MFD'[ab]] == 0.03010", "'95%CI' == '[-0.03481, 0.09500]'"),
             parse=T,size=3.5)+
    theme_regular()    
    
    )

###### merge plots for Fig. S7 ######
library(cowplot)

Fig.S7 = ggdraw() +
  draw_plot(Fig.S7_mlgfd.a_mpd.a_all, 0.02, 0.5, .48, .46) +
  draw_plot(Fig.S7_mlgfd.a_mfunc_d.a_all, 0.5, 0.5, .48, .46) +
  draw_plot(Fig.S7_mnd.a_mpd.a_all, 0.02, 0.08, .48, .46) +
  draw_plot(Fig.S7_mnd.a_mfunc_d.a_all, 0.5, 0.08, .48, .46) +
  draw_plot_label(c("(a)", "(b)", "(c)", "(d)"),
                  c(0.015, 0.495, 0.015, 0.495),
                  c(0.97, 0.97, 0.55, 0.55)
                  ,size = 12
  )
emf('results/figures_sameages_top50/Fig.S7.emf',
    width = 20, height = 20, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.S7
dev.off() #turn off device and finalize file

#### Relationships among mean differences ####
##### predictive curve for mnd ~ mpd_all #####
lincombs.data.mnd.mpd_all = data.frame(mpd_all=seq(0.0001,max(dat_suc_sp$mpd_all),length=100),
                                     mfunc_d_all=mean(dat_suc_sp$mfunc_d_all))

lincombs.matrix.mnd.mpd_all=model.matrix(~mpd_all+mfunc_d_all,
                                       data=lincombs.data.mnd.mpd_all)
lincombs.matrix.mnd.mpd_all=as.data.frame(lincombs.matrix.mnd.mpd_all)
lincombs.mnd.mpd_all=inla.make.lincombs(lincombs.matrix.mnd.mpd_all)

inla.model_lincombs.mnd.mpd_all = pglmm(mnd ~ mpd_all+mfunc_d_all+(1|species) + 
                                        (1|f_p) + (1|field), data = dat_suc_sp,
                                      family = "gaussian", cov_ranef = list(species = tree),
                                      bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                  config = TRUE),
                                                           quantiles=c(0.025,0.5,0.975),
                                                           lincomb=lincombs.mnd.mpd_all,
                                                           control.predictor=list(compute=T)),
                                      bayes = T)


inla.model_lincombs.mnd.mpd_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mnd.mpd_all$predicted.value=inla.model_lincombs.mnd.mpd_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mnd.mpd_all$lower=inla.model_lincombs.mnd.mpd_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mnd.mpd_all$upper=inla.model_lincombs.mnd.mpd_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mnd.mpd_all

save(lincombs.data.mnd.mpd_all, 
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mnd.mpd_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mnd.mpd_all.rdata')
(Fig.S8_mnd_mpd_all = ggplot(data=lincombs.data.mnd.mpd_all, aes(x=mpd_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='grey1',alpha=0.2)+
    geom_line(color='grey1', linetype = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #scale_y_continuous(breaks = seq(0, 0.8, length.out = 5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mpd_all, y=mnd),
               shape=1,
               alpha = 0.1,
               position=position_jitter(height=0.001))+
    labs(x= 'MPD',
         y='MND')+
    annotate(geom="text",
             x=c((0.0001+max(dat_suc_sp$mpd_all))/2,
                 (0.0001+max(dat_suc_sp$mpd_all))/2),
             y=c(0.0,-0.1),
             label=c("italic(β)['MPD'] == 0.00021", "'95%CI' == '[-0.00016, 0.00057]'"),
             parse=T,size=3.5)+
    theme_regular()
  )

##### predictive curve for mnd ~ mfunc_d_all #####
lincombs.data.mnd.mfunc_d_all = data.frame(mfunc_d_all=seq(0.0001,
                                                                 max(dat_suc_sp$mfunc_d_all),
                                                                 length=100),
                                               mpd_all=mean(dat_suc_sp$mpd_all))

lincombs.matrix.mnd.mfunc_d_all=model.matrix(~mpd_all+mfunc_d_all,
                                             data=lincombs.data.mnd.mfunc_d_all)
lincombs.matrix.mnd.mfunc_d_all=as.data.frame(lincombs.matrix.mnd.mfunc_d_all)
lincombs.mnd.mfunc_d_all=inla.make.lincombs(lincombs.matrix.mnd.mfunc_d_all)

inla.model_lincombs.mnd.mfunc_d_all = pglmm(mnd ~ mpd_all+mfunc_d_all+(1|species) + 
                                              (1|f_p) + (1|field), data = dat_suc_sp,
                                            family = "gaussian", cov_ranef = list(species = tree),
                                            bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                        config = TRUE),
                                                                 quantiles=c(0.025,0.5,0.975),
                                                                 lincomb=lincombs.mnd.mfunc_d_all,
                                                                 control.predictor=list(compute=T)),
                                            bayes = T)

inla.model_lincombs.mnd.mfunc_d_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mnd.mfunc_d_all$predicted.value=inla.model_lincombs.mnd.mfunc_d_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mnd.mfunc_d_all$lower=inla.model_lincombs.mnd.mfunc_d_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mnd.mfunc_d_all$upper=inla.model_lincombs.mnd.mfunc_d_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mnd.mfunc_d_all

save(lincombs.data.mnd.mfunc_d_all, 
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mnd.mfunc_d_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mnd.mfunc_d_all.rdata')
(Fig.S8_mnd_mfunc_d_all = ggplot(data=lincombs.data.mnd.mfunc_d_all, 
                                     aes(x=mfunc_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='grey1',alpha=0.2)+
    geom_line(color='grey1', linetype = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #scale_y_continuous(breaks = seq(0, 0.8, length.out = 5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mfunc_d_all, y=mnd),
               shape=1,
               alpha = 0.1,
               position=position_jitter(height=0.001))+
    labs(x= 'MFD',
         y=' ')+
    annotate(geom="text",
             x=c((0.0001+max(dat_suc_sp$mfunc_d_all))/2,
                 (0.0001+max(dat_suc_sp$mfunc_d_all))/2),
             y=c(-0.1,-0.2),
             label=c("italic(β)['MFD'] == -0.08068", "'95%CI' == '[-0.23304, 0.06976]'"),
             parse=T,size=3.5)+
    theme_regular()
  )


##### predictive curve for mlgfd ~ mpd_all #####
lincombs.data.mlgfd.mpd_all = data.frame(mpd_all=seq(0.0001,max(dat_suc_sp$mpd_all),length=100),
                                           mfunc_d_all=mean(dat_suc_sp$mfunc_d_all))

lincombs.matrix.mlgfd.mpd_all=model.matrix(~mpd_all+mfunc_d_all,
                                             data=lincombs.data.mlgfd.mpd_all)
lincombs.matrix.mlgfd.mpd_all=as.data.frame(lincombs.matrix.mlgfd.mpd_all)
lincombs.mlgfd.mpd_all=inla.make.lincombs(lincombs.matrix.mlgfd.mpd_all)

inla.model_lincombs.mlgfd.mpd_all = pglmm(mlgfd ~ mpd_all+mfunc_d_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mlgfd.mpd_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mlgfd.mpd_all.rdata')
(Fig.S8_mlgfd_mpd_all = ggplot(data=lincombs.data.mlgfd.mpd_all, aes(x=mpd_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='grey1',alpha=0.2)+
    geom_line(color='grey1', linetype = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #scale_y_continuous(breaks = seq(0, 0.8, length.out = 5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mpd_all, y=mlgfd),
               shape=1,
               alpha = 0.1,
               position=position_jitter(height=0.001))+
    labs(x= ' ',
         y='MRFD')+
    annotate(geom="text",
             x=c((0.0001+max(dat_suc_sp$mpd_all))/2,
                 (0.0001+max(dat_suc_sp$mpd_all))/2),
             y=c(0.5,0.4),
             label=c("italic(β)['MPD'] == -0.00027", "'95%CI' == '[-0.00068, 0.00013]'"),
             parse=T,size=3.5)+
    theme_regular()
  )

##### predictive curve for mlgfd ~ mfunc_d_all #####
lincombs.data.mlgfd.mfunc_d_all = data.frame(mfunc_d_all=seq(0.0001,
                                                                 max(dat_suc_sp$mfunc_d_all),
                                                                 length=100),
                                               mpd_all=mean(dat_suc_sp$mpd_all))

lincombs.matrix.mlgfd.mfunc_d_all=model.matrix(~mpd_all+mfunc_d_all,
                                                 data=lincombs.data.mlgfd.mfunc_d_all)
lincombs.matrix.mlgfd.mfunc_d_all=as.data.frame(lincombs.matrix.mlgfd.mfunc_d_all)
lincombs.mlgfd.mfunc_d_all=inla.make.lincombs(lincombs.matrix.mlgfd.mfunc_d_all)

inla.model_lincombs.mlgfd.mfunc_d_all = pglmm(mlgfd ~ mpd_all+mfunc_d_all+(1|species) + 
                                                  (1|f_p) + (1|field), data = dat_suc_sp,
                                                family = "gaussian", cov_ranef = list(species = tree),
                                                bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                            config = TRUE),
                                                                     quantiles=c(0.025,0.5,0.975),
                                                                     lincomb=lincombs.mlgfd.mfunc_d_all,
                                                                     control.predictor=list(compute=T)),
                                                bayes = T)

inla.model_lincombs.mlgfd.mfunc_d_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mlgfd.mfunc_d_all$predicted.value=inla.model_lincombs.mlgfd.mfunc_d_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mlgfd.mfunc_d_all$lower=inla.model_lincombs.mlgfd.mfunc_d_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mlgfd.mfunc_d_all$upper=inla.model_lincombs.mlgfd.mfunc_d_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mlgfd.mfunc_d_all

save(lincombs.data.mlgfd.mfunc_d_all, 
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mlgfd.mfunc_d_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mlgfd.mfunc_d_all.rdata')
(Fig.S8_mlgfd_mfunc_d_all = ggplot(data=lincombs.data.mlgfd.mfunc_d_all, 
                                     aes(x=mfunc_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='grey1',alpha=0.2)+
    geom_line(color='grey1', linetype = 1)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #scale_y_continuous(breaks = seq(0, 0.8, length.out = 5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mfunc_d_all, y=mlgfd),
               shape=1,
               alpha = 0.1,
               position=position_jitter(height=0.001))+
    labs(x= ' ',
         y=' ')+
    annotate(geom="text",
             x=c((0.0001+max(dat_suc_sp$mfunc_d_all))/2,
                 (0.0001+max(dat_suc_sp$mfunc_d_all))/2),
             y=c(0.5,0.4),
             label=c("italic(β)['MFD'] == 0.25460", "'95%CI' == '[0.09231, 0.41621]'"),
             parse=T,size=3.5)+
    theme_regular()
  )

###### merge plots for Fig. S8 ######
library(cowplot)

Fig.S8 = ggdraw() +
  draw_plot(Fig.S8_mlgfd_mpd_all, 0.02, 0.5, .48, .46) +
  draw_plot(Fig.S8_mlgfd_mfunc_d_all, 0.5, 0.5, .48, .46) +
  draw_plot(Fig.S8_mnd_mpd_all, 0.02, 0.08, .48, .46) +
  draw_plot(Fig.S8_mnd_mfunc_d_all, 0.5, 0.08, .48, .46) +
  draw_plot_label(c("(a)", "(b)", "(c)", "(d)"),
                  c(0.015, 0.495, 0.015, 0.495),
                  c(0.97, 0.97, 0.55, 0.55)
                  ,size = 12
  )
emf('results/figures_sameages_top50/Fig.S8.emf',
    width = 20, height = 20, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.S8
dev.off() #turn off device and finalize file

#### Relationships among mean nearest differences ####
##### predictive curve for mnnd ~ mntd_all #####
lincombs.data.mnnd.mntd_all = data.frame(mntd_all=seq(0.0001,max(dat_suc_sp$mntd_all),length=100),
                                       mnfunc_d_all=mean(dat_suc_sp$mnfunc_d_all))

lincombs.matrix.mnnd.mntd_all=model.matrix(~mntd_all+mnfunc_d_all,
                                         data=lincombs.data.mnnd.mntd_all)
lincombs.matrix.mnnd.mntd_all=as.data.frame(lincombs.matrix.mnnd.mntd_all)
lincombs.mnnd.mntd_all=inla.make.lincombs(lincombs.matrix.mnnd.mntd_all)

inla.model_lincombs.mnnd.mntd_all = pglmm(mnnd ~ mntd_all+mnfunc_d_all+(1|species) + 
                                          (1|f_p) + (1|field), data = dat_suc_sp,
                                        family = "gaussian", cov_ranef = list(species = tree),
                                        bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                    config = TRUE),
                                                             quantiles=c(0.025,0.5,0.975),
                                                             lincomb=lincombs.mnnd.mntd_all,
                                                             control.predictor=list(compute=T)),
                                        bayes = T)


inla.model_lincombs.mnnd.mntd_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mnnd.mntd_all$predicted.value=inla.model_lincombs.mnnd.mntd_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mnnd.mntd_all$lower=inla.model_lincombs.mnnd.mntd_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mnnd.mntd_all$upper=inla.model_lincombs.mnnd.mntd_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mnnd.mntd_all

save(lincombs.data.mnnd.mntd_all, 
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mnnd.mntd_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mnnd.mntd_all.rdata')
(Fig.S9_mnnd_mntd_all = ggplot(data=lincombs.data.mnnd.mntd_all, aes(x=mntd_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='grey1',alpha=0.2)+
    geom_line(color='grey1', linetype = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #scale_y_continuous(breaks = seq(0, 0.8, length.out = 5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mntd_all, y=mnnd),
               shape=1,
               alpha = 0.1,
               position=position_jitter(height=0.001))+
    labs(x= 'MNTD',
         y='MNND')+
    annotate(geom="text",
             x=c((0.0001+max(dat_suc_sp$mntd_all))/2,
                 (0.0001+max(dat_suc_sp$mntd_all))/2),
             y=c(-1,-1.2),
             label=c("italic(β)['MNTD'] == '-0.00008'", "'95%CI' == '[-0.00019, 0.00002]'"),
             parse=T,size=3.5)+
    theme_regular()
    )

##### predictive curve for mnnd ~ mnfunc_d_all #####
lincombs.data.mnnd.mnfunc_d_all = data.frame(mnfunc_d_all=seq(0.0001,
                                                           max(dat_suc_sp$mnfunc_d_all),
                                                           length=100),
                                           mntd_all=mean(dat_suc_sp$mntd_all))

lincombs.matrix.mnnd.mnfunc_d_all=model.matrix(~mntd_all+mnfunc_d_all,
                                             data=lincombs.data.mnnd.mnfunc_d_all)
lincombs.matrix.mnnd.mnfunc_d_all=as.data.frame(lincombs.matrix.mnnd.mnfunc_d_all)
lincombs.mnnd.mnfunc_d_all=inla.make.lincombs(lincombs.matrix.mnnd.mnfunc_d_all)

inla.model_lincombs.mnnd.mnfunc_d_all = pglmm(mnnd ~ mntd_all+mnfunc_d_all+(1|species) + 
                                              (1|f_p) + (1|field), data = dat_suc_sp,
                                            family = "gaussian", cov_ranef = list(species = tree),
                                            bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                        config = TRUE),
                                                                 quantiles=c(0.025,0.5,0.975),
                                                                 lincomb=lincombs.mnnd.mnfunc_d_all,
                                                                 control.predictor=list(compute=T)),
                                            bayes = T)

inla.model_lincombs.mnnd.mnfunc_d_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mnnd.mnfunc_d_all$predicted.value=inla.model_lincombs.mnnd.mnfunc_d_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mnnd.mnfunc_d_all$lower=inla.model_lincombs.mnnd.mnfunc_d_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mnnd.mnfunc_d_all$upper=inla.model_lincombs.mnnd.mnfunc_d_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mnnd.mnfunc_d_all

save(lincombs.data.mnnd.mnfunc_d_all, 
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mnnd.mnfunc_d_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mnnd.mnfunc_d_all.rdata')
(Fig.S9_mnnd_mnfunc_d_all = ggplot(data=lincombs.data.mnnd.mnfunc_d_all, 
                                 aes(x=mnfunc_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='grey1',alpha=0.2)+
    geom_line(color='grey1', linetype = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #scale_y_continuous(breaks = seq(0, 0.8, length.out = 5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnfunc_d_all, y=mnnd),
               shape=1,
               alpha = 0.1,
               position=position_jitter(height=0.001))+
    labs(x= 'MNFD',
         y=' ')+
    annotate(geom="text",
             x=c((0.0001+max(dat_suc_sp$mnfunc_d_all))/2,
                 (0.0001+max(dat_suc_sp$mnfunc_d_all))/2),
             y=c(-1,-1.2),
             label=c("italic(β)['MNFD'] == -0.10503", "'95%CI' == '[-0.25055, 0.04019]'"),
             parse=T,size=3.5)+
    theme_regular()
  )


##### predictive curve for mnlgfd ~ mntd_all #####
lincombs.data.mnlgfd.mntd_all = data.frame(mntd_all=seq(0.0001,max(dat_suc_sp$mntd_all),length=100),
                                         mnfunc_d_all=mean(dat_suc_sp$mnfunc_d_all))

lincombs.matrix.mnlgfd.mntd_all=model.matrix(~mntd_all+mnfunc_d_all,
                                           data=lincombs.data.mnlgfd.mntd_all)
lincombs.matrix.mnlgfd.mntd_all=as.data.frame(lincombs.matrix.mnlgfd.mntd_all)
lincombs.mnlgfd.mntd_all=inla.make.lincombs(lincombs.matrix.mnlgfd.mntd_all)

inla.model_lincombs.mnlgfd.mntd_all = pglmm(mnlgfd ~ mntd_all+mnfunc_d_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mnlgfd.mntd_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mnlgfd.mntd_all.rdata')
(Fig.S9_mnlgfd_mntd_all = ggplot(data=lincombs.data.mnlgfd.mntd_all, aes(x=mntd_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='grey1',alpha=0.2)+
    geom_line(color='grey1', linetype = 1)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #scale_y_continuous(breaks = seq(0, 0.8, length.out = 5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mntd_all, y=mnlgfd),
               shape=1,
               alpha = 0.1,
               position=position_jitter(height=0.001))+
    labs(x= ' ',
         y='MNRFD')+
    annotate(geom="text",
             x=c((0.0001+max(dat_suc_sp$mntd_all))/2,
                 (0.0001+max(dat_suc_sp$mntd_all))/2),
             y=c(0.5,0.4),
             label=c("italic(β)['MNTD'] == '0.00009'", "'95%CI' == '[0.00003, 0.00015]'"),
             parse=T,size=3.5)+
    theme_regular()
  )

##### predictive curve for mnlgfd ~ mnfunc_d_all #####
lincombs.data.mnlgfd.mnfunc_d_all = data.frame(mnfunc_d_all=seq(0.0001,
                                                             max(dat_suc_sp$mnfunc_d_all),
                                                             length=100),
                                             mntd_all=mean(dat_suc_sp$mntd_all))

lincombs.matrix.mnlgfd.mnfunc_d_all=model.matrix(~mntd_all+mnfunc_d_all,
                                               data=lincombs.data.mnlgfd.mnfunc_d_all)
lincombs.matrix.mnlgfd.mnfunc_d_all=as.data.frame(lincombs.matrix.mnlgfd.mnfunc_d_all)
lincombs.mnlgfd.mnfunc_d_all=inla.make.lincombs(lincombs.matrix.mnlgfd.mnfunc_d_all)

inla.model_lincombs.mnlgfd.mnfunc_d_all = pglmm(mnlgfd ~ mntd_all+mnfunc_d_all+(1|species) + 
                                                (1|f_p) + (1|field), data = dat_suc_sp,
                                              family = "gaussian", cov_ranef = list(species = tree),
                                              bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                          config = TRUE),
                                                                   quantiles=c(0.025,0.5,0.975),
                                                                   lincomb=lincombs.mnlgfd.mnfunc_d_all,
                                                                   control.predictor=list(compute=T)),
                                              bayes = T)

inla.model_lincombs.mnlgfd.mnfunc_d_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mnlgfd.mnfunc_d_all$predicted.value=inla.model_lincombs.mnlgfd.mnfunc_d_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mnlgfd.mnfunc_d_all$lower=inla.model_lincombs.mnlgfd.mnfunc_d_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mnlgfd.mnfunc_d_all$upper=inla.model_lincombs.mnlgfd.mnfunc_d_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mnlgfd.mnfunc_d_all

save(lincombs.data.mnlgfd.mnfunc_d_all, 
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mnlgfd.mnfunc_d_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_sameages_top50_data/lincombs.data.mnlgfd.mnfunc_d_all.rdata')
(Fig.S9_mnlgfd_mnfunc_d_all = ggplot(data=lincombs.data.mnlgfd.mnfunc_d_all, 
                                   aes(x=mnfunc_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='grey1',alpha=0.2)+
    geom_line(color='grey1', linetype = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #scale_y_continuous(breaks = seq(0, 0.8, length.out = 5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnfunc_d_all, y=mnlgfd),
               shape=1,
               alpha = 0.1,
               position=position_jitter(height=0.001))+
    labs(x= ' ',
         y=' ')+
    annotate(geom="text",
             x=c((0.0001+max(dat_suc_sp$mnfunc_d_all))/2,
                 (0.0001+max(dat_suc_sp$mnfunc_d_all))/2),
             y=c(0.5,0.4),
             label=c("italic(β)['MNFD'] == '0.04997'", "'95%CI' == '[-0.04118, 0.14106]'"),
             parse=T,size=3.5)+
    theme_regular()
  )


###### merge plots for Fig. S9 ######
library(cowplot)

Fig.S9 = ggdraw() +
  draw_plot(Fig.S9_mnlgfd_mntd_all, 0.02, 0.5, .48, .46) +
  draw_plot(Fig.S9_mnlgfd_mnfunc_d_all, 0.5, 0.5, .48, .46) +
  draw_plot(Fig.S9_mnnd_mntd_all, 0.02, 0.08, .48, .46) +
  draw_plot(Fig.S9_mnnd_mnfunc_d_all, 0.5, 0.08, .48, .46) +
  draw_plot_label(c("(a)", "(b)", "(c)", "(d)"),
                  c(0.015, 0.495, 0.015, 0.495),
                  c(0.97, 0.97, 0.55, 0.55)
                  ,size = 12
  )
emf('results/figures_sameages_top50/Fig.S9.emf',
    width = 20, height = 20, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.S9
dev.off() #turn off device and finalize file


