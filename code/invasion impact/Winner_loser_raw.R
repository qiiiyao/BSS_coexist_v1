# R code for the analyses in: 
# Author: Helge Bruelheide, 23.05.2022
# with contributions of Francesco Maria Sabatini 19.05.2022

library(reshape2)    # for casting matrices
library(data.table)  # for aggregating data
#library(ggplot2)     # for graphics
#library(DescTools)   # Gini coefficient and Lorenz curve
#library(BSDA)        # sign test
#library(Hmisc)       # error bars in Fig. 4
library(vegan)       # Shannon diversity  
library(matrixStats) # to calculate rowRanks
#library(lmerTest)    # for mixed effects models (lmer)
#library(parameters)  # for credible intervals for lmer


### Step. 1: Load data
##### Data loading
load('code/data preparation/transformed data/sp_cover_f1.rdata')
load('code/data preparation/transformed data/sp_cover_f2.rdata')
load('code/data preparation/transformed data/sp_cover_f3.rdata')

trait = read.csv('data/original data/traits332.csv',
                 header = T)

### Step. 2: Calculation of yearly change of real time data 
species.change.fp.realtime_raw = data.frame(FIELD=NULL, PLOT=NULL,
                                        from.n=NULL, to.n=NULL,from=NULL, to=NULL,species=NULL,
                                        absolute.change=NULL, relative.change=NULL, relative.rank.change=NULL,
                                        log.repsonse.change=NULL, absolute.change.colonizer=NULL,
                                        absolute.change.extinct=NULL)

# loop for all fields
for (i in c(1:length(unique(sp_cover_f1$Field)))){
  #i = 1
  print(i)
  # select the species cover data from a given project
  species.long = sp_cover_f1[sp_cover_f1$Field == i,]
  field = unique(species.long$Field)
  species.long = arrange(species.long, species.long$Plot,
                         species.long$fake_age)
  # turn the long format into a RELEVE_NR by TaxonName matrix
  species = species.long[,10:ncol(species.long)]
  
  # RS_PLOT holds a unique (within the site) code of the resurveyed plot; 
  # it is used to pair observations from different times recorded in 
  # the same plot; gives a unique identifier for the resurveyed plot or 
  # set of plots in time if combined with RS_PROJECT (see RS_PROJECT_PLOT). 
  # Several plots in the same year might have the same RS_PLOT code if they have 
  # to be summarised for temporal comparisons.
  plot.list = sort(unique(species.long$Plot))
  # produce a presence/absence version of the RELEVE_NR by TaxonName matrix
  species_pa = species.long[,10:ncol(species.long)]
  species_pa[species_pa > 0] = 1
  
  # loop for all RS_PLOT
  for (j in 1:length(plot.list)){
    #j = 1
    target.year.field = sort(unique(species.long$fake_age))
    if (length(target.year.field)>1){
      # Check whether the plot series actually has more than one data
      # this might happen when a previous plot was recorded several times
      # but no new record was made
      
      # produce a subset of the plot by species matrix
      species.target.plot = species.long[species.long$Plot == plot.list[j],10:ncol(species.long)]
      target.year = sort(unique(species.long[species.long$Plot == plot.list[j],])$fake_age)
      # Remove empty cols, i.e. species that do not occur in the subset
      species.target.plot = species.target.plot[,colSums(species.target.plot)!=0, drop=F]
      # produce a presence/absence version of the species.target.plot matrix
      species.target.plot.pa = species.target.plot 
      species.target.plot.pa[species.target.plot.pa>0] = 1
      
      # take the mean of replicated plots
      species.target.plot = aggregate(species.target.plot, by=list(target.year),FUN=mean)
      rownames(species.target.plot) = species.target.plot$Group.1
      species.target.plot = species.target.plot[,-1, drop=F] 
      # remove the first column from the aggregation function
      species.target.plot = species.target.plot[order(as.numeric(rownames(species.target.plot)),decreasing=T),,drop=F]
      # reorder the plot records from the latest year being the first row
      
      # calculate differences in rank abundance curves
      # according to Avolio, M.L., Carroll, I.T., Collins, S.L., Houseman, G.R., 
      # Hallett, L.M., Isbell, F., Koerner, S.E., Komatsu, K.J., Smith, M.D., 
      # Wilcox, K.R., 2019. A comprehensive approach to analyzing community 
      # dynamics using rank abundance curves 10: e02881. 10.1002/ecs2.2881
      species.target.plot2 = sweep(species.target.plot, 1,
                                   apply(species.target.plot,1,FUN=sum),
                                   FUN="/")
      species.target.plot2 = species.target.plot2[order(as.numeric(rownames(species.target.plot2)),
                                                        decreasing=T),
                                                  ,drop=F]
      # species.target.plot2 holds relative abundance values for each plot
      # the species are ranked per row using the matrixStats package
      species.target.plot.ranks2 = as.matrix(rowRanks(as.matrix(-species.target.plot2,
                                                                ties.method="average")))
      species.target.plot.ranks2[which(rowSums(species.target.plot) == 0),] = ncol(species.target.plot)
      dimnames(species.target.plot.ranks2)[[1]] = row.names(species.target.plot2)
      dimnames(species.target.plot.ranks2)[[2]] = colnames(species.target.plot2)
      # turn ranks into relative ranks
      species.target.plot.ranks2 = sweep(species.target.plot.ranks2,
                                         1,
                                         apply(species.target.plot.ranks2,1,FUN=max), FUN="/")
      
      # calculate delta curve (curve.diff) according to Avolio et al. (2019)
      unique.relative.ranks = c(0,unique(sort(species.target.plot.ranks2)))
      species.target.plot.ranks3 = apply(species.target.plot.ranks2,2,
                                         FUN=cut,
                                         breaks=unique.relative.ranks, labels=F)
      dimnames(species.target.plot.ranks3)[[1]] = rownames(species.target.plot2)
      curve.diff = 0
      y1sum = 0
      y2sum = 0
      for (k in 2: length(unique.relative.ranks)){
        y1 = species.target.plot2[1,species.target.plot.ranks3[1,]==k-1]
        y2 = species.target.plot2[2,species.target.plot.ranks3[2,]==k-1]
        y1sum = y1sum + ifelse(length(y1)>0,sum(y1, na.rm = T),0)
        y2sum = y2sum + ifelse(length(y2)>0,sum(y2, na.rm = T),0)
        curve.diff = curve.diff+y1sum-y2sum
      }
      curve.diff
      
      # calculate rank change (rank.change) according to Avolio et al. (2019)
      species.target.plot.ranks4 = species.target.plot.ranks2[1,]-species.target.plot.ranks2[2,]
      rank.change = sum(abs(species.target.plot.ranks4))/ncol(species.target.plot.ranks2)
      
      # calculate rank change separately for pos. and neg. change (rank.change.sign)
      rank.change.neg = mean(species.target.plot.ranks4[species.target.plot.ranks4<0])
      rank.change.pos = mean(species.target.plot.ranks4[species.target.plot.ranks4>0])
      
      rank.change.neg.sum = sum(species.target.plot.ranks4[species.target.plot.ranks4<0])
      rank.change.pos.sum = sum(species.target.plot.ranks4[species.target.plot.ranks4>0])
      
      # Number of years in the time series
      n = table(target.year)[order(as.numeric(names(table(target.year))),decreasing=T)]
      # Compare subsequent changes in a time series
      for (k in 1:(length(n)-1)){
        # diff.absolute.change is the difference in absolute cover values
        diff.absolute.change = species.target.plot[k,,drop=F] - species.target.plot[k+1,,drop=F]
        #diff.absolute.change[species.target.plot[k,]==0 & species.target.plot[k+1,]==0] = NA
        # diff.relative.change is the difference in absolute cover values
        diff.relative.change = species.target.plot2[k,,drop=F] - species.target.plot2[k+1,,drop=F]
        #diff.relative.change[species.target.plot2[k,]==0 & species.target.plot2[k+1,]==0] = NA
        # diff.relative.rank.change is the difference in relative ranks
        diff.relative.rank.change = species.target.plot.ranks2[k,] - species.target.plot.ranks2[k+1,]
        #diff.relative.rank.change[species.target.plot.ranks2[k,]==0 & species.target.plot.ranks2[k+1,]==0] = NA
        # diff.absolute.change.colonizer is diff.absolute.change only for 
        # new colonizers in the interval
        diff.absolute.change.colonizer = diff.absolute.change
        # only keep records that are new and did not occur in Year 1
        diff.absolute.change.colonizer[species.target.plot[k+1,]!=0] = NA
        # diff.absolute.change.extinct is diff.absolute.change only for 
        # species that went extinct in the interval
        diff.absolute.change.extinct = diff.absolute.change
        # only keep records that went extinct and did not occur in Year 2
        diff.absolute.change.extinct[species.target.plot[k,]!=0] = NA
        # collect all metrics for resurvey ID x species x time interval combinations
        
        ## FMS: Switching to LIST increases the speed of this loop substantially (since it avoids R overwriting a vector with 10^5 rows at every epoch)
        species.change.fp.realtime_raw = rbind(species.change.fp.realtime_raw,
                                           data.frame(
                                             FIELD=field,
                                             PLOT=plot.list[j],
                                             from.n=as.numeric(n[k+1]), to.n=as.numeric(n[k]), from=as.numeric(names(n)[k+1]),
                                             to=as.numeric(names(n)[k]),
                                             species=names(diff.relative.change),
                                             absolute.change=as.numeric(diff.absolute.change),
                                             relative.change=as.numeric(diff.relative.change),
                                             relative.rank.change=as.numeric(diff.relative.rank.change),
                                             absolute.change.colonizer=as.numeric(diff.absolute.change.colonizer),
                                             absolute.change.extinct=as.numeric(diff.absolute.change.extinct)))
        
      }
    }
  }
}


# remove empty rows from the resurvey ID x species x time interval combinations
str(species.change.fp.realtime_raw) #579449 obs. of  13 variables:
#species.change.fp.realtime_raw2 = species.change.fp.realtime_raw[!is.na(species.change.fp.realtime_raw$absolute.change),]
# add some variables for analysis
species.change.fp.realtime_raw$F_P = paste(species.change.fp.realtime_raw$FIELD,
                                       species.change.fp.realtime_raw$PLOT,
                                       sep = '_')
species.change.fp.realtime_raw = species.change.fp.realtime_raw %>% relocate(F_P, .after = PLOT)

save(species.change.fp.realtime_raw,
     file = 'code/invasion impact/species.change.fp.realtime_raw.rdata')


### Step. 3: Real time: Calculation of change at the species level at all time and moving windows

load('code/invasion impact/species.change.fp.realtime_raw.rdata')

## All time
# making use of the data.table function to aggregate changes in cover 

species.change.fp.realtime_raw2 = data.table(species.change.fp.realtime_raw)
change.alltime.fp.realtime_raw = species.change.fp.realtime_raw2[,{n_pos = length(absolute.change[absolute.change>0]);
n_neg = length(absolute.change[absolute.change<0]);
n_all = ifelse(n_pos+n_neg>0,n_pos+n_neg,1);
est = binom.test(n_pos,n_all);
list(n=length(absolute.change),
     pos=length(absolute.change[absolute.change>0]),
     equal=length(absolute.change[absolute.change==0]),
     neg=length(absolute.change[absolute.change<0]),
     est.binom=est$estimate,
     conf.binom.minus=est$conf.int[1],
     conf.binom.plus=est$conf.int[2],
     p.values.binom=est$p.value,
     mean.absolute.change=mean(absolute.change, na.rm=T),
     mean.relative.change=mean(relative.change, na.rm=T))}, by=.(F_P, species)]
change.alltime.fp.realtime_raw = change.alltime.fp.realtime_raw %>% 
  filter(species %in% trait[trait$Origin == 'Native',]$Species)
change.alltime.fp.realtime_raw$gain_loss = NA
change.alltime.fp.realtime_raw[change.alltime.fp.realtime_raw$mean.relative.change < 0,]$gain_loss = 0
change.alltime.fp.realtime_raw[change.alltime.fp.realtime_raw$mean.relative.change > 0,]$gain_loss = 1
change.alltime.fp.realtime_raw$f_p_species = paste(change.alltime.fp.realtime_raw$F_P,
                                                   change.alltime.fp.realtime_raw$species,
                                               sep = '_')
change.alltime.fp.realtime_raw = change.alltime.fp.realtime_raw %>% relocate(f_p_species,
                                                                     .after = species)
save(change.alltime.fp.realtime_raw,
     file = 'code/invasion impact/change_nati_alltime_raw.rdata')

#### Abundance change of native species at different time moving windows
fake_age_f2 = sort(unique(sp_cover_f2$fake_age))

change.moving.fp.realtime.l = foreach(i = 1:14) %dopar% {
  #i = 1
  seq_t = seq(fake_age_f2[i], fake_age_f2[i+12],1)
  species.change.fp.realtime_raw2_t = species.change.fp.realtime_raw2 %>% filter(
    to %in% seq_t
  ) 
  species.change.fp.realtime_raw2_t = data.table(species.change.fp.realtime_raw2_t)
  
  change.moving.fp.t = species.change.fp.realtime_raw2_t[,{n_pos = length(absolute.change[absolute.change>0]);
  n_neg = length(absolute.change[absolute.change<0]);
  n_all = ifelse(n_pos+n_neg>0,n_pos+n_neg,1);
  est = binom.test(n_pos,n_all);
  list(n=length(absolute.change),
       pos=length(absolute.change[absolute.change>0]),
       equal=length(absolute.change[absolute.change==0]),
       neg=length(absolute.change[absolute.change<0]),
       est.binom=est$estimate,
       conf.binom.minus=est$conf.int[1],
       conf.binom.plus=est$conf.int[2],
       p.values.binom=est$p.value,
       mean.absolute.change=mean(absolute.change, na.rm=T),
       mean.relative.change=mean(relative.change, na.rm=T))}, by=.(F_P, species)]
  change.moving.fp.t = change.moving.fp.t %>%
    filter(species %in% trait[trait$Origin == 'Native',]$Species)
  change.moving.fp.t$gain_loss = NA
  change.moving.fp.t[change.moving.fp.t$mean.relative.change < 0,]$gain_loss = 0
  change.moving.fp.t[change.moving.fp.t$mean.relative.change > 0,]$gain_loss = 1
  return(change.moving.fp.t)
  
}
save(change.moving.fp.realtime.l,
     file = 'code/invasion impact/change.moving.fp.realtime.l.rdata')


### Step. 4: Calculation of yearly change of fitted time data 
species.change.fp.fittedtime = data.frame(FIELD=NULL, PLOT=NULL,
                                          from.n=NULL, to.n=NULL,from=NULL, to=NULL,species=NULL,
                                          absolute.change=NULL, relative.change=NULL, relative.rank.change=NULL,
                                          log.repsonse.change=NULL, absolute.change.colonizer=NULL,
                                          absolute.change.extinct=NULL)

# loop for all fields
for (i in c(1:length(unique(sp_cover_f3$Field)))){
  #i = 1
  print(i)
  # select the species cover data from a given project
  species.long = sp_cover_f3[sp_cover_f3$Field == i,]
  field = unique(species.long$Field)
  species.long = arrange(species.long, species.long$Plot,
                         species.long$fake_age)
  # turn the long format into a RELEVE_NR by TaxonName matrix
  species = species.long[,10:ncol(species.long)]
  
  # RS_PLOT holds a unique (within the site) code of the resurveyed plot; 
  # it is used to pair observations from different times recorded in 
  # the same plot; gives a unique identifier for the resurveyed plot or 
  # set of plots in time if combined with RS_PROJECT (see RS_PROJECT_PLOT). 
  # Several plots in the same year might have the same RS_PLOT code if they have 
  # to be summarised for temporal comparisons.
  plot.list = sort(unique(species.long$Plot))
  # produce a presence/absence version of the RELEVE_NR by TaxonName matrix
  species_pa = species.long[,10:ncol(species.long)]
  species_pa[species_pa > 0] = 1
  
  # loop for all RS_PLOT
  for (j in 1:length(plot.list)){
    #j = 1
    target.year.field = sort(unique(species.long$fake_age))
    if (length(target.year.field)>1){
      # Check whether the plot series actually has more than one data
      # this might happen when a previous plot was recorded several times
      # but no new record was made
      
      # produce a subset of the plot by species matrix
      species.target.plot = species.long[species.long$Plot == plot.list[j],10:ncol(species.long)]
      target.year = sort(unique(species.long[species.long$Plot == plot.list[j],])$fake_age)
      # Remove empty cols, i.e. species that do not occur in the subset
      species.target.plot = species.target.plot[,colSums(species.target.plot)!=0, drop=F]
      # produce a presence/absence version of the species.target.plot matrix
      species.target.plot.pa = species.target.plot 
      species.target.plot.pa[species.target.plot.pa>0] = 1
      
      # take the mean of replicated plots
      species.target.plot = aggregate(species.target.plot, by=list(target.year),FUN=mean)
      rownames(species.target.plot) = species.target.plot$Group.1
      species.target.plot = species.target.plot[,-1, drop=F] 
      # remove the first column from the aggregation function
      species.target.plot = species.target.plot[order(as.numeric(rownames(species.target.plot)),decreasing=T),,drop=F]
      # reorder the plot records from the latest year being the first row
      
      # calculate differences in rank abundance curves
      # according to Avolio, M.L., Carroll, I.T., Collins, S.L., Houseman, G.R., 
      # Hallett, L.M., Isbell, F., Koerner, S.E., Komatsu, K.J., Smith, M.D., 
      # Wilcox, K.R., 2019. A comprehensive approach to analyzing community 
      # dynamics using rank abundance curves 10: e02881. 10.1002/ecs2.2881
      species.target.plot2 = sweep(species.target.plot, 1,
                                   apply(species.target.plot,1,FUN=sum),
                                   FUN="/")
      species.target.plot2 = species.target.plot2[order(as.numeric(rownames(species.target.plot2)),
                                                        decreasing=T),
                                                  ,drop=F]
      # species.target.plot2 holds relative abundance values for each plot
      # the species are ranked per row using the matrixStats package
      species.target.plot.ranks2 = as.matrix(rowRanks(as.matrix(-species.target.plot2,
                                                                ties.method="average")))
      species.target.plot.ranks2[which(rowSums(species.target.plot) == 0),] = ncol(species.target.plot)
      dimnames(species.target.plot.ranks2)[[1]] = row.names(species.target.plot2)
      dimnames(species.target.plot.ranks2)[[2]] = colnames(species.target.plot2)
      # turn ranks into relative ranks
      species.target.plot.ranks2 = sweep(species.target.plot.ranks2,
                                         1,
                                         apply(species.target.plot.ranks2,1,FUN=max), FUN="/")
      
      # calculate delta curve (curve.diff) according to Avolio et al. (2019)
      unique.relative.ranks = c(0,unique(sort(species.target.plot.ranks2)))
      species.target.plot.ranks3 = apply(species.target.plot.ranks2,2,FUN=cut,breaks=unique.relative.ranks, labels=F)
      dimnames(species.target.plot.ranks3)[[1]] = rownames(species.target.plot2)
      curve.diff = 0
      y1sum = 0
      y2sum = 0
      for (k in 2: length(unique.relative.ranks)){
        y1 = species.target.plot2[1,species.target.plot.ranks3[1,]==k-1]
        y2 = species.target.plot2[2,species.target.plot.ranks3[2,]==k-1]
        y1sum = y1sum + ifelse(length(y1)>0,sum(y1),0)
        y2sum = y2sum + ifelse(length(y2)>0,sum(y2),0)
        curve.diff = curve.diff+y1sum-y2sum
      }
      curve.diff
      
      # calculate rank change (rank.change) according to Avolio et al. (2019)
      species.target.plot.ranks4 = species.target.plot.ranks2[1,]-species.target.plot.ranks2[2,]
      rank.change = sum(abs(species.target.plot.ranks4))/ncol(species.target.plot.ranks2)
      
      # calculate rank change separately for pos. and neg. change (rank.change.sign)
      rank.change.neg = mean(species.target.plot.ranks4[species.target.plot.ranks4<0])
      rank.change.pos = mean(species.target.plot.ranks4[species.target.plot.ranks4>0])
      
      rank.change.neg.sum = sum(species.target.plot.ranks4[species.target.plot.ranks4<0])
      rank.change.pos.sum = sum(species.target.plot.ranks4[species.target.plot.ranks4>0])
      
      # Number of years in the time series
      n = table(target.year)[order(as.numeric(names(table(target.year))),decreasing=T)]
      # Compare subsequent changes in a time series
      for (k in 1:(length(n)-1)){
        # diff.absolute.change is the difference in absolute cover values
        diff.absolute.change = species.target.plot[k,,drop=F] - species.target.plot[k+1,,drop=F]
        #diff.absolute.change[species.target.plot[k,]==0 & species.target.plot[k+1,]==0] = NA
        # diff.relative.change is the difference in absolute cover values
        diff.relative.change = species.target.plot2[k,,drop=F] - species.target.plot2[k+1,,drop=F]
        #diff.relative.change[species.target.plot2[k,]==0 & species.target.plot2[k+1,]==0] = NA
        # diff.relative.rank.change is the difference in relative ranks
        diff.relative.rank.change = species.target.plot.ranks2[k,] - species.target.plot.ranks2[k+1,]
        #diff.relative.rank.change[species.target.plot.ranks2[k,]==0 & species.target.plot.ranks2[k+1,]==0] = NA
        # diff.absolute.change.colonizer is diff.absolute.change only for 
        # new colonizers in the interval
        diff.absolute.change.colonizer = diff.absolute.change
        # only keep records that are new and did not occur in Year 1
        diff.absolute.change.colonizer[species.target.plot[k+1,]!=0] = NA
        # diff.absolute.change.extinct is diff.absolute.change only for 
        # species that went extinct in the interval
        diff.absolute.change.extinct = diff.absolute.change
        # only keep records that went extinct and did not occur in Year 2
        diff.absolute.change.extinct[species.target.plot[k,]!=0] = NA
        # collect all metrics for resurvey ID x species x time interval combinations
        
        ## FMS: Switching to LIST increases the speed of this loop substantially (since it avoids R overwriting a vector with 10^5 rows at every epoch)
        species.change.fp.fittedtime = rbind(species.change.fp.fittedtime,
                                             data.frame(
                                               FIELD=field,
                                               PLOT=plot.list[j],
                                               from.n=as.numeric(n[k+1]), to.n=as.numeric(n[k]), from=as.numeric(names(n)[k+1]),
                                               to=as.numeric(names(n)[k]),
                                               species=names(diff.relative.change),
                                               absolute.change=as.numeric(diff.absolute.change),
                                               relative.change=as.numeric(diff.relative.change),
                                               relative.rank.change=as.numeric(diff.relative.rank.change),
                                               absolute.change.colonizer=as.numeric(diff.absolute.change.colonizer),
                                               absolute.change.extinct=as.numeric(diff.absolute.change.extinct)))
        
      }
    }
  }
}


# remove empty rows from the resurvey ID x species x time interval combinations
str(species.change.fp.fittedtime) #769579 obs. of  12 variables
#species.change.fp.fittedtime2 = species.change.fp.fittedtime[!is.na(species.change.fp.fittedtime$absolute.change),]
# add some variables for analysis
species.change.fp.fittedtime$F_P = paste(species.change.fp.fittedtime$FIELD,
                                         species.change.fp.fittedtime$PLOT,
                                         sep = '_')

save(species.change.fp.fittedtime,
     file = 'code/invasion impact/species.change.fp.fittedtime.rdata')

### Step. 5: Fitted time: Calculation of change at the species level at all time and moving windows

load('code/invasion impact/species.change.fp.fittedtime.rdata')

## All time
# making use of the data.table function to aggregate changes in cover 
species.change.fp.fittedtime2 = data.table(species.change.fp.fittedtime)
change.alltime.fp.fittedtime = species.change.fp.fittedtime2[,{n_pos = length(absolute.change[absolute.change>0]);
n_neg = length(absolute.change[absolute.change<0]);
n_all = ifelse(n_pos+n_neg>0,n_pos+n_neg,1);
est = binom.test(n_pos,n_all);
list(n=length(absolute.change),
     pos=length(absolute.change[absolute.change>0]),
     equal=length(absolute.change[absolute.change==0]),
     neg=length(absolute.change[absolute.change<0]),
     est.binom=est$estimate,
     conf.binom.minus=est$conf.int[1],
     conf.binom.plus=est$conf.int[2],
     p.values.binom=est$p.value,
     mean.absolute.change=mean(absolute.change, na.rm=T),
     mean.relative.change=mean(relative.change, na.rm=T))}, by=.(F_P, species)]
change.alltime.fp.fittedtime = change.alltime.fp.fittedtime %>% 
  filter(species %in% trait[trait$Origin == 'Native',]$Species)
change.alltime.fp.fittedtime$gain_loss = NA
change.alltime.fp.fittedtime[change.alltime.fp.fittedtime$mean.relative.change < 0,]$gain_loss = 0
change.alltime.fp.fittedtime[change.alltime.fp.fittedtime$mean.relative.change > 0,]$gain_loss = 1
#save(change.alltime.fp.fittedtime, file = 'code/invasion impact/change_nati_alltime.rdata')

#### Merge data from real time and fitted time
change.alltime.fp.fittedtime$f_p_species = paste(change.alltime.fp.fittedtime$F_P,
                                                 change.alltime.fp.fittedtime$species,
                                                 sep = '_')
change.alltime.fp.fittedtime = change.alltime.fp.fittedtime %>% relocate(f_p_species,
                                                                         .after = species)
change.alltime.fp.realtime$f_p_species = paste(change.alltime.fp.realtime$F_P,
                                               change.alltime.fp.realtime$species,
                                               sep = '_')
change.alltime.fp.realtime = change.alltime.fp.realtime %>% relocate(f_p_species,
                                                                     .after = species)
change.alltime.fp.real.fitted = change.alltime.fp.realtime[,c(1:3,12:14)] %>%
  left_join(change.alltime.fp.fittedtime[,c(1:3,12:14)],
            by = join_by(f_p_species))
colnames(change.alltime.fp.real.fitted)[c(4:6, 9:11)] = c("mean.absolute.change.realtime",
                                                          "mean.relative.change.realtime",
                                                          "gain_loss.realtime",
                                                          "mean.absolute.change.fittedtime",
                                                          "mean.relative.change.fittedtime",
                                                          "gain_loss.fittedtime")
save(change.alltime.fp.real.fitted,
     file = 'code/invasion impact/change.alltime.fp.real.fitted.rdata')
#### Abundance change of native species at different time moving windows
fake_age_f2 = sort(unique(sp_cover_f2$fake_age))

change.moving.fp.fittedtime.l = foreach(i = 1:14) %dopar% {
  #i = 1
  seq_t = seq(fake_age_f2[i], fake_age_f2[i+12],1)
  species.change.fp.fittedtime2_t = species.change.fp.fittedtime2 %>% filter(
    to %in% seq_t
  ) 
  species.change.fp.fittedtime2_t = data.table(species.change.fp.fittedtime2_t)
  
  change.moving.fp.t = species.change.fp.fittedtime2_t[,{n_pos = length(absolute.change[absolute.change>0]);
  n_neg = length(absolute.change[absolute.change<0]);
  n_all = ifelse(n_pos+n_neg>0,n_pos+n_neg,1);
  est = binom.test(n_pos,n_all);
  list(n=length(absolute.change),
       pos=length(absolute.change[absolute.change>0]),
       equal=length(absolute.change[absolute.change==0]),
       neg=length(absolute.change[absolute.change<0]),
       est.binom=est$estimate,
       conf.binom.minus=est$conf.int[1],
       conf.binom.plus=est$conf.int[2],
       p.values.binom=est$p.value,
       mean.absolute.change=mean(absolute.change, na.rm=T),
       mean.relative.change=mean(relative.change, na.rm=T))}, by=.(F_P, species)]
  change.moving.fp.t = change.moving.fp.t %>%
    filter(species %in% trait[trait$Origin == 'Native',]$Species)
  change.moving.fp.t$gain_loss = NA
  change.moving.fp.t[change.moving.fp.t$mean.relative.change < 0,]$gain_loss = 0
  change.moving.fp.t[change.moving.fp.t$mean.relative.change > 0,]$gain_loss = 1
  return(change.moving.fp.t)
  
}
save(change.moving.fp.fittedtime.l,
     file = 'code/invasion impact/change.moving.fp.fittedtime.l.rdata')


