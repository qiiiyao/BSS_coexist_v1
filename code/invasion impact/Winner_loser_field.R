# R code for the analyses in: 
# Author: Helge Bruelheide, 23.05.2022
# with contributions of Francesco Maria Sabatini 19.05.2022
rm(list = ls())
library(reshape2)    # for casting matrices
library(data.table)  # for aggregating data
library(dplyr)
#library(ggplot2)     # for graphics
#library(DescTools)   # Gini coefficient and Lorenz curve
#library(BSDA)        # sign test
#library(Hmisc)       # error bars in Fig. 4
library(vegan)       # Shannon diversity  
library(matrixStats) # to calculate rowRanks
library(doParallel)
#library(lmerTest)    # for mixed effects models (lmer)
#library(parameters)  # for credible intervals for lmer


### Step. 1: Load data
##### Data loading
load('code/data preparation/transformed data/sp_cover_f1_field.rdata')
load('code/data preparation/transformed data/sp_cover_f2_field.rdata')
load('code/data preparation/transformed data/sp_cover_f3_field.rdata')

# Just retain the fake age from 1 to 51
sp_cover_f1_field = sp_cover_f1_field %>% filter(fake_age %in% seq(1,51,1))

trait = read.csv('data/original data/traits332.csv',
                 header = T)

### Step. 2: Calculation of yearly change of real time data 
species.change.field.realtime = data.frame(FIELD=NULL,
                                          from.n=NULL, to.n=NULL,
                                          from=NULL, to=NULL,species=NULL,
                                          absolute.change=NULL,
                                          relative.change=NULL,
                                          relative.rank.change=NULL,
                                          log.repsonse.change=NULL, absolute.change.colonizer=NULL,
                                          absolute.change.extinct=NULL)

# loop for all fields
for (i in c(1:length(unique(sp_cover_f1_field$Field)))){
  #i = 1
  print(i)
  field.list = sort(unique(sp_cover_f1_field$Field))
  # select the species cover data from a given project
  target.year.field = sort(unique(sp_cover_f1_field$fake_age))
  if (length(target.year.field)>1){
      # Check whether the plot series actually has more than one data
      # this might happen when a previous plot was recorded several times
      # but no new record was made
      
      # produce a subset of the plot by species matrix
      species.target.field = sp_cover_f1_field[sp_cover_f1_field$Field== field.list[i],
                                               8:ncol(sp_cover_f1_field)]
      target.year = sort(unique(sp_cover_f1_field[sp_cover_f1_field$Field == field.list[i],])$fake_age)
      # Remove empty cols, i.e. species that do not occur in the subset
      species.target.field = species.target.field[,colSums(species.target.field)!=0, drop=F]
      # produce a presence/absence version of the species.target.field matrix
      species.target.field.pa = species.target.field 
      species.target.field.pa[species.target.field.pa>0] = 1
      
      # take the mean of replicated plots
      species.target.field = aggregate(species.target.field, by=list(target.year),FUN=mean)
      rownames(species.target.field) = species.target.field$Group.1
      species.target.field = species.target.field[,-1, drop=F] 
      # remove the first column from the aggregation function
      species.target.field = species.target.field[order(as.numeric(rownames(species.target.field)),
                                                        decreasing=T),,drop=F]
      # reorder the plot records from the latest year being the first row
      
      # calculate differences in rank abundance curves
      # according to Avolio, M.L., Carroll, I.T., Collins, S.L., Houseman, G.R., 
      # Hallett, L.M., Isbell, F., Koerner, S.E., Komatsu, K.J., Smith, M.D., 
      # Wilcox, K.R., 2019. A comprehensive approach to analyzing community 
      # dynamics using rank abundance curves 10: e02881. 10.1002/ecs2.2881
      species.target.field2 = sweep(species.target.field, 1,
                                   apply(species.target.field,1,FUN=sum),
                                   FUN="/")
      species.target.field2 = species.target.field2[order(as.numeric(rownames(species.target.field2)),
                                                        decreasing=T),
                                                  ,drop=F]
      # species.target.field2 holds relative abundance values for each plot
      # the species are ranked per row using the matrixStats package
      species.target.field.ranks2 = as.matrix(rowRanks(as.matrix(-species.target.field2,
                                                                ties.method="average")))
      species.target.field.ranks2[which(rowSums(species.target.field) == 0),] = ncol(species.target.field)
      dimnames(species.target.field.ranks2)[[1]] = row.names(species.target.field2)
      dimnames(species.target.field.ranks2)[[2]] = colnames(species.target.field2)
      # turn ranks into relative ranks
      species.target.field.ranks2 = sweep(species.target.field.ranks2,
                                         1,
                                         apply(species.target.field.ranks2,1,FUN=max), FUN="/")
      
      # calculate delta curve (curve.diff) according to Avolio et al. (2019)
      unique.relative.ranks = c(0,unique(sort(species.target.field.ranks2)))
      species.target.field.ranks3 = apply(species.target.field.ranks2,2,FUN=cut,breaks=unique.relative.ranks, labels=F)
      dimnames(species.target.field.ranks3)[[1]] = rownames(species.target.field2)
      curve.diff = 0
      y1sum = 0
      y2sum = 0
      for (k in 2: length(unique.relative.ranks)){
        y1 = species.target.field2[1,species.target.field.ranks3[1,]==k-1]
        y2 = species.target.field2[2,species.target.field.ranks3[2,]==k-1]
        y1sum = y1sum + ifelse(length(y1)>0,sum(y1),0)
        y2sum = y2sum + ifelse(length(y2)>0,sum(y2),0)
        curve.diff = curve.diff+y1sum-y2sum
      }
      curve.diff
      
      # calculate rank change (rank.change) according to Avolio et al. (2019)
      species.target.field.ranks4 = species.target.field.ranks2[1,]-species.target.field.ranks2[2,]
      rank.change = sum(abs(species.target.field.ranks4))/ncol(species.target.field.ranks2)
      
      # calculate rank change separately for pos. and neg. change (rank.change.sign)
      rank.change.neg = mean(species.target.field.ranks4[species.target.field.ranks4<0])
      rank.change.pos = mean(species.target.field.ranks4[species.target.field.ranks4>0])
      
      rank.change.neg.sum = sum(species.target.field.ranks4[species.target.field.ranks4<0])
      rank.change.pos.sum = sum(species.target.field.ranks4[species.target.field.ranks4>0])
      
      # Number of years in the time series
      n = table(target.year)[order(as.numeric(names(table(target.year))),decreasing=T)]
      # Compare subsequent changes in a time series
      for (k in 1:(length(n)-1)){
        # diff.absolute.change is the difference in absolute cover values
        diff.absolute.change = species.target.field[k,,drop=F] - species.target.field[k+1,,drop=F]
        #diff.absolute.change[species.target.field[k,]==0 & species.target.field[k+1,]==0] = NA
        # diff.relative.change is the difference in absolute cover values
        diff.relative.change = species.target.field2[k,,drop=F] - species.target.field2[k+1,,drop=F]
        #diff.relative.change[species.target.field2[k,]==0 & species.target.field2[k+1,]==0] = NA
        # diff.relative.rank.change is the difference in relative ranks
        diff.relative.rank.change = species.target.field.ranks2[k,] - species.target.field.ranks2[k+1,]
        #diff.relative.rank.change[species.target.field.ranks2[k,]==0 & species.target.field.ranks2[k+1,]==0] = NA
        # diff.absolute.change.colonizer is diff.absolute.change only for 
        # new colonizers in the interval
        diff.absolute.change.colonizer = diff.absolute.change
        # only keep records that are new and did not occur in Year 1
        diff.absolute.change.colonizer[species.target.field[k+1,]!=0] = NA
        # diff.absolute.change.extinct is diff.absolute.change only for 
        # species that went extinct in the interval
        diff.absolute.change.extinct = diff.absolute.change
        # only keep records that went extinct and did not occur in Year 2
        diff.absolute.change.extinct[species.target.field[k,]!=0] = NA
        # collect all metrics for resurvey ID x species x time interval combinations
        
        ## FMS: Switching to LIST increases the speed of this loop substantially (since it avoids R overwriting a vector with 10^5 rows at every epoch)
        species.change.field.realtime = rbind(species.change.field.realtime,
                                             data.frame(
                                               FIELD=field.list[i],
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



# remove empty rows from the resurvey ID x species x time interval combinations
str(species.change.field.realtime) #579449 obs. of  13 variables:

save(species.change.field.realtime,
     file = 'code/invasion impact/species.change.field.realtime.rdata')

### Step. 3: Real time: Calculation of change at the species level at all time and moving windows

load('code/invasion impact/species.change.field.realtime.rdata')

## All time
# making use of the data.table function to aggregate changes in cover 

species.change.field.realtime2 = data.table(species.change.field.realtime)
change.alltime.field.realtime = species.change.field.realtime2[,{n_pos = length(absolute.change[absolute.change>0]);
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
     mean.relative.change=mean(relative.change, na.rm=T))}, by=.(FIELD, species)]
change.alltime.field.realtime = change.alltime.field.realtime %>% 
  filter(species %in% trait[trait$Origin == 'Native',]$Species)
change.alltime.field.realtime$gain_loss = NA
change.alltime.field.realtime[change.alltime.field.realtime$mean.relative.change < 0,]$gain_loss = 0
change.alltime.field.realtime[change.alltime.field.realtime$mean.relative.change > 0,]$gain_loss = 1
#save(change.alltime.field.realtime, file = 'code/invasion impact/change_nati_alltime.rdata')

#### Abundance change of native species at different time moving windows
fake_age_f2 = sort(unique(sp_cover_f2_field$fake_age))

change.moving.field.realtime.l = foreach(i = 1:14) %dopar% {
  #i = 1
  seq_t = seq(fake_age_f2[i], fake_age_f2[i+12],1)
  species.change.field.realtime2_t = species.change.field.realtime2 %>% filter(
    to %in% seq_t
  ) 
  species.change.field.realtime2_t = data.table(species.change.field.realtime2_t)
  
  change.moving.field.t = species.change.field.realtime2_t[,{n_pos = length(absolute.change[absolute.change>0]);
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
       mean.relative.change=mean(relative.change, na.rm=T))}, by=.(FIELD, species)]
  change.moving.field.t = change.moving.field.t %>%
    filter(species %in% trait[trait$Origin == 'Native',]$Species)
  change.moving.field.t$gain_loss = NA
  change.moving.field.t[change.moving.field.t$mean.relative.change < 0,]$gain_loss = 0
  change.moving.field.t[change.moving.field.t$mean.relative.change > 0,]$gain_loss = 1
  return(change.moving.field.t)
  
}
save(change.moving.field.realtime.l,
     file = 'code/invasion impact/change.moving.field.realtime.l.rdata')

### Step. 4: Calculation of yearly change of real time data 
species.change.field.fittedtime = data.frame(FIELD=NULL,
                                           from.n=NULL, to.n=NULL,
                                           from=NULL, to=NULL,species=NULL,
                                           absolute.change=NULL,
                                           relative.change=NULL,
                                           relative.rank.change=NULL,
                                           log.repsonse.change=NULL, absolute.change.colonizer=NULL,
                                           absolute.change.extinct=NULL)

# loop for all fields
for (i in c(1:length(unique(sp_cover_f3_field$Field)))){
  #i = 1
  print(i)
  field.list = sort(unique(sp_cover_f3_field$Field))
  # select the species cover data from a given project
  target.year.field = sort(unique(sp_cover_f3_field$fake_age))
  if (length(target.year.field)>1){
    # Check whether the plot series actually has more than one data
    # this might happen when a previous plot was recorded several times
    # but no new record was made
    
    # produce a subset of the plot by species matrix
    species.target.field = sp_cover_f3_field[sp_cover_f3_field$Field== field.list[i],
                                             8:ncol(sp_cover_f3_field)]
    target.year = sort(unique(sp_cover_f3_field[sp_cover_f3_field$Field == field.list[i],])$fake_age)
    # Remove empty cols, i.e. species that do not occur in the subset
    species.target.field = species.target.field[,colSums(species.target.field)!=0, drop=F]
    # produce a presence/absence version of the species.target.field matrix
    species.target.field.pa = species.target.field 
    species.target.field.pa[species.target.field.pa>0] = 1
    
    # take the mean of replicated plots
    species.target.field = aggregate(species.target.field, by=list(target.year),FUN=mean)
    rownames(species.target.field) = species.target.field$Group.1
    species.target.field = species.target.field[,-1, drop=F] 
    # remove the first column from the aggregation function
    species.target.field = species.target.field[order(as.numeric(rownames(species.target.field)),
                                                      decreasing=T),,drop=F]
    # reorder the plot records from the latest year being the first row
    
    # calculate differences in rank abundance curves
    # according to Avolio, M.L., Carroll, I.T., Collins, S.L., Houseman, G.R., 
    # Hallett, L.M., Isbell, F., Koerner, S.E., Komatsu, K.J., Smith, M.D., 
    # Wilcox, K.R., 2019. A comprehensive approach to analyzing community 
    # dynamics using rank abundance curves 10: e02881. 10.1002/ecs2.2881
    species.target.field2 = sweep(species.target.field, 1,
                                  apply(species.target.field,1,FUN=sum),
                                  FUN="/")
    species.target.field2 = species.target.field2[order(as.numeric(rownames(species.target.field2)),
                                                        decreasing=T),
                                                  ,drop=F]
    # species.target.field2 holds relative abundance values for each plot
    # the species are ranked per row using the matrixStats package
    species.target.field.ranks2 = as.matrix(rowRanks(as.matrix(-species.target.field2,
                                                               ties.method="average")))
    species.target.field.ranks2[which(rowSums(species.target.field) == 0),] = ncol(species.target.field)
    dimnames(species.target.field.ranks2)[[1]] = row.names(species.target.field2)
    dimnames(species.target.field.ranks2)[[2]] = colnames(species.target.field2)
    # turn ranks into relative ranks
    species.target.field.ranks2 = sweep(species.target.field.ranks2,
                                        1,
                                        apply(species.target.field.ranks2,1,FUN=max), FUN="/")
    
    # calculate delta curve (curve.diff) according to Avolio et al. (2019)
    unique.relative.ranks = c(0,unique(sort(species.target.field.ranks2)))
    species.target.field.ranks3 = apply(species.target.field.ranks2,2,FUN=cut,breaks=unique.relative.ranks, labels=F)
    dimnames(species.target.field.ranks3)[[1]] = rownames(species.target.field2)
    curve.diff = 0
    y1sum = 0
    y2sum = 0
    for (k in 2: length(unique.relative.ranks)){
      y1 = species.target.field2[1,species.target.field.ranks3[1,]==k-1]
      y2 = species.target.field2[2,species.target.field.ranks3[2,]==k-1]
      y1sum = y1sum + ifelse(length(y1)>0,sum(y1),0)
      y2sum = y2sum + ifelse(length(y2)>0,sum(y2),0)
      curve.diff = curve.diff+y1sum-y2sum
    }
    curve.diff
    
    # calculate rank change (rank.change) according to Avolio et al. (2019)
    species.target.field.ranks4 = species.target.field.ranks2[1,]-species.target.field.ranks2[2,]
    rank.change = sum(abs(species.target.field.ranks4))/ncol(species.target.field.ranks2)
    
    # calculate rank change separately for pos. and neg. change (rank.change.sign)
    rank.change.neg = mean(species.target.field.ranks4[species.target.field.ranks4<0])
    rank.change.pos = mean(species.target.field.ranks4[species.target.field.ranks4>0])
    
    rank.change.neg.sum = sum(species.target.field.ranks4[species.target.field.ranks4<0])
    rank.change.pos.sum = sum(species.target.field.ranks4[species.target.field.ranks4>0])
    
    # Number of years in the time series
    n = table(target.year)[order(as.numeric(names(table(target.year))),decreasing=T)]
    # Compare subsequent changes in a time series
    for (k in 1:(length(n)-1)){
      # diff.absolute.change is the difference in absolute cover values
      diff.absolute.change = species.target.field[k,,drop=F] - species.target.field[k+1,,drop=F]
      #diff.absolute.change[species.target.field[k,]==0 & species.target.field[k+1,]==0] = NA
      # diff.relative.change is the difference in absolute cover values
      diff.relative.change = species.target.field2[k,,drop=F] - species.target.field2[k+1,,drop=F]
      #diff.relative.change[species.target.field2[k,]==0 & species.target.field2[k+1,]==0] = NA
      # diff.relative.rank.change is the difference in relative ranks
      diff.relative.rank.change = species.target.field.ranks2[k,] - species.target.field.ranks2[k+1,]
      #diff.relative.rank.change[species.target.field.ranks2[k,]==0 & species.target.field.ranks2[k+1,]==0] = NA
      # diff.absolute.change.colonizer is diff.absolute.change only for 
      # new colonizers in the interval
      diff.absolute.change.colonizer = diff.absolute.change
      # only keep records that are new and did not occur in Year 1
      diff.absolute.change.colonizer[species.target.field[k+1,]!=0] = NA
      # diff.absolute.change.extinct is diff.absolute.change only for 
      # species that went extinct in the interval
      diff.absolute.change.extinct = diff.absolute.change
      # only keep records that went extinct and did not occur in Year 2
      diff.absolute.change.extinct[species.target.field[k,]!=0] = NA
      # collect all metrics for resurvey ID x species x time interval combinations
      
      ## FMS: Switching to LIST increases the speed of this loop substantially (since it avoids R overwriting a vector with 10^5 rows at every epoch)
      species.change.field.fittedtime = rbind(species.change.field.fittedtime,
                                            data.frame(
                                              FIELD=field.list[i],
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



# remove empty rows from the resurvey ID x species x time interval combinations
str(species.change.field.fittedtime) #579449 obs. of  13 variables:

save(species.change.field.fittedtime,
     file = 'code/invasion impact/species.change.field.fittedtime.rdata')

### Step. 5: Real time: Calculation of change at the species level at all time and moving windows

load('code/invasion impact/species.change.field.fittedtime.rdata')

## All time
# making use of the data.table function to aggregate changes in cover 

species.change.field.fittedtime2 = data.table(species.change.field.fittedtime)
change.alltime.field.fittedtime = species.change.field.fittedtime2[,{n_pos = length(absolute.change[absolute.change>0]);
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
     mean.relative.change=mean(relative.change, na.rm=T))}, by=.(FIELD, species)]
change.alltime.field.fittedtime = change.alltime.field.fittedtime %>% 
  filter(species %in% trait[trait$Origin == 'Native',]$Species)
change.alltime.field.fittedtime$gain_loss = NA
change.alltime.field.fittedtime[change.alltime.field.fittedtime$mean.relative.change < 0,]$gain_loss = 0
change.alltime.field.fittedtime[change.alltime.field.fittedtime$mean.relative.change > 0,]$gain_loss = 1
#save(change.alltime.field.fittedtime, file = 'code/invasion impact/change_nati_alltime.rdata')

#### Abundance change of native species at different time moving windows
fake_age_f2 = sort(unique(sp_cover_f2_field$fake_age))

change.moving.field.fittedtime.l = foreach(i = 1:14) %dopar% {
  #i = 1
  seq_t = seq(fake_age_f2[i], fake_age_f2[i+12],1)
  species.change.field.fittedtime2_t = species.change.field.fittedtime2 %>% filter(
    to %in% seq_t
  ) 
  species.change.field.fittedtime2_t = data.table(species.change.field.fittedtime2_t)
  
  change.moving.field.t = species.change.field.fittedtime2_t[,{n_pos = length(absolute.change[absolute.change>0]);
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
       mean.relative.change=mean(relative.change, na.rm=T))}, by=.(FIELD, species)]
  change.moving.field.t = change.moving.field.t %>%
    filter(species %in% trait[trait$Origin == 'Native',]$Species)
  change.moving.field.t$gain_loss = NA
  change.moving.field.t[change.moving.field.t$mean.relative.change < 0,]$gain_loss = 0
  change.moving.field.t[change.moving.field.t$mean.relative.change > 0,]$gain_loss = 1
  return(change.moving.field.t)
  
}
save(change.moving.field.fittedtime.l,
     file = 'code/invasion impact/change.moving.field.fittedtime.l.rdata')

