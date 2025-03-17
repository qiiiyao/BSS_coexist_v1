library(dplyr)

load("code/data preparation/invasion_orginaldata.RData")

cover_ab_all_nati = list()
cover_ab_all_intro = list()
cover_ab_all_estab = list()
cover_ab_all_domin = list()

for (i in 1:length(cover_all_ab)) {
  #i = 1
  trans_filter_nati = re_cover_ab_nati_l[[i]]
  trans_filter_intro = re_cover_ab_intro_l[[i]]
  trans_filter_estab = re_cover_ab_estab_l[[i]]
  trans_filter_domin = re_cover_ab_domin_l[[i]]
  trans_all = cover_all_ab[[i]]
  
  trans_all_intro = trans_all %>% 
                       select(colnames(trans_filter_intro))
  trans_all_estab = trans_all %>%
                       select(colnames(trans_filter_estab))
  trans_all_domin = trans_all %>%
                       select(colnames(trans_filter_domin))
  trans_all_nati = trans_all %>%
                       select(colnames(trans_filter_nati))
  
  cover_ab_all_nati[[i]] = trans_all_nati
  cover_ab_all_intro[[i]] = trans_all_intro
  cover_ab_all_estab[[i]] = trans_all_estab
  cover_ab_all_domin[[i]] = trans_all_domin
    }

### check the relationship between alien species abundance and natives richness
source('D:/R projects/BSS/code/invasion impact/measure_impact_func.R')

#### community level merge estab species to one species
cover_ab_estab_for_impact = list()



for (i in 1:length(cover_ab_all_domin)) {
  i = 51
  alien = cover_ab_all_estab[[i]]
  intro = cover_ab_all_intro[[i]]

  if (ncol(alien) > 7) {
  nati = cover_ab_all_nati[[i]]
  sp_name = c(colnames(nati)[8:ncol(nati)],
              colnames(alien)[8:ncol(alien)])
  plot_all = cbind(nati[,1:ncol(nati)],
                   alien[,8:ncol(alien)])
  colnames(plot_all)[8:ncol(plot_all)] = sp_name
  
  resi_richness = apply(plot_all[,8:ncol(plot_all)], 1, function(x){
  length(x[x>0])
  })

  sp_col = c()
  for (j in 1:(ncol(alien)-7)) {
    j = 1
    focal_alien = alien[,c(2, 7+j)]
    
    
    ab_rich = as.data.frame(cbind(focal_alien, resi_richness))
    colnames(ab_rich)[2] = 'focal_alien'
    ab_rich_no0 = ab_rich %>% filter(focal_alien != 0)
    standard_seq = focal_alien$Absolute_year
    sp_name_focal = colnames(focal_alien)[2]
    ddd = standard_seq %in% ab_rich_no0$Absolute_year
    
    sep_t = get_longest_consec_rep(ddd) ## get the longest consecutive sequence in
    ## in time series data
    
    year = standard_seq[sep_t]
    
    ab_rich_longest = ab_rich_no0[ab_rich_no0$Absolute_year %in% year,]
    
    one_t1 = plot_all %>% 
             filter(Absolute_year %in% year[c(1)]) %>% 
             select(colnames(plot_all)[8:ncol(plot_all)])
    one_t2 = plot_all %>% 
             filter(Absolute_year %in% year[c(length(year)-1)]) %>% 
             select(colnames(plot_all)[8:ncol(plot_all)])
    
    #get list of communities for 2019 as well.
    DSV.change_1 = one_t2 %>%
                   select(sp_name_focal) - one_t1 %>%
                                           select(sp_name_focal)
    DSV.change = data.frame(site.plot = rownames(one_t1),
                           DSV.changes = unlist(DSV.change_1),
                           DSV.t2 = unlist(one_t2 %>%
                             select(sp_name_focal)),
                           DSV.t1 = unlist(one_t1 %>%
                             select(sp_name_focal)))
    out_t1<-list()
    
    for (i in 1:nrow(one_t1)){
      #i = 1
      tmp<-as.numeric(one_t1[i,])
      names(tmp)<-colnames(one_t1)
      tmp<-tmp[names(tmp)!= sp_name_focal]#get rid of VIRO
      tmp<-tmp[tmp>0]
      out_t1[[i]]<-data.frame(Abun=tmp[order(tmp,decreasing=TRUE)],Rank=1:length(tmp))
    }
    
    names(out_t1)<-rownames(one_t1)
    rich_t1<-sapply(out_t1,nrow) 
    DSV.change$Rich._t1<-rich_t1
    
    out_t2<-list()
    
    for (i in 1:nrow(one_t2)){
      #i = 1
      tmp<-as.numeric(one_t2[i,])
      names(tmp)<-colnames(one_t2)
      tmp<-tmp[names(tmp)!= sp_name_focal]#get rid of VIRO
      tmp<-tmp[tmp>0]
      out_t2[[i]]<-data.frame(Abun=tmp[order(tmp,decreasing=TRUE)],Rank=1:length(tmp))
    }
    
    names(out_t2) <- rownames(one_t2)
    rich_t2<-sapply(out_t2,nrow) 
    DSV.change$Rich._t2 <- rich_t2
    
    #focus on extinctions
    
    exts<-NULL
    ranks<-NULL
    
    ## comapre the t and t+1, the extinction species number and mean rank
    for (i in 1:length(out_t1)){
      #i = 1
      a<-out_t1[[i]]
      b<-out_t2[[i]]
      tmp.e<-is.na(match(rownames(a),rownames(b)))
      exts[i]<-sum(tmp.e)
      ranks[i]<-mean(a$Rank[tmp.e])
    }
    ranks[is.nan(ranks)]<-NA
    
    DSV.change$obs.ext<-exts
    DSV.change$obs.rank.ext<-ranks
    
    ##do ses on each community
    
    mod.out<-list()
    
    #calculating only for plots with increasing DSV
    for (i in 1:length(out_t1)){
      tmp<-out_t1[[i]]
      if (DSV.change$DSV.change[i] <= 0) mod.out[[i]]<-NA
      if (DSV.change$DSV.change[i] > 0){
        mod.out[[i]]<-ext.sum.stat(tmp,DSV.change$DSV.change[i], 999,
                                   DSV.change$obs.ext[i],DSV.change$obs.rank.ext[i])
      }
      
    }
    
    names(mod.out)<-names(out_t1)
    
    
    focal_alien = alien[,c(2, 7+j)]


    ab_rich = as.data.frame(cbind(focal_alien, nati_richness))
    colnames(ab_rich)[2] = 'focal_alien'
    ab_rich_no0 = ab_rich %>% filter(focal_alien != 0)
    standard_seq = focal_alien$Absolute_year
      
    ddd = standard_seq %in% ab_rich_no0$Absolute_year

    sep_t = get_longest_consec_rep(ddd) ## get the longest consecutive sequence in
                                        ## in time series data
    
    year = standard_seq[sep_t]
    
    ab_rich_longest = ab_rich_no0[ab_rich_no0$Absolute_year %in% year,]
    
    mod_ab_rich = lm(nati_richness~focal_alien, ab_rich_longest)
    cof = as.data.frame(summary(mod_ab_rich)$coefficients)
    if(!is.na(cof$`Pr(>|t|)`[2]) & cof$`Pr(>|t|)`[2] < 0.1 & cof$Estimate[2] < 0) {
        sp_col[j] = j
          }

  }
        sp_col = sp_col[!is.na(sp_col)]
        if(length(sp_col) == 0) {
        alien_focal = alien[,c(1:7)]
        } else {
        sp_col = 7+sp_col
        alien_focal = alien[,c(1:7, sp_col)]}

        cover_ab_estab_for_impact[[i]] = alien_focal
        
  } else {
        cover_ab_estab_for_impact[[i]] = NULL
    }
  
}




#### community level merge estab species to one species
cover_ab_estab_for_impact = list()

for (i in 1:length(cover_ab_all_domin)) {
  #i = 51
  if (ncol(cover_ab_all_estab[[i]]) > 7) {
    if (ncol(cover_ab_all_estab[[i]]) > 8) {
    alien = rowSums(cover_ab_all_estab[[i]][,8:ncol(cover_ab_all_estab[[i]])])
    } else {
      alien = cover_ab_all_estab[[i]][,8:ncol(cover_ab_all_estab[[i]])] 
    }
    nati = cover_ab_all_nati[[i]]
    f_p_name = unique(nati$f_p)
    
    sp_name = c(colnames(nati)[8:ncol(nati)])
    plot_all = cbind(nati[,1:ncol(nati)], alien)
    colnames(plot_all)[8:ncol(plot_all)] = c(sp_name, 'estab')
    
    nati_richness_1 = apply(plot_all[,8:(ncol(plot_all)-1)], 1, function(x){
      length(x[x>0])
    })
    plot_all$nati_richness = nati_richness_1
    
    plot_all = plot_all %>% filter(nati_richness > 0)
    
      year_t2 = (plot_all %>% filter(estab == max(plot_all$estab)))$Absolute_year
      if (length(year_t2) > 1) {
      year_t2 = sample(year_t2, 1)}
      year_t1 = min(plot_all$Absolute_year)
      
      one_t1 = plot_all %>% 
        filter(Absolute_year %in% year_t1) %>% 
        select(colnames(plot_all)[8:ncol(plot_all)])
      one_t2 = plot_all %>% 
        filter(Absolute_year %in% year_t2) %>% 
        select(colnames(plot_all)[8:ncol(plot_all)])
      
      #get list of communities for 2019 as well.
      DSV.change_1 = one_t2 %>%
        select('estab') - one_t1 %>%
        select('estab')
      DSV.change = data.frame(DSV.changes = unlist(DSV.change_1),
                              DSV.t2 = unlist(one_t2 %>%
                                                select('estab')),
                              DSV.t1 = unlist(one_t1 %>%
                                                select('estab')))
      
      out_t1 = data.frame()
      tmp<-as.numeric(one_t1[1,])
      names(tmp)<-colnames(one_t1)
      tmp<-tmp[names(tmp) != 'estab']#get rid of VIRO
      tmp<-tmp[names(tmp) != 'nati_richness']#get rid of VIRO
      tmp<-tmp[tmp>0]
      out_t1<-data.frame(Abun=tmp[order(tmp,decreasing=TRUE)],Rank=1:length(tmp))
      rich_t1<-nrow(out_t1)
      DSV.change$Rich._t1<-rich_t1
      
      out_t2 = data.frame()
      tmp<-as.numeric(one_t2[1,])
      names(tmp)<-colnames(one_t2)
      tmp<-tmp[names(tmp) != 'estab']#get rid of VIRO
      tmp<-tmp[names(tmp) != 'nati_richness']#get rid of VIRO
      tmp<-tmp[tmp>0]
      out_t2<-data.frame(Abun=tmp[order(tmp,decreasing=TRUE)],Rank=1:length(tmp))
      rich_t2<-nrow(out_t2)
      DSV.change$Rich._t2 <- rich_t2
      
      #focus on extinctions
      
      exts<-NULL
      ranks<-NULL
      
      ## comapre the t and t+1, the extinction species number and mean rank
        #i = 1
      a<-out_t1
      b<-out_t2
      tmp.e<-is.na(match(rownames(a),rownames(b)))
      exts<-sum(tmp.e)
      ranks<-mean(a$Rank[tmp.e])

      ranks[is.nan(ranks)]<-NA
      
      DSV.change$obs.ext<-exts
      DSV.change$obs.rank.ext<-ranks
      
      ##do ses on each community
      #calculating only for plots with increasing DSV

      tmp<-out_t1
      if (DSV.change$DSV.changes <= 0) mod.out[[i]]<-NA
      if (DSV.change$DSV.changes > 0){
          mod.out <- ext.sum.stat(tmp,DSV.change$DSV.changes, 999,
                                     DSV.change$obs.ext, DSV.change$obs.rank.ext)
        }
      mod.out$neg.rank.z.value = -1*(mod.out$rank.z.value)
      mod.out$sis.score = (mod.out$neg.rank.z.value+mod.out$ext.z.value)/2
      mod.out$f_p = f_p_name
      cover_ab_estab_for_impact[[i]] = mod.out
  } else {
    cover_ab_estab_for_impact[[i]] = NULL
  }
  
}
library(data.table)
cover_ab_estab_for_impact_dat = rbindlist(cover_ab_estab_for_impact)
write.table(cover_ab_estab_for_impact_dat, 'cover_ab_estab_for_impact_dat.txt')
