###spatio-temporal data
###data sorting
library(spacetime)
library(sp)
library(foreign)
library(dplyr)

##### Data loading
sp_cover = read.csv('D:/BSS study/data/original data/BSS_community_332.csv',
                    header = T)
field = read.csv('D:/BSS study/data/original data/FIELDS.csv',
                 header = T)
trait = read.csv('D:/BSS study/data/original data/traits332.csv',
                 header = T)

sp_cover$f_p = paste(sp_cover$Field, sp_cover$Plot, sep = '_')
sp_cover = sp_cover %>% relocate(f_p, .after = Plot)
trait$species_id = paste('sp', trait$species_id, sep = '_')
sp_real = vector()
sp = colnames(sp_cover)[8:ncol(sp_cover)]

for (i in 1:length(sp)) {
  sp_real[i] = trait[trait$species_id == sp[i],]$Species
}

colnames(sp_cover)[8:ncol(sp_cover)] = sp_real # change the species name of data

sp_cover_fp = split(sp_cover, sp_cover$f_p)

## species persistence probability
cover_end_all = list()
for (k in 1:length(sp_cover_fp)) {
  
  #k = 1  
  cover = sp_cover_fp[[k]]
  year = unique(cover$Year)
  cover_end = cover[cover$Year == max(year),]
  cover_end[,8:ncol(cover_end)][cover_end[,8:ncol(cover_end)] > 0] = 1
  cover_end_all[[k]] = cover_end
  
}
library(data.table)
cover_end_all = rbindlist(cover_end_all)

# prepare for fit 

load("D:/BSS study/code/invasion_orginaldata.RData")

source('D:/BSS study/code/data collation_function.r')

# separate different time to fit
# About 20 years
re_cover_prop_estab_t1_l = list()
re_cover_prop_domin_t1_l = list()
re_cover_prop_intro_t1_l = list()
re_cover_prop_nati_t1_l = list()

re_cover_prop_estab_t2_l = list()
re_cover_prop_domin_t2_l = list()
re_cover_prop_intro_t2_l = list()
re_cover_prop_nati_t2_l = list()

re_cover_prop_estab_t3_l = list()
re_cover_prop_domin_t3_l = list()
re_cover_prop_intro_t3_l = list()
re_cover_prop_nati_t3_l = list()

re_cover_ab_estab_t1_l = list()
re_cover_ab_domin_t1_l = list()
re_cover_ab_intro_t1_l = list()
re_cover_ab_nati_t1_l = list()

re_cover_ab_estab_t2_l = list()
re_cover_ab_domin_t2_l = list()
re_cover_ab_intro_t2_l = list()
re_cover_ab_nati_t2_l = list()

re_cover_ab_estab_t3_l = list()
re_cover_ab_domin_t3_l = list()
re_cover_ab_intro_t3_l = list()
re_cover_ab_nati_t3_l = list()

for (z in 1:length(re_cover_ab_estab_l)) {
  #z = 145
  trans_ab_plot = re_cover_ab_nati_l[[z]]
  age = trans_ab_plot$Age
  y_end = max(age)
  y_star = min(age)
  
  if (y_end %in% c(53:70)) {
    d = length(unique(age))/3
    d1 = floor(d)
    d2 = d-floor(d)
    if (d2 == 1/3) {
      t1 = seq(y_star, y_star + d1*2, by = 2)
      t2 = seq(y_star + d1*2 + 2, y_star + d1*2 + d1*2, by = 2)
      t3 = seq(y_star + d1*2 + 2 + d1*2, y_star + d1*2 + d1*2 + d1*2, by = 2)
    } else if (d2 == 2/3) {
      t1 = seq(y_star, y_star + d1*2, by = 2)
      t2 = seq(y_star + d1*2 + 2, y_star + d1*2 + (d1+1)*2, by = 2)
      t3 = seq( y_star + d1*2 + (d1+1)*2 + 2, y_star + d1*2 + (d1+1)*2 + d1*2, by = 2)
      
    } else {
      t1 = seq(y_star, y_star + (d1-1)*2, by = 2)
      t2 = seq(y_star + (d1-1)*2 + 2, y_star + (d1-1)*2 + d1*2, by = 2)
      t3 = seq(y_star + (d1-1)*2 + d1*2 + 2, y_star + (d1-1)*2 + d1*2+ d1*2, by = 2)
    }
  } else if (y_end %in% c(35:52)){
    d = length(unique(age))/2
    d1 = floor(d)
    d2 = d-floor(d)
    if (d2 == 0.5) {
      t1 = seq(y_star, y_star + (d1)*2, by = 2)
      t2 = seq(y_star + (d1)*2 + 2, y_star + (d1)*2 + d1*2, by = 2)
      t3 = 0
    } else {
      t1 = seq(y_star, y_star + (d1-1)*2, by = 2)
      t2 = seq(y_star + (d1-1)*2 + 2, y_star + (d1-1)*2 + d1*2, by = 2)
      t3 = 0
    }
  } else {print('have not enough years!')}
  
  re_cover_prop_estab = re_cover_prop_estab_l[[z]]
  re_cover_prop_domin = re_cover_prop_domin_l[[z]]
  re_cover_prop_intro = re_cover_prop_intro_l[[z]]
  re_cover_prop_nati = re_cover_prop_nati_l[[z]]
  
  re_cover_ab_estab = re_cover_ab_estab_l[[z]]
  re_cover_ab_domin = re_cover_ab_domin_l[[z]]
  re_cover_ab_intro = re_cover_ab_intro_l[[z]]
  re_cover_ab_nati = re_cover_ab_nati_l[[z]]
  
  re_cover_prop_estab_t1 = re_cover_prop_estab %>% filter(Age %in% t1)
  re_cover_prop_domin_t1 = re_cover_prop_domin %>% filter(Age %in% t1)
  re_cover_prop_intro_t1 = re_cover_prop_intro %>% filter(Age %in% t1)
  re_cover_prop_nati_t1 = re_cover_prop_nati %>% filter(Age %in% t1)
  
  re_cover_prop_estab_t2 = re_cover_prop_estab %>% filter(Age %in% t2)
  re_cover_prop_domin_t2 = re_cover_prop_domin %>% filter(Age %in% t2)
  re_cover_prop_intro_t2 = re_cover_prop_intro %>% filter(Age %in% t2)
  re_cover_prop_nati_t2 = re_cover_prop_nati %>% filter(Age %in% t2)
  
  re_cover_ab_estab_t1 = re_cover_ab_estab %>% filter(Age %in% t1)
  re_cover_ab_domin_t1 = re_cover_ab_domin %>% filter(Age %in% t1)
  re_cover_ab_intro_t1 = re_cover_ab_intro %>% filter(Age %in% t1)
  re_cover_ab_nati_t1 = re_cover_ab_nati %>% filter(Age %in% t1)
  
  re_cover_ab_estab_t2 = re_cover_ab_estab %>% filter(Age %in% t2)
  re_cover_ab_domin_t2 = re_cover_ab_domin %>% filter(Age %in% t2)
  re_cover_ab_intro_t2 = re_cover_ab_intro %>% filter(Age %in% t2)
  re_cover_ab_nati_t2 = re_cover_ab_nati %>% filter(Age %in% t2)
  
  if (sum(t3) != 0) {
    re_cover_prop_estab_t3 = re_cover_prop_estab %>% filter(Age %in% t3)
    re_cover_prop_domin_t3 = re_cover_prop_domin %>% filter(Age %in% t3)
    re_cover_prop_intro_t3 = re_cover_prop_intro %>% filter(Age %in% t3)
    re_cover_prop_nati_t3 = re_cover_prop_nati %>% filter(Age %in% t3)
    
    re_cover_ab_estab_t3 = re_cover_ab_estab %>% filter(Age %in% t3)
    re_cover_ab_domin_t3 = re_cover_ab_domin %>% filter(Age %in% t3)
    re_cover_ab_intro_t3 = re_cover_ab_intro %>% filter(Age %in% t3)
    re_cover_ab_nati_t3 = re_cover_ab_nati %>% filter(Age %in% t3)
  } else {
    re_cover_prop_estab_t3 = NULL
    re_cover_prop_domin_t3 = NULL
    re_cover_prop_intro_t3 = NULL
    re_cover_prop_nati_t3 = NULL
    
    re_cover_ab_estab_t3 = NULL
    re_cover_ab_domin_t3 = NULL
    re_cover_ab_intro_t3 = NULL
    re_cover_ab_nati_t3 = NULL    
  }
  
  re_cover_prop_estab_t1_l[[z]] = re_cover_prop_estab_t1
  re_cover_prop_domin_t1_l[[z]] = re_cover_prop_domin_t1
  re_cover_prop_intro_t1_l[[z]] = re_cover_prop_intro_t1
  re_cover_prop_nati_t1_l[[z]] = re_cover_prop_nati_t1
  
  re_cover_prop_estab_t2_l[[z]] = re_cover_prop_estab_t2
  re_cover_prop_domin_t2_l[[z]] = re_cover_prop_domin_t2
  re_cover_prop_intro_t2_l[[z]] = re_cover_prop_intro_t2
  re_cover_prop_nati_t2_l[[z]] = re_cover_prop_nati_t2
  
  re_cover_prop_estab_t3_l[[z]] = re_cover_prop_estab_t3
  re_cover_prop_domin_t3_l[[z]] = re_cover_prop_domin_t3
  re_cover_prop_intro_t3_l[[z]] = re_cover_prop_intro_t3
  re_cover_prop_nati_t3_l[[z]] = re_cover_prop_nati_t3
  
  re_cover_ab_estab_t1_l[[z]] = re_cover_ab_estab_t1
  re_cover_ab_domin_t1_l[[z]] = re_cover_ab_domin_t1
  re_cover_ab_intro_t1_l[[z]] = re_cover_ab_intro_t1
  re_cover_ab_nati_t1_l[[z]] = re_cover_ab_nati_t1
  
  re_cover_ab_estab_t2_l[[z]] = re_cover_ab_estab_t2
  re_cover_ab_domin_t2_l[[z]] = re_cover_ab_domin_t2
  re_cover_ab_intro_t2_l[[z]] = re_cover_ab_intro_t2
  re_cover_ab_nati_t2_l[[z]] = re_cover_ab_nati_t2
  
  re_cover_ab_estab_t3_l[[z]] = re_cover_ab_estab_t3
  re_cover_ab_domin_t3_l[[z]] = re_cover_ab_domin_t3
  re_cover_ab_intro_t3_l[[z]] = re_cover_ab_intro_t3
  re_cover_ab_nati_t3_l[[z]] = re_cover_ab_nati_t3
}

### For t1
re_cover_ab_estab_t1_l_real = list()
re_cover_ab_domin_t1_l_real = list()
re_cover_ab_intro_t1_l_real = list()

re_cover_prop_estab_t1_l_real = list()
re_cover_prop_domin_t1_l_real = list()
re_cover_prop_intro_t1_l_real = list()

stage_t1 = data.frame()

for (z in 1:length(unique(sp_cover$f_p))) {
  #z = 1
  cover = sp_cover_fp[[z]]
  age = unique(cover$Age)
 
  intro_ab_t1 = re_cover_ab_intro_t1_l[[z]]
  estab_ab_t1 = re_cover_ab_estab_t1_l[[z]]
  nati_ab_t1 = re_cover_ab_nati_t1_l[[z]]
    
  intro_prop_t1 = re_cover_prop_intro_t1_l[[z]]
  estab_prop_t1 = re_cover_prop_estab_t1_l[[z]]
  nati_prop_t1 = re_cover_prop_nati_t1_l[[z]]

  re_cover_ab_t1 = cbind(nati_ab_t1[,1:7], intro_ab_t1[,8:ncol(intro_ab_t1)], estab_ab_t1[,8:ncol(estab_ab_t1)])
  re_cover_prop_t1 =  cbind(nati_prop_t1[,1:7], intro_prop_t1[,8:ncol(intro_prop_t1)], estab_prop_t1[,8:ncol(estab_prop_t1)])
  
  t_1_range = c(min(nati_ab_t1$Age):max(nati_ab_t1$Age))
  
  t1 = unique(nati_ab_t1$Age)

  cover_t1 = cover %>% filter(Age %in% t_1_range)
  cover_ab_t1 = cbind(cover_t1[,1:7], cover_t1[,8:(ncol(cover_t1))][colSums(cover_t1[,8:(ncol(cover_t1))]) != 0])
  cover_prop_t1 = cbind(cover_ab_t1[,1:7], t(apply(cover_ab_t1[,8:(ncol(cover_ab_t1))], 1, function(x){x/sum(x)})))
  
  cover_prop_t1_inva = cbind(cover_prop_t1[,1:7],
                          cover_prop_t1 %>% select(intersect(colnames(cover_prop_t1)[8:ncol(cover_prop_t1)], unlist(trait %>% filter(Origin == 'Exotic') %>% select(Species)))))
  cover_prop_t1_nati = cbind(cover_prop_t1[,1:7],
                          cover_prop_t1 %>% select(intersect(colnames(cover_prop_t1)[8:ncol(cover_prop_t1)], unlist(trait %>% filter(Origin == 'Native') %>% select(Species)))))
  
  inva_sp = colnames(cover_prop_t1_inva)[8:ncol(cover_prop_t1_inva)]
  estab_sp = vector() ## species at establishment stage (occur more than 10 successive years)
  domin_sp = vector() ## species at dominant stage (occur more than 10 successive years and mean cover is beyond any natives)
  for (j in 1:length(unique(cover_prop_t1_inva$Absolute_year))) {
    #j = 13
    y_now = unique(cover_prop_t1_inva$Absolute_year)[j]
    x = cover_prop_t1_inva %>% filter(Absolute_year %in% (y_now:(y_now+9)))
    x = arrange(x, x$Absolute_year)
    year_dif = max(x$Absolute_year)-min(x$Absolute_year)
    if (year_dif %in% c(6:8)){
      x = cover_prop_t1_inva %>% filter(Absolute_year %in% (y_now:(y_now+10)))
      x = arrange(x, x$Absolute_year)
      year_dif = max(x$Absolute_year)-min(x$Absolute_year)
      if (year_dif == 9) {
      y = data.frame(freq = apply(x[,8:ncol(x)], 2, function(y){length(which(y != 0))})) %>% filter(freq > (nrow(x)-1))
      
      } else if (year_dif > 9) {
        y = data.frame(freq = apply(x[,8:ncol(x)], 2, function(y){length(which(y != 0))})) %>% filter(freq > (nrow(x)-2))  
      } 
      }else if (year_dif == 9) {
      y = data.frame(freq = apply(x[,8:ncol(x)], 2, function(y){length(which(y != 0))})) %>% filter(freq > (nrow(x)-1))
    } else if (year_dif > 9) {
      y = data.frame(freq = apply(x[,8:ncol(x)], 2, function(y){length(which(y != 0))})) %>% filter(freq > (nrow(x)-2))  
    }
    
    if(nrow(y) == 0){print('no')} else{ 
      y_1 = cbind(x[,1:7],x %>% select(row.names(y)))
      x_nati = cover_prop_t1_nati %>% filter(Absolute_year %in% y_1$Absolute_year)
      max_nati = max(apply(x_nati[,8:ncol(x_nati)], 2, function(x){mean(x)}))
      y_2 = data.frame(domin = apply(data.frame(y_1[,8:ncol(y_1)]), 2, function(y){mean(y) > max_nati})) %>% filter(domin == T)
      
      if(ncol(y_1) == 8 & nrow(y_2) == 1){
        row.names(y_2) = colnames(y_1)[8]} else if (nrow(y_2) >= 1) {row.names(y_2) = row.names(y_2)}
      
      domin_sp = c(domin_sp, row.names(y_2))
      estab_sp = c(estab_sp, row.names(y))
    }
  }
  estab_sp = unique(estab_sp)
  domin_sp = unique(domin_sp)
  intro_sp = setdiff(inva_sp, estab_sp)
  
  if (length(estab_sp) != 0) {
  re_cover_prop_t1_estab = cbind(re_cover_prop_t1[,1:7], re_cover_prop_t1 %>% select(any_of(estab_sp)))
  } else {re_cover_prop_t1_estab = NULL}
  if (length(domin_sp) != 0) {
  re_cover_prop_t1_domin = cbind(re_cover_prop_t1[,1:7], re_cover_prop_t1 %>% select(any_of(domin_sp)))
  } else {re_cover_prop_t1_domin = NULL}
  re_cover_prop_t1_intro = cbind(re_cover_prop_t1[,1:7], re_cover_prop_t1 %>% select(any_of(intro_sp)))
  
  if (length(estab_sp) != 0) {
  re_cover_ab_t1_estab = cbind(re_cover_ab_t1[,1:7], re_cover_ab_t1 %>% select(any_of(estab_sp)))
  } else {re_cover_ab_t1_estab = NULL}
  if (length(domin_sp) != 0) {
  re_cover_ab_t1_domin = cbind(re_cover_ab_t1[,1:7], re_cover_ab_t1 %>% select(any_of(domin_sp)))
  } else {re_cover_ab_t1_domin = NULL}
  re_cover_ab_t1_intro = cbind(re_cover_ab_t1[,1:7], re_cover_ab_t1 %>% select(any_of(intro_sp)))
  
  re_cover_ab_estab_t1_l_real[[z]] = re_cover_ab_t1_estab
  re_cover_ab_domin_t1_l_real[[z]] = re_cover_ab_t1_domin
  re_cover_ab_intro_t1_l_real[[z]] = re_cover_ab_t1_intro
  
  re_cover_prop_estab_t1_l_real[[z]] = re_cover_prop_t1_estab
  re_cover_prop_domin_t1_l_real[[z]] = re_cover_prop_t1_domin
  re_cover_prop_intro_t1_l_real[[z]] = re_cover_prop_t1_intro
  
  stage_t1_1 = data.frame(stage = c('intro', 'estab', 'domin'),
                          species = c(paste(intro_sp, collapse = ','),
                                      paste(estab_sp, collapse = ','),
                                      paste(domin_sp, collapse = ',')))
  stage_t1_1$f_p = unique(re_cover_ab_t1_intro$f_p)
  stage_t1_1$time = 't1'
  stage_t1 = rbind(stage_t1, stage_t1_1)

}


### For t2
re_cover_ab_estab_t2_l_real = list()
re_cover_ab_domin_t2_l_real = list()
re_cover_ab_intro_t2_l_real = list()

re_cover_prop_estab_t2_l_real = list()
re_cover_prop_domin_t2_l_real = list()
re_cover_prop_intro_t2_l_real = list()

stage_t2 = data.frame()

for (z in 1:length(unique(sp_cover$f_p))) {
  #z = 1
  cover = sp_cover_fp[[z]]
  age = unique(cover$Age)
  
  intro_ab_t2 = re_cover_ab_intro_t2_l[[z]]
  estab_ab_t2 = re_cover_ab_estab_t2_l[[z]]
  nati_ab_t2 = re_cover_ab_nati_t2_l[[z]]
  
  intro_prop_t2 = re_cover_prop_intro_t2_l[[z]]
  estab_prop_t2 = re_cover_prop_estab_t2_l[[z]]
  nati_prop_t2 = re_cover_prop_nati_t2_l[[z]]
  
  re_cover_ab_t2 = cbind(nati_ab_t2[,1:7], intro_ab_t2[,8:ncol(intro_ab_t2)], estab_ab_t2[,8:ncol(estab_ab_t2)])
  re_cover_prop_t2 =  cbind(nati_prop_t2[,1:7], intro_prop_t2[,8:ncol(intro_prop_t2)], estab_prop_t2[,8:ncol(estab_prop_t2)])
  
  t_1_range = c(min(nati_ab_t2$Age):max(nati_ab_t2$Age))
  
  t2 = unique(nati_ab_t2$Age)
  
  cover_t2 = cover %>% filter(Age %in% t_1_range)
  cover_ab_t2 = cbind(cover_t2[,1:7], cover_t2[,8:(ncol(cover_t2))][colSums(cover_t2[,8:(ncol(cover_t2))]) != 0])
  cover_prop_t2 = cbind(cover_ab_t2[,1:7], t(apply(cover_ab_t2[,8:(ncol(cover_ab_t2))], 1, function(x){x/sum(x)})))
  
  cover_prop_t2_inva = cbind(cover_prop_t2[,1:7],
                             cover_prop_t2 %>% select(intersect(colnames(cover_prop_t2)[8:ncol(cover_prop_t2)], unlist(trait %>% filter(Origin == 'Exotic') %>% select(Species)))))
  cover_prop_t2_nati = cbind(cover_prop_t2[,1:7],
                             cover_prop_t2 %>% select(intersect(colnames(cover_prop_t2)[8:ncol(cover_prop_t2)], unlist(trait %>% filter(Origin == 'Native') %>% select(Species)))))
  
  inva_sp = colnames(cover_prop_t2_inva)[8:ncol(cover_prop_t2_inva)]
  estab_sp = vector() ## species at establishment stage (occur more than 10 successive years)
  domin_sp = vector() ## species at dominant stage (occur more than 10 successive years and mean cover is beyond any natives)
  for (j in 1:length(unique(cover_prop_t2_inva$Absolute_year))) {
    #j = 13
    y_now = unique(cover_prop_t2_inva$Absolute_year)[j]
    x = cover_prop_t2_inva %>% filter(Absolute_year %in% (y_now:(y_now+9)))
    x = arrange(x, x$Absolute_year)
    year_dif = max(x$Absolute_year)-min(x$Absolute_year)
    if (year_dif %in% c(6:8)){
      x = cover_prop_t2_inva %>% filter(Absolute_year %in% (y_now:(y_now+10)))
      x = arrange(x, x$Absolute_year)
      year_dif = max(x$Absolute_year)-min(x$Absolute_year)
      if (year_dif == 9) {
        y = data.frame(freq = apply(x[,8:ncol(x)], 2, function(y){length(which(y != 0))})) %>% filter(freq > (nrow(x)-1))
        
      } else if (year_dif > 9) {
        y = data.frame(freq = apply(x[,8:ncol(x)], 2, function(y){length(which(y != 0))})) %>% filter(freq > (nrow(x)-2))  
      } 
    }else if (year_dif == 9) {
      y = data.frame(freq = apply(x[,8:ncol(x)], 2, function(y){length(which(y != 0))})) %>% filter(freq > (nrow(x)-1))
    } else if (year_dif > 9) {
      y = data.frame(freq = apply(x[,8:ncol(x)], 2, function(y){length(which(y != 0))})) %>% filter(freq > (nrow(x)-2))  
    }
    
    if(nrow(y) == 0){print('no')} else{ 
      y_1 = cbind(x[,1:7],x %>% select(row.names(y)))
      x_nati = cover_prop_t2_nati %>% filter(Absolute_year %in% y_1$Absolute_year)
      max_nati = max(apply(x_nati[,8:ncol(x_nati)], 2, function(x){mean(x)}))
      y_2 = data.frame(domin = apply(data.frame(y_1[,8:ncol(y_1)]), 2, function(y){mean(y) > max_nati})) %>% filter(domin == T)
      
      if(ncol(y_1) == 8 & nrow(y_2) == 1){
        row.names(y_2) = colnames(y_1)[8]} else if (nrow(y_2) >= 1) {row.names(y_2) = row.names(y_2)}
      
      domin_sp = c(domin_sp, row.names(y_2))
      estab_sp = c(estab_sp, row.names(y))
    }
  }
  estab_sp = unique(estab_sp)
  domin_sp = unique(domin_sp)
  intro_sp = setdiff(inva_sp, estab_sp)
  
  if (length(estab_sp) != 0) {
    re_cover_prop_t2_estab = cbind(re_cover_prop_t2[,1:7], re_cover_prop_t2 %>% select(any_of(estab_sp)))
  } else {re_cover_prop_t2_estab = NULL}
  if (length(domin_sp) != 0) {
    re_cover_prop_t2_domin = cbind(re_cover_prop_t2[,1:7], re_cover_prop_t2 %>% select(any_of(domin_sp)))
  } else {re_cover_prop_t2_domin = NULL}
  re_cover_prop_t2_intro = cbind(re_cover_prop_t2[,1:7], re_cover_prop_t2 %>% select(any_of(intro_sp)))
  
  if (length(estab_sp) != 0) {
    re_cover_ab_t2_estab = cbind(re_cover_ab_t2[,1:7], re_cover_ab_t2 %>% select(any_of(estab_sp)))
  } else {re_cover_ab_t2_estab = NULL}
  if (length(domin_sp) != 0) {
    re_cover_ab_t2_domin = cbind(re_cover_ab_t2[,1:7], re_cover_ab_t2 %>% select(any_of(domin_sp)))
  } else {re_cover_ab_t2_domin = NULL}
  re_cover_ab_t2_intro = cbind(re_cover_ab_t2[,1:7], re_cover_ab_t2 %>% select(any_of(intro_sp)))
  
  re_cover_ab_estab_t2_l_real[[z]] = re_cover_ab_t2_estab
  re_cover_ab_domin_t2_l_real[[z]] = re_cover_ab_t2_domin
  re_cover_ab_intro_t2_l_real[[z]] = re_cover_ab_t2_intro
  
  re_cover_prop_estab_t2_l_real[[z]] = re_cover_prop_t2_estab
  re_cover_prop_domin_t2_l_real[[z]] = re_cover_prop_t2_domin
  re_cover_prop_intro_t2_l_real[[z]] = re_cover_prop_t2_intro
  
  stage_t2_1 = data.frame(stage = c('intro', 'estab', 'domin'),
                          species = c(paste(intro_sp, collapse = ','),
                                      paste(estab_sp, collapse = ','),
                                      paste(domin_sp, collapse = ',')))
  stage_t2_1$f_p = unique(re_cover_ab_t2_intro$f_p)
  stage_t2_1$time = 't2'
  stage_t2 = rbind(stage_t2, stage_t2_1)
  
}


### For t3
re_cover_ab_estab_t3_l_real = list()
re_cover_ab_domin_t3_l_real = list()
re_cover_ab_intro_t3_l_real = list()

re_cover_prop_estab_t3_l_real = list()
re_cover_prop_domin_t3_l_real = list()
re_cover_prop_intro_t3_l_real = list()

stage_t3 = data.frame()

for (z in 1:length(unique(sp_cover$f_p))) {
  #z = 96
  cover = sp_cover_fp[[z]]
  age = unique(cover$Age)
  
  intro_ab_t3 = re_cover_ab_intro_t3_l[[z]]
  estab_ab_t3 = re_cover_ab_estab_t3_l[[z]]
  nati_ab_t3 = re_cover_ab_nati_t3_l[[z]]
  if (!is.null(intro_ab_t3)) {
  intro_prop_t3 = re_cover_prop_intro_t3_l[[z]]
  estab_prop_t3 = re_cover_prop_estab_t3_l[[z]]
  nati_prop_t3 = re_cover_prop_nati_t3_l[[z]]
  
  re_cover_ab_t3 = cbind(nati_ab_t3[,1:7], intro_ab_t3[,8:ncol(intro_ab_t3)], estab_ab_t3[,8:ncol(estab_ab_t3)])
  re_cover_prop_t3 =  cbind(nati_prop_t3[,1:7], intro_prop_t3[,8:ncol(intro_prop_t3)], estab_prop_t3[,8:ncol(estab_prop_t3)])
  
  t_1_range = c(min(nati_ab_t3$Age):max(nati_ab_t3$Age))
  
  t3 = unique(nati_ab_t3$Age)
  
  cover_t3 = cover %>% filter(Age %in% t_1_range)
  cover_ab_t3 = cbind(cover_t3[,1:7], cover_t3[,8:(ncol(cover_t3))][colSums(cover_t3[,8:(ncol(cover_t3))]) != 0])
  cover_prop_t3 = cbind(cover_ab_t3[,1:7], t(apply(cover_ab_t3[,8:(ncol(cover_ab_t3))], 1, function(x){x/sum(x)})))
  
  cover_prop_t3_inva = cbind(cover_prop_t3[,1:7],
                             cover_prop_t3 %>% select(intersect(colnames(cover_prop_t3)[8:ncol(cover_prop_t3)], unlist(trait %>% filter(Origin == 'Exotic') %>% select(Species)))))
  cover_prop_t3_nati = cbind(cover_prop_t3[,1:7],
                             cover_prop_t3 %>% select(intersect(colnames(cover_prop_t3)[8:ncol(cover_prop_t3)], unlist(trait %>% filter(Origin == 'Native') %>% select(Species)))))
  
  inva_sp = colnames(cover_prop_t3_inva)[8:ncol(cover_prop_t3_inva)]
  estab_sp = vector() ## species at establishment stage (occur more than 10 successive years)
  domin_sp = vector() ## species at dominant stage (occur more than 10 successive years and mean cover is beyond any natives)
  for (j in 1:length(unique(cover_prop_t3_inva$Absolute_year))) {
    #j = 13
    y_now = unique(cover_prop_t3_inva$Absolute_year)[j]
    x = cover_prop_t3_inva %>% filter(Absolute_year %in% (y_now:(y_now+9)))
    x = arrange(x, x$Absolute_year)
    year_dif = max(x$Absolute_year)-min(x$Absolute_year)
    if (year_dif %in% c(6:8)){
      x = cover_prop_t3_inva %>% filter(Absolute_year %in% (y_now:(y_now+10)))
      x = arrange(x, x$Absolute_year)
      year_dif = max(x$Absolute_year)-min(x$Absolute_year)
      if (year_dif == 9) {
        y = data.frame(freq = apply(x[,8:ncol(x)], 2, function(y){length(which(y != 0))})) %>% filter(freq > (nrow(x)-1))
        
      } else if (year_dif > 9) {
        y = data.frame(freq = apply(x[,8:ncol(x)], 2, function(y){length(which(y != 0))})) %>% filter(freq > (nrow(x)-2))  
      } 
    }else if (year_dif == 9) {
      y = data.frame(freq = apply(x[,8:ncol(x)], 2, function(y){length(which(y != 0))})) %>% filter(freq > (nrow(x)-1))
    } else if (year_dif > 9) {
      y = data.frame(freq = apply(x[,8:ncol(x)], 2, function(y){length(which(y != 0))})) %>% filter(freq > (nrow(x)-2))  
    }
    
    if(nrow(y) == 0){print('no')} else{ 
      y_1 = cbind(x[,1:7],x %>% select(row.names(y)))
      x_nati = cover_prop_t3_nati %>% filter(Absolute_year %in% y_1$Absolute_year)
      max_nati = max(apply(x_nati[,8:ncol(x_nati)], 2, function(x){mean(x)}))
      y_2 = data.frame(domin = apply(data.frame(y_1[,8:ncol(y_1)]), 2, function(y){mean(y) > max_nati})) %>% filter(domin == T)
      
      if(ncol(y_1) == 8 & nrow(y_2) == 1){
        row.names(y_2) = colnames(y_1)[8]} else if (nrow(y_2) >= 1) {row.names(y_2) = row.names(y_2)}
      
      domin_sp = c(domin_sp, row.names(y_2))
      estab_sp = c(estab_sp, row.names(y))
    }
  }
  estab_sp = unique(estab_sp)
  domin_sp = unique(domin_sp)
  intro_sp = setdiff(inva_sp, estab_sp)
  
  if (length(estab_sp) != 0) {
    re_cover_prop_t3_estab = cbind(re_cover_prop_t3[,1:7], re_cover_prop_t3 %>% select(any_of(estab_sp)))
  } else {re_cover_prop_t3_estab = NULL}
  if (length(domin_sp) != 0) {
    re_cover_prop_t3_domin = cbind(re_cover_prop_t3[,1:7], re_cover_prop_t3 %>% select(any_of(domin_sp)))
  } else {re_cover_prop_t3_domin = NULL}
  re_cover_prop_t3_intro = cbind(re_cover_prop_t3[,1:7], re_cover_prop_t3 %>% select(any_of(intro_sp)))
  
  if (length(estab_sp) != 0) {
    re_cover_ab_t3_estab = cbind(re_cover_ab_t3[,1:7], re_cover_ab_t3 %>% select(any_of(estab_sp)))
  } else {re_cover_ab_t3_estab = NULL}
  if (length(domin_sp) != 0) {
    re_cover_ab_t3_domin = cbind(re_cover_ab_t3[,1:7], re_cover_ab_t3 %>% select(any_of(domin_sp)))
  } else {re_cover_ab_t3_domin = NULL}
  re_cover_ab_t3_intro = cbind(re_cover_ab_t3[,1:7], re_cover_ab_t3 %>% select(any_of(intro_sp)))
  
  re_cover_ab_estab_t3_l_real[[z]] = re_cover_ab_t3_estab
  re_cover_ab_domin_t3_l_real[[z]] = re_cover_ab_t3_domin
  re_cover_ab_intro_t3_l_real[[z]] = re_cover_ab_t3_intro
  
  re_cover_prop_estab_t3_l_real[[z]] = re_cover_prop_t3_estab
  re_cover_prop_domin_t3_l_real[[z]] = re_cover_prop_t3_domin
  re_cover_prop_intro_t3_l_real[[z]] = re_cover_prop_t3_intro
  
  stage_t3_1 = data.frame(stage = c('intro', 'estab', 'domin'),
                          species = c(paste(intro_sp, collapse = ','),
                                      paste(estab_sp, collapse = ','),
                                      paste(domin_sp, collapse = ',')))
  stage_t3_1$f_p = unique(re_cover_ab_t3_intro$f_p)
  stage_t3_1$time = 't3'
  stage_t3 = rbind(stage_t3, stage_t3_1)}
  
}

save.image(file = "D:/BSS study/code/data_for_timesepfit.RData")

### Filter different time fit data ###
re_cover_prop_estab_t1_l_f = list()
re_cover_prop_domin_t1_l_f = list()
re_cover_prop_intro_t1_l_f = list()
re_cover_prop_nati_t1_l_f = list()

re_cover_prop_estab_t2_l_f = list()
re_cover_prop_domin_t2_l_f = list()
re_cover_prop_intro_t2_l_f = list()
re_cover_prop_nati_t2_l_f = list()

re_cover_prop_estab_t3_l_f = list()
re_cover_prop_domin_t3_l_f = list()
re_cover_prop_intro_t3_l_f = list()
re_cover_prop_nati_t3_l_f = list()

### for prop
# for t1
for (i in 1:length(re_cover_prop_estab_l)){
  #i = 313
  x = re_cover_prop_estab_t1_l[[i]]
  if (ncol(x) < 9){
    x = cbind(x, fake_sp = 0, fake_sp2 = 0)
  }
  y = x[,8:ncol(x)]
  z = filter_comm_mean(y, per = 0.3, rare = F)
  if(is.null(z)){re_cover_prop_estab_t1_l_f[[i]] = NULL
  sr_estab = 0
  } else {
    #  colnames(z)[ncol(z)] = 'estab_rare'
    f = cbind(x[,1:7], z)
    re_cover_prop_estab_t1_l_f[[i]] = f
    sr_estab = ncol(f)-7}
  
  x_1 = re_cover_prop_intro_t1_l[[i]]
  if (ncol(x_1) < 9){
    x_1 = cbind(x_1, fake_sp = 0, fake_sp2 = 0)
  }
  y_1 = x_1[,8:ncol(x_1)]
  z_1 = filter_comm_mean(y_1, per = 0.15, rare = F)
  if(is.null(z_1)){
    re_cover_prop_intro_t1_l_f[[i]] = NULL
    sr_intro = 0} else {
      #   colnames(z_1)[ncol(z_1)] = 'intro_rare'
      f_1 = cbind(x_1[,1:7], z_1)
      re_cover_prop_intro_t1_l_f[[i]] = f_1
      sr_intro = ncol(f_1)-7}
  
  x_2 = re_cover_prop_nati_t1_l[[i]]
  if (ncol(x_2) < 9){
    x_2 = cbind(x_2, fake_sp = 0, fake_sp2 = 0)
  }
  y_2 = x_2[,8:ncol(x_2)]
  z_2 = filter_comm_mean(x = y_2, per = 0.3, rare = F)
  if(is.null(z_2)){re_cover_prop_nati_t1_l_f[[i]] = NULL
  sr_nati = 0
  } else {
    #  colnames(z_2)[ncol(z_2)] = 'nati_rare'
    f_2 = cbind(x_2[,1:7], z_2)
    re_cover_prop_nati_t1_l_f[[i]] = f_2
    sr_nati = ncol(f_2)-7}
  
  x_3 = re_cover_prop_domin_t1_l[[i]]
  if (ncol(x_3) < 9){
    x_3 = cbind(x_3, fake_sp = 0, fake_sp2 = 0)
  }
  y_3 = x_3[,8:ncol(x_3)]
  z_3 = filter_comm_mean(x = y_3, per = 0.3, rare = F)
  if(is.null(z_3)){re_cover_prop_domin_t1_l_f[[i]] = NULL
  sr_domin = 0
  } else {
    #  colnames(z_3)[ncol(z_3)] = 'domin_rare'
    f_3 = cbind(x_3[,1:7], z_3)
    re_cover_prop_domin_t1_l_f[[i]] = f_3
    sr_domin = ncol(f_3)-7}
  
}

# for t2
for (i in 1:length(re_cover_prop_estab_l)){
  #i = 313
  x = re_cover_prop_estab_t2_l[[i]]
  if (ncol(x) < 9){
    x = cbind(x, fake_sp = 0, fake_sp2 = 0)
  }
  y = x[,8:ncol(x)]
  z = filter_comm_mean(y, per = 0.3, rare = F)
  if(is.null(z)){re_cover_prop_estab_t2_l_f[[i]] = NULL
  sr_estab = 0
  } else {
    #  colnames(z)[ncol(z)] = 'estab_rare'
    f = cbind(x[,1:7], z)
    re_cover_prop_estab_t2_l_f[[i]] = f
    sr_estab = ncol(f)-7}
  
  x_1 = re_cover_prop_intro_t2_l[[i]]
  if (ncol(x_1) < 9){
    x_1 = cbind(x_1, fake_sp = 0, fake_sp2 = 0)
  }
  y_1 = x_1[,8:ncol(x_1)]
  z_1 = filter_comm_mean(y_1, per = 0.15, rare = F)
  if(is.null(z_1)){
    re_cover_prop_intro_t2_l_f[[i]] = NULL
    sr_intro = 0} else {
      #   colnames(z_1)[ncol(z_1)] = 'intro_rare'
      f_1 = cbind(x_1[,1:7], z_1)
      re_cover_prop_intro_t2_l_f[[i]] = f_1
      sr_intro = ncol(f_1)-7}
  
  x_2 = re_cover_prop_nati_t2_l[[i]]
  if (ncol(x_2) < 9){
    x_2 = cbind(x_2, fake_sp = 0, fake_sp2 = 0)
  }
  y_2 = x_2[,8:ncol(x_2)]
  z_2 = filter_comm_mean(x = y_2, per = 0.3, rare = F)
  if(is.null(z_2)){re_cover_prop_nati_t2_l_f[[i]] = NULL
  sr_nati = 0
  } else {
    #  colnames(z_2)[ncol(z_2)] = 'nati_rare'
    f_2 = cbind(x_2[,1:7], z_2)
    re_cover_prop_nati_t2_l_f[[i]] = f_2
    sr_nati = ncol(f_2)-7}
  
  x_3 = re_cover_prop_domin_t2_l[[i]]
  if (ncol(x_3) < 9){
    x_3 = cbind(x_3, fake_sp = 0, fake_sp2 = 0)
  }
  y_3 = x_3[,8:ncol(x_3)]
  z_3 = filter_comm_mean(x = y_3, per = 0.3, rare = F)
  if(is.null(z_3)){re_cover_prop_domin_t2_l_f[[i]] = NULL
  sr_domin = 0
  } else {
    #  colnames(z_3)[ncol(z_3)] = 'domin_rare'
    f_3 = cbind(x_3[,1:7], z_3)
    re_cover_prop_domin_t2_l_f[[i]] = f_3
    sr_domin = ncol(f_3)-7}
  
}

# for t3
for (i in 1:length(re_cover_prop_estab_l)){
  #i = 313
  x1 = re_cover_prop_estab_t3_l[[1]]
  x = re_cover_prop_estab_t3_l[[i]]
  if (is.null(x)) {
    x = cbind(x1[,1:7], fake_sp = 0, fake_sp2 = 0)
    
  } else if (ncol(x) < 9){
    x = cbind(x, fake_sp = 0, fake_sp2 = 0)
  } 
  y = x[,8:ncol(x)]
  z = filter_comm_mean(y, per = 0.3, rare = F)
  if(is.null(z)){re_cover_prop_estab_t3_l_f[[i]] = NULL
  sr_estab = 0
  } else {
    #  colnames(z)[ncol(z)] = 'estab_rare'
    f = cbind(x[,1:7], z)
    re_cover_prop_estab_t3_l_f[[i]] = f
    sr_estab = ncol(f)-7}
  
  x_1 = re_cover_prop_intro_t3_l[[i]]
  if (is.null(x_1)) {
    x_1 = cbind(x1[,1:7], fake_sp = 0, fake_sp2 = 0)
    
  } else if (ncol(x_1) < 9){
    x_1 = cbind(x_1, fake_sp = 0, fake_sp2 = 0)
  }
  y_1 = x_1[,8:ncol(x_1)]
  z_1 = filter_comm_mean(y_1, per = 0.15, rare = F)
  if(is.null(z_1)){
    re_cover_prop_intro_t3_l_f[[i]] = NULL
    sr_intro = 0} else {
      #   colnames(z_1)[ncol(z_1)] = 'intro_rare'
      f_1 = cbind(x_1[,1:7], z_1)
      re_cover_prop_intro_t3_l_f[[i]] = f_1
      sr_intro = ncol(f_1)-7}
  
  x_2 = re_cover_prop_nati_t3_l[[i]]
  if (is.null(x_2)) {
    x_2 = cbind(x1[,1:7], fake_sp = 0, fake_sp2 = 0)
    
  } else if (ncol(x_2) < 9){
    x_2 = cbind(x_2, fake_sp = 0, fake_sp2 = 0)
  }
  y_2 = x_2[,8:ncol(x_2)]
  z_2 = filter_comm_mean(x = y_2, per = 0.3, rare = F)
  if(is.null(z_2)){re_cover_prop_nati_t3_l_f[[i]] = NULL
  sr_nati = 0
  } else {
    #  colnames(z_2)[ncol(z_2)] = 'nati_rare'
    f_2 = cbind(x_2[,1:7], z_2)
    re_cover_prop_nati_t3_l_f[[i]] = f_2
    sr_nati = ncol(f_2)-7}
  
  x_3 = re_cover_prop_domin_t3_l[[i]]
  if (is.null(x_3)) {
    x_3 = cbind(x1[,1:7], fake_sp = 0, fake_sp2 = 0)
    
  } else if (ncol(x_3) < 9){
    x_3 = cbind(x_3, fake_sp = 0, fake_sp2 = 0)
  }
  y_3 = x_3[,8:ncol(x_3)]
  z_3 = filter_comm_mean(x = y_3, per = 0.3, rare = F)
  if(is.null(z_3)){re_cover_prop_domin_t3_l_f[[i]] = NULL
  sr_domin = 0
  } else {
    #  colnames(z_3)[ncol(z_3)] = 'domin_rare'
    f_3 = cbind(x_3[,1:7], z_3)
    re_cover_prop_domin_t3_l_f[[i]] = f_3
    sr_domin = ncol(f_3)-7}
  
}

### for ab
# for t1
re_cover_ab_estab_t1_l_f = list()
re_cover_ab_domin_t1_l_f = list()
re_cover_ab_intro_t1_l_f = list()
re_cover_ab_nati_t1_l_f = list()

for (i in 1:length(re_cover_ab_estab_l)){
  #i = 313
  x = re_cover_ab_estab_t1_l[[i]]
  if (ncol(x) < 9){
    x = cbind(x, fake_sp = 0, fake_sp2 = 0)
  }
  y = x[,8:ncol(x)]
  z = filter_comm_mean(y, per = 0.3, rare = F)
  if(is.null(z)){re_cover_ab_estab_t1_l_f[[i]] = NULL
  sr_estab = 0
  } else {
    #  colnames(z)[ncol(z)] = 'estab_rare'
    f = cbind(x[,1:7], z)
    re_cover_ab_estab_t1_l_f[[i]] = f
    sr_estab = ncol(f)-7}
  
  x_1 = re_cover_ab_intro_t1_l[[i]]
  if (ncol(x_1) < 9){
    x_1 = cbind(x_1, fake_sp = 0, fake_sp2 = 0)
  }
  y_1 = x_1[,8:ncol(x_1)]
  z_1 = filter_comm_mean(y_1, per = 0.15, rare = F)
  if(is.null(z_1)){
    re_cover_ab_intro_t1_l_f[[i]] = NULL
    sr_intro = 0} else {
      #   colnames(z_1)[ncol(z_1)] = 'intro_rare'
      f_1 = cbind(x_1[,1:7], z_1)
      re_cover_ab_intro_t1_l_f[[i]] = f_1
      sr_intro = ncol(f_1)-7}
  
  x_2 = re_cover_ab_nati_t1_l[[i]]
  if (ncol(x_2) < 9){
    x_2 = cbind(x_2, fake_sp = 0, fake_sp2 = 0)
  }
  y_2 = x_2[,8:ncol(x_2)]
  z_2 = filter_comm_mean(x = y_2, per = 0.3, rare = F)
  if(is.null(z_2)){re_cover_ab_nati_t1_l_f[[i]] = NULL
  sr_nati = 0
  } else {
    #  colnames(z_2)[ncol(z_2)] = 'nati_rare'
    f_2 = cbind(x_2[,1:7], z_2)
    re_cover_ab_nati_t1_l_f[[i]] = f_2
    sr_nati = ncol(f_2)-7}
  
  x_3 = re_cover_ab_domin_t1_l[[i]]
  if (ncol(x_3) < 9){
    x_3 = cbind(x_3, fake_sp = 0, fake_sp2 = 0)
  }
  y_3 = x_3[,8:ncol(x_3)]
  z_3 = filter_comm_mean(x = y_3, per = 0.3, rare = F)
  if(is.null(z_3)){re_cover_ab_domin_t1_l_f[[i]] = NULL
  sr_domin = 0
  } else {
    #  colnames(z_3)[ncol(z_3)] = 'domin_rare'
    f_3 = cbind(x_3[,1:7], z_3)
    re_cover_ab_domin_t1_l_f[[i]] = f_3
    sr_domin = ncol(f_3)-7}
  
}

# for t2
re_cover_ab_estab_t2_l_f = list()
re_cover_ab_domin_t2_l_f = list()
re_cover_ab_intro_t2_l_f = list()
re_cover_ab_nati_t2_l_f = list()

for (i in 1:length(re_cover_ab_estab_l)){
  #i = 313
  x = re_cover_ab_estab_t2_l[[i]]
  if (ncol(x) < 9){
    x = cbind(x, fake_sp = 0, fake_sp2 = 0)
  }
  y = x[,8:ncol(x)]
  z = filter_comm_mean(y, per = 0.3, rare = F)
  if(is.null(z)){re_cover_ab_estab_t2_l_f[[i]] = NULL
  sr_estab = 0
  } else {
    #  colnames(z)[ncol(z)] = 'estab_rare'
    f = cbind(x[,1:7], z)
    re_cover_ab_estab_t2_l_f[[i]] = f
    sr_estab = ncol(f)-7}
  
  x_1 = re_cover_ab_intro_t2_l[[i]]
  if (ncol(x_1) < 9){
    x_1 = cbind(x_1, fake_sp = 0, fake_sp2 = 0)
  }
  y_1 = x_1[,8:ncol(x_1)]
  z_1 = filter_comm_mean(y_1, per = 0.15, rare = F)
  if(is.null(z_1)){
    re_cover_ab_intro_t2_l_f[[i]] = NULL
    sr_intro = 0} else {
      #   colnames(z_1)[ncol(z_1)] = 'intro_rare'
      f_1 = cbind(x_1[,1:7], z_1)
      re_cover_ab_intro_t2_l_f[[i]] = f_1
      sr_intro = ncol(f_1)-7}
  
  x_2 = re_cover_ab_nati_t2_l[[i]]
  if (ncol(x_2) < 9){
    x_2 = cbind(x_2, fake_sp = 0, fake_sp2 = 0)
  }
  y_2 = x_2[,8:ncol(x_2)]
  z_2 = filter_comm_mean(x = y_2, per = 0.3, rare = F)
  if(is.null(z_2)){re_cover_ab_nati_t2_l_f[[i]] = NULL
  sr_nati = 0
  } else {
    #  colnames(z_2)[ncol(z_2)] = 'nati_rare'
    f_2 = cbind(x_2[,1:7], z_2)
    re_cover_ab_nati_t2_l_f[[i]] = f_2
    sr_nati = ncol(f_2)-7}
  
  x_3 = re_cover_ab_domin_t2_l[[i]]
  if (ncol(x_3) < 9){
    x_3 = cbind(x_3, fake_sp = 0, fake_sp2 = 0)
  }
  y_3 = x_3[,8:ncol(x_3)]
  z_3 = filter_comm_mean(x = y_3, per = 0.3, rare = F)
  if(is.null(z_3)){re_cover_ab_domin_t2_l_f[[i]] = NULL
  sr_domin = 0
  } else {
    #  colnames(z_3)[ncol(z_3)] = 'domin_rare'
    f_3 = cbind(x_3[,1:7], z_3)
    re_cover_ab_domin_t2_l_f[[i]] = f_3
    sr_domin = ncol(f_3)-7}
  
}

# for t3
re_cover_ab_estab_t3_l_f = list()
re_cover_ab_domin_t3_l_f = list()
re_cover_ab_intro_t3_l_f = list()
re_cover_ab_nati_t3_l_f = list()

for (i in 1:length(re_cover_ab_estab_l)){
  #i = 313
  x1 = re_cover_ab_estab_t3_l[[1]]
  x = re_cover_ab_estab_t3_l[[i]]
  if (is.null(x)) {
    x = cbind(x1[,1:7], fake_sp = 0, fake_sp2 = 0)
    
  } else if (ncol(x) < 9){
    x = cbind(x, fake_sp = 0, fake_sp2 = 0)
  } 
  y = x[,8:ncol(x)]
  z = filter_comm_mean(y, per = 0.3, rare = F)
  if(is.null(z)){re_cover_ab_estab_t3_l_f[[i]] = NULL
  sr_estab = 0
  } else {
    #  colnames(z)[ncol(z)] = 'estab_rare'
    f = cbind(x[,1:7], z)
    re_cover_ab_estab_t3_l_f[[i]] = f
    sr_estab = ncol(f)-7}
  
  x_1 = re_cover_ab_intro_t3_l[[i]]
  if (is.null(x_1)) {
    x_1 = cbind(x1[,1:7], fake_sp = 0, fake_sp2 = 0)
    
  } else if (ncol(x_1) < 9){
    x_1 = cbind(x_1, fake_sp = 0, fake_sp2 = 0)
  }
  y_1 = x_1[,8:ncol(x_1)]
  z_1 = filter_comm_mean(y_1, per = 0.15, rare = F)
  if(is.null(z_1)){
    re_cover_ab_intro_t3_l_f[[i]] = NULL
    sr_intro = 0} else {
      #   colnames(z_1)[ncol(z_1)] = 'intro_rare'
      f_1 = cbind(x_1[,1:7], z_1)
      re_cover_ab_intro_t3_l_f[[i]] = f_1
      sr_intro = ncol(f_1)-7}
  
  x_2 = re_cover_ab_nati_t3_l[[i]]
  if (is.null(x_2)) {
    x_2 = cbind(x1[,1:7], fake_sp = 0, fake_sp2 = 0)
    
  } else if (ncol(x_2) < 9){
    x_2 = cbind(x_2, fake_sp = 0, fake_sp2 = 0)
  }
  y_2 = x_2[,8:ncol(x_2)]
  z_2 = filter_comm_mean(x = y_2, per = 0.3, rare = F)
  if(is.null(z_2)){re_cover_ab_nati_t3_l_f[[i]] = NULL
  sr_nati = 0
  } else {
    #  colnames(z_2)[ncol(z_2)] = 'nati_rare'
    f_2 = cbind(x_2[,1:7], z_2)
    re_cover_ab_nati_t3_l_f[[i]] = f_2
    sr_nati = ncol(f_2)-7}
  
  x_3 = re_cover_ab_domin_t3_l[[i]]
  if (is.null(x_3)) {
    x_3 = cbind(x1[,1:7], fake_sp = 0, fake_sp2 = 0)
    
  } else if (ncol(x_3) < 9){
    x_3 = cbind(x_3, fake_sp = 0, fake_sp2 = 0)
  }
  y_3 = x_3[,8:ncol(x_3)]
  z_3 = filter_comm_mean(x = y_3, per = 0.3, rare = F)
  if(is.null(z_3)){re_cover_ab_domin_t3_l_f[[i]] = NULL
  sr_domin = 0
  } else {
    #  colnames(z_3)[ncol(z_3)] = 'domin_rare'
    f_3 = cbind(x_3[,1:7], z_3)
    re_cover_ab_domin_t3_l_f[[i]] = f_3
    sr_domin = ncol(f_3)-7}
  
}

# Species proportion
library(data.table)

load("D:/BSS study/code/data_for_timesepfit.RData")

cover_mprop_all_t1 = list()
cover_mprop_all_t2 = list()
cover_mprop_all_t3 = list()

## for t1
for (z in 1:length(unique(sp_cover$f_p))) {
  #z = 1
  cover = sp_cover_fp[[z]]
  cover_t1 = re_cover_ab_nati_t1_l[[z]]
  age = cover_t1$Age
  cover_ab = cbind(cover[,1:7], cover[,8:(ncol(cover))][colSums(cover[,8:(ncol(cover))]) != 0])
  cover_prop = cbind(cover_ab[,1:7], t(apply(cover_ab[,8:(ncol(cover_ab))], 1, function(x){x/sum(x)})))
  cover_prop_t1 = cover_prop %>% filter(Age %in% age)
  cover_mprop_t1 = cbind(cover_prop_t1[1,c(1, 5:7)], t(data.frame(colMeans(cover_prop_t1[,8:ncol(cover_prop_t1)]))))
  cover_mprop_all_t1[[z]] = cover_mprop_t1
}
cover_mprop_all_t1 = rbindlist(cover_mprop_all_t1, fill = T)

## for t2
for (z in 1:length(unique(sp_cover$f_p))) {
  #z = 1
  cover = sp_cover_fp[[z]]
  cover_t2 = re_cover_ab_nati_t2_l[[z]]
  age = cover_t2$Age
  cover_ab = cbind(cover[,1:7], cover[,8:(ncol(cover))][colSums(cover[,8:(ncol(cover))]) != 0])
  cover_prop = cbind(cover_ab[,1:7], t(apply(cover_ab[,8:(ncol(cover_ab))], 1, function(x){x/sum(x)})))
  cover_prop_t2 = cover_prop %>% filter(Age %in% age)
  cover_mprop_t2 = cbind(cover_prop_t2[1,c(1, 5:7)], t(data.frame(colMeans(cover_prop_t2[,8:ncol(cover_prop_t2)]))))
  cover_mprop_all_t2[[z]] = cover_mprop_t2
}
cover_mprop_all_t2 = rbindlist(cover_mprop_all_t2, fill = T)

## for t3
for (z in 1:length(unique(sp_cover$f_p))) {
  #z = 1
  cover = sp_cover_fp[[z]]
  cover_t3 = re_cover_ab_nati_t3_l[[z]]
  age = cover_t3$Age
  cover_ab = cbind(cover[,1:7], cover[,8:(ncol(cover))][colSums(cover[,8:(ncol(cover))]) != 0])
  cover_prop = cbind(cover_ab[,1:7], t(apply(cover_ab[,8:(ncol(cover_ab))], 1, function(x){x/sum(x)})))
  cover_prop_t3 = cover_prop %>% filter(Age %in% age)
  cover_mprop_t3 = cbind(cover_prop_t3[1,c(1, 5:7)], t(data.frame(colMeans(cover_prop_t3[,8:ncol(cover_prop_t3)]))))
  cover_mprop_all_t3[[z]] = cover_mprop_t3
}
cover_mprop_all_t3 = rbindlist(cover_mprop_all_t3, fill = T)
