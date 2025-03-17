###spatio-temporal data
###data sorting
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
cover_fp_c = list()
cover_fp_prop = list()
cover_prop_estab_l = list()
cover_prop_domin_l = list()
cover_prop_intro_l = list()
cover_prop_nati_l = list()

re_cover_fp_c = list()
re_cover_fp_prop = list()
re_cover_prop_estab_l = list()
re_cover_prop_domin_l = list()
re_cover_prop_intro_l = list()
re_cover_prop_nati_l = list()

re_cover_ab_estab_l = list()
re_cover_ab_domin_l = list()
re_cover_ab_intro_l = list()
re_cover_ab_nati_l = list()

for (z in 1:length(unique(sp_cover$f_p))) {
    #z = 1
    cover = sp_cover_fp[[z]]
    age = unique(cover$Age)
    if (max(age) %in% seq(2,1000,by=2) & min(age) %in% seq(2,1000,by=2)){
    re_age = seq(min(age), max(age), by = 2)} else if (max(age) %in% seq(2,1000,by=2) & min(age) %in% seq(1,1000,by=2)){
    re_age = seq(min(age)+1, max(age), by = 2) } else if (max(age) %in% seq(1,1000,by=2) & min(age) %in% seq(2,1000,by=2)) {
      re_age = seq(min(age)+1, max(age), by = 2) } else if (max(age) %in% seq(1,1000,by=2) & min(age) %in% seq(1,1000,by=2)) {
        re_age = seq(min(age), max(age), by = 2)}
    
    age_dif = setdiff(re_age, intersect(re_age, age))
    if (length(age_dif) == 0){re_cover = cover %>% filter(Age %in% re_age)
    } else { 
      re_cover_1 = cover %>% add_row(ID = NA, Absolute_year = NA,Year = NA, Age = age_dif[1], cover[1,5:7],
                                   (filter(cover, age == age_dif[1]-1)[,8:ncol(cover)]+
                                      filter(cover, age == age_dif[1]+1)[,8:ncol(cover)])/2)
      if (length(age_dif) > 1){
      for (i in 2:length(age_dif)) {
      add.row = data.frame(ID = NA, Absolute_year = NA, Year = NA, Age = age_dif[i], cover[1,5:7],
                   (filter(cover, age == age_dif[i]-1)[,8:ncol(cover)]+
                      filter(cover, age == age_dif[i]+1)[,8:ncol(cover)])/2)
      re_cover_1 = rbind(re_cover_1, add.row)
      }
      }
      re_cover = re_cover_1 %>% filter(Age %in% re_age)
    }
    re_cover = arrange(re_cover, re_cover$Age)
    if (length(setdiff(re_age, re_cover$Age)) == 0){print(paste('Filter the equal age suceesfully', z, '!',sep = ''))}
    
    re_cover_ab = cbind(re_cover[,1:7], re_cover[,8:(ncol(re_cover))][colSums(re_cover[,8:(ncol(re_cover))]) != 0])
    re_cover_prop = cbind(re_cover_ab[,1:7], t(apply(re_cover_ab[,8:(ncol(re_cover_ab))], 1, function(x){x/sum(x)})))
    re_cover_prop_inva = cbind(re_cover_prop[,1:7],
                            re_cover_prop %>% select(intersect(colnames(re_cover_prop)[8:ncol(re_cover_prop)], unlist(trait %>% filter(Origin == 'Exotic') %>% select(Species)))))
    re_cover_prop_nati = cbind(re_cover_prop[,1:7],
                            re_cover_prop %>% select(intersect(colnames(re_cover_prop)[8:ncol(re_cover_prop)], unlist(trait %>% filter(Origin == 'Native') %>% select(Species)))))
    
    
    cover_ab = cbind(cover[,1:7], cover[,8:(ncol(cover))][colSums(cover[,8:(ncol(cover))]) != 0])
    cover_prop = cbind(cover_ab[,1:7], t(apply(cover_ab[,8:(ncol(cover_ab))], 1, function(x){x/sum(x)})))

    cover_prop_inva = cbind(cover_prop[,1:7],
                            cover_prop %>% select(intersect(colnames(cover_prop)[8:ncol(cover_prop)],
                                                            unlist(trait %>% filter(Origin == 'Exotic') %>% select(Species)))))
    cover_prop_nati = cbind(cover_prop[,1:7],
                            cover_prop %>% select(intersect(colnames(cover_prop)[8:ncol(cover_prop)],
                                                            unlist(trait %>% filter(Origin == 'Native') %>% select(Species)))))
    
    re_cover_ab_inva = cbind(re_cover_ab[,1:7],
                             re_cover_ab %>% select(intersect(colnames(re_cover_ab)[8:ncol(re_cover_ab)],
                                                              unlist(trait %>% filter(Origin == 'Exotic') %>% select(Species)))))
    re_cover_ab_nati = cbind(re_cover_ab[,1:7],
                             re_cover_ab %>% select(intersect(colnames(re_cover_ab)[8:ncol(re_cover_ab)],
                                                              unlist(trait %>% filter(Origin == 'Native') %>% select(Species)))))
    
    #svg(filename = paste('D:/BSS study/Figures/', unique(sp_cover$f_p)[i], '.svg', sep = ''), width = 8, height = 5)
    #op = par(mfrow = c(1,2))
    #plot(cover_prop_inva$Age, rowSums(cover_prop_inva[,8:ncol(cover_prop_inva)]),main = "", xlab = "Age", ylab = 'Inva',
    #     las = 1, xpd = NA, cex.lab = 1.5)
    #plot(cover_prop_nati$Age, rowSums(cover_prop_nati[,8:ncol(cover_prop_nati)]), main = "", xlab = "Age", ylab = 'Nati',
    #     las = 1, xpd = NA, cex.lab = 1.5)
    #dev.off()
    estab_sp = vector() ## species at establishment stage (occur more than 10 successive years)
    domin_sp = vector() ## species at dominant stage (occur more than 10 successive years and mean cover is beyond any natives)
    for (j in 1:length(unique(cover_prop_inva$Absolute_year))) {
      #j = 2
      y_now = unique(cover_prop_inva$Absolute_year)[j]
      x = cover_prop_inva %>% filter(Absolute_year %in% (y_now:(y_now+9)))
      x = arrange(x, x$Absolute_year)
      year_dif = max(x$Absolute_year)-min(x$Absolute_year)
      if (year_dif %in% c(6:8)){
        x = cover_prop_inva %>% filter(Absolute_year %in% (y_now:(y_now+10)))
        x = arrange(x, x$Absolute_year)
        year_dif = max(x$Absolute_year)-min(x$Absolute_year)
        y = data.frame(freq = apply(x[,8:ncol(x)], 2, function(y){length(which(y != 0))})) %>% filter(freq > (nrow(x)-1))
                                                                                                
      } else if (year_dif == 9) {
      y = data.frame(freq = apply(x[,8:ncol(x)], 2, function(y){length(which(y != 0))})) %>% filter(freq > (nrow(x)-1))
      } else if (year_dif > 9) {
      y = data.frame(freq = apply(x[,8:ncol(x)], 2, function(y){length(which(y != 0))})) %>% filter(freq > (nrow(x)-2))  
      }
      
      if(nrow(y) == 0){print('no')} else{ 
      y_1 = cbind(x[,1:7],x %>% select(row.names(y)))
      x_nati = cover_prop_nati %>% filter(Absolute_year %in% y_1$Absolute_year)
      max_nati = max(apply(x_nati[,8:ncol(x_nati)], 2, function(x){mean(x)}))
      y_2 = data.frame(domin = apply(data.frame(y_1[,8:ncol(y_1)]), 2, function(y){mean(y) > max_nati})) %>% filter(domin == T)
      
      if(ncol(y_1) == 8 & nrow(y_2) == 1){
        row.names(y_2) = colnames(y_1)[8]} else if (nrow(y_2) >= 1) {row.names(y_2) = row.names(y_2)}
      
      domin_sp = c(domin_sp, row.names(y_2))
      estab_sp = c(estab_sp, row.names(y))
      }
    }
    unique(estab_sp)
    unique(domin_sp)
    cover_prop_estab = cbind(cover_prop[,1:7], cover_prop %>% select(unique(estab_sp)))
    cover_prop_domin = cbind(cover_prop[,1:7], cover_prop %>% select(unique(domin_sp)))
    cover_prop_intro = cbind(cover_prop[,1:7], cover_prop %>% select(setdiff(colnames(cover_prop_inva)[8:ncol(cover_prop_inva)], unique(estab_sp))))
    
    re_cover_prop_estab = cbind(re_cover_prop[,1:7], re_cover_prop %>% select(unique(estab_sp)))
    re_cover_prop_domin = cbind(re_cover_prop[,1:7], re_cover_prop %>% select(unique(domin_sp)))
    re_cover_prop_intro = cbind(re_cover_prop[,1:7], re_cover_prop %>% select(setdiff(colnames(re_cover_prop_inva)[8:ncol(re_cover_prop_inva)], unique(estab_sp))))
    
    re_cover_ab_estab = cbind(re_cover_ab[,1:7], re_cover_ab %>% select(unique(estab_sp)))
    re_cover_ab_domin = cbind(re_cover_ab[,1:7], re_cover_ab %>% select(unique(domin_sp)))
    re_cover_ab_intro = cbind(re_cover_ab[,1:7], re_cover_ab %>% select(setdiff(colnames(re_cover_ab_inva)[8:ncol(re_cover_ab_inva)], unique(estab_sp))))
    
    
    cover_fp_c[[z]] = cover_ab
    cover_fp_prop[[z]] = cover_prop
    cover_prop_estab_l[[z]] = cover_prop_estab
    cover_prop_domin_l[[z]] = cover_prop_domin
    cover_prop_intro_l[[z]] = cover_prop_intro
    cover_prop_nati_l[[z]] = cover_prop_nati
    
    re_cover_fp_c[[z]] = re_cover_ab
    re_cover_fp_prop[[z]] = re_cover_prop
    re_cover_prop_estab_l[[z]] = re_cover_prop_estab
    re_cover_prop_domin_l[[z]] = re_cover_prop_domin
    re_cover_prop_intro_l[[z]] = re_cover_prop_intro
    re_cover_prop_nati_l[[z]] = re_cover_prop_nati
    
    re_cover_ab_estab_l[[z]] = re_cover_ab_estab
    re_cover_ab_domin_l[[z]] = re_cover_ab_domin
    re_cover_ab_intro_l[[z]] = re_cover_ab_intro
    re_cover_ab_nati_l[[z]] = re_cover_ab_nati
    
}

save.image(file = "D:/BSS study/code/invasion_orginaldata.RData")

load("code/invasion_orginaldata.RData")

source('D:/BSS study/code/data collation_function.r')

re_cover_prop_estab_l_f = list()
re_cover_prop_intro_l_f = list()
re_cover_prop_nati_l_f = list()
re_cover_prop_domin_l_f = list()

for (i in 1:length(re_cover_prop_estab_l)){
  #i = 313
  x = re_cover_prop_estab_l[[i]]
  if (ncol(x) < 9){
    x = cbind(x, fake_sp = 0, fake_sp2 = 0)
  }
  y = x[,8:ncol(x)]
  z = filter_comm_mean(y, per = 0.3, rare = F)
  if(is.null(z)){re_cover_prop_estab_l_f[[i]] = NULL
  sr_estab = 0
  } else {
    #  colnames(z)[ncol(z)] = 'estab_rare'
    f = cbind(x[,1:7], z)
    re_cover_prop_estab_l_f[[i]] = f
    sr_estab = ncol(f)-7}
  
  x_1 = re_cover_prop_intro_l[[i]]
  if (ncol(x_1) < 9){
    x_1 = cbind(x_1, fake_sp = 0, fake_sp2 = 0)
  }
  y_1 = x_1[,8:ncol(x_1)]
  z_1 = filter_comm_mean(y_1, per = 0.15, rare = F)
  if(is.null(z_1)){
    re_cover_prop_intro_l_f[[i]] = NULL
    sr_intro = 0} else {
      #   colnames(z_1)[ncol(z_1)] = 'intro_rare'
      f_1 = cbind(x_1[,1:7], z_1)
      re_cover_prop_intro_l_f[[i]] = f_1
      sr_intro = ncol(f_1)-7}
  
  x_2 = re_cover_prop_nati_l[[i]]
  if (ncol(x_2) < 9){
    x_2 = cbind(x_2, fake_sp = 0, fake_sp2 = 0)
  }
  y_2 = x_2[,8:ncol(x_2)]
  z_2 = filter_comm_mean(x = y_2, per = 0.3, rare = F)
  if(is.null(z_2)){re_cover_prop_nati_l_f[[i]] = NULL
  sr_nati = 0
  } else {
    #  colnames(z_2)[ncol(z_2)] = 'nati_rare'
    f_2 = cbind(x_2[,1:7], z_2)
    re_cover_prop_nati_l_f[[i]] = f_2
    sr_nati = ncol(f_2)-7}
  
  x_3 = re_cover_prop_domin_l[[i]]
  if (ncol(x_3) < 9){
    x_3 = cbind(x_3, fake_sp = 0, fake_sp2 = 0)
  }
  y_3 = x_3[,8:ncol(x_3)]
  z_3 = filter_comm_mean(x = y_3, per = 0.3, rare = F)
  if(is.null(z_3)){re_cover_prop_domin_l_f[[i]] = NULL
  sr_domin = 0
  } else {
    #  colnames(z_3)[ncol(z_3)] = 'domin_rare'
    f_3 = cbind(x_3[,1:7], z_3)
    re_cover_prop_domin_l_f[[i]] = f_3
    sr_domin = ncol(f_3)-7}
  
}
 
re_cover_ab_estab_l_f = list()
re_cover_ab_intro_l_f = list()
re_cover_ab_nati_l_f = list()
re_cover_ab_domin_l_f = list()
sr_fit = data.frame()
for (i in 1:length(re_cover_ab_estab_l)){
  #i = 313
  x = re_cover_ab_estab_l[[i]]
  if (ncol(x) < 9){
    x = cbind(x, fake_sp = 0, fake_sp2 = 0)
  }
  y = x[,8:ncol(x)]
  z = filter_comm_mean(y, per = 0.3, rare = F)
  if(is.null(z)){re_cover_ab_estab_l_f[[i]] = NULL
  sr_estab = 0
  } else {
  #  colnames(z)[ncol(z)] = 'estab_rare'
    f = cbind(x[,1:7], z)
    re_cover_ab_estab_l_f[[i]] = f
    sr_estab = ncol(f)-7}
  
  
  x_1 = re_cover_ab_intro_l[[i]]
  if (ncol(x_1) < 9){
    x_1 = cbind(x_1, fake_sp = 0, fake_sp2 = 0)
  }
  y_1 = x_1[,8:ncol(x_1)]
  z_1 = filter_comm_mean(y_1, per = 0.15, rare = F)
  if(is.null(z_1)){
    re_cover_ab_intro_l_f[[i]] = NULL
    sr_intro = 0} else {
   #   colnames(z_1)[ncol(z_1)] = 'intro_rare'
      f_1 = cbind(x_1[,1:7], z_1)
      re_cover_ab_intro_l_f[[i]] = f_1
      sr_intro = ncol(f_1)-7}
  
  
  x_2 = re_cover_ab_nati_l[[i]]
  if (ncol(x_2) < 9){
    x_2 = cbind(x_2, fake_sp = 0, fake_sp2 = 0)
  }
  y_2 = x_2[,8:ncol(x_2)]
  z_2 = filter_comm_mean(x = y_2, per = 0.3, rare = F)
  if(is.null(z_2)){re_cover_ab_nati_l_f[[i]] = NULL
  sr_nati = 0
  } else {
  #  colnames(z_2)[ncol(z_2)] = 'nati_rare'
    f_2 = cbind(x_2[,1:7], z_2)
    re_cover_ab_nati_l_f[[i]] = f_2
    sr_nati = ncol(f_2)-7}
    
    x_3 = re_cover_ab_domin_l[[i]]
    if (ncol(x_3) < 9){
      x_3 = cbind(x_3, fake_sp = 0, fake_sp2 = 0)
    }
    y_3 = x_3[,8:ncol(x_3)]
    z_3 = filter_comm_mean(x = y_3, per = 0.3, rare = F)
    if(is.null(z_3)){re_cover_ab_domin_l_f[[i]] = NULL
    sr_domin = 0
    } else {
    #  colnames(z_3)[ncol(z_3)] = 'domin_rare'
    f_3 = cbind(x_3[,1:7], z_3)
    re_cover_ab_domin_l_f[[i]] = f_3
    sr_domin = ncol(f_3)-7}
  sr_1 = data.frame(unique(x$f_p), sr_all = sr_nati + sr_intro + sr_estab,
                    sr_nati = sr_nati,
                    sr_intro = sr_intro, sr_estab = sr_estab,
                    sr_domin = sr_domin)
  sr_fit = rbind(sr_1, sr_fit)
    
}

# Species proportion
library(data.table)

cover_mprop_all = list()

for (z in 1:length(unique(sp_cover$f_p))) {
  # z = 1
  cover = sp_cover_fp[[z]]
  cover_ab = cbind(cover[,1:7], cover[,8:(ncol(cover))][colSums(cover[,8:(ncol(cover))]) != 0])
  cover_prop = cbind(cover_ab[,1:7], t(apply(cover_ab[,8:(ncol(cover_ab))], 1, function(x){x/sum(x)})))
  cover_mprop = cbind(cover_prop[1,c(1, 5:7)], t(data.frame(colMeans(cover_prop[,8:ncol(cover_prop)]))))
  cover_mprop_all[[z]] = cover_mprop
}

cover_mprop_all = rbindlist(cover_mprop_all, fill = T)

for (z in 1:length(unique(sp_cover$f_p))) {
  #z = 218
  ra_Nati = re_cover_prop_nati_l_f[[z]]
  if (length(ra_Nati) == 0) {
    sp_Nati = NULL  
  } else {
    sp_Nati = colnames(ra_Nati)[8:ncol(ra_Nati)]
  }
  
  ra_Estab = re_cover_prop_estab_l_f[[z]]
  if (length(ra_Estab) == 0) {
    sp_Estab = NULL  
    ra_Estab = cbind(ra_Nati[,1:7], fake_sp1 = NA, fake_sp2 = NA)
  } else {
    sp_Estab = colnames(ra_Estab)[8:ncol(ra_Estab)]
  }
  
  ra_Intro = re_cover_prop_intro_l_f[[z]]
  if (length(ra_Intro) == 0) {
    sp_Intro = NULL  
    ra_Intro = cbind(ra_Nati[,1:7], fake_sp1 = NA, fake_sp2 = NA)
  } else {
    sp_Intro = colnames(ra_Intro)[8:ncol(ra_Intro)]
  }
  
 
  sp_name = c(sp_Estab, sp_Intro, sp_Nati)
  ra = cbind(ra_Estab, ra_Intro[,8:ncol(ra_Intro)],
             ra_Nati[,8:ncol(ra_Nati)])
  ra = ra %>% select(colnames(ra)[which(colnames(ra)!= 'fake_sp1' & colnames(ra) != 'fake_sp2')])
  colnames(ra)[8:ncol(ra)] = sp_name
  
  cover_mprop = cbind(ra[1,c(1, 5:7)], t(data.frame(colMeans(ra[,8:ncol(ra)]))))
  cover_mprop_all[[z]] = cover_mprop
}

cover_mprop_all = rbindlist(cover_mprop_all, fill = T)


save.image(file = "D:/BSS study/code/invasion_stage.RData")
