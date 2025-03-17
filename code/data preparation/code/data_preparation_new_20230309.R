###spatio-temporal data
###data sorting
library(dplyr)

##### Data loading
sp_cover = read.csv('data/original data/BSS_community_332.csv',
                    header = T)
field = read.csv('data/original data/FIELDS.csv',
                 header = T)
trait = read.csv('data/original data/traits332.csv',
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
sp_cover_year = split(sp_cover, sp_cover$Absolute_year)
sp_cover_age = split(sp_cover, sp_cover$Age)

sapply(sp_cover_fp, function(x) {sort(unique(x$Absolute_year))})
sapply(sp_cover_year, function(x) {sort(unique(x$f_p))})

##### Follow Li et al., (2015), condense the two year (half of field survey) to
##### one year and Follow Yin et al., (2022),, our study restricted to 3-51 age
sp_cover_f1 = sp_cover %>% filter(Absolute_year %in% c(1958:2017))
sp_cover_f1$fake_year = sp_cover_f1$Absolute_year
sp_cover_f1$fake_age = sp_cover_f1$Age
full_f_p = unique(sp_cover_f1$f_p)
for (i in seq(1959, 2017, 2)){
  #i = 1991
  f_p_year_2 = unique(sp_cover_f1[sp_cover_f1$Absolute_year == i-1,]$f_p)
  diff_f_p_year_2 = setdiff(full_f_p, f_p_year_2)
  row_diff = nrow(sp_cover_f1[sp_cover_f1$Absolute_year == i & 
                                sp_cover_f1$f_p %in% diff_f_p_year_2,
  ])
  if (row_diff == 0) {print(paste('cannot make it up', 'for',
                                  i-1, sep = ' '))}
  else {
    sp_cover_f1[sp_cover_f1$Absolute_year == i & sp_cover_f1$f_p %in% 
                  diff_f_p_year_2,]$fake_year = i-1
    age = unique(sp_cover_f1[sp_cover_f1$Absolute_year == i & sp_cover_f1$f_p %in% 
                               diff_f_p_year_2,]$Age)
    sp_cover_f1[sp_cover_f1$Absolute_year == i & sp_cover_f1$f_p %in%
                  diff_f_p_year_2,]$fake_age = rep(age-1,
                                                   each = row_diff/length(age-1))
  }
}

sp_cover_f2 = sp_cover_f1 %>% filter(fake_age %in% seq(3, 52, 2))
sp_cover_f2_y = split(sp_cover_f2, sp_cover_f2$fake_age)
sapply(sp_cover_f2_y, function(x) {sort(unique(x$f_p))})

sp_cover_f2_l = split(sp_cover_f2, sp_cover_f2$f_p)

# prepare for fit 
cover_all_ab = list()
cover_all_prop = list()
cover_prop_estab_l = list()
cover_prop_domin_l = list()
cover_prop_intro_l = list()
cover_prop_nati_l = list()

re_cover_all_prop = list()
re_cover_prop_estab_l = list()
re_cover_prop_domin_l = list()
re_cover_prop_intro_l = list()
re_cover_prop_nati_l = list()

re_cover_all_ab = list()
re_cover_ab_estab_l = list()
re_cover_ab_domin_l = list()
re_cover_ab_intro_l = list()
re_cover_ab_nati_l = list()

for (z in 1:length(unique(sp_cover$f_p))) {
  #z = 7
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
  
  ### filter the different invasion stage species
  estab_sp = vector() ## species at establishment stage (occur more than 10 successive years)
  domin_sp = vector() ## species at dominant stage (occur more than 10 successive years and mean cover is beyond any natives)
  
  for (j in 1:(ncol(cover_prop_inva)-7)) {
    #j = 15
    
    focal_alien = cover_prop_inva[,c(2, 7+j)]
    sp_name = colnames(focal_alien)[2]
    colnames(focal_alien)[2] = 'focal_sp'
    ab_rich_no0 = focal_alien %>% filter(focal_sp != 0)
    
    standard_seq = unique(focal_alien$Absolute_year)
    
    ddd = standard_seq %in% ab_rich_no0$Absolute_year
    sep_t_l = get_consec_rep(ddd)
    for (k in 1:length(sep_t_l)) {
      #k = 1
      year = standard_seq[sep_t_l[[k]]]
      ab_rich_rep = ab_rich_no0[ab_rich_no0$Absolute_year %in% year,]
      range_t = range(ab_rich_rep$Absolute_year)
      range_t1 = range_t[2]-range_t[1]
      if (range_t1 >= 10) {
        estab_sp[j] = sp_name
        x_nati = cover_prop_nati %>% filter(Absolute_year %in% ab_rich_rep$Absolute_year)
        max_nati = max(apply(x_nati[,8:ncol(x_nati)], 2, function(x){mean(x)}))
        mena_alien = mean(ab_rich_rep$focal_sp)
        
        if (mena_alien > max_nati) {
          domin_sp[j] = sp_name
        }
      }
      
    }
  }
  
  estab_sp = estab_sp[!is.na(estab_sp)]
  domin_sp = domin_sp[!is.na(domin_sp)]
  
  if (length(domin_sp) == 0) {
    cover_prop_domin = cbind(cover_prop[,1:7])
    re_cover_prop_domin = cbind(re_cover_prop[,1:7])
    re_cover_ab_domin = cbind(re_cover_ab[,1:7])
  } else {
    cover_prop_domin = cbind(cover_prop[,1:7],
                             cover_prop %>% select(domin_sp))
    re_cover_prop_domin = cbind(re_cover_prop[,1:7],
                                re_cover_prop %>% select(domin_sp))
    re_cover_ab_domin = cbind(re_cover_ab[,1:7],
                              re_cover_ab %>% select(domin_sp))
  }
  
  if (length(estab_sp) == 0) {
    cover_prop_estab = cbind(cover_prop[,1:7])
    re_cover_prop_estab = cbind(re_cover_prop[,1:7])
    re_cover_ab_estab = cbind(re_cover_ab[,1:7])
  } else {
    cover_prop_estab = cbind(cover_prop[,1:7],
                             cover_prop %>% select(estab_sp))
    re_cover_prop_estab = cbind(re_cover_prop[,1:7],
                                re_cover_prop %>% select(estab_sp))
    re_cover_ab_estab = cbind(re_cover_ab[,1:7],
                              re_cover_ab %>% select(estab_sp))
  }
  
  cover_prop_intro = cbind(cover_prop[,1:7],
                           cover_prop %>%
                             select(setdiff(colnames(cover_prop_inva)[8:ncol(cover_prop_inva)],
                                            estab_sp)))
  re_cover_prop_intro = cbind(re_cover_prop[,1:7],
                              re_cover_prop %>%
                                select(setdiff(colnames(re_cover_prop_inva)[8:ncol(re_cover_prop_inva)], estab_sp)))
  re_cover_ab_intro = cbind(re_cover_ab[,1:7],
                            re_cover_ab %>%
                              select(setdiff(colnames(re_cover_ab_inva)[8:ncol(re_cover_ab_inva)], estab_sp)))
  
  ## raw absoulte cover data
  cover_all_ab[[z]] = cover_ab
  
  ## raw propotion cover data
  cover_all_prop[[z]] = cover_prop
  cover_prop_estab_l[[z]] = cover_prop_estab
  cover_prop_domin_l[[z]] = cover_prop_domin
  cover_prop_intro_l[[z]] = cover_prop_intro
  cover_prop_nati_l[[z]] = cover_prop_nati
  
  ## fit propotion cover data
  re_cover_all_prop[[z]] = re_cover_prop
  re_cover_prop_estab_l[[z]] = re_cover_prop_estab
  re_cover_prop_domin_l[[z]] = re_cover_prop_domin
  re_cover_prop_intro_l[[z]] = re_cover_prop_intro
  re_cover_prop_nati_l[[z]] = re_cover_prop_nati
  
  ## fit absoulte cover data
  re_cover_all_ab[[z]] = re_cover_ab
  re_cover_ab_estab_l[[z]] = re_cover_ab_estab
  re_cover_ab_domin_l[[z]] = re_cover_ab_domin
  re_cover_ab_intro_l[[z]] = re_cover_ab_intro
  re_cover_ab_nati_l[[z]] = re_cover_ab_nati
  
}

save.image(file = "code/invasion_orginaldata.RData")

load("code/invasion_orginaldata.RData")

source('code/data collation_function.r')
