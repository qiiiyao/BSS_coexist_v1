fit_prepa_top = function(w, t = NULL, posi_sp_start = 10) {
  #w = sp_cover_f2_top50_ages1_35_merge_all
  #t = 2:14
  #t = NULL
  #w = sp_cover_fp[[1]]
  #nati_per = 0.3
  #intro_per = 0.3 
  #estab_per = 0.3
  #domin_per = 0.3
  ### filter the different invasion stage species
  no_sp_num = posi_sp_start-1
  if (is.null(t)) { 
    re_cover_ab = cbind(w[,1:no_sp_num],
                        w[,posi_sp_start:(ncol(w))][colSums(w[,
                                                   posi_sp_start:(ncol(w))]) != 0])
  } else { 
    t_seq = seq_t_all[t]
    w_t = w[w$fake_age %in% t_seq,]
    re_cover_ab = cbind(w_t[,1:no_sp_num],
                        w_t[,posi_sp_start:(ncol(w_t))][colSums(w_t[,
                                                         posi_sp_start:(ncol(w_t))]) != 0])
  }
  re_cover_ab = as.data.frame(re_cover_ab) ## ungrouped grouped dataframe 
  re_cover_ab_inva = cbind(re_cover_ab[,1:no_sp_num],
                           re_cover_ab %>%
                             select(all_of(intersect(colnames(re_cover_ab)[posi_sp_start:ncol(re_cover_ab)],
                                              unlist(trait %>% filter(Origin == 'Exotic') %>% select(Species))))))
  re_cover_ab_nati = cbind(re_cover_ab[,1:no_sp_num],
                           re_cover_ab %>%
                             select(all_of(intersect(colnames(re_cover_ab)[posi_sp_start:ncol(re_cover_ab)],
                                              unlist(trait %>% filter(Origin == 'Native') %>% select(Species))))))
  
  re_cover_prop = cbind(re_cover_ab[,1:no_sp_num],
                        t(apply(re_cover_ab[,posi_sp_start:(ncol(re_cover_ab))], 1, function(x){x/sum(x)})))
  re_cover_prop_inva = cbind(re_cover_prop[,1:no_sp_num],
                             re_cover_prop %>% 
                               select(all_of(intersect(colnames(re_cover_prop)[posi_sp_start:ncol(re_cover_prop)],
                                                unlist(trait %>% filter(Origin == 'Exotic') %>% select(Species))))))
  re_cover_prop_nati = cbind(re_cover_prop[,1:no_sp_num],
                             re_cover_prop %>%
                               select(all_of(intersect(colnames(re_cover_prop)[posi_sp_start:ncol(re_cover_prop)],
                                                unlist(trait %>% filter(Origin == 'Native') %>% select(Species))))))
  
  estab_sp = vector() ## species at establishment stage (occur more than
  #10 successive years)
  domin_sp = vector() ## species at dominant stage (occur more than 10
  #successive years and mean cover is beyond any natives)
  
  for (j in 1:(ncol(re_cover_prop_inva)-no_sp_num)) {
    #j = 15
    
    focal_alien = re_cover_prop_inva[,c(2, no_sp_num+j)]
    sp_name = colnames(focal_alien)[2]
    colnames(focal_alien)[2] = 'focal_sp'
    ab_rich_no0 = focal_alien %>% filter(focal_sp != 0)
    
    standard_seq = unique(focal_alien$Absolute_year)
    
    ddd = standard_seq %in% ab_rich_no0$Absolute_year
    sep_t_l = get_consec_rep(ddd)
    for (k in 1:length(sep_t_l)) {
      #k = 1
      year = standard_seq[sep_t_l[[k]]]
      ab_rich_rep = as.data.frame(ab_rich_no0[ab_rich_no0$Absolute_year %in% year,])
      range_t = range(ab_rich_rep$Absolute_year)
      range_t1 = range_t[2]-range_t[1]
      if (range_t1 >= 10) {
        estab_sp[j] = sp_name
        x_nati = as.data.frame(re_cover_prop_nati %>% filter(Absolute_year %in% ab_rich_rep$Absolute_year))
        max_nati = max(apply(x_nati[,posi_sp_start:ncol(x_nati)], 2, function(x){mean(x)}))
        mean_alien = mean(ab_rich_rep$focal_sp)
        
        if (mean_alien > max_nati) {
          domin_sp[j] = sp_name
        }
      }
      
    }
  }
  
  #### Get the stage information
  estab_sp = estab_sp[!is.na(estab_sp)]
  domin_sp = domin_sp[!is.na(domin_sp)]
  estab_spname = paste(estab_sp, collapse = ', ')
  domin_spname = paste(domin_sp, collapse = ', ')
  ####
  
  
  if (length(domin_sp) == 0) {
    re_cover_prop_domin = cbind(re_cover_prop[,1:no_sp_num])
    re_cover_ab_domin = cbind(re_cover_ab[,1:no_sp_num])
  } else {
    re_cover_prop_domin = cbind(re_cover_prop[,1:no_sp_num],
                                re_cover_prop %>% select(all_of(domin_sp)))
    re_cover_ab_domin = cbind(re_cover_ab[,1:no_sp_num],
                              re_cover_ab %>% select(all_of(domin_sp)))
  }
  
  if (length(estab_sp) == 0) {
    re_cover_prop_estab = cbind(re_cover_prop[,1:no_sp_num])
    re_cover_ab_estab = cbind(re_cover_ab[,1:no_sp_num])
  } else {
    re_cover_prop_estab = cbind(re_cover_prop[,1:no_sp_num],
                                re_cover_prop %>% select(all_of(estab_sp)))
    re_cover_ab_estab = cbind(re_cover_ab[,1:no_sp_num],
                              re_cover_ab %>% select(all_of(estab_sp)))
  }
  re_cover_prop_intro = cbind(re_cover_prop[,1:no_sp_num],
                              re_cover_prop %>%
                                select(all_of(setdiff(colnames(re_cover_prop_inva)[posi_sp_start:ncol(re_cover_prop_inva)], estab_sp))))
  re_cover_ab_intro = cbind(re_cover_ab[,1:no_sp_num],
                            re_cover_ab %>%
                              select(all_of(setdiff(colnames(re_cover_ab_inva)[posi_sp_start:ncol(re_cover_ab_inva)], estab_sp))))
  
  intro_spname = paste(setdiff(colnames(re_cover_prop_inva)[posi_sp_start:ncol(re_cover_prop_inva)],
                               estab_sp),
                       collapse = ', ')
  
  stage_sp_name = cbind(f_p = unique(re_cover_prop_nati$f_p),
                        intro_sp = intro_spname,
                        estab_sp = estab_spname, 
                        domin_sp = domin_spname)
  
  
  ### Create the data list for stan fit 
  sp = colnames(re_cover_ab)[posi_sp_start:ncol(re_cover_ab)]
  re_cover_ab = arrange(re_cover_ab, re_cover_ab$fake_age)
  re_cover_ab[re_cover_ab==0] = 1e-06
  Pred = re_cover_ab[2:(nrow(re_cover_ab)),1:ncol(re_cover_ab)]
  Pred_t2 = re_cover_ab[1:(nrow(re_cover_ab)-1),1:ncol(re_cover_ab)]
  Yrs = nrow(Pred)
  yid = as.numeric(as.factor(Pred$fake_age))
  sp_matrix = re_cover_ab[1:(nrow(re_cover_ab)-1),posi_sp_start:ncol(re_cover_ab)]
  t_diff = Pred$Age-Pred_t2$Age
  
  datalist = list()
  for (l in 1:length(sp)) {
    
    Y = unlist(Pred %>% select(all_of(sp[l])))
    I = unlist(sp_matrix %>% select(all_of(sp[l])))
    G = Y/I
    
    datalist[[l]] = list(N = Yrs, P = ncol(sp_matrix), ftime=Yrs, time=yid,
                         time_diff = t_diff,
                         G = G, X=sp_matrix) }
  
  ### Create a list to aggregate filter results
  filter_inf = list(stage_sp_name = stage_sp_name,
                    re_cover_ab = re_cover_ab,
                    datalist
  )
  
  return(filter_inf)
  
}



fit_prepa_top_prop = function(w, t = NULL, posi_sp_start = 10) {
  #w = sp_cover_f2_top50_ages1_35_merge_all
  #t = 2:14
  #t = NULL
  #w = sp_cover_fp[[1]]
  #nati_per = 0.3
  #intro_per = 0.3 
  #estab_per = 0.3
  #domin_per = 0.3
  ### filter the different invasion stage species
  no_sp_num = posi_sp_start-1
  if (is.null(t)) { 
    re_cover_ab = cbind(w[,1:no_sp_num],
                        w[,posi_sp_start:(ncol(w))][colSums(w[,
                                                              posi_sp_start:(ncol(w))]) != 0])
  } else { 
    t_seq = seq_t_all[t]
    w_t = w[w$fake_age %in% t_seq,]
    re_cover_ab = cbind(w_t[,1:no_sp_num],
                        w_t[,posi_sp_start:(ncol(w_t))][colSums(w_t[,
                                                                    posi_sp_start:(ncol(w_t))]) != 0])
  }
  re_cover_ab = as.data.frame(re_cover_ab) ## ungrouped grouped dataframe 
  re_cover_ab_inva = cbind(re_cover_ab[,1:no_sp_num],
                           re_cover_ab %>%
                             select(all_of(intersect(colnames(re_cover_ab)[posi_sp_start:ncol(re_cover_ab)],
                                                     unlist(trait %>% filter(Origin == 'Exotic') %>% select(Species))))))
  re_cover_ab_nati = cbind(re_cover_ab[,1:no_sp_num],
                           re_cover_ab %>%
                             select(all_of(intersect(colnames(re_cover_ab)[posi_sp_start:ncol(re_cover_ab)],
                                                     unlist(trait %>% filter(Origin == 'Native') %>% select(Species))))))
  
  re_cover_prop = cbind(re_cover_ab[,1:no_sp_num],
                        t(apply(re_cover_ab[,posi_sp_start:(ncol(re_cover_ab))], 1, function(x){x/sum(x)})))
  re_cover_prop_inva = cbind(re_cover_prop[,1:no_sp_num],
                             re_cover_prop %>% 
                               select(all_of(intersect(colnames(re_cover_prop)[posi_sp_start:ncol(re_cover_prop)],
                                                       unlist(trait %>% filter(Origin == 'Exotic') %>% select(Species))))))
  re_cover_prop_nati = cbind(re_cover_prop[,1:no_sp_num],
                             re_cover_prop %>%
                               select(all_of(intersect(colnames(re_cover_prop)[posi_sp_start:ncol(re_cover_prop)],
                                                       unlist(trait %>% filter(Origin == 'Native') %>% select(Species))))))
  
  estab_sp = vector() ## species at establishment stage (occur more than
  #10 successive years)
  domin_sp = vector() ## species at dominant stage (occur more than 10
  #successive years and mean cover is beyond any natives)
  
  for (j in 1:(ncol(re_cover_prop_inva)-no_sp_num)) {
    #j = 15
    
    focal_alien = re_cover_prop_inva[,c(2, no_sp_num+j)]
    sp_name = colnames(focal_alien)[2]
    colnames(focal_alien)[2] = 'focal_sp'
    ab_rich_no0 = focal_alien %>% filter(focal_sp != 0)
    
    standard_seq = unique(focal_alien$Absolute_year)
    
    ddd = standard_seq %in% ab_rich_no0$Absolute_year
    sep_t_l = get_consec_rep(ddd)
    for (k in 1:length(sep_t_l)) {
      #k = 1
      year = standard_seq[sep_t_l[[k]]]
      ab_rich_rep = as.data.frame(ab_rich_no0[ab_rich_no0$Absolute_year %in% year,])
      range_t = range(ab_rich_rep$Absolute_year)
      range_t1 = range_t[2]-range_t[1]
      if (range_t1 >= 10) {
        estab_sp[j] = sp_name
        x_nati = as.data.frame(re_cover_prop_nati %>% filter(Absolute_year %in% ab_rich_rep$Absolute_year))
        max_nati = max(apply(x_nati[,posi_sp_start:ncol(x_nati)], 2, function(x){mean(x)}))
        mean_alien = mean(ab_rich_rep$focal_sp)
        
        if (mean_alien > max_nati) {
          domin_sp[j] = sp_name
        }
      }
      
    }
  }
  
  #### Get the stage information
  estab_sp = estab_sp[!is.na(estab_sp)]
  domin_sp = domin_sp[!is.na(domin_sp)]
  estab_spname = paste(estab_sp, collapse = ', ')
  domin_spname = paste(domin_sp, collapse = ', ')
  ####
  
  
  if (length(domin_sp) == 0) {
    re_cover_prop_domin = cbind(re_cover_prop[,1:no_sp_num])
    re_cover_ab_domin = cbind(re_cover_ab[,1:no_sp_num])
  } else {
    re_cover_prop_domin = cbind(re_cover_prop[,1:no_sp_num],
                                re_cover_prop %>% select(all_of(domin_sp)))
    re_cover_ab_domin = cbind(re_cover_ab[,1:no_sp_num],
                              re_cover_ab %>% select(all_of(domin_sp)))
  }
  
  if (length(estab_sp) == 0) {
    re_cover_prop_estab = cbind(re_cover_prop[,1:no_sp_num])
    re_cover_ab_estab = cbind(re_cover_ab[,1:no_sp_num])
  } else {
    re_cover_prop_estab = cbind(re_cover_prop[,1:no_sp_num],
                                re_cover_prop %>% select(all_of(estab_sp)))
    re_cover_ab_estab = cbind(re_cover_ab[,1:no_sp_num],
                              re_cover_ab %>% select(all_of(estab_sp)))
  }
  re_cover_prop_intro = cbind(re_cover_prop[,1:no_sp_num],
                              re_cover_prop %>%
                                select(all_of(setdiff(colnames(re_cover_prop_inva)[posi_sp_start:ncol(re_cover_prop_inva)], estab_sp))))
  re_cover_ab_intro = cbind(re_cover_ab[,1:no_sp_num],
                            re_cover_ab %>%
                              select(all_of(setdiff(colnames(re_cover_ab_inva)[posi_sp_start:ncol(re_cover_ab_inva)], estab_sp))))
  
  intro_spname = paste(setdiff(colnames(re_cover_prop_inva)[posi_sp_start:ncol(re_cover_prop_inva)],
                               estab_sp),
                       collapse = ', ')
  
  stage_sp_name = cbind(f_p = unique(re_cover_prop_nati$f_p),
                        intro_sp = intro_spname,
                        estab_sp = estab_spname, 
                        domin_sp = domin_spname)
  
  
  ### Create the data list for stan fit 
  sp = colnames(re_cover_ab)[posi_sp_start:ncol(re_cover_ab)]
  re_cover_ab = arrange(re_cover_ab, re_cover_ab$fake_age)
  re_cover_ab[re_cover_ab==0] = 1e-06
  re_cover_ab_bio = re_cover_ab[,posi_sp_start:ncol(re_cover_ab)]
  re_cover_prop_bio = re_cover_ab_bio/apply(re_cover_ab_bio, 1, sum)
  re_cover_prop = cbind(re_cover_ab[,1:(posi_sp_start-1)], re_cover_prop_bio)
  Pred_prop = re_cover_prop[2:(nrow(re_cover_prop)),1:ncol(re_cover_prop)]
  Pred_ab = re_cover_ab[2:(nrow(re_cover_ab)),1:ncol(re_cover_ab)]
  Pred_ab_t2 = re_cover_ab[1:(nrow(re_cover_ab)-1),1:ncol(re_cover_ab)]
  Yrs = nrow(Pred_ab)
  yid = as.numeric(as.factor(Pred_ab$fake_age))
  sp_matrix_prop = re_cover_prop[1:(nrow(re_cover_prop)-1),
                                 posi_sp_start:ncol(re_cover_prop)]
  sp_matrix_ab = re_cover_ab[1:(nrow(re_cover_ab)-1),
                             posi_sp_start:ncol(re_cover_ab)]
  t_diff = Pred_ab$fake_age-Pred_ab_t2$fake_age
  
  datalist_prop = list()
  datalist_ab = list()
  for (l in 1:length(sp)) {
    
    Y_prop = unlist(Pred_prop %>% select(all_of(sp[l])))
    I_prop = unlist(sp_matrix_prop %>% select(all_of(sp[l])))
    G_prop = Y_prop/I_prop
    datalist_prop[[l]] = list(N = Yrs, P = ncol(sp_matrix_prop), ftime=Yrs, time=yid,
                              time_diff = t_diff,
                              G = G_prop, X=sp_matrix_prop)
    
    Y_ab = unlist(Pred_ab %>% select(all_of(sp[l])))
    I_ab = unlist(sp_matrix_ab %>% select(all_of(sp[l])))
    G_ab = Y_ab/I_ab
    datalist_ab[[l]] = list(N = Yrs, P = ncol(sp_matrix_ab), ftime=Yrs, time=yid,
                            time_diff = t_diff,
                            G = G_ab, X = sp_matrix_ab)}
  
  ### Create a list to aggregate filter results
  filter_inf = list(stage_sp_name = stage_sp_name,
                    re_cover_ab = re_cover_ab,
                    fit_data_ab = datalist_ab,
                    fit_data_prop = datalist_prop
  )
  
  return(filter_inf)
  
}


fit_prepa_top_merge = function(w, t = NULL, posi_sp_start = 10) {
  #w = sp_cover_f1_ages1_35_field_top50_l[[1]]
  #t = 2:14
  #t = NULL
  #posi_sp_start = 3
  #w = sp_cover_fp[[1]]
  #nati_per = 0.3
  #intro_per = 0.3 
  #estab_per = 0.3
  #domin_per = 0.3
  ### filter the different invasion stage species
  no_sp_num = posi_sp_start-1
  if (is.null(t)) { 
    re_cover_ab = cbind(w[,1:no_sp_num],
                        w[,posi_sp_start:(ncol(w))][colSums(w[,
                                                              posi_sp_start:(ncol(w))]) != 0])
  } else { 
    t_seq = seq_t_all[t]
    w_t = w[w$fake_age %in% t_seq,]
    re_cover_ab = cbind(w_t[,1:no_sp_num],
                        w_t[,posi_sp_start:(ncol(w_t))][colSums(w_t[,
                                                                    posi_sp_start:(ncol(w_t))]) != 0])
  }
  re_cover_ab = as.data.frame(re_cover_ab) ## ungrouped grouped dataframe 
  fake_age = as.data.frame(re_cover_ab[,1:no_sp_num])
  #colnames(fake_age) = 'fake_age'
  re_cover_ab_inva = cbind(fake_age,
                           re_cover_ab %>%
                             select(all_of(intersect(colnames(re_cover_ab)[posi_sp_start:ncol(re_cover_ab)],
                                                     unlist(trait %>% filter(Origin == 'Exotic') %>% select(Species))))))
  re_cover_ab_nati = cbind(fake_age,
                           re_cover_ab %>%
                             select(all_of(intersect(colnames(re_cover_ab)[posi_sp_start:ncol(re_cover_ab)],
                                                     unlist(trait %>% filter(Origin == 'Native') %>% select(Species))))))
  
  re_cover_prop = cbind(fake_age,
                        t(apply(re_cover_ab[,posi_sp_start:(ncol(re_cover_ab))], 1, function(x){x/sum(x)})))
  re_cover_prop_inva = cbind(fake_age,
                             re_cover_prop %>% 
                               select(all_of(intersect(colnames(re_cover_prop)[posi_sp_start:ncol(re_cover_prop)],
                                                       unlist(trait %>% filter(Origin == 'Exotic') %>% select(Species))))))
  re_cover_prop_nati = cbind(fake_age,
                             re_cover_prop %>%
                               select(all_of(intersect(colnames(re_cover_prop)[posi_sp_start:ncol(re_cover_prop)],
                                                       unlist(trait %>% filter(Origin == 'Native') %>% select(Species))))))
  
  estab_sp = vector() ## species at establishment stage (occur more than
  #10 successive years)
  domin_sp = vector() ## species at dominant stage (occur more than 10
  #successive years and mean cover is beyond any natives)
  
  for (j in 1:(ncol(re_cover_prop_inva)-no_sp_num)) {
    #j = 15
    
    focal_alien = re_cover_prop_inva[,c(1, no_sp_num+j)]
    sp_name = colnames(focal_alien)[2]
    colnames(focal_alien)[2] = 'focal_sp'
    ab_rich_no0 = focal_alien %>% filter(focal_sp != 0)
    
    standard_seq = unique(focal_alien$fake_age)
    
    ddd = standard_seq %in% ab_rich_no0$fake_age
    sep_t_l = get_consec_rep(ddd)
    for (k in 1:length(sep_t_l)) {
      #k = 1
      year = standard_seq[sep_t_l[[k]]]
      ab_rich_rep = as.data.frame(ab_rich_no0[ab_rich_no0$fake_age %in% year,])
      range_t = range(ab_rich_rep$fake_age)
      range_t1 = range_t[2]-range_t[1]
      if (range_t1 >= 8) {
        estab_sp[j] = sp_name
        x_nati = as.data.frame(re_cover_prop_nati %>% filter(fake_age %in% ab_rich_rep$fake_age))
        max_nati = max(apply(x_nati[,posi_sp_start:ncol(x_nati)], 2, function(x){mean(x)}))
        mean_alien = mean(ab_rich_rep$focal_sp)
        
        if (mean_alien > max_nati) {
          domin_sp[j] = sp_name
        }
      }
      
    }
  }
  
  #### Get the stage information
  estab_sp = estab_sp[!is.na(estab_sp)]
  domin_sp = domin_sp[!is.na(domin_sp)]
  estab_spname = paste(estab_sp, collapse = ', ')
  domin_spname = paste(domin_sp, collapse = ', ')
  ####
  
  
  if (length(domin_sp) == 0) {
    re_cover_prop_domin = cbind(re_cover_prop[,1:no_sp_num])
    re_cover_ab_domin = cbind(re_cover_ab[,1:no_sp_num])
  } else {
    re_cover_prop_domin = cbind(re_cover_prop[,1:no_sp_num],
                                re_cover_prop %>% select(all_of(domin_sp)))
    re_cover_ab_domin = cbind(re_cover_ab[,1:no_sp_num],
                              re_cover_ab %>% select(all_of(domin_sp)))
  }
  
  if (length(estab_sp) == 0) {
    re_cover_prop_estab = cbind(re_cover_prop[,1:no_sp_num])
    re_cover_ab_estab = cbind(re_cover_ab[,1:no_sp_num])
  } else {
    re_cover_prop_estab = cbind(re_cover_prop[,1:no_sp_num],
                                re_cover_prop %>% select(all_of(estab_sp)))
    re_cover_ab_estab = cbind(re_cover_ab[,1:no_sp_num],
                              re_cover_ab %>% select(all_of(estab_sp)))
  }
  re_cover_prop_intro = cbind(re_cover_prop[,1:no_sp_num],
                              re_cover_prop %>%
                                select(all_of(setdiff(colnames(re_cover_prop_inva)[posi_sp_start:ncol(re_cover_prop_inva)], estab_sp))))
  re_cover_ab_intro = cbind(re_cover_ab[,1:no_sp_num],
                            re_cover_ab %>%
                              select(all_of(setdiff(colnames(re_cover_ab_inva)[posi_sp_start:ncol(re_cover_ab_inva)], estab_sp))))
  
  intro_spname = paste(setdiff(colnames(re_cover_prop_inva)[posi_sp_start:ncol(re_cover_prop_inva)],
                               estab_sp),
                       collapse = ', ')
  
  stage_sp_name = cbind(f_p = unique(re_cover_prop_nati$f_p),
                        intro_sp = intro_spname,
                        estab_sp = estab_spname, 
                        domin_sp = domin_spname)
  
  
  ### Create the data list for stan fit 
  sp = colnames(re_cover_ab)[posi_sp_start:ncol(re_cover_ab)]
  re_cover_ab = arrange(re_cover_ab, re_cover_ab$fake_age)
  re_cover_ab[re_cover_ab==0] = 1e-06
  re_cover_ab_bio = re_cover_ab[,posi_sp_start:ncol(re_cover_ab)]
  re_cover_prop_bio = re_cover_ab_bio/apply(re_cover_ab_bio, 1, sum)
  re_cover_prop = cbind(re_cover_ab[,1:(posi_sp_start-1)], re_cover_prop_bio)
  Pred_prop = re_cover_prop[2:(nrow(re_cover_prop)),1:ncol(re_cover_prop)]
  Pred_ab = re_cover_ab[2:(nrow(re_cover_ab)),1:ncol(re_cover_ab)]
  Pred_ab_t2 = re_cover_ab[1:(nrow(re_cover_ab)-1),1:ncol(re_cover_ab)]
  Yrs = nrow(Pred_ab)
  yid = as.numeric(as.factor(Pred_ab$fake_age))
  sp_matrix_prop = re_cover_prop[1:(nrow(re_cover_prop)-1),
                            posi_sp_start:ncol(re_cover_prop)]
  sp_matrix_ab = re_cover_ab[1:(nrow(re_cover_ab)-1),
                            posi_sp_start:ncol(re_cover_ab)]
  t_diff = Pred_ab$fake_age-Pred_ab_t2$fake_age
  
  datalist_prop = list()
  datalist_ab = list()
  for (l in 1:length(sp)) {
    
    Y_prop = unlist(Pred_prop %>% select(all_of(sp[l])))
    I_prop = unlist(sp_matrix_prop %>% select(all_of(sp[l])))
    G_prop = Y_prop/I_prop
    datalist_prop[[l]] = list(N = Yrs, P = ncol(sp_matrix_prop), ftime=Yrs, time=yid,
                         time_diff = t_diff,
                         G = G_prop, X=sp_matrix_prop)
    
    Y_ab = unlist(Pred_ab %>% select(all_of(sp[l])))
    I_ab = unlist(sp_matrix_ab %>% select(all_of(sp[l])))
    G_ab = Y_ab/I_ab
    datalist_ab[[l]] = list(N = Yrs, P = ncol(sp_matrix_ab), ftime=Yrs, time=yid,
                              time_diff = t_diff,
                              G = G_ab, X = sp_matrix_ab)}
  
    ### Create a list to aggregate filter results
    filter_inf = list(stage_sp_name = stage_sp_name,
                      re_cover_ab = re_cover_ab,
                      fit_data_ab = datalist_ab,
                      fit_data_prop = datalist_prop
  )
  
  return(filter_inf)
  
}
