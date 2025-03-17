fit_prepa_neighbor_group = function(w, t = NULL) {
  #w = sp_cover_f2_top50_early_suc_fp[[1]]
  #t = NULL
  #w = sp_cover_fp[[1]]
  #nati_per = 0.3
  #intro_per = 0.3 
  #estab_per = 0.3
  #domin_per = 0.3
  ### filter the different invasion stage species
  if (is.null(t)) { 
    re_cover_ab = cbind(w[,1:9],
                        w[,10:(ncol(w))][colSums(w[,
                                                   10:(ncol(w))]) != 0])
  } else { 
    t_seq = seq_t_all[t]
    w_t = w[w$fake_age %in% t_seq,]
    re_cover_ab = cbind(w_t[,1:9],
                        w_t[,10:(ncol(w_t))][colSums(w_t[,
                                                         10:(ncol(w_t))]) != 0])
  }
  re_cover_ab = as.data.frame(re_cover_ab) ## ungrouped grouped dataframe 
  re_cover_ab_inva = cbind(re_cover_ab[,1:9],
                           re_cover_ab %>%
                             select(all_of(intersect(colnames(re_cover_ab)[10:ncol(re_cover_ab)],
                                                     unlist(trait %>% filter(Origin == 'Exotic') %>% select(Species))))))
  re_cover_ab_nati = cbind(re_cover_ab[,1:9],
                           re_cover_ab %>%
                             select(all_of(intersect(colnames(re_cover_ab)[10:ncol(re_cover_ab)],
                                                     unlist(trait %>% filter(Origin == 'Native') %>% select(Species))))))
  
  re_cover_prop = cbind(re_cover_ab[,1:9],
                        t(apply(re_cover_ab[,10:(ncol(re_cover_ab))], 1, function(x){x/sum(x)})))
  re_cover_prop_inva = cbind(re_cover_prop[,1:9],
                             re_cover_prop %>% 
                               select(all_of(intersect(colnames(re_cover_prop)[10:ncol(re_cover_prop)],
                                                       unlist(trait %>% filter(Origin == 'Exotic') %>% select(Species))))))
  re_cover_prop_nati = cbind(re_cover_prop[,1:9],
                             re_cover_prop %>%
                               select(all_of(intersect(colnames(re_cover_prop)[10:ncol(re_cover_prop)],
                                                       unlist(trait %>% filter(Origin == 'Native') %>% select(Species))))))
  
  estab_sp = vector() ## species at establishment stage (occur more than
  #10 successive years)
  domin_sp = vector() ## species at dominant stage (occur more than 10
  #successive years and mean cover is beyond any natives)
  
  for (j in 1:(ncol(re_cover_prop_inva)-9)) {
    #j = 15
    
    focal_alien = re_cover_prop_inva[,c(2, 9+j)]
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
        max_nati = max(apply(x_nati[,10:ncol(x_nati)], 2, function(x){mean(x)}))
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
    re_cover_prop_domin = cbind(re_cover_prop[,1:9])
    re_cover_ab_domin = cbind(re_cover_ab[,1:9])
  } else {
    re_cover_prop_domin = cbind(re_cover_prop[,1:9],
                                re_cover_prop %>% select(all_of(domin_sp)))
    re_cover_ab_domin = cbind(re_cover_ab[,1:9],
                              re_cover_ab %>% select(all_of(domin_sp)))
  }
  
  if (length(estab_sp) == 0) {
    re_cover_prop_estab = cbind(re_cover_prop[,1:9])
    re_cover_ab_estab = cbind(re_cover_ab[,1:9])
  } else {
    re_cover_prop_estab = cbind(re_cover_prop[,1:9],
                                re_cover_prop %>% select(all_of(estab_sp)))
    re_cover_ab_estab = cbind(re_cover_ab[,1:9],
                              re_cover_ab %>% select(all_of(estab_sp)))
  }
  re_cover_prop_intro = cbind(re_cover_prop[,1:9],
                              re_cover_prop %>%
                                select(all_of(setdiff(colnames(re_cover_prop_inva)[10:ncol(re_cover_prop_inva)], estab_sp))))
  re_cover_ab_intro = cbind(re_cover_ab[,1:9],
                            re_cover_ab %>%
                              select(all_of(setdiff(colnames(re_cover_ab_inva)[10:ncol(re_cover_ab_inva)], estab_sp))))
  
  intro_spname = paste(setdiff(colnames(re_cover_prop_inva)[10:ncol(re_cover_prop_inva)],
                               estab_sp),
                       collapse = ', ')
  
  stage_sp_name = cbind(f_p = unique(re_cover_prop_nati$f_p),
                        intro_sp = intro_spname,
                        estab_sp = estab_spname, 
                        domin_sp = domin_spname)
  
  pure_intro = setdiff(colnames(re_cover_prop_inva)[10:ncol(re_cover_prop_inva)],
                       estab_sp)
  pure_estab = setdiff(estab_sp, domin_sp)
  pure_domin = domin_sp
  pure_native = colnames(re_cover_ab_nati)[10:ncol(re_cover_ab_nati)]
  sp_new = c(pure_native, pure_intro, pure_estab,
             pure_domin)  
  
  ### Create the data list for stan fit 
  re_cover_ab = as.data.frame(cbind(re_cover_ab[,1:9], 
                      re_cover_ab %>% select(all_of(sp_new)))) ## change the species order
  sp = colnames(re_cover_ab)[10:ncol(re_cover_ab)]
  re_cover_ab = arrange(re_cover_ab, re_cover_ab$fake_age)
  re_cover_ab[re_cover_ab==0] = 1e-06
  Pred = re_cover_ab[2:(nrow(re_cover_ab)),1:ncol(re_cover_ab)]
  Pred_t2 = re_cover_ab[1:(nrow(re_cover_ab)-1),1:ncol(re_cover_ab)]
  Yrs = nrow(Pred)
  yid = as.numeric(as.factor(Pred$fake_age))

  t_diff = Pred$Age-Pred_t2$Age
  
  datalist = list()
  for (l in 1:length(sp)) {
    #l = 1
    focal_sp = sp[l]
    focal_sp_posi = which(sp == focal_sp)
    pure_intro_fit = setdiff(pure_intro, focal_sp)
    if(length(pure_intro_fit) < 1) {pure_intro_fit = character(0)}
    pure_estab_fit = setdiff(pure_estab, focal_sp)
    if(length(pure_estab_fit) < 1) {pure_estab_fit = character(0)}
    pure_domin_fit = setdiff(pure_domin, focal_sp)
    if(length(pure_domin_fit) < 1) {pure_domin_fit = character(0)}
    pure_native_fit = setdiff(pure_native, focal_sp)
    
    sp_matrix = re_cover_ab[1:(nrow(re_cover_ab)-1),10:ncol(re_cover_ab)]
    sp_matrix_fit = cbind(sp_matrix[,1], 
                          natives = rowMeans(sp_matrix %>% select(all_of(pure_native_fit))),
                          intros = rowMeans(sp_matrix %>% select(all_of(pure_intro_fit))),
                          estabs = rowMeans(sp_matrix %>% select(all_of(pure_estab_fit))),
                          domins = rowMeans(sp_matrix %>% select(all_of(pure_domin_fit)))
                          )
    colnames(sp_matrix_fit)[1] = focal_sp
    sp_matrix_fit = sp_matrix_fit[,colSums(sp_matrix_fit!=0)!=0]
    sp_matrix_fit = sp_matrix_fit[, colSums(is.na(sp_matrix_fit)) == 0]
    
    Y = unlist(Pred %>% select(all_of(focal_sp)))
    I = sp_matrix_fit[,focal_sp]
    G = log(Y/I)/t_diff
    
    datalist[[l]] = list(N = Yrs, P = ncol(sp_matrix_fit), ftime=Yrs, time=yid,
                         G = G, X=sp_matrix_fit,
                         pure_native_fit = pure_native_fit,
                         pure_intro_fit = pure_intro_fit,
                         pure_estab_fit = pure_estab_fit,
                         pure_domin_fit = pure_domin_fit) 
  }
  
  stage_sp_infor = list(stage_sp_name = stage_sp_name,
                          pure_native = pure_native,
                          pure_intro = pure_intro,
                          pure_estab = pure_estab,
                          pure_domin = pure_domin,
                          sp = sp)
  ### Create a list to aggregate filter results
  filter_inf = list(stage_sp_infor = stage_sp_infor,
                    re_cover_ab = re_cover_ab,
                    datalist
  )
  
  return(filter_inf)
  
}