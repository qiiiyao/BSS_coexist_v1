fit_prepa_original_field = function(w, t = NULL, seq_t_all, nati_per, intro_per, estab_per,
                              domin_per) {
  #w = sp_cover_f3_fp[[36]]
  #t = 2:14
 # t = NULL
  #w = sp_cover_field[[10]]
  #nati_per = 0.2
  #intro_per = 0.2 
  #estab_per = 0.2
  #domin_per = 0.2
  ### filter the different invasion stage species
  if (is.null(t)) { 
    re_cover_ab = cbind(w[,1:4],
                        w[,5:(ncol(w))][colSums(w[,
                                                   5:(ncol(w))]) != 0])
  } else { 
    t_seq = seq_t_all[t]
    w_t = w[w$Age %in% t_seq,]
    re_cover_ab = cbind(w_t[,1:4],
                        w_t[,5:(ncol(w_t))][colSums(w_t[,
                                                         5:(ncol(w_t))]) != 0])
  }
  re_cover_ab = as.data.frame(re_cover_ab) ## ungrouped grouped dataframe 
  re_cover_ab_inva = cbind(re_cover_ab[,1:4],
                           re_cover_ab %>%
                             select(all_of(intersect(colnames(re_cover_ab)[5:ncol(re_cover_ab)],
                                                     unlist(trait %>% filter(Origin == 'Exotic') %>% select(Species))))))
  re_cover_ab_nati = cbind(re_cover_ab[,1:4],
                           re_cover_ab %>%
                             select(all_of(intersect(colnames(re_cover_ab)[5:ncol(re_cover_ab)],
                                                     unlist(trait %>% filter(Origin == 'Native') %>% select(Species))))))
  
  re_cover_prop = cbind(re_cover_ab[,1:4],
                        t(apply(re_cover_ab[,5:(ncol(re_cover_ab))], 1, function(x){x/sum(x)})))
  re_cover_prop_inva = cbind(re_cover_prop[,1:4],
                             re_cover_prop %>% 
                               select(all_of(intersect(colnames(re_cover_prop)[5:ncol(re_cover_prop)],
                                                       unlist(trait %>% filter(Origin == 'Exotic') %>% select(Species))))))
  re_cover_prop_nati = cbind(re_cover_prop[,1:4],
                             re_cover_prop %>%
                               select(all_of(intersect(colnames(re_cover_prop)[5:ncol(re_cover_prop)],
                                                       unlist(trait %>% filter(Origin == 'Native') %>% select(Species))))))
  
  estab_sp = vector() ## species at establishment stage (occur more than
  #10 successive years)
  domin_sp = vector() ## species at dominant stage (occur more than 10
  #successive years and mean cover is beyond any natives)
  
  for (j in 1:(ncol(re_cover_prop_inva)-4)) {
    #j = 15
    
    focal_alien = re_cover_prop_inva[,c(1, 4+j)]
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
        max_nati = max(apply(x_nati[,5:ncol(x_nati)], 2, function(x){mean(x)}))
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
    re_cover_prop_domin = cbind(re_cover_prop[,1:4])
    re_cover_ab_domin = cbind(re_cover_ab[,1:4])
  } else {
    re_cover_prop_domin = cbind(re_cover_prop[,1:4],
                                re_cover_prop %>% select(all_of(domin_sp)))
    re_cover_ab_domin = cbind(re_cover_ab[,1:4],
                              re_cover_ab %>% select(all_of(domin_sp)))
  }
  
  if (length(estab_sp) == 0) {
    re_cover_prop_estab = cbind(re_cover_prop[,1:4])
    re_cover_ab_estab = cbind(re_cover_ab[,1:4])
  } else {
    re_cover_prop_estab = cbind(re_cover_prop[,1:4],
                                re_cover_prop %>% select(all_of(estab_sp)))
    re_cover_ab_estab = cbind(re_cover_ab[,1:4],
                              re_cover_ab %>% select(all_of(estab_sp)))
  }
  re_cover_prop_intro = cbind(re_cover_prop[,1:4],
                              re_cover_prop %>%
                                select(all_of(setdiff(colnames(re_cover_prop_inva)[5:ncol(re_cover_prop_inva)], estab_sp))))
  re_cover_ab_intro = cbind(re_cover_ab[,1:4],
                            re_cover_ab %>%
                              select(all_of(setdiff(colnames(re_cover_ab_inva)[5:ncol(re_cover_ab_inva)], estab_sp))))
  
  intro_spname = paste(setdiff(colnames(re_cover_prop_inva)[5:ncol(re_cover_prop_inva)],
                               estab_sp),
                       collapse = ', ')
  
  stage_sp_name = cbind(f_p = unique(re_cover_prop_nati$f_p),
                        intro_sp = intro_spname,
                        estab_sp = estab_spname, 
                        domin_sp = domin_spname)
  
  x = re_cover_ab_estab
  if (ncol(x) < 11){
    x = cbind(x, fake_sp = 0, fake_sp2 = 0)
  }
  y = x[,5:ncol(x)]
  z = filter_comm_mean(y, per = estab_per, rare = F)
  if(is.null(z)){re_cover_ab_estab_f = NULL
  sr_estab = 0
  sp_estab_f = NULL
  } else {
    #  colnames(z)[ncol(z)] = 'estab_rare'
    f = cbind(x[,1:4], z)
    re_cover_ab_estab_f = f
    sr_estab = ncol(f)-4
    sp_estab_f = colnames(re_cover_ab_estab_f)[5:ncol(re_cover_ab_estab_f)]}
  
  x_1 = re_cover_ab_intro
  if (ncol(x_1) < 11){
    x_1 = cbind(x_1, fake_sp = 0, fake_sp2 = 0)
  }
  y_1 = x_1[,5:ncol(x_1)]
  z_1 = filter_comm_mean(y_1, per = intro_per, rare = F)
  if(is.null(z_1)){
    re_cover_ab_intro_f = NULL
    sr_intro = 0
    sp_intro_f = NULL} else {
      #   colnames(z_1)[ncol(z_1)] = 'intro_rare'
      f_1 = cbind(x_1[,1:4], z_1)
      re_cover_ab_intro_f = f_1
      sr_intro = ncol(f_1)-4
      sp_intro_f = colnames(re_cover_ab_intro_f)[5:ncol(re_cover_ab_intro_f)]}
  
  
  x_2 = re_cover_ab_nati
  if (ncol(x_2) < 11){
    x_2 = cbind(x_2, fake_sp = 0, fake_sp2 = 0)
  }
  y_2 = x_2[,5:ncol(x_2)]
  z_2 = filter_comm_mean(x = y_2, per = nati_per, rare = F)
  if(is.null(z_2)){re_cover_ab_nati_f = NULL
  sr_nati = 0
  sp_nati_f = NULL
  } else {
    #  colnames(z_2)[ncol(z_2)] = 'nati_rare'
    f_2 = cbind(x_2[,1:4], z_2)
    re_cover_ab_nati_f = f_2
    sr_nati = ncol(f_2)-4
    sp_nati_f = colnames(re_cover_ab_nati_f)[5:ncol(re_cover_ab_nati_f)]}
  
  x_3 = re_cover_ab_domin
  if (ncol(x_3) < 11){
    x_3 = cbind(x_3, fake_sp = 0, fake_sp2 = 0)
  }
  y_3 = x_3[,5:ncol(x_3)]
  z_3 = filter_comm_mean(x = y_3, per = domin_per, rare = F)
  if(is.null(z_3)){re_cover_ab_domin_f = NULL
  sr_domin = 0
  sp_domin_f = NULL
  } else {
    #  colnames(z_3)[ncol(z_3)] = 'domin_rare'
    f_3 = cbind(x_3[,1:4], z_3)
    re_cover_ab_domin_f = f_3
    sr_domin = ncol(f_3)-4
    sp_domin_f = colnames(re_cover_ab_domin_f)[5:ncol(re_cover_ab_domin_f)]}
  
  stage_sp_name_f = cbind(field = unique(re_cover_prop_nati$Field),
                          nati_sp = paste(sp_nati_f, collapse = ', '),
                          intro_sp = paste(sp_intro_f, collapse = ', '),
                          estab_sp = paste(sp_estab_f, collapse = ', '), 
                          domin_sp = paste(sp_domin_f, collapse = ', '))
  
  
  if (is.null(re_cover_ab_intro_f)) {
    col_intro = 0
  } else {
    col_intro = ncol(as.data.frame(re_cover_ab_intro_f[,5:ncol(re_cover_ab_intro_f)]))
  }
  if (is.null(re_cover_ab_estab_f)) {
    col_estab = 0
  } else {
    col_estab = ncol(as.data.frame(re_cover_ab_estab_f[,5:ncol(re_cover_ab_estab_f)]))
  }
  
  if (col_intro>0 & col_estab>0) {
    re_cover_ab_f = cbind(re_cover_ab_nati_f, 
                          re_cover_ab_intro_f[,5:ncol(re_cover_ab_intro_f)],
                          re_cover_ab_estab_f[,5:ncol(re_cover_ab_estab_f)])
  } else if (col_intro>0) {
    re_cover_ab_f = cbind(re_cover_ab_nati_f, 
                          re_cover_ab_intro_f[,5:ncol(re_cover_ab_intro_f)])
  } else if (col_estab>0) {
    re_cover_ab_f = cbind(re_cover_ab_nati_f, 
                          as.data.frame(re_cover_ab_estab_f[,5:ncol(re_cover_ab_estab_f)]))
  } else {
    re_cover_ab_f = cbind(re_cover_ab_nati_f)
  }
  
  ### Create the data list for stan fit 
  colnames(re_cover_ab_f)[5:ncol(re_cover_ab_f)] = c(sp_nati_f, sp_intro_f, sp_estab_f)
  sp = colnames(re_cover_ab_f)[5:ncol(re_cover_ab_f)]
  re_cover_ab_f = arrange(re_cover_ab_f, re_cover_ab_f$Age)
  re_cover_ab_f[re_cover_ab_f==0] = 1e-2
  Pred = re_cover_ab_f[2:(nrow(re_cover_ab_f)),1:ncol(re_cover_ab_f)]
  Pred_t2 = re_cover_ab_f[1:(nrow(re_cover_ab_f)-1),1:ncol(re_cover_ab_f)]
  Yrs = nrow(Pred)
  yid = as.numeric(as.factor(sort(Pred$Age)))
  sp_matrix = re_cover_ab_f[1:(nrow(re_cover_ab_f)-1),5:ncol(re_cover_ab_f)]
  t_diff = Pred$Age-Pred_t2$Age
  
  datalist = list()
  for (l in 1:length(sp)) {
    #l = 1
    
    Y = unlist(Pred %>% select(all_of(sp[l])))
    I = unlist(sp_matrix %>% select(all_of(sp[l])))
    G = log(Y/I)/t_diff
    
    datalist[[l]] = list(N = Yrs, P = ncol(sp_matrix), ftime=Yrs, time=yid,
                         G = G, X=sp_matrix) }
  
  ### Create a list to aggregate filter results
  filter_inf = list(stage_sp_name = stage_sp_name,
                    stage_sp_name_f = stage_sp_name_f,
                    re_cover_ab = re_cover_ab,
                    re_cover_ab_f = re_cover_ab_f,
                    datalist
  )
  
  return(filter_inf)
  
}