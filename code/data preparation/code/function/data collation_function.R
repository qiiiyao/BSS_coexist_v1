filter_comm_sum = function(x, per, rare = FALSE){
  if (is.null(x)) {print('no species!')
  } else if (ncol(x) < 2){z = x
                          return(z)} else {
    
  y = data.frame(freq = apply(x, 2, function(y){length(which(y != 0))})/nrow(x)) %>% filter(freq > per)
  re_sp = row.names(y)
  rare_sp = setdiff(colnames(x), re_sp)
  
  if (length(re_sp) > 0 & rare == T & length(rare_sp) > 0){
    rare = rowSums(x %>% select(rare_sp))
    z = cbind(x %>% select(all_of(re_sp)), rare = rare)
    return(z)
  } else if (rare == T & length(rare_sp) > 0) {
    rare = rowSums(x %>% select(rare_sp))
    z = data.frame(rare_sum = rare)
    return(z)
  } else { print('No species meet your conditions!')
    z = NULL
    return(z)}

}
}

filter_comm_mean = function(x, per, rare = FALSE){
  if (is.null(x)) {print('no species!')
  } else if (ncol(x) < 2){z = x
                          return(z)} else {
    
    y = data.frame(freq = apply(x, 2, function(y){length(which(y != 0))})/nrow(x)) %>% filter(freq > per)
    re_sp = row.names(y)
    rare_sp = setdiff(colnames(x), re_sp)
    
    if (length(re_sp) > 0 & rare == T & length(rare_sp) > 0){
      rare = rowMeans(x %>% select(rare_sp))
      z = cbind(x %>% select(all_of(re_sp)), rare = rare)
      return(z)
    } else if (length(re_sp) > 0) {
      z = x %>% select(all_of(re_sp))
      return(z)
    } else if (rare == T & length(rare_sp) > 0) {
      rare = rowMeans(x %>% select(rare_sp))
      z = data.frame(rare_mean = rare)
      return(z)
    } else { print('No species meet your conditions!')
             z = NULL
             return(z)}
  }
  
}