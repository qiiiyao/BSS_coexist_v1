get_longest_consec_rep = function(x) {
  #x = ddd
  if (!is.vector(x) && !is.list(x)) 
    stop("'x' must be a vector of an atomic type")
  n <- length(x)
  if (n == 0L) 
    return(structure(list(lengths = integer(), values = x), 
                     class = "rle"))
  y <- x[-1L] != x[-n]
  i <- c(which(y | is.na(y)), n)
  lengths = diff(c(0L, i))
  
  rep_c_l = list()
  for (j in 1:length(i)) {
    rep_c = x[(i[j]-lengths[j]+1):i[j]]
    seq_c = c((i[j]-lengths[j]+1):i[j])
    if(unique(rep_c == T)){
      rep_c_l[[j]] = seq_c
    }
  }
  rep_c_l = rep_c_l[sapply(rep_c_l, function(x){!is.null(x)})]
  len_rep_c = sapply(rep_c_l, function(x){length(x)})
  if (sum(len_rep_c) == length(len_rep_c)*len_rep_c[1]){
    focal_seq = unlist(sample(rep_c_l, 1))
  } else {
    focal_seq = rep_c_l[[which.max(len_rep_c)]]
  }
  
  return(focal_seq)
  
}




get_consec_rep = function(x) {
  #x = y
  if (!is.vector(x) && !is.list(x)) 
    stop("'x' must be a vector of an atomic type")
  n <- length(x)
  if (n == 0L) 
    return(structure(list(lengths = integer(), values = x), 
                     class = "rle"))
  y <- x[-1L] != x[-n]
  i <- c(which(y | is.na(y)), n)
  lengths = diff(c(0L, i))
  
  rep_c_l = list()
  for (j in 1:length(i)) {
    #j = 8
    rep_c = x[(i[j]-lengths[j]+1):i[j]]
    seq_c = c((i[j]-lengths[j]+1):i[j])
    if(unique(rep_c == T)){
      rep_c_l[[j]] = seq_c
    }
  }
  rep_c_l = rep_c_l[sapply(rep_c_l, function(x){!is.null(x)})]
  
  return(rep_c_l)
  
}
