#function to estimate loss of species. c.ab is a vector of abundances from community i and inv.ab is the abundance of invasive species

exp.ext<-function(c.ab,inv.ab){
  if (sum(c.ab<=0)>0) print("Removing species with abundance = 0")
  c.ab<-c.ab[c.ab>0]
  tmp<-c.ab-(inv.ab*(c.ab/sum(c.ab)))
  return(sum(tmp<=1))
}

# function to randomly remove resident abundance in units of 1,
# rounded down to nearest whole number. c.ab is a vector of abundances
# from community i and inv.ab is the abundance of invasive species. 
# Times is the number of randomizations and hist = TRUE returns a histogram 
# along with summary statistics
ext.rnd<-function(c.ab,inv.ab,times = 999){
  
  if (sum(c.ab<=0)>0) print("Removing species with abundance = 0")
  c.ab<-c.ab[c.ab>0]
  c.ab<-c.ab[order(c.ab,decreasing=TRUE)]
  inv.ab<-floor(inv.ab)
  rank<-1:length(c.ab)
  
  ave.rank<-NULL
  num.ext<-NULL
  for (j in 1:times){
    #j = 1
    tmp<-data.frame(c.ab,rank)
    
    for (i in 1:inv.ab){
      ## random decrease a species by 1
      pos<-sample(nrow(tmp),1,prob=(c.ab/sum(c.ab)))
      tmp[pos,1]<-tmp[pos,1]-1
    }
    
    tmp$c.ab<-floor(tmp$c.ab)
    num.ext[j]<-sum(tmp$c.ab<=1) # the number of extinction species: index k
    # average rank of extinction species: index R
    if (sum(tmp$c.ab<=1)>0)	ave.rank[j]<-mean(tmp$rank[tmp$c.ab<=1])
    if (sum(tmp$c.ab<=1)<=0) ave.rank[j]<-length(c.ab)+1
  }
  
  out<-data.frame(num.ext,ave.rank)
  return(out)
  
}


### summary stat function. takes output from previous function (ext.rnd)
# and user adds the number of extinctions for a community 
# and a value or vector of the ranks of the extinct species
ext.sum.stat<-function(c.ab,inv.ab,times = 999,obs.num.ext,obs.ext.rank){
  ext.rnd.dat<-ext.rnd(c.ab,inv.ab,times)
  if (sd(ext.rnd.dat[,1]) == 0) {
    sd_ext_sim = 0.00001
  } else {
    sd_ext_sim = sd(ext.rnd.dat[,1])
  }
  ext.rnd.mean<-mean(ext.rnd.dat[,1])
  ext.z.value<-(obs.num.ext-ext.rnd.mean)/sd_ext_sim
  
  ext.rnd.95<-quantile(ext.rnd.dat[,1],probs=c(0.025,0.975))
  rank.mean<-mean(ext.rnd.dat[,2])
  
  if (sd(ext.rnd.dat[,2]) == 0) {
    sd_rank_sim = 0.00001
  } else {
    sd_rank_sim = sd(ext.rnd.dat[,2])
  }
  rank.z.value<-(mean(obs.ext.rank)-rank.mean)/sd_rank_sim
  rank.95<-quantile(ext.rnd.dat[,2],probs=c(0.025,0.975))
  ext.p.value<-rank(c(obs.num.ext,ext.rnd.dat[,1]))[1]/(length(ext.rnd.dat[,1])+1)
  if (ext.p.value>0.5) ext.p.value<- 1-ext.p.value
  if (ext.p.value<=0.5) ext.p.value<- ext.p.value
  
  rank.p.value=rank(c(mean(obs.ext.rank),ext.rnd.dat[,2]))[1]/(length(ext.rnd.dat[,2])+1)
  if (rank.p.value>0.5) rank.p.value<- 1-rank.p.value
  if (rank.p.value<=0.5) rank.p.value<- rank.p.value
  
  
  out<-data.frame(obs.num.ext,
                  ext.rnd.mean,
                  ext.rnd.lower=ext.rnd.95[1],
                  ext.rnd.upper=ext.rnd.95[2],
                  ext.z.value,
                  ext.p.value,
                  obs.ext.rank=mean(obs.ext.rank),
                  rank.mean,
                  rank.lower=rank.95[1],
                  rank.upper=rank.95[2],
                  rank.z.value,
                  rank.p.value,
                  row.names=1)
  
  return(out)
}
