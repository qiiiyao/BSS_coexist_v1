getwd()
library(codyn)
library(dplyr)
library(vegan)
library(tidyr)

#read community data
bsplotdata<-read.csv("BSS_community_332.csv", header=T)
species_unique<-read.csv("species332.csv")
colnames(bsplotdata)[7:338]<-species_unique$Species

spdata1 <- filter(bsplotdata[,c(1:338)],Absolute_year%%2==0 & Field %in%"1")  # %%返回余数
spdata6 <- filter(bsplotdata[,c(1:338)],Absolute_year%%2==1 & Field %in%"6")
spdata7 <- filter(bsplotdata[,c(1:338)],Absolute_year%%2==0 & Field %in%"7")
spdata8 <- filter(bsplotdata[,c(1:338)],Absolute_year%%2==1 & Field %in%"8")
spdata9 <- filter(bsplotdata[,c(1:338)],Absolute_year%%2==0 & Field %in%"9" & Absolute_year>1967)
spdata10<- filter(bsplotdata[,c(1:338)],Absolute_year%%2==1 & Field %in%"10")
spdata4 <- filter(bsplotdata[,c(1:338)],Absolute_year%%2==0 & Field %in%"4" & Absolute_year>1967)
spdata5 <- filter(bsplotdata[,c(1:338)],Absolute_year%%2==1 & Field %in%"5")
spdata2 <- filter(bsplotdata[,c(1:338)],Absolute_year%%2==0 & Field %in%"2")
spdata3 <- filter(bsplotdata[,c(1:338)],Absolute_year%%2==1 & Field %in%"3")

###---------------------------------------- ---------------------------------------------------
###-------------------------------------------------------------------------------------------
for (l in c(3:7)) {       #set interval l=4, time=10 year
  #creat a NULL dataframe
  fielddata<- data.frame()  
  #field=k
  for (k in 1:10) {
    #filter 5 years data and creat a new dataframe 'spdata_5_i"
    initial_year<-unique(get(paste('spdata',k,sep=''))$Absolute_year)
    for (i in initial_year[1:(length(initial_year)-l)]) {
      
      spdata_5<-filter(get(paste('spdata',k,sep='')), 
                       Absolute_year >= i & 
                         Absolute_year <i+2*(l+1))
      spdata_5_sp<-spdata_5%>%select(-c(1:6))
      #translate absolute cover into relative_cover data
      spdata1_5_sp_re<-cbind(spdata_5 %>% select(c(1:6)),decostand(spdata_5_sp, method = "total"))
      
      
      ##1.gamma stability
      field_cover<-spdata_5%>%
        mutate(coversum = rowSums(.[7:338], na.rm = TRUE))
      field_cover_year<-field_cover%>%
        group_by(Absolute_year)%>%
        summarise_at(vars(coversum),funs(sum))
      cover_metacommunity_ini<-field_cover_year$coversum[1]
      cover_metacommunity_mean<-mean(field_cover_year$coversum)
      gamma_stab<-cover_metacommunity_mean/sd(field_cover_year$coversum)
      
      ##2.alpha stability (Averaged local community stability) 
      #plot temporal relative cover
      re_cover_temporalmean<- field_cover%>%
        group_by(Plot)%>%
        summarise_at(vars(coversum),funs(sum(., na.rm = TRUE)))%>% 
        mutate(re_cover=coversum/sum(coversum))
      alpha_cv_unweighted<- field_cover%>%
        select(Absolute_year,Plot,coversum)%>%
        group_by(Plot)%>%
        summarise_at(vars(coversum),funs(mean(., na.rm = TRUE),sd(.,na.rm = TRUE)))%>% 
        mutate(alpha_cv_unweighted=sd/mean)
      alpha_stab<-1/sum(re_cover_temporalmean %>% select(re_cover)*
                          alpha_cv_unweighted %>% select(alpha_cv_unweighted))
      ##3.species stability
      #methods1
      sp_data_unweighted  <- spdata_5 %>%
        gather("species","cover", -c(1:6))%>% 
        group_by(Plot,species)%>%
        summarise_at(vars(cover), funs(sd(., na.rm = TRUE)/mean(., na.rm = TRUE)))
      sp_data_unweighted <-as.data.frame(sp_data_unweighted)
      
      sp_re_cover_data <- spdata1_5_sp_re %>%    
        gather("species","re_cover", -c(1:6))%>% 
        group_by(Plot,species)%>%
        summarise_at(vars(re_cover), funs(mean(., na.rm = TRUE)))
      sp_re_cover_data<-as.data.frame(sp_re_cover_data)
      
      sp_data<-cbind(sp_data_unweighted , sp_re_cover_data$re_cover) 
      colnames(sp_data)[4] <- "stability"
      
      sp_stab<-as.numeric(sp_data%>%mutate(re_stability=stability*cover)%>% 
                            group_by(Plot)%>%
                            summarise_at(vars(re_stability), funs(sum(., na.rm = TRUE)))%>%
                            mutate(sp_stability=1/re_stability)%>%
                            summarise(mean(sp_stability)))
      #methods2
      sp_sd_unweighted <- spdata_5 %>%
        gather("species","cover", -c(1:6))%>% 
        group_by(Plot,species)%>%
        summarise_at(vars(cover), funs(sd(., na.rm = TRUE)))%>% 
        group_by(Plot)%>%
        summarise_at(vars(cover), funs(sum(., na.rm = TRUE)))
      cover_mean<-field_cover%>%
        group_by(Plot)%>%
        summarise_at(vars(coversum),funs(mean))
      species_stab<-mean(cover_mean$coversum/sp_sd_unweighted$cover) 
      ##4.species asynchrony 
      
      #alpha stability divide species stability
      sp_asynchrony<-alpha_stab/species_stab
      species_asynchrony<-alpha_stab/sp_stab
      ##5.spatial asynchrony     
      spatial_asynchrony<-gamma_stab/alpha_stab
      ##6.population asynchrony 
      Field_sp_data<-spdata_5 %>%   
        gather("species","cover",-c(1:6)) 
      sp_re_cover_amongfield<-Field_sp_data%>% 
        group_by(species)%>% 
        summarise_at(vars(cover), funs(sum))%>% 
        mutate(re_cover=cover/sum(cover))
      
      #Loreau 
      pop_synchrony_Loreau_unweighted <- codyn::synchrony(df            = Field_sp_data, 
                                                          time.var      = "Absolute_year", 
                                                          replicate.var = "species",
                                                          species.var   = "Plot", 
                                                          abundance.var = "cover", 
                                                          metric        = "Loreau")
      pop_synchrony<-sum(pop_synchrony_Loreau_unweighted  %>% select(synchrony)*
                           sp_re_cover_amongfield %>% select(re_cover),na.rm = T)
      
      ini_year<-i
      field<-k
      onefield<-cbind(field,ini_year,
                      cover_metacommunity_ini,cover_metacommunity_mean,
                      gamma_stab,alpha_stab,sp_stab,species_stab,
                      sp_asynchrony,spatial_asynchrony,pop_synchrony)  
      fielddata<-rbind(fielddata,onefield)
      
    }
  }

  #Even-numbered years remain unchanged, odd-numbered years minus one
  
  fielddata1<-fielddata
  for (j in 1:length(fielddata1$ini_year)) {
    if(fielddata1$ini_year[j]%%2>0){
      fielddata1$ini_year[j]<-fielddata1$ini_year[j]-1
    }else{
      fielddata1$ini_year[j]<-fielddata1$ini_year[j]
    }
  }
#Output  
  write.csv(fielddata1,paste('moving size=',l+1,'.csv',sep=''))
}
