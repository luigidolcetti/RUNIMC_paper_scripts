myRaster<-RUNIMC:::retrieve.RsCollection(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/rasterStacks"
)

geo_list_BO<-sf::st_read(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/analysis/BODEN/BP_SOM.sqlite"
)

geo_list_RU<-sf::st_read(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/analysis/RUNIMC1/RU_SOM.sqlite"
)

uids<-unique(geo_list_BO$uid)
names(uids)<-uids
BO_mat<-lapply(uids,function(u){
  out<-exactextractr::exact_extract(myRaster[[u]][[1]],geo_list_BO[geo_list_BO$uid==u,],include_cell=T,fun=NULL,include_cols=c('som','uid'))
  out<-do.call(rbind.data.frame,out)[,c('som','uid','cell')]
  return(out)
})
BO_mat<-do.call(rbind.data.frame,BO_mat)

RU_mat<-lapply(uids,function(u){
  out<-exactextractr::exact_extract(myRaster[[u]][[1]],geo_list_RU[geo_list_RU$uid==u,],include_cell=T,fun=NULL,include_cols=c('som','uid'))
  out<-do.call(rbind.data.frame,out)[,c('som','uid','cell')]
  return(out)
})

RU_mat<-do.call(rbind.data.frame,RU_mat)


tot_mat<-dplyr::full_join(BO_mat,RU_mat,by=c('uid','cell'))

colnames(tot_mat)<-c("som.BO","uid","cell","som.RU")
tot_mat$som.BO[is.na(tot_mat$som.BO)]<-'00'
tot_mat$som.RU[is.na(tot_mat$som.RU)]<-'00'
tot_mat$som.BO<-factor(tot_mat$som.BO)
tot_mat$som.RU<-factor(tot_mat$som.RU)
tot_table<-lapply(uids,function(u){table(tot_mat[tot_mat$uid==u,c('som.BO','som.RU')])})

concordance<-lapply(uids,function(u){
  sR<-rowSums(tot_table[[u]])
  sC<-colSums(tot_table[[u]])
  dg<-diag(tot_table[[u]])
  
  dg/(sR+sC-dg)*100
})
                   
concordance<-do.call(rbind,concordance)

apply(concordance,2,function(x){mean(na.omit(x))})      
apply(concordance,2,function(x){sd(na.omit(x))})      


heatmap(tot_table,Colv = NA,Rowv = NA,scale = 'row')

sum(diag(tot_table)/rowSums(tot_table))
sum(diag(tot_table)/colSums(tot_table))
