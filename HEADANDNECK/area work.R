geo_list_BO<-sf::st_read(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/analysis/BODEN/ex.sqlite"
)

geo_list_RU<-sf::st_read(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/analysis/RUNIMC1/topLayers.sqlite"
)



trns<-scales::modulus_trans(0)

colToUse<-colnames(geo_list_BO)[3:32]

# geo_list_BO_trans<-sf::st_drop_geometry(geo_list_BO)[,colToUse]
geo_list_BO[,colToUse]<-trns$transform(sf::st_drop_geometry(geo_list_BO)[,colToUse])

# geo_list_RU_trans<-sf::st_drop_geometry(geo_list_RU)[,colToUse]
geo_list_RU[,colToUse]<-trns$transform(sf::st_drop_geometry(geo_list_RU)[,colToUse])

fsom<-readRDS(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001//som"
)


SOM_Label<-FlowSOM::GetClusters(fsom)

SOM_Label<-factor(as.character(formatC(SOM_Label,flag='0',format = 'd',digits=1)))


geo_list_BO_annot<-dplyr::bind_cols(geo_list_BO,data.frame(SOM=SOM_Label[1:nrow(geo_list_BO)]))


geo_list_RU_annot<-dplyr::bind_cols(geo_list_RU,data.frame(SOM=factor(as.character(SOM_Label[(1+nrow(geo_list_BO)):(nrow(geo_list_BO)+nrow(geo_list_RU))]))))

geo_list_BO_annot<-dplyr::bind_cols(
  geo_list_BO_annot,
  area = sf::st_area(geo_list_BO_annot)
)

geo_list_RU_annot$area<-sf::st_area(geo_list_RU_annot)

plot.new()
par(mfrow=c(5,5),mar=c(1,1,1,1))
for (i in levels(SOM_Label)){
BO_sub<-geo_list_BO_annot[geo_list_BO_annot$SOM==i,]
RU_sub<-geo_list_RU_annot[geo_list_RU_annot$SOM==i,]

chnl<-colnames(BO_sub)[1:30]
chnl_sec<-chnl[29]
BO_x<-unlist(sf::st_drop_geometry(BO_sub[,'area']))
RU_x<-unlist(sf::st_drop_geometry(RU_sub[,'area']))
x_min<-0
x_max<-quantile(c(unlist(BO_x),unlist(RU_x)),0.99)
BO_x<-BO_x[BO_x<x_max]
RU_x<-RU_x[RU_x<x_max]
brks<-seq(x_min,x_max,length.out=100)
BO_h<-hist(BO_x,breaks = brks,plot = F)
RU_h<-hist(RU_x,breaks = brks,plot = F)
y_min=0
y_max=max(c(BO_h$counts,RU_h$counts))
plot(NA,xlim=c(x_min,x_max),ylim=c(0,y_max))
points(BO_h$mids,BO_h$counts,type='l',col='blue')
points(RU_h$mids,RU_h$counts,type='l',col='red')
}
aggregate(area~SOM,geo_list_BO_annot,sum)
cbind(aggregate(area~SOM,geo_list_BO_annot,sum)
,aggregate(area~SOM,geo_list_RU_annot,sum)
)


BO_tot<-aggregate(area~uid,geo_list_BO_annot,sum)
RU_tot<-aggregate(area~uid,geo_list_BO_annot,sum)
BO_agg<-aggregate(area~uid+SOM,geo_list_BO_annot,sum)
RU_agg<-aggregate(area~uid+SOM,geo_list_RU_annot,sum)

uids<-unique(geo_list_BO_annot$uid)
names(uids)<-uids
BO_agg<-lapply(uids,function(u){
  out<-BO_agg[BO_agg$uid==u,c('SOM','area')]
  out$area<-out$area/BO_tot$area[BO_tot$uid==u]*100
  return(out)
})

RU_agg<-lapply(uids,function(u){
  out<-RU_agg[RU_agg$uid==u,c('SOM','area')]
  out$area<-out$area/RU_tot$area[RU_tot$uid==u]*100
  return(out)
})

tot_add<-lapply(uids,function(u){
 
  out<-dplyr::full_join(BO_agg[[u]],RU_agg[[u]],by=c('SOM'))
  names(out)[c(2,3)]<-c('BO','RU')
  out$BO[is.na(out$BO)]<-0
  out$RU[is.na(out$RU)]<-0
  for (i in levels(out$SOM)[!(levels(out$SOM) %in% out$SOM)]){
    out<-rbind.data.frame(out,data.frame(SOM=i,BO=0,RU=0))
  }
  out<-out[order(out$SOM),]
  return(out)
})

tot_som<-do.call(rbind.data.frame,tot_add)
pairwise.wilcox.test(tot_som$,)
