geo_list_BO<-sf::st_read(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/analysis/BODEN/ex.sqlite"
)

geo_list_RU<-sf::st_read(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/analysis/RUNIMC1/topLayers.sqlite"
)


trns<-scales::modulus_trans(0)

colToUse<-colnames(geo_list_BO)[3:32]

geo_list_BO_trans<-sf::st_drop_geometry(geo_list_BO)[,colToUse]
geo_list_BO_trans<-trns$transform(geo_list_BO_trans)

geo_list_RU_trans<-sf::st_drop_geometry(geo_list_RU)[,colToUse]
geo_list_RU_trans<-trns$transform(geo_list_RU_trans)

geo_list_total<-as.matrix(rbind(geo_list_BO_trans,geo_list_RU_trans))

fsom<-readRDS(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001//som"
)


SOM_Label<-FlowSOM::GetClusters(fsom)

SOM_Label<-factor(as.character(formatC(SOM_Label,flag='0',format = 'd',digits=1)))


geo_list_BO_annot<-dplyr::bind_cols(geo_list_BO_trans,data.frame(SOM=SOM_Label[1:nrow(geo_list_BO)]))


geo_list_RU_annot<-dplyr::bind_cols(geo_list_RU_trans,data.frame(SOM=factor(as.character(SOM_Label[(1+nrow(geo_list_BO)):(nrow(geo_list_BO)+nrow(geo_list_RU))]))))

BO_sub<-geo_list_BO_annot[geo_list_BO_annot$SOM=='02',]
RU_sub<-geo_list_RU_annot[geo_list_RU_annot$SOM=='02',]

chnl<-colnames(BO_sub)[1:30]
chnl_sec<-chnl[29]
BO_x<-BO_sub[,chnl_sec]
RU_x<-RU_sub[,chnl_sec]
x_min<-0
x_max<-max(c(BO_x,RU_x))
brks<-seq(x_min,x_max,length.out=100)
BO_h<-hist(BO_x,breaks = brks,plot = F)
RU_h<-hist(RU_x,breaks = brks,plot = F)
y_min=0
y_max=max(c(BO_h$counts,RU_h$counts))
plot(NA,xlim=c(x_min,x_max),ylim=c(0,y_max))
points(BO_h$mids,BO_h$counts,type='l',col='blue')
points(RU_h$mids,RU_h$counts,type='l',col='red')
