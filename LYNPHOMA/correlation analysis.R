
geo_list_BO<-sf::st_read(
  "C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/analysis/BODEN/ex.sqlite"
)

geo_list_RU<-sf::st_read(
  "C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/analysis/RUNIMC1/topLayers.sqlite"
)


trns<-scales::modulus_trans(0)

colToUse<-colnames(geo_list_BO)[3:17]

geo_list_BO_trans<-sf::st_drop_geometry(geo_list_BO)[,colToUse]
geo_list_BO_trans<-trns$transform(geo_list_BO_trans)

geo_list_RU_trans<-sf::st_drop_geometry(geo_list_RU)[,colToUse]
geo_list_RU_trans<-trns$transform(geo_list_RU_trans)

colnames(geo_list_BO_trans)<-unlist(lapply(strsplit(colnames(geo_list_BO_trans),'.',fixed = T),'[',2))

geo_list_BO_trans<-cbind(uid=geo_list_BO$uid[drop=T],geo_list_BO_trans,area=sf::st_area(geo_list_BO))

colnames(geo_list_RU_trans)<-unlist(lapply(strsplit(colnames(geo_list_RU_trans),'.',fixed = T),'[',2))

geo_list_RU_trans<-cbind(uid=geo_list_RU$uid[drop=T],geo_list_RU_trans,area=sf::st_area(geo_list_RU))

RU_CD4_CD8<-lapply(unique(geo_list_RU_trans$uid),function(u){
  cor.test(geo_list_RU_trans$cd4[geo_list_RU_trans$uid==u],geo_list_RU_trans$cd8[geo_list_RU_trans$uid==u])
})


BO_CD4_CD8<-lapply(unique(geo_list_RU_trans$uid),function(u){
  cor.test(geo_list_BO_trans$cd4[geo_list_BO_trans$uid==u],geo_list_BO_trans$cd8[geo_list_BO_trans$uid==u])
})

RU_CD3_CD20<-lapply(unique(geo_list_RU_trans$uid),function(u){
  cor.test(geo_list_RU_trans$cd3[geo_list_RU_trans$uid==u],geo_list_RU_trans$cd20[geo_list_RU_trans$uid==u])
})


BO_CD3_CD20<-lapply(unique(geo_list_RU_trans$uid),function(u){
  cor.test(geo_list_BO_trans$cd3[geo_list_BO_trans$uid==u],geo_list_BO_trans$cd20[geo_list_BO_trans$uid==u])
})


RU_CD3_CD68<-lapply(unique(geo_list_RU_trans$uid),function(u){
  cor.test(geo_list_RU_trans$cd3[geo_list_RU_trans$uid==u],geo_list_RU_trans$cd68[geo_list_RU_trans$uid==u])
})


BO_CD3_CD68<-lapply(unique(geo_list_RU_trans$uid),function(u){
  cor.test(geo_list_BO_trans$cd3[geo_list_BO_trans$uid==u],geo_list_BO_trans$cd68[geo_list_BO_trans$uid==u])
})


ss<-sample(1:10000,100)
cor.test(geo_list_BO_trans$cd3[ss],geo_list_BO_trans$cd68[ss])

cor.test(geo_list_RU_trans$cd3[ss],geo_list_RU_trans$cd68[ss])



CD3_CD68<-lapply(unique(geo_list_RU_trans$uid),function(u){
  BO<-geo_list_BO_trans[geo_list_BO_trans$uid==u & geo_list_BO_trans$cd68>2.5,]
  RU<-geo_list_RU_trans[geo_list_RU_trans$uid==u & geo_list_RU_trans$cd68>2.5,]
  # BO<-BO[sample(1:nrow(BO),1000),]
  # RU<-RU[sample(1:nrow(RU),1000),]
  cocor::cocor(~cd3+cd68|cd3+cd68,list(BO,RU))
})

CD3_CD68<-lapply(unique(geo_list_RU_trans$uid),function(u){
  BO<-geo_list_BO_trans[geo_list_BO_trans$uid==u,]
  RU<-geo_list_RU_trans[geo_list_RU_trans$uid==u,]
  # BO<-BO[sample(1:nrow(BO),1000),]
  # RU<-RU[sample(1:nrow(RU),1000),]
  cocor::cocor(~cd68+cd20|cd68+cd20,list(BO,RU))
})

TEMP<-mixtools::boot.comp(geo_list_RU_trans[sample(1:nrow(geo_list_RU_trans),500),'cd68',drop=T],max.comp = 1,mix.type = 'normalmix')
TEMP1<-mixtools::boot.comp(geo_list_BO_trans[sample(1:nrow(geo_list_BO_trans),500),'cd68',drop=T],max.comp = 1,mix.type = 'normalmix')


TEMP<-mixtools::normalmixEM(geo_list_BO_trans$cd16,k = 2,arbmean = T,arbvar = F,fast = T)
plot(TEMP,2)

TEMP<-mixtools::normalmixEM(geo_list_RU_trans$cd16,k = 2,arbmean = T,arbvar = F,fast = T)
plot(TEMP,2)

TEMP<-mclust::mclustBIC(geo_list_BO_trans[sample(1:nrow(geo_list_BO_trans),100),'cd8',drop=T],modelNames = 'E')
plot(TEMP)

TEMP<-mclust::mclustBIC(geo_list_RU_trans[sample(1:nrow(geo_list_RU_trans),100),'cd68',drop=T],modelNames = 'E')
plot(TEMP)


TEMP_RU<-geo_list_RU_trans[geo_list_RU_trans$uid==unique(geo_list_RU_trans$uid)[1],]
TEMP_max<-which.max(TEMP_RU$)

