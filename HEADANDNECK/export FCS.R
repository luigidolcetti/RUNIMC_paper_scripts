targetFolder<-"D:/Luigi/FCS.HeadAndNeak"

geo_list_BO<-sf::st_read(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/analysis/BODEN/ex.sqlite"
)

geo_list_RU<-sf::st_read(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/analysis/RUNIMC1/topLayers.sqlite"
)


trns<-scales::modulus_trans(0)

colToUse<-colnames(geo_list_BO)[3:32]

geo_list_BO_trans<-sf::st_drop_geometry(geo_list_BO)[,c('uid',colToUse)]
geo_list_BO_trans[,colToUse]<-trns$transform(geo_list_BO_trans[,colToUse])

geo_list_RU_trans<-sf::st_drop_geometry(geo_list_RU)[,c('uid',colToUse)]
geo_list_RU_trans[,colToUse]<-trns$transform(geo_list_RU_trans[,colToUse])


for (u in unique(geo_list_BO_trans$uid)){
  TEMP<-as.matrix(geo_list_BO_trans[geo_list_BO_trans$uid==u,colToUse])
  TEMP<-flowCore::flowFrame(TEMP)
  flowCore::write.FCS(TEMP,
                    file.path(targetFolder,paste0('BO_',u,'.fcs')))
  TEMP<-as.matrix(geo_list_RU_trans[geo_list_RU_trans$uid==u,colToUse])
  TEMP<-flowCore::flowFrame(TEMP)
  flowCore::write.FCS(TEMP,
                     file.path(targetFolder,paste0('RU_',u,'.fcs')))
  
  
}
