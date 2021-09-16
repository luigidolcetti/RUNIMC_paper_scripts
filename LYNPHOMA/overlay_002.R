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

fsom<-readRDS(
  file = "C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/COMPARISONS/5x5_som"  )

fsomMFI<-FlowSOM::GetMFIs(fsom)

fsomMFI<-apply(fsomMFI,2,function(x){(x-min(x))/(max(x)-min(x))})

colnames(fsomMFI)<-unlist(lapply(strsplit(colnames(fsomMFI),'.',fixed = T),'[',2))

rownames(fsomMFI)<-as.character(formatC(1:nrow(fsomMFI),flag='0',digits = 1,format = 'd'))

SOM_Label<-FlowSOM::GetClusters(fsom)

SOM_Label<-factor(as.character(formatC(SOM_Label,flag='0',format = 'd',digits=1)),levels = rownames(fsomMFI))

geo_list_BO_annot<-dplyr::bind_cols(geo_list_BO[,c('uid')],
                                    geo_list_BO_trans,
                                    data.frame(area = sf::st_area(geo_list_BO)),
                                    data.frame(SOM=SOM_Label[1:nrow(geo_list_BO)]))

geo_list_RU_annot<-dplyr::bind_cols(geo_list_RU[,c('uid')],
                                    geo_list_RU_trans,
                                    data.frame(area = sf::st_area(geo_list_RU)),
                                    data.frame(SOM=SOM_Label[(1+nrow(geo_list_BO)):(nrow(geo_list_BO)+nrow(geo_list_RU))]))

rst<-RUNIMC:::retrieve.RsCollection(
"C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/rasterStacks"
)

rstcl<-RUNIMC:::retrieve.RsCollection(
  "C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/analysis/RUNIMC1/test/classification/rasterStacks"
)

uids<-names(rst)
raster::plot(RUNIMC::quantNorm(rst[[uids[1]]]$x142nd.asma.nd142di.,0.95),xlim=c(300,400),ylim=c(300,400),col=grey(seq(1,0,-0.1)))
plot(geo_list_BO_annot[geo_list_BO_annot$uid==uids[1],'GEOMETRY'],add=T,border='blue')

raster::plot(RUNIMC::quantNorm(rst[[uids[1]]]$x142nd.asma.nd142di.,0.95),xlim=c(300,400),ylim=c(300,400),col=grey(seq(1,0,-0.1)))
plot(geo_list_RU_annot[geo_list_RU_annot$uid==uids[1],'GEOMETRY'],add=T,border='blue')


raster::plot(rstcl[[uids[1]]]$sko_CD8,xlim=c(300,400),ylim=c(300,400),col=grey(seq(1,0,-0.1)))
plot(geo_list_RU_annot[geo_list_RU_annot$uid==uids[1] & geo_list_RU_annot$area>20,'GEOMETRY'],add=T,border='blue')
