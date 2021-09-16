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


importanceX<-c(0.1,0.1,0.3,
               0.8,0.8,0.8,
               0.1,0.8,0.8,
               0.5,0.8,0.5,
               0.8,0.8,1,
               0.8,0.8,1,
               1,0.1,0.8,
               0.8,0.8,0.3,
               1,0.8,0.8,
               0.8,0.1,0.1)



dimx=7
dimy=7
meta=10


BO_fsom<-FlowSOM::ReadInput(as.matrix(geo_list_BO_trans),
                         compensate = F,
                         transform = F,
                         scale = T,
                         silent = F)
BO_fsom<-FlowSOM::BuildSOM(fsom = BO_fsom,
                        colsToUse = colToUse,
                        silent = F,
                        xdim=dimx,
                        ydim=dimy,
                        # rlen=10,
                        # mst = 30,
                        # distf = 1,
                        importance=importanceX)
BO_fsom<-FlowSOM::BuildMST(BO_fsom)

saveRDS(BO_fsom,
        file = "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/two_soms/BO_som"  )

BO_fsomMFI<-FlowSOM::GetMFIs(BO_fsom)

BO_fsomMFI<-apply(BO_fsomMFI,2,function(x){(x-min(x))/(max(x)-min(x))})

colnames(BO_fsomMFI)<-c('CD38','SMA','VIM',
                     'CD14','CD33','CD16',
                     'PANK','CD11b','CD15',
                     'IgD','CD45','CD24',
                     'CD11c','FoxP3','CD4',
                     'E-CAD','CD68','CD20',
                     'CD8','Arg1','CD45ra',
                     'CD74','CD127','COLL',
                     'CD3','CD27','CD45ro',
                     'CD25','DNA1','DNA2')

rownames(BO_fsomMFI)<-as.character(formatC(1:nrow(BO_fsomMFI),flag='0',digits = 1,format = 'd'))

pdf(file="C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/two_soms/BO_heatmap.pdf",
    width = 8,
    height = 5
)
pheatmap::pheatmap(t(BO_fsomMFI),scale = 'none')

dev.off()

BO_SOM_Label<-FlowSOM::GetClusters(BO_fsom)

BO_SOM_Label<-factor(as.character(formatC(BO_SOM_Label,flag='0',format = 'd',digits=1)),levels = rownames(BO_fsomMFI))


geo_list_BO_annot<-dplyr::bind_cols(geo_list_BO,data.frame(SOM=BO_SOM_Label[1:nrow(geo_list_BO)]))



####### SOM RU #######


RU_fsom<-FlowSOM::ReadInput(as.matrix(geo_list_RU_trans),
                            compensate = F,
                            transform = F,
                            scale = T,
                            silent = F)
RU_fsom<-FlowSOM::BuildSOM(fsom = RU_fsom,
                           colsToUse = colToUse,
                           silent = F,
                           xdim=dimx,
                           ydim=dimy,
                           # rlen=10,
                           # mst = 30,
                           # distf = 1,
                           importance=importanceX)
RU_fsom<-FlowSOM::BuildMST(RU_fsom)

saveRDS(RU_fsom,
        file = "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/two_soms/RU_som"  )

RU_fsomMFI<-FlowSOM::GetMFIs(RU_fsom)

RU_fsomMFI<-apply(RU_fsomMFI,2,function(x){(x-min(x))/(max(x)-min(x))})

colnames(RU_fsomMFI)<-c('CD38','SMA','VIM',
                        'CD14','CD33','CD16',
                        'PANK','CD11b','CD15',
                        'IgD','CD45','CD24',
                        'CD11c','FoxP3','CD4',
                        'E-CAD','CD68','CD20',
                        'CD8','Arg1','CD45ra',
                        'CD74','CD127','COLL',
                        'CD3','CD27','CD45ro',
                        'CD25','DNA1','DNA2')

rownames(RU_fsomMFI)<-as.character(formatC(1:nrow(RU_fsomMFI),flag='0',digits = 1,format = 'd'))

pdf(file="C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/two_soms/RU_heatmap.pdf",
    width = 8,
    height = 5
)
pheatmap::pheatmap(t(RU_fsomMFI),scale = 'none')

dev.off()

RU_SOM_Label<-FlowSOM::GetClusters(RU_fsom)

RU_SOM_Label<-factor(as.character(formatC(RU_SOM_Label,flag='0',format = 'd',digits=1)),levels = rownames(RU_fsomMFI))


geo_list_RU_annot<-dplyr::bind_cols(geo_list_BO,data.frame(SOM=RU_SOM_Label[1:nrow(geo_list_BO)]))

TEMP<-apply(BO_fsomMFI,1,function(x){
  apply(RU_fsomMFI,1,function(y){
    dist(rbind(x,y))
  })
})

duplicated(apply(TEMP,1,function(x){rownames(TEMP)[which.min(x)]}))
      
sum(dist(RU_fsomMFI))
sum(dist(BO_fsomMFI))
