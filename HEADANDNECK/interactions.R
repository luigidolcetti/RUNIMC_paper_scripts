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



dimx=5
dimy=5
meta=10


fsom<-FlowSOM::ReadInput(geo_list_total,
                         compensate = F,
                         transform = F,
                         scale = T,
                         silent = F)
fsom<-FlowSOM::BuildSOM(fsom = fsom,
                        colsToUse = colToUse,
                        silent = F,
                        xdim=dimx,
                        ydim=dimy,
                        # rlen=10,
                        # mst = 30,
                        # distf = 1,
                        importance=importanceX)
fsom<-FlowSOM::BuildMST(fsom)

# saveRDS(fsom,
#         file = "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001//som"  )

fsomMFI<-FlowSOM::GetMFIs(fsom)

fsomMFI<-apply(fsomMFI,2,function(x){(x-min(x))/(max(x)-min(x))})

colnames(fsomMFI)<-c('CD38','SMA','VIM',
                     'CD14','CD33','CD16',
                     'PANK','CD11b','CD15',
                     'IgD','CD45','CD24',
                     'CD11c','FoxP3','CD4',
                     'E-CAD','CD68','CD20',
                     'CD8','Arg1','CD45ra',
                     'CD74','CD127','COLL',
                     'CD3','CD27','CD45ro',
                     'CD25','DNA1','DNA2')

rownames(fsomMFI)<-as.character(formatC(1:nrow(fsomMFI),flag='0',digits = 1,format = 'd'))

# pdf(file="C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/heatmap.pdf",
#     width = 8,
#     height = 5
# )
# pheatmap::pheatmap(t(fsomMFI),scale = 'none')
# 
# dev.off()

SOM_Label<-FlowSOM::GetClusters(fsom)

SOM_Label<-factor(as.character(formatC(SOM_Label,flag='0',format = 'd',digits=1)),levels = rownames(fsomMFI))


geo_list_BO_annot<-dplyr::bind_cols(geo_list_BO,data.frame(SOM=SOM_Label[1:nrow(geo_list_BO)]))


geo_list_RU_annot<-dplyr::bind_cols(geo_list_RU,data.frame(SOM=factor(as.character(SOM_Label[(1+nrow(geo_list_BO)):(nrow(geo_list_BO)+nrow(geo_list_RU))]))))

uids<-unique(geo_list_BO$uid)

homo<-lapply(uids[1],function(u){
  BO_sub<-geo_list_BO_annot[geo_list_BO_annot$uid==u,]
  BO_int<-sf::st_intersects(BO_sub)
  
  BO_int_rap<-lapply(1:length(BO_int),function(x){
    partners<-BO_int[[x]]
    partners<-partners[partners!=x]
    data.frame(source=as.character(unlist(BO_sub$SOM[x,drop=T])),
         partner=as.character(unlist(BO_sub$SOM[BO_int[[x]],drop=T])),stringsAsFactors = F)})
  
  # BO_homo<-unlist(lapply(BO_int_rap,function(x){sum(x$source==x$partner)/length(x$partner)}))
  
  RU_sub<-geo_list_RU_annot[geo_list_RU_annot$uid==u,]
  RU_int<-sf::st_intersects(RU_sub)
  
  RU_int_rap<-lapply(1:length(RU_int),function(x){
    partners<-RU_int[[x]]
    partners<-partners[partners!=x]
    data.frame(source=as.character(unlist(RU_sub$SOM[x,drop=T])),
         partner=as.character(unlist(RU_sub$SOM[RU_int[[x]],drop=T])),stringsAsFactors = F)})
  
  # RU_homo<-unlist(lapply(RU_int_rap,function(x){sum(x$source==x$partner)/length(x$partner)}))
 
  out_BO<-do.call(rbind,BO_int_rap)
  out_BO<-out_BO[order(out_BO[,1],out_BO[,2]),]
  out_RU<-do.call(rbind,RU_int_rap)
  out_RU<-out_RU[order(out_RU[,1],out_RU[,2]),]
  return(list(BO=out_BO,RU=out_RU))
})

hTemp<-hclust(dist(fsomMFI[,c('CD4','CD8','CD68','CD14','E-CAD','PANK','CD20')],method = 'euclidean'))

BO<-table(homo[[1]]$BO)
plus<-hTemp$labels[!(hTemp$labels %in% colnames(BO))]
minus<-hTemp$labels[(hTemp$labels %in% colnames(BO))]
extra<-array(0,
            dim = c(length(plus),length(hTemp$labels)),
            dimnames=list(source=plus,partner=hTemp$labels))
BO<-rbind(BO,extra[,minus])
BO<-cbind(BO,t(extra))
BO<-BO[order(rownames(BO)),order(colnames(BO))]

RU<-table(homo[[1]]$RU)
plus<-hTemp$labels[!(hTemp$labels %in% colnames(RU))]
minus<-hTemp$labels[(hTemp$labels %in% colnames(RU))]
extra<-array(0,
             dim = c(length(plus),length(hTemp$labels)),
             dimnames=list(source=plus,partner=hTemp$labels))
RU<-rbind(RU,extra[,minus])
RU<-cbind(RU,t(extra))
RU<-RU[order(rownames(RU)),order(colnames(RU))]

# 
# TOT<-array(NA,dim=c(ncol(BO),nrow(BO),2))
# TOT[,,1]<-BO
# TOT[,,2]<-RU
# mantelhaen.test(TOT)
# 
# vcd::woolf_test(TOT)
# 
# rcompanion::groupwiseCMH(TOT)
# 
# 
# 
# par(mfrow=c(1,2))
# plot.new()
pheatmap::pheatmap(scales::boxcox_trans(0.1)$transform(BO),
                   scale ='none' ,
                   cluster_rows = hTemp,
                  cluster_cols = hTemp,
                  
                   Rowv = as.dendrogram(hTemp),
                   Colv =  as.dendrogram(hTemp))
pheatmap::pheatmap(scales::boxcox_trans(0.1)$transform(RU),
                   scale ='none' ,
                   cluster_rows = hTemp,
                   cluster_cols = hTemp,
                   
                   Rowv = as.dendrogram(hTemp),
                   Colv =  as.dendrogram(hTemp))
heatmap(RU,scale ='row' ,Rowv = as.dendrogram(hTemp),Colv =  as.dendrogram(hTemp))




chisq.test(table(homo[[1]]$BO),table(homo[[1]]$RU))


hh<-do.call(rbind,homo)

t.test(unlist(hh[,1]),unlist(hh[,2]),paired = T)

mean(unlist(hh[,1]),)
mean(unlist(hh[,2]),)