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

SOM_Label<-factor(as.character(formatC(SOM_Label,flag='0',format = 'd',digits=1)))


geo_list_BO_annot<-dplyr::bind_cols(geo_list_BO,data.frame(SOM=factor(SOM_Label[1:nrow(geo_list_BO)])))


geo_list_RU_annot<-dplyr::bind_cols(geo_list_RU,data.frame(SOM=factor(as.character(SOM_Label[(1+nrow(geo_list_BO)):(nrow(geo_list_BO)+nrow(geo_list_RU))]))))


uids<-unique(geo_list_BO_annot$uid)
rst<-RUNIMC:::retrieve.RsCollection(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/rasterStacks"
)

BO_mat<-lapply(setNames(uids,uids)[1],function(u,
                             frst=rst,
                             fsql=geo_list_BO_annot){
  ncl<-raster::ncol(frst[[u]])
  nrw<-raster::nrow(frst[[u]])
 
  out<-rep('00',ncl*nrw)
  newCoord<-exactextractr::exact_extract(frst[[u]][[1]],
                                         fsql[fsql$uid==u,],
                                         include_cell = T,
                                         include_cols = 'SOM')
  newCoord<-do.call(rbind.data.frame,newCoord)
  out[newCoord$cell]<-as.character(newCoord$SOM)
  out<-factor(out,levels = c('00',levels(SOM_Label)))
  return(out)
})


RU_mat<-lapply(setNames(uids,uids)[1],function(u,
                                            frst=rst,
                                            fsql=geo_list_RU_annot){
  ncl<-raster::ncol(frst[[u]])
  nrw<-raster::nrow(frst[[u]])
  
  out<-rep('00',ncl*nrw)
  newCoord<-exactextractr::exact_extract(frst[[u]][[1]],
                                         fsql[fsql$uid==u,],
                                         include_cell = T,
                                         include_cols = 'SOM')
  newCoord<-do.call(rbind.data.frame,newCoord)
  out[newCoord$cell]<-as.character(newCoord$SOM)
  out<-factor(out,levels = c('00',levels(SOM_Label)))
  return(out)
})

tab_mat<-lapply(setNames(uids,uids)[1], function(u){
  table(BO_mat[[u]],RU_mat[[u]],dnn = c('BO','RU'))
})

sum(diag(tab_mat[[1]]))/sum(tab_mat[[1]])
fsomMFI[c(15,5),]

hTemp<-hclust(dist(fsomMFI[,c('CD4','CD8','CD68','CD14','E-CAD','PANK','CD20')],method = 'euclidean'))
plot(hclust(dist(fsomMFI[,c('CD4','CD8','CD68','CD14','E-CAD','PANK','CD20')],method = 'euclidean')))
fsomMFI[c(16,11,7,13),]
fsomMFI[c(03,05,15,04,10),]
fsomMFI[c(13,19,24),]
fsomMFI[order(fsomMFI[,'CD4']),]

uu=4

BO_Area<-data.frame(SOM=geo_list_BO_annot[geo_list_BO_annot$uid==uids[uu],'SOM',drop=T],
                    Area=sf::st_area(geo_list_BO_annot[geo_list_BO_annot$uid==uids[uu],]))
BO_Area<-stats::aggregate(Area~SOM,BO_Area,sum,drop=F)
BO_Area$Area<-BO_Area$Area/sum(BO_Area$Area)*1000
BO_Area<-rbind.data.frame(BO_Area,
                          data.frame(SOM=levels(BO_Area$SOM)[!(levels(BO_Area$SOM) %in% unique(BO_Area$SOM))],
                                     Area = rep(0,length(which(!(levels(BO_Area$SOM) %in% unique(BO_Area$SOM)))))))
BO_Area<-BO_Area[order(BO_Area$SOM),]
BO_Area<-BO_Area[hTemp$order,]

RU_Area<-data.frame(SOM=geo_list_RU_annot[geo_list_RU_annot$uid==uids[uu],'SOM',drop=T],
                    Area=sf::st_area(geo_list_RU_annot[geo_list_RU_annot$uid==uids[uu],]))
RU_Area<-stats::aggregate(Area~SOM,RU_Area,sum,drop=F)
RU_Area$Area<-RU_Area$Area/sum(RU_Area$Area)*1000
RU_Area<-rbind.data.frame(RU_Area,
                          data.frame(SOM=levels(RU_Area$SOM)[!(levels(RU_Area$SOM) %in% unique(RU_Area$SOM))],
                                     Area = rep(0,length(which(!(levels(RU_Area$SOM) %in% unique(RU_Area$SOM)))))))
RU_Area<-RU_Area[order(RU_Area$SOM),]
RU_Area<-RU_Area[hTemp$order,]

out<-rbind(BO_Area$Area,RU_Area$Area)
colnames(out)<-RU_Area$SOM






BO_Area<-data.frame(SOM=geo_list_BO_annot[,'SOM',drop=T],
                    uid=geo_list_BO_annot[,'uid',drop=T],
                    Area=sf::st_area(geo_list_BO_annot))
BO_total<-stats::aggregate(Area~uid,BO_Area,sum,drop=F)
BO_Area<-stats::aggregate(Area~SOM+uid,BO_Area,sum,drop=F)
for (u in uids){
  BO_Area$Area[BO_Area$uid==u]<-BO_Area$Area[BO_Area$uid==u]/BO_total$Area[BO_total$uid==u]
}
BO_Area$Area[is.na(BO_Area$Area)]<-0
BO_Area<-cbind.data.frame(Method='BO',BO_Area)
BO_Area$SOM<-as.factor(BO_Area$SOM)
BO_Area$Method<-as.factor(BO_Area$Method)
BO_Area$uid<-as.factor(BO_Area$uid)

RU_Area<-data.frame(SOM=geo_list_RU_annot[,'SOM',drop=T],
                    uid=geo_list_RU_annot[,'uid',drop=T],
                    Area=sf::st_area(geo_list_RU_annot))
RU_total<-stats::aggregate(Area~uid,RU_Area,sum,drop=F)
RU_Area<-stats::aggregate(Area~SOM+uid,RU_Area,sum,drop=F)
for (u in uids){
  RU_Area$Area[RU_Area$uid==u]<-RU_Area$Area[RU_Area$uid==u]/RU_total$Area[RU_total$uid==u]
}
RU_Area$Area[is.na(RU_Area$Area)]<-0
RU_Area<-cbind.data.frame(Method='RU',RU_Area)
RU_Area$SOM<-as.factor(RU_Area$SOM)
RU_Area$Method<-as.factor(RU_Area$Method)
RU_Area$uid<-as.factor(RU_Area$uid)

TOT_Area<-rbind.data.frame(BO_Area,RU_Area)

tt<-vapply(setNames(levels(TOT_Area$SOM),levels(TOT_Area$SOM)),function(x){
  out<-t.test(TOT_Area$Area[TOT_Area$SOM==x & TOT_Area$Method=='BO'],
         TOT_Area$Area[TOT_Area$SOM==x & TOT_Area$Method=='RU'],
         paired=T)
  out$p.value
},c(1))

tt<-p.adjust(tt,'bonferroni')
formatC(tt,digits = 2,mode = 'double',drop0trailing = T)









cs<-lapply(uids,function(x){
  browser()
  inp<-as.table(rbind(BO_Area$Area[BO_Area$uid==x],RU_Area$Area[BO_Area$uid==x]))
  dimnames(inp)<-list(method=c('BO','RU'),SOM=levels(BO_Area$SOM))
  chisq.test(inp)
})

cs<-lapply(setNames(uids,uids),function(x){
  BO<-table(geo_list_BO_annot$SOM[geo_list_BO_annot$uid==x])
  BO<-BO/sum(BO)*100
  RU<-table(geo_list_RU_annot$SOM[geo_list_RU_annot$uid==x])
  RU<-RU/sum(RU)*100
  
  rbind.data.frame(data.frame(val=BO,method='BO',SOM=names(BO)),
                   data.frame(val=RU,method='RU',SOM=names(RU)))
})

cs<-do.call(rbind,cs)
an<-aov(val.Freq~SOM:method,cs)
out<-TukeyHSD(an)
out

  