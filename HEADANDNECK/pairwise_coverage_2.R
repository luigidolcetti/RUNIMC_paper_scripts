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

subsetMakers<-c('CD4','CD8','CD68','E-CAD','PANK','CD20')
hTemp<-hclust(dist(fsomMFI[,subsetMakers],method = 'euclidean'))

ranK<-lapply(hTemp$order,function(x){
  names(fsomMFI[x,order(fsomMFI[x,],decreasing = T)])
})
ranK<-do.call(cbind,ranK)
colnames(ranK)<-levels(SOM_Label)[hTemp$order]


BO_Area<-data.frame(SOM=geo_list_BO_annot[,'SOM',drop=T],
                    uid=geo_list_BO_annot[,'uid',drop=T],
                    Area=sf::st_area(geo_list_BO_annot))
BO_total<-stats::aggregate(Area~uid,BO_Area,sum,drop=F)
BO_Area<-stats::aggregate(Area~SOM+uid,BO_Area,sum,drop=F)
# for (u in uids){
#   BO_Area$Area[BO_Area$uid==u]<-BO_Area$Area[BO_Area$uid==u]/BO_total$Area[BO_total$uid==u]
# }
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
# for (u in uids){
#   RU_Area$Area[RU_Area$uid==u]<-RU_Area$Area[RU_Area$uid==u]/RU_total$Area[RU_total$uid==u]
# }
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
tt1<-formatC(tt,digits = 2,mode = 'double',drop0trailing = T)
tt1<-gtools::stars.pval(tt)

plot.new()
par(mfrow=c(2,1),mar=c(0,0,0,0))
plot(hclust(dist(fsomMFI[,subsetMakers],method = 'euclidean')),
     xaxt='n',
     yaxt='n')
plot(NA,xlim=c(1,length(levels(SOM_Label))),
     ylim=c(0,nrow(ranK)+1),
     xaxt='n',
     yaxt='n')
for (cx in seq_along(hTemp$order)){
text (x=cx,y=1,tt1[hTemp$order[cx]],cex=0.9)
    for (rx in 1:nrow(ranK)){
    text(x=cx,y=nrow(ranK)-rx+2,ranK[rx,cx],cex=0.6)
  }
}

plot.new()
par(mfrow=c(dimx,dimy),mar=c(1,2,1,1))
for (i in levels(SOM_Label)){
  
  perc_lim<-c(min(TOT_Area$Area[TOT_Area$SOM==i]),max(TOT_Area$Area[TOT_Area$SOM==i]))
  plot(NA,xlim=c(0.5,2.5),ylim=perc_lim,main=paste(i,tt1[i],sep=' '),xaxt='n',yaxt='n')
  axis(2,at = pretty(perc_lim),labels = scales::percent(pretty(perc_lim)))
  points(x=rep(1,length(uids)),
         y=TOT_Area$Area[TOT_Area$SOM==i & TOT_Area$Method=='BO'],
         pch=16)
  points(x=rep(2,length(uids)),
         y=TOT_Area$Area[TOT_Area$SOM==i & TOT_Area$Method=='RU'],
         pch=17)
  lines(matrix(c(0.8,1.2,
                mean(TOT_Area$Area[TOT_Area$SOM==i & TOT_Area$Method=='BO']),
                mean(TOT_Area$Area[TOT_Area$SOM==i & TOT_Area$Method=='BO'])),
               ncol=2,
               byrow = F),
        lwd=2)
  lines(matrix(c(1.8,2.2,
                 mean(TOT_Area$Area[TOT_Area$SOM==i & TOT_Area$Method=='RU']),
                 mean(TOT_Area$Area[TOT_Area$SOM==i & TOT_Area$Method=='RU'])),
               ncol=2,
               byrow = F),
        lwd=2)
  for (u in uids){
    lines(matrix(c(1,2,
      TOT_Area$Area[TOT_Area$SOM==i & TOT_Area$Method=='BO' &TOT_Area$uid==u],
      y=TOT_Area$Area[TOT_Area$SOM==i & TOT_Area$Method=='RU' &TOT_Area$uid==u]),
      ncol = 2, byrow = F
    ))
  }
}

plot.new()
par(mfrow=c(dimx,dimy*2),mar=c(0,0,1,0))
for (u in uids[10]){
for (i in levels(SOM_Label)){
  
  BO_sub<-geo_list_BO_annot[geo_list_BO_annot$uid==u
                            & geo_list_BO_annot$SOM==i,'GEOMETRY']
  RU_sub<-geo_list_RU_annot[geo_list_RU_annot$uid==u
                            & geo_list_RU_annot$SOM==i,'GEOMETRY']
  BO_tot<-sf::st_union(geo_list_BO_annot[geo_list_BO_annot$uid==u,'GEOMETRY'])
  RU_tot<-sf::st_union(geo_list_RU_annot[geo_list_RU_annot$uid==u,'GEOMETRY'])
  
  
  plot(BO_tot,main= paste('Ref: ',i,tt1[i]),col='gray80',border=NA,lwd=1,add=F,cex=0.5)
  if (nrow(BO_sub)>0){
    plot(BO_sub,col='gray30',border='black',lwd=0.1,add=T,cex=0.5)
  }
  
  plot(RU_tot,main= paste('Alt: ',i,tt1[i]),col='gray80',border=NA,lwd=1,add=F,cex=0.5)
  if (nrow(RU_sub)>0){
    plot(RU_sub,col='gray30',border='black',lwd=0.1,add=T,cex=0.5)
  }
}
}
