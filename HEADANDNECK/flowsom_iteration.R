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

mrkr<-c('CD38','SMA','VIM',
        'CD14','CD33','CD16',
        'PANK','CD11b','CD15',
        'IgD','CD45','CD24',
        'CD11c','FoxP3','CD4',
        'E-CAD','CD68','CD20',
        'CD8','Arg1','CD45ra',
        'CD74','CD127','COLL',
        'CD3','CD27','CD45ro',
        'CD25','DNA1','DNA2')

LLL<-expand.grid(LETTERS,LETTERS)
LLL<-apply(LLL,1,paste0,collapse='')[1:length(mrkr)]
names(LLL)<-mrkr
names(mrkr)<-LLL[1:length(mrkr)]


out<-pbapply::pblapply(1:1000, function(iii){
  
  
  fsom<-FlowSOM::ReadInput(geo_list_total,
                           compensate = F,
                           transform = F,
                           scale = T,
                           silent = T)
  fsom<-FlowSOM::BuildSOM(fsom = fsom,
                          colsToUse = colToUse,
                          silent = T,
                          xdim=dimx,
                          ydim=dimy,
                          # rlen=10,
                          # mst = 30,
                          # distf = 1,
                          importance=importanceX)
  fsom<-FlowSOM::BuildMST(fsom,silent = T,tSNE = F)
  
  
  fsomMFI<-FlowSOM::GetClusterMFIs(fsom)
  
  fsomMFI<-apply(fsomMFI,2,function(x){(x-min(x))/(max(x)-min(x))})
  
  colnames(fsomMFI)<-mrkr
  
  rownames(fsomMFI)<-as.character(formatC(1:nrow(fsomMFI),flag='0',digits = 1,format = 'd'))
  
  SOM_Label<-FlowSOM::GetClusters(fsom)
  
  SOM_Label<-factor(as.character(formatC(SOM_Label,flag='0',format = 'd',digits=1)))
  
  geo_list_BO_annot<-dplyr::bind_cols(geo_list_BO,data.frame(SOM=factor(SOM_Label[1:nrow(geo_list_BO)])))
  
  geo_list_RU_annot<-dplyr::bind_cols(geo_list_RU,data.frame(SOM=factor(as.character(SOM_Label[(1+nrow(geo_list_BO)):(nrow(geo_list_BO)+nrow(geo_list_RU))]))))
  
  uids<-unique(geo_list_BO_annot$uid)
  
  subsetMakers<-c('CD4','CD8','CD68','E-CAD','PANK','CD20')
  hTemp<-hclust(dist(fsomMFI[,subsetMakers],method = 'euclidean'))
  
  ranK<-lapply(hTemp$order,function(x){
    names(fsomMFI[x,order(fsomMFI[x,mrkr[1:28]],decreasing = T)])
  })
  ranK<-do.call(cbind,ranK)
  colnames(ranK)<-levels(SOM_Label)[hTemp$order]
  ranKLLL<-matrix(unname(LLL[ranK]),nrow=nrow(ranK),ncol=ncol(ranK),dimnames = dimnames(ranK))
  ranKLLL<-apply(ranKLLL,2,paste0,collapse='.')

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
  tt1<-formatC(tt,digits = 2,mode = 'double',drop0trailing = T)
  
  out<-list(POP=ranKLLL,
            pval=tt)
  return(out)
})


# ALL_out<-lapply(out,'[[','POP')
# 
# ALL_out<-do.call(c,ALL_out)
# 
# ALL_out<-substr(ALL_out,1,15)
# All_TAB<-table(ALL_out)
# All_TAB<-All_TAB[order(All_TAB,decreasing = T)]
# All_TAB[1:10]
# max(All_TAB)
# length(All_TAB)
# plot(All_TAB[1:100])
# 
# DescTools::StrDist(names(All_TAB)[1],names(All_TAB)[2])
# 
# TEMP<-names(All_TAB)[4]
# TEMP<-strsplit(TEMP,'.',fixed = T)
# mrkr[TEMP[[1]]]
# 
# All_names<-names(All_TAB)
# 
# TEMP<-cbind.data.frame(expand.grid(All_names[1:10],All_names[1:10]),dist=NA)
# 
# TEMP$dist<-apply(TEMP,1,function(x) phrDist(x[1],x[2],split='.'))
# phrDist<-function(x,y,split){
#   xs<-strsplit(x,split,fixed = T)[[1]]
#   ys<-strsplit(y,split,fixed = T)[[1]]
#   ds<-unlist(lapply(1:length(xs),function(i){
#     ws<-which(ys==xs[i])
#     if (length(ws)==0) {
#       out<-length(ys)
#     } else {
#       out<-abs(i-ws)
#     }
#     return(out)
#   }))
#   return(sum(ds))
# }

# TEMP<-TEMP[TEMP$dist!=0,]
# TEMP<-TEMP[order(TEMP$dist),]
# TEMP<-TEMP[!duplicated(TEMP[,1]),]
# 
# TEMP_net<-igraph::graph(edges = as.vector(t(TEMP[,1:2]))) 
# plot(TEMP_net)
# 

# 
# TEMP<-cbind.data.frame(expand.grid(ALL_out[[1]],ALL_out[[2]]),dist=NA)
# 
# TEMP$dist<-apply(TEMP,1,function(x) phrDist(x[1],x[2],split='.'))
# 
# 
# TEMP<-tidyr::pivot_wider(data = TEMP,
#                          names_from = 'Var2',
#                          values_from = 'dist')
# TEMP1<-as.matrix(TEMP[,-1])
# rownames(TEMP1)<-TEMP[,1]
# colnames(TEMP1)<-NULL
# rankCol<-apply(TEMP1,2,rank)
# rankRow<-apply(TEMP1,1,rank)
# rankTot<-rankCol*rankRow
# rankService<-rankTot
# rankVector<-vector('numeric',ncol(rankTot))
# for (i in 1:nrow(rankTot)){
#   currentRow<-rankService[i,]
#   ww<-which.min(currentRow)
#   rankVector[i]<-ww
#   rankService[,ww]<-rep(+Inf,nrow(rankService))
# }
# rankTot<-rankTot[,rankVector]
# mtch<-cbind(ALL_out[[1]],ALL_out[[2]][rankVector])[20,]
# mrkr[strsplit(mtch[1],'.',fixed = T)[[1]]]
# mrkr[strsplit(mtch[2],'.',fixed = T)[[1]]]
# 
# 
# 
# ####### pictogram #####
# ####### 
# 
# 
# P_Out<-lapply(out,function(x){
#   x$POP[x$pval<0.05]
# })
# 
# P_Out<-do.call(c,P_Out)
# P_Out<-lapply(P_Out,function(x){strsplit(x,'.',fixed=T)[[1]]})
# P_Out<-do.call(rbind,P_Out)
# P_Out<-apply(P_Out,2,table)
# 
# 
# 
# 
# txt<-'CIAO'
# Cairo::Cairo(file='/dev/null')
#              # units = 'px',
#              # height = 100,width = 100*nchar(txt),bg = 'green')
# par(mar=c(0,0,0,0),bty='n',family='sans')
# plot(NA,xlim=c(0,1),ylim=c(0,nchar(txt)))
# cx<-par('pin')[1]/strwidth(txt,'inches')
# text(0.5,nchar(txt)/2,txt,cex=cx)
# i = Cairo:::.image(dev.cur())
# r = Cairo:::.ptr.to.raw(i$ref, 0, i$width * i$height * 4)
# dev.off()
# dim(r) = c(4,i$width, i$height) # RGBA planes
# r1<-aperm(r,c(3,2,1))
# plot(NA,xlim=c(0,100),ylim=c(0,100))
# rasterImage(as.raster(r1),0,0,20,20,interpolate = T)
# dev.off()



#####  areaPlot ####
#####  

P_Out<-lapply(out,function(x){
  x$POP[x$pval<0.0001]
})

P_Out<-unlist(P_Out)
P_Out<-lapply(P_Out,function(x){strsplit(x,'.',fixed=T)[[1]]})
P_Out<-do.call(rbind,P_Out)
P_Out<-apply(P_Out,2,table)
P_Out<-lapply(P_Out,function(x)x[order(x,decreasing = T)])
P_Out<-lapply(P_Out,function(x)x/sum(x))
P_levels<-unique(names(unlist(P_Out)))
P_colors<-setNames(rainbow(length(P_levels)),P_levels)

plot(NA,xlim=c(1,length(P_Out)),ylim=c(0,1))
for(c in 1:length(P_Out)){
  P_extra<-c(P_Out[[c]],0)
  for (m in 1:length(P_Out[[c]])){
    rect(c-0.5,
         sum(P_extra[m:length(P_extra)]),
         c+0.5,
         sum(P_extra[(m+1):length(P_extra)]),
         col = P_colors[names(P_Out[[c]][m])])
    text(c, sum(P_extra[m:length(P_extra)]),labels = mrkr[names(P_Out[[c]][m])],
         c(0.5,1),cex = 0.8)
  }
}
  
  
  
