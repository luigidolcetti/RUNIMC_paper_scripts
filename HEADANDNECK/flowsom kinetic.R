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

uids<-unique(geo_list_BO$uid)
meta=10
a=0
bioDiv<-pbapply::pblapply(1:20,function(dimxy){
  dimx=dimxy
  dimy=dimxy
    out<-pbapply::pblapply(1:10, function(rpt){
      
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
                              importance=importanceX)
      fsom<-FlowSOM::BuildMST(fsom,silent = T)
      
      fsomMFI<-FlowSOM::GetClusterMFIs(fsom)
      
      fsomMFI<-apply(fsomMFI,2,function(x){(x-min(x))/(max(x)-min(x))})
      
      SOM_Label<-FlowSOM::GetClusters(fsom)
      
      SOM_Label<-factor(as.character(formatC(SOM_Label,flag='0',format = 'd',digits=1)))
      
      
      geo_list_BO_annot<-dplyr::bind_cols(geo_list_BO,data.frame(SOM=SOM_Label[1:nrow(geo_list_BO)]))
      
      
      geo_list_RU_annot<-dplyr::bind_cols(geo_list_RU,data.frame(SOM=factor(as.character(SOM_Label[(1+nrow(geo_list_BO)):(nrow(geo_list_BO)+nrow(geo_list_RU))]))))
      
     
      RU<-vapply(uids,function(x){
        RU_TEMP<-table(geo_list_RU_annot$SOM[geo_list_RU_annot$uid==x,drop=T])
        RU_TEMP<-RU_TEMP/sum(RU_TEMP)
        if (a==1) { out<-exp(-sum(RU_TEMP*log(RU_TEMP)))} else {
          out<-sum(RU_TEMP^a)^(1/(1-a))
        }
      },c(1))
      
      
      BO<-vapply(uids,function(x){
        BO_TEMP<-table(geo_list_BO_annot$SOM[geo_list_BO_annot$uid==x,drop=T])
        BO_TEMP<-BO_TEMP/sum(BO_TEMP)
        if (a==1) { out<-exp(-sum(BO_TEMP*log(BO_TEMP)))} else {
          out<-sum(BO_TEMP^a)^(1/(1-a))
        }
      },c(1))
      
      wt<-wilcox.test(RU,BO,paired = T)
      tt<-t.test(RU,BO,paired = T)
      
      rt<-list(BO=BO,
               RU=RU,
               mBO=mean(BO),
               mRU=mean(RU),
               r=log(mean(BO)/mean(RU),2),
               wt=wt,
               tt=tt)
      return(rt)
    })

})


rr<-unlist(lapply(bioDiv,function(x){x[[1]]$tt$p.value}),recursive = T)
           

rr<-lapply(1:length(bioDiv),function (x){
  out<-lapply(bioDiv[[x]],function(y){
    c(x,y$mBO,y$mRU)
  })
  do.call(rbind,out)
}) 
rr<-do.call(rbind,rr)

plot(NA,xlim=c(0,max(rr[,1])),ylim=c(0,max(rr[,2:3])))

points(rr[,c(1,2)],col='red',cex=0.1)
points(rr[,c(1,3)],col='blue',cex=0.1)
