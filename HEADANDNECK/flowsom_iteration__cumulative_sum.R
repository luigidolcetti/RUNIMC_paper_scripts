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

geo_list_total<-as.matrix(rbind(geo_list_BO_trans))


cl<-parallel::makeCluster(24)
parallel::clusterExport(cl,
                        varlist = list(
                        'geo_list_total',
                        'colToUse'))

multiSom<-parallel::parLapply(1:20,function(ms){
  
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
  
  
  out<-pbapply::pblapply(1:100, function(iii){
    
    fsom<-FlowSOM::ReadInput(geo_list_total,
                             compensate = F,
                             transform = F,
                             scale = F,
                             silent = T)
    fsom<-FlowSOM::BuildSOM(fsom = fsom,
                            colsToUse = colToUse[-c(3,24,29,30)],
                            silent = T,
                            xdim=dimx,
                            ydim=dimy)
    fsom<-FlowSOM::BuildMST(fsom,silent = T,tSNE = F)
    
    
    SOM_Label<-FlowSOM::GetClusters(fsom)
    
    SOM_Label<-factor(as.character(formatC(SOM_Label,flag='0',format = 'd',digits=1)))
    
    return(SOM_Label)
  })
  
  out<-do.call(cbind,out)
  
  out_subset<-out[sample(1:nrow(out),1000),]
  bigMatrix<-matrix(0,ncol=nrow(out_subset),nrow = nrow(out_subset))
  totBig<-length(bigMatrix)
  pb<-progress::progress_bar$new(total = totBig)
  
  for (cc in 1:ncol(out_subset)){
    for (rr in 1:nrow(out_subset)){
      pointerClass<-out_subset[rr,cc]
      pointerOut<-which(out_subset[,cc]==pointerClass)
      bigMatrix[rr,pointerOut]<-bigMatrix[rr,pointerOut]+1
      pb$tick()
    }
  }
  
  table(bigMatrix[lower.tri(bigMatrix,diag = F)])
},cl = cl)

parallel::stopCluster(cl)

multiSom<-do.call(rbind,multiSom)

MS_BO<-multiSom


saveRDS(MS_BO,
        "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/flowsom_cum_sum/MS_BO.R")

geo_list_total<-as.matrix(rbind(geo_list_RU_trans))


cl<-parallel::makeCluster(24)
parallel::clusterExport(cl,
                        varlist = list(
                          'geo_list_total',
                          'colToUse'))

multiSom<-parallel::parLapply(1:20,function(ms){
  
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
  
  
  out<-pbapply::pblapply(1:100, function(iii){
    
    fsom<-FlowSOM::ReadInput(geo_list_total,
                             compensate = F,
                             transform = F,
                             scale = F,
                             silent = T)
    fsom<-FlowSOM::BuildSOM(fsom = fsom,
                            colsToUse = colToUse[-c(3,24,29,30)],
                            silent = T,
                            xdim=dimx,
                            ydim=dimy)
    fsom<-FlowSOM::BuildMST(fsom,silent = T,tSNE = F)
    
    
    SOM_Label<-FlowSOM::GetClusters(fsom)
    
    SOM_Label<-factor(as.character(formatC(SOM_Label,flag='0',format = 'd',digits=1)))
    
    return(SOM_Label)
  })
  
  out<-do.call(cbind,out)
  
  out_subset<-out[sample(1:nrow(out),1000),]
  bigMatrix<-matrix(0,ncol=nrow(out_subset),nrow = nrow(out_subset))
  totBig<-length(bigMatrix)
  pb<-progress::progress_bar$new(total = totBig)
  
  for (cc in 1:ncol(out_subset)){
    for (rr in 1:nrow(out_subset)){
      pointerClass<-out_subset[rr,cc]
      pointerOut<-which(out_subset[,cc]==pointerClass)
      bigMatrix[rr,pointerOut]<-bigMatrix[rr,pointerOut]+1
      pb$tick()
    }
  }
  
  table(bigMatrix[lower.tri(bigMatrix,diag = F)])
},cl = cl)

parallel::stopCluster(cl)

multiSom<-do.call(rbind,multiSom)

MS_RU<-multiSom

saveRDS(MS_RU,
        "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/flowsom_cum_sum/MS_RU.R")

MS_BO_CUMSUM<-apply(MS_BO,1,function(x){
  out<-cumsum(x)
  out<-(out-min(out))/(max(out)-min(out))
})
MS_RU_CUMSUM<-apply(MS_RU,1,function(x){
  out<-cumsum(x)
  out<-(out-min(out))/(max(out)-min(out))
})

MS_BO_CI<-apply(MS_BO_CUMSUM,1,Rmisc::CI)
MS_RU_CI<-apply(MS_RU_CUMSUM,1,Rmisc::CI)
postscript("C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/flowsom_cum_sum/flow_cumsum.eps",
           onefile = F,
           width = 5,
           height = 5,
           horizontal = F,
           bg='white',
           paper = 'special')
plot(NA,xlim=c(1,101),ylim=c(0,1))
lines(MS_BO_CI[1,],col='red',lty=2)
lines(MS_BO_CI[2,],col='red')
lines(MS_BO_CI[3,],col='red',lty=2)
lines(MS_RU_CI[1,],col='blue',lty=2)
lines(MS_RU_CI[2,],col='blue')
lines(MS_RU_CI[3,],col='blue',lty=2)

lines(x = matrix(c(1,101,0,1),ncol=2,byrow = F),col='black')

dev.off()
