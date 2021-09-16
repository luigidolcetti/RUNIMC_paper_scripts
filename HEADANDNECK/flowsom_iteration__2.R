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
  
  SOM_Lab<-vector('list',2)
  SOM_Lab[['BO']]=factor(SOM_Label[1:nrow(geo_list_BO)])

  SOM_Lab[['RU']]=factor(as.character(SOM_Label[(1+nrow(geo_list_BO)):(nrow(geo_list_BO)+nrow(geo_list_RU))]))
  
  return(list(MFI=fsomMFI,
              LAB=SOM_Lab))
})


TEMP<-lapply(out[-1],function(x){
  TEMP<-as.matrix(pdist::pdist(x$MFI,out[[1]]$MFI))
  dimnames(TEMP)<-list(rownames(out[[1]]$MFI),rownames(out[[1]]$MFI))
  re<-apply(TEMP,1,which.min)
  re<-setNames(levels(out[[1]]$LAB$BO)[re],levels(out[[1]]$LAB$BO))
})

TEMP1<-lapply(1:length(TEMP),function(i){
  newLAB<-list(BO=character(0),RU=character(0))
  newLAB$BO<-unname(TEMP[[i]][out[[i+1]]$LAB$BO])
  newLAB$BO<-factor(newLAB$BO)
  newLAB$RU<-unname(TEMP[[i]][out[[i+1]]$LAB$RU])
  newLAB$RU<-factor(newLAB$RU)
  return(newLAB)
})

comp_BO<-do.call(cbind,lapply(TEMP1,'[[','BO'))
comp_BO<-apply(comp_BO,1,function(x){
  y<-table(x)
  length(y)
  # y<-y/sum(y)
  # out<--sum(y*log(y)/log(length(y)))
})

comp_RU<-do.call(cbind,lapply(TEMP1,'[[','RU'))
comp_RU<-apply(comp_RU,1,function(x){
  y<-table(x)
  length(y)
  # y<-y/sum(y)
  # out<--sum(y*log(y)/log(length(y)))
})

mean(comp_BO)
mean(comp_RU)
boxplot(comp_BO,comp_RU)

plot(density(comp_BO))
lines(density(comp_RU),col='red')
