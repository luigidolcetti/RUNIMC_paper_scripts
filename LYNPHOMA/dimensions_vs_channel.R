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

####### with rescale #######

geo_list_BO_annot<-dplyr::bind_cols(geo_list_BO[,c('uid')],
                                    geo_list_BO_trans,
                                    data.frame(area = sf::st_area(geo_list_BO)),
                                    data.frame(SOM=SOM_Label[1:nrow(geo_list_BO)]))

colnames(geo_list_BO_annot)[2:16]<-vapply(strsplit(colnames(geo_list_BO_annot)[2:16],'.',fixed = T),function(x)x[2],FUN.VALUE = character(1))

geo_list_RU_annot<-dplyr::bind_cols(geo_list_RU[,c('uid')],
                                    geo_list_RU_trans,
                                    data.frame(area = sf::st_area(geo_list_RU)),
                                    data.frame(SOM=SOM_Label[(1+nrow(geo_list_BO)):(nrow(geo_list_BO)+nrow(geo_list_RU))]))

colnames(geo_list_RU_annot)[2:16]<-vapply(strsplit(colnames(geo_list_RU_annot)[2:16],'.',fixed = T),function(x)x[2],FUN.VALUE = character(1))

uids<-unique(geo_list_BO$uid)

geo_list_BO_annot[,2:16]<-scale(geo_list_BO_annot[,2:16,drop=T])

BO_stats<-lapply(levels(SOM_Label),function(x,y=geo_list_BO_annot){
  colToConsider<-colnames(y)[c(3,5,6,8,9,10,11,12)]
  TEMP<-y[y$SOM==x,colToConsider,drop=T]
  meanCol<-colMeans(TEMP)
  sdCol<-apply(TEMP,2,sd)
  leader<-names(meanCol)[which.max(meanCol)]
  return(list(leader=leader,
              means = meanCol,
              sd = sdCol))
})

BO_leaders<-do.call(cbind,lapply(BO_stats,function(x){
  nm<-x$leader
  return(matrix(c(x$means[nm],x$sd[nm]),ncol = 1,nrow=2,dimnames = list(c('mn','sd'),nm)))
}))

BO_leaders_stats<-do.call(cbind,lapply(unique(colnames(BO_leaders)),function(x,y=BO_leaders){
  TEMP<-y[,x,drop=F]
  TEMP<-TEMP[,which.max(TEMP['mn',]),drop=F]
}))

BO_purity<-lapply(setNames(1:length(levels(SOM_Label)),levels(SOM_Label)),function(x,y=geo_list_BO_annot,z=BO_leaders,zz=BO_leaders_stats){
  leader<-colnames(z)[x]
  cl<-colnames(zz)[colnames(zz)!=leader]
  out<-vapply(setNames(cl,cl),function(ll){
    TEMP<-y[y$SOM==levels(SOM_Label)[x],ll,drop=T]
    mn<-zz['mn',ll]
    sd<-zz['sd',ll]
    trsh<-mn-3*sd
    length(TEMP[TEMP<trsh])/length(TEMP)*100
  },c(1))
})


geo_list_RU_annot[,2:16]<-scale(geo_list_RU_annot[,2:16,drop=T])

RU_stats<-lapply(levels(SOM_Label),function(x,y=geo_list_RU_annot){
  colToConsider<-colnames(y)[c(3,5,6,8,9,10,11,12)]
  TEMP<-y[y$SOM==x,colToConsider,drop=T]
  meanCol<-colMeans(TEMP)
  sdCol<-apply(TEMP,2,sd)
  leader<-names(meanCol)[which.max(meanCol)]
  return(list(leader=leader,
              means = meanCol,
              sd = sdCol))
})

RU_leaders<-do.call(cbind,lapply(RU_stats,function(x){
  nm<-x$leader
  return(matrix(c(x$means[nm],x$sd[nm]),ncol = 1,nrow=2,dimnames = list(c('mn','sd'),nm)))
}))

RU_leaders_stats<-do.call(cbind,lapply(unique(colnames(RU_leaders)),function(x,y=RU_leaders){
  TEMP<-y[,x,drop=F]
  TEMP<-TEMP[,which.max(TEMP['mn',]),drop=F]
}))

RU_purity<-lapply(setNames(1:length(levels(SOM_Label)),levels(SOM_Label)),function(x,y=geo_list_RU_annot,z=RU_leaders,zz=RU_leaders_stats){
  leader<-colnames(z)[x]
  cl<-colnames(zz)[colnames(zz)!=leader]
  out<-vapply(setNames(cl,cl),function(ll){
    TEMP<-y[y$SOM==levels(SOM_Label)[x],ll,drop=T]
    mn<-zz['mn',ll]
    sd<-zz['sd',ll]
    trsh<-mn-3*sd
    length(TEMP[TEMP<trsh])/length(TEMP)*100
  },c(1))
})

final<-do.call(rbind,lapply(levels(SOM_Label),function(x){
  c(mean(BO_purity[[x]]),mean(RU_purity[[x]]))
}))
colMeans(final)
t.test(final[,1],final[,2],paired = T)


######### with 0-1 range #####


geo_list_BO_annot<-dplyr::bind_cols(geo_list_BO[,c('uid')],
                                    geo_list_BO_trans,
                                    data.frame(area = sf::st_area(geo_list_BO)),
                                    data.frame(SOM=SOM_Label[1:nrow(geo_list_BO)]))

colnames(geo_list_BO_annot)[2:16]<-vapply(strsplit(colnames(geo_list_BO_annot)[2:16],'.',fixed = T),function(x)x[2],FUN.VALUE = character(1))

geo_list_RU_annot<-dplyr::bind_cols(geo_list_RU[,c('uid')],
                                    geo_list_RU_trans,
                                    data.frame(area = sf::st_area(geo_list_RU)),
                                    data.frame(SOM=SOM_Label[(1+nrow(geo_list_BO)):(nrow(geo_list_BO)+nrow(geo_list_RU))]))

colnames(geo_list_RU_annot)[2:16]<-vapply(strsplit(colnames(geo_list_RU_annot)[2:16],'.',fixed = T),function(x)x[2],FUN.VALUE = character(1))

uids<-unique(geo_list_BO$uid)

geo_list_BO_annot[,2:16]<-apply(geo_list_BO_annot[,2:16,drop=T],2,function(x){(x-min(x))/(max(x)-min(x))})

geo_list_BO_annot<-geo_list_BO_annot[geo_list_BO_annot$area>0 & geo_list_BO_annot$area<50,]

BO_purity<-lapply(levels(SOM_Label),function(x,y=geo_list_BO_annot){
  colToConsider<-colnames(y)[c(3,5,6,8,9,10,11,12)]
  TEMP<-y[y$SOM==x,colToConsider,drop=T]
  out<-apply(TEMP,1,function(x){
    x<-x/sum(x)
    x<-x[x>0]
    out<--sum(x*log(x,2))
    
  })
  return(mean(out))
})


geo_list_RU_annot[,2:16]<-apply(geo_list_RU_annot[,2:16,drop=T],2,function(x){(x-min(x))/(max(x)-min(x))})

geo_list_RU_annot<-geo_list_RU_annot[geo_list_RU_annot$area>0 & geo_list_RU_annot$area<50,]

RU_purity<-lapply(levels(SOM_Label),function(x,y=geo_list_RU_annot){
  colToConsider<-colnames(y)[c(3,5,6,8,9,10,11,12)]
  TEMP<-y[y$SOM==x,colToConsider,drop=T]
  out<-apply(TEMP,1,function(x){
    x<-x/sum(x)
    x<-x[x>0]
    out<--sum(x*log(x,2))
  })
  return(mean(out))
})

t.test(unlist(BO_purity),unlist(RU_purity),paired = T)






plot(geo_list_BO_annot[geo_list_BO_annot$SOM=='22',c('histone','area'),drop=T],cex=0.1,xlim=c(0,4),ylim=c(0,200))

plot(geo_list_RU_annot[geo_list_RU_annot$SOM=='22',c('histone','area'),drop=T],cex=0.1,xlim=c(0,4),ylim=c(0,200))



plot(geo_list_BO_annot[geo_list_BO_annot$uid==uids[3] & geo_list_BO_annot$SOM=='22','GEOMETRY'])
plot(geo_list_RU_annot[geo_list_RU_annot$uid==uids[3] & geo_list_RU_annot$SOM=='22','area'],border=NA,xlim=c(100,300),ylim=c(0,200))

TEMP<-RUNIMC:::retrieve.RsCollection(
  "C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/rasterStacks"
)

raster::plot(RUNIMC::quantNorm(TEMP[[3]]$x142nd.asma.nd142di.,0.9),col=gray(level = seq(0,1,0.1)),xlim=c(200,400),ylim=c(0,200))
plot(geo_list_RU_annot[geo_list_RU_annot$uid==uids[3] & geo_list_RU_annot$SOM=='19','GEOMETRY'],col=NA,border='red',xlim=c(100,300),ylim=c(0,200),add=T)
