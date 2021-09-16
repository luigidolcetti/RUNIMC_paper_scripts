targetFolder<-"C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/SUBCELLULAR"

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

saveRDS(fsom,
        file = file.path(targetFolder,'som.R'  ))

fsomMFI<-FlowSOM::GetClusterMFIs(fsom)

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

SOM_Label<-FlowSOM::GetClusters(fsom)

SOM_Label<-factor(as.character(formatC(SOM_Label,flag='0',format = 'd',digits=1)),levels = rownames(fsomMFI))

hTemp<-hclust(dist(fsomMFI[,c('CD4','CD8','CD68','CD14','E-CAD','PANK','CD20')],method = 'euclidean'))

postscript(file=file.path(targetFolder,'SOM_heatmap.eps'),
           onefile = F,
           horizontal = F,
           paper = 'special',
           width = 8,
           height = 8,
           bg = 'white',
           pointsize = 18
)
pheatmap::pheatmap(t(fsomMFI),
                   scale = 'none',
                   cellwidth = 15,
                   cellheight = 15,
                   cluster_cols = hTemp,
                   angle_col = 0)

dev.off()

geo_list_BO_annot<-dplyr::bind_cols(geo_list_BO,data.frame(SOM=SOM_Label[1:nrow(geo_list_BO)]))

geo_list_RU_annot<-dplyr::bind_cols(geo_list_RU,data.frame(SOM=factor(as.character(SOM_Label[(1+nrow(geo_list_BO)):(nrow(geo_list_BO)+nrow(geo_list_RU))]))))


prePostTable<-table(sf::st_drop_geometry(geo_list_RU_annot)[,c('primary_id','SOM')])

prePostTable<-apply(prePostTable,2,function(x)x/sum(x))

postscript(file=file.path(targetFolder,'PrePost_heatmap.eps'),
           onefile = F,
           horizontal = F,
           paper = 'special',
           width = 8,
           height = 3,
           bg = 'white',
           pointsize = 18
)

pheatmap::pheatmap(prePostTable,
                   scale = 'none',
                   cellwidth = 15,
                   cellheight = 15,
                   angle_col = 0)

dev.off()

uids<-unique(geo_list_BO$uid)


for (u in uids){
  postscript(file = paste0(targetFolder,'/',u,"_BO_map.eps"),
             height = 30,
             width = 30,
             onefile = F,
             paper = 'special',
             horizontal = F,
             bg = 'white',
             pointsize = 18)
  par(mar=c(1,1,1,1))
  plot(geo_list_BO_annot[geo_list_BO_annot$uid==u,]['SOM'],
       key.pos = 1,
       key.length = 1,
       key.width = 0.05,
       lwd=0.1)
  dev.off()
  postscript(file = paste0(targetFolder,'/',u,"_RU_map.eps"),
             height = 30,
             width = 30,
             onefile = F,
             paper = 'special',
             horizontal = F,
             bg = 'white',
             pointsize = 18)
  plot(geo_list_RU_annot[geo_list_RU_annot$uid==u,]['SOM'],
       key.pos = 1,
       key.length = 1,
       key.width = 0.05,
       lwd=0.1)
  dev.off()
}


rst<-RUNIMC:::retrieve.RsCollection(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/rasterStacks"
)

out_p.val<-list()
out_r.val<-list()
for (u in 1:length(uids)){
  trns<-scales::modulus_trans(0)
  
  RU_test<-exactextractr::exact_extract(x = rst[[uids[u]]],
                                        y= geo_list_RU_annot[geo_list_RU_annot$uid==uids[u],],
                                        include_xy=T,
                                        include_cols='SOM')
  
  whichColumn<-colnames(RU_test[[1]])[2:31]
  setNames(whichColumn,whichColumn)
  
  RU_testCr<-lapply(RU_test,function(x){
    xc<-mean(x$x)
    yc<-mean(x$y)
    out<-apply(x,1,function(y){
      dist(matrix(c(xc,yc,y['x'],y['y']),ncol = 2,byrow = T))
    })
    out<-out/max(out)
    out<-cbind.data.frame(x,dist=out)
    cr<-lapply(whichColumn,function(wc){
      cor(out[,'dist'],trns$transform(out[,wc]))
    })
    # cor(out[,'dist'],out[,'value'])
    cr<-do.call(cbind,cr)
    colnames(cr)<-whichColumn
    out<-data.frame(SOM=out[1,'SOM'],cor=cr)
  })
  
  RU_testCr<-do.call(rbind.data.frame,RU_testCr)
  RU_mean<-aggregate(.~SOM,RU_testCr,mean)
  RU_sd<-aggregate(.~SOM,RU_testCr,sd)
  RU_tt<-aggregate(.~SOM,RU_testCr,function(x){
    if (length(x)>3) {
      out<-t.test(x,mu=0)
      out<-out$p.value} else out<-NA
      return(out)
  })
  RU_tt<-apply(RU_tt[,-1],2,p.adjust,method='bonferroni',n=length(levels(SOM_Label))*length(whichColumn))
  RU_tt<-cbind.data.frame(RU_mean[,1],gtools::stars.pval(RU_tt))
  colnames(RU_tt)<-colnames(RU_mean)
  
  BO_test<-exactextractr::exact_extract(x = rst[[uids[u]]],
                                        y= geo_list_BO_annot[geo_list_BO_annot$uid==uids[u],],
                                        include_xy=T,
                                        include_cols='SOM')
  whichColumn<-colnames(BO_test[[1]])[2:31]
  setNames(whichColumn,whichColumn)
  BO_testCr<-lapply(BO_test,function(x){
    xc<-mean(x$x)
    yc<-mean(x$y)
    out<-apply(x,1,function(y){
      dist(matrix(c(xc,yc,y['x'],y['y']),ncol = 2,byrow = T))
    })
    out<-out/max(out)
    out<-cbind.data.frame(x,dist=out)
    cr<-lapply(whichColumn,function(wc){
      cor(out[,'dist'],trns$transform(out[,wc]))
    })
    # cor(out[,'dist'],out[,'value'])
    cr<-do.call(cbind,cr)
    colnames(cr)<-whichColumn
    out<-data.frame(SOM=out[1,'SOM'],cor=cr)
  })
  
  BO_testCr<-do.call(rbind.data.frame,BO_testCr)
  BO_mean<-aggregate(.~SOM,BO_testCr,mean)
  BO_sd<-aggregate(.~SOM,BO_testCr,sd)
  BO_tt<-aggregate(.~SOM,BO_testCr,function(x){
    if (length(x)>3) {
      out<-t.test(x,mu=0)
      out<-out$p.value} else out<-NA
      return(out)
  })
  BO_tt<-apply(BO_tt[,-1],2,p.adjust,method='bonferroni',n=length(levels(SOM_Label))*length(whichColumn))
  BO_tt<-cbind.data.frame(BO_mean[,1],gtools::stars.pval(BO_tt))
  colnames(BO_tt)<-colnames(BO_mean)
  
  BO_n<-table(geo_list_BO_annot$SOM[geo_list_BO_annot$uid==uids[u]])
  BO_n<-round((BO_n-min(BO_n))/(max(BO_n)-min(BO_n))*(length(whichColumn)/2-2))+1
  RU_n<-table(geo_list_RU_annot$SOM[geo_list_RU_annot$uid==uids[u]])
  RU_n<-round((RU_n-min(RU_n))/(max(RU_n)-min(RU_n))*(length(whichColumn)/2-2))+1
  
  bg_matrix<-matrix('white',
                    ncol = length(whichColumn),
                    nrow = length(levels(SOM_Label)),
                    dimnames = list(levels(SOM_Label),
                                    whichColumn))
  for (i in levels(SOM_Label)){
    bg_matrix[i,1:BO_n[i]]<-'olivedrab'
    bg_matrix[i,length(whichColumn)-(1:RU_n[i])+1]<-'orange'
  }
  
  postscript(file = paste0(targetFolder,'/',uids[u],"_cor_map.eps"),
             height = 8,
             width = 20,
             onefile = F,
             paper = 'special',
             horizontal = F,
             bg = 'white',
             pointsize = 18)
  
  par(mfrow=c(length(levels(SOM_Label))+1,length(whichColumn)+1),
      mar=c(0,0,0,0))
  # par(mfrow=c(3,3),
  #     mar=c(0,0,0,0))
  plot(NA,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',bty='n')
  for (mr in whichColumn){
    plot(NA,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',bty='n')
    text(0.5,0.5,strsplit(mr,'.',fixed = T)[[1]][2])
  }
  
  for (sm in levels(SOM_Label)){
    plot(NA,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',bty='n')
    text(0.5,0.5,sm)
    for (mr in whichColumn){
      plot(NA,xlim=c(0.5,2.5),ylim=c(-1,1),xaxt='n',yaxt='n',bty='n')
      rect(0.5,-1,2.5,1,col = bg_matrix[sm,mr],border = NA)
      BOy<-BO_mean[BO_mean$SOM==sm,paste0('cor.',mr)]
      BOv<-BO_sd[BO_sd$SOM==sm,paste0('cor.',mr)]
      BOt<-BO_tt[BO_tt$SOM==sm,paste0('cor.',mr)]
      
      RUy<-RU_mean[RU_mean$SOM==sm,paste0('cor.',mr)]
      RUv<-RU_sd[RU_sd$SOM==sm,paste0('cor.',mr)]
      RUt<-RU_tt[RU_tt$SOM==sm,paste0('cor.',mr)]
      
      clr<-ifelse(length(BOy)==0,'gray',ifelse(BOy>0,'red','blue'))
      clr<-ifelse(BOt %in% c(' ',''),colorspace::lighten(clr,0.75),clr)
      rect(1-0.4,-0.5,1+0.4,0.5,col = clr,border=NA)
      clr<-ifelse(length(RUy)==0,'gray',ifelse(RUy>0,'red','blue'))
      clr<-ifelse(RUt %in% c(' ',''),colorspace::lighten(clr,0.75),clr)
      rect(2-0.4,-0.5,2+0.4,0.5,col = clr,border=NA)
      abline(v=c(0.5,2.5))
      
      if (length(BOy)!=0){
        clr<-ifelse(BOy==0,'gray',ifelse(BOy>0,'red','blue'))
        clr<-'white'
        text(1,0,BO_tt[BO_tt$SOM==sm,paste0('cor.',mr)],col=clr,font=2)
      } else {
        text(1,0,'nd')
      }
      
      if(length(RUy)!=0){
        clr<-ifelse(BOy==0,'gray',ifelse(BOy>0,'red','blue'))
        clr<-'white'
        text(2,0,RU_tt[RU_tt$SOM==sm,paste0('cor.',mr)],col=clr,font=2)
      } else {
        text(2,0,'nd')
      }
    }
  }
  dev.off()
  out_p.val[['BO']][[uids[u]]]<-BO_tt
  out_p.val[['RU']][[uids[u]]]<-RU_tt
  out_r.val[['BO']][[uids[u]]]<-BO_mean
  out_r.val[['RU']][[uids[u]]]<-RU_mean
}


TEMP_BO<-do.call(rbind.data.frame,out_r.val$BO)

TEMP_RU<-do.call(rbind.data.frame,out_r.val$RU)

out.p<-matrix(NA,
              nrow=length(levels(SOM_Label)),
              ncol=ncol(TEMP_BO)-1,
              dimnames = list(SOM=levels(SOM_Label),
                              marker=whichColumn))
out.m_BO<-matrix(NA,
                 nrow=length(levels(SOM_Label)),
                 ncol=ncol(TEMP_BO)-1,
                 dimnames = list(SOM=levels(SOM_Label),
                                 marker=whichColumn))
out.m_RU<-matrix(NA,
                 nrow=length(levels(SOM_Label)),
                 ncol=ncol(TEMP_BO)-1,
                 dimnames = list(SOM=levels(SOM_Label),
                                 marker=whichColumn))
out.s_BO<-matrix(NA,
                 nrow=length(levels(SOM_Label)),
                 ncol=ncol(TEMP_BO)-1,
                 dimnames = list(SOM=levels(SOM_Label),
                                 marker=whichColumn))
out.s_RU<-matrix(NA,
                 nrow=length(levels(SOM_Label)),
                 ncol=ncol(TEMP_BO)-1,
                 dimnames = list(SOM=levels(SOM_Label),
                                 marker=whichColumn))

for (rr in levels(SOM_Label)){
  for (cc in whichColumn){
    BOB<-TEMP_BO[TEMP_BO$SOM==rr,paste0('cor.',cc)]
    RUB<-TEMP_RU[TEMP_RU$SOM==rr,paste0('cor.',cc)]
    out.p[rr,cc]<-t.test(BOB,RUB)[['p.value']]
    out.m_BO[rr,cc]<-mean(BOB)
    out.m_RU[rr,cc]<-mean(RUB)
    out.s_BO[rr,cc]<-sd(BOB)
    out.s_RU[rr,cc]<-sd(RUB)
  }
}

out.p<-matrix(gtools::stars.pval(p.adjust(out.p,method = 'bonferroni')),
              nrow=length(levels(SOM_Label)),
              ncol=ncol(TEMP_BO)-1,
              dimnames = list(SOM=levels(SOM_Label),
                              marker=whichColumn))

out.m_BO['11','x160gd.cd68.gd160di.']
out.m_RU['11','x160gd.cd68.gd160di.']
lapply(out_p.val[['BO']],function(x)x['11','cor.x160gd.cd68.gd160di.'])
lapply(out_p.val[['RU']],function(x)x['11','cor.x160gd.cd68.gd160di.'])



################# alternative method ##########
################# 
################# 

RU_dist<-lapply(uids[1],function(u){
  out<-exactextractr::exact_extract(x = rst[[u]],
                                    y= geo_list_RU_annot[geo_list_RU_annot$uid==u,],
                                    include_xy=T,
                                    include_cols='SOM')
  out<-lapply(out,function(x){
    xc<-mean(x$x)
    yc<-mean(x$y)
    out<-apply(x,1,function(y){
      dist(matrix(c(xc,yc,y['x'],y['y']),ncol = 2,byrow = T))})
    out<-out/max(out)
    out<-cbind.data.frame(x,dist=out)})
  
  out<-do.call(rbind.data.frame,out)
  return(out)
})
RU_dist<-do.call(rbind.data.frame,RU_dist)

BO_dist<-lapply(uids[1],function(u){
  out<-exactextractr::exact_extract(x = rst[[u]],
                                    y= geo_list_BO_annot[geo_list_BO_annot$uid==u,],
                                    include_xy=T,
                                    include_cols='SOM')
  out<-lapply(out,function(x){
    xc<-mean(x$x)
    yc<-mean(x$y)
    out<-apply(x,1,function(y){
      dist(matrix(c(xc,yc,y['x'],y['y']),ncol = 2,byrow = T))})
    out<-out/max(out)
    out<-cbind.data.frame(x,dist=out)})
  
  out<-do.call(rbind.data.frame,out)
  return(out)
})
BO_dist<-do.call(rbind.data.frame,BO_dist)

List_dist<-list(BO_dist[BO_dist$SOM=='12',],RU_dist[RU_dist$SOM=='12',])

cc<-cocor::cocor(~dist+x160gd.cd68.gd160di.|dist+x160gd.cd68.gd160di.,List_dist)
cor.test(List_dist[[1]]$dist,List_dist[[1]]$x160gd.cd68.gd160di.)

RU_dist<-do.call(rbind.data.frame,RU_dist)

RU_cor.P<-lapply(setNames(levels(SOM_Label),levels(SOM_Label)),function(s){
  TEMP<-RU_dist[RU_dist$SOM==s,]
  cr<-lapply(setNames(whichColumn,whichColumn),function(wc){
    out<-try(cor.test(TEMP[,'dist'],trns$transform(TEMP[,wc])))
    if (inherits(out,'try-error')) out<-NA else out<-out[['p.value']]
    return (out)
  })
  do.call(cbind,cr)
})
RU_cor.P<-do.call(rbind,RU_cor.P)
rownames(RU_cor.P)<-levels(SOM_Label)
RU_cor<-matrix(p.adjust(RU_cor.P,'bonferroni'),
               nrow = nrow(RU_cor.P),
               ncol=ncol(RU_cor.P),
               dimnames = dimnames(RU_cor.P))
RU_cor.P<-gtools::stars.pval(RU_cor.P)

RU_cor.C<-lapply(setNames(levels(SOM_Label),levels(SOM_Label)),function(s){
  TEMP<-RU_dist[RU_dist$SOM==s,]
  cr<-lapply(setNames(whichColumn,whichColumn),function(wc){
    out<-try(cor.test(TEMP[,'dist'],trns$transform(TEMP[,wc])))
    if (inherits(out,'try-error')) out<-NA else out<-out[['estimate']]
    return (out)
  })
  do.call(cbind,cr)
})
RU_cor.C<-do.call(rbind,RU_cor.C)
rownames(RU_cor.C)<-levels(SOM_Label)


################# other option#####
################# 

RU_dist<-do.call(rbind.data.frame,lapply(uids,function(u){
  out<-exactextractr::exact_extract(x = rst[[u]],
                                    y= geo_list_RU_annot[geo_list_RU_annot$uid==u,],
                                    include_xy=T,
                                    include_cols='SOM')
  out<-lapply(out,function(x){
    xc<-mean(x$x)
    yc<-mean(x$y)
    out<-apply(x,1,function(y){
      dist(matrix(c(xc,yc,y['x'],y['y']),ncol = 2,byrow = T))})
    out<-out/max(out)
    out<-cbind.data.frame(x,dist=out)})
  
  out<-do.call(rbind.data.frame,out)
  return(out)
}))


BO_dist<-do.call(rbind.data.frame,lapply(uids,function(u){
  out<-exactextractr::exact_extract(x = rst[[u]],
                                    y= geo_list_BO_annot[geo_list_BO_annot$uid==u,],
                                    include_xy=T,
                                    include_cols='SOM')
  out<-lapply(out,function(x){
    xc<-mean(x$x)
    yc<-mean(x$y)
    out<-apply(x,1,function(y){
      dist(matrix(c(xc,yc,y['x'],y['y']),ncol = 2,byrow = T))})
    out<-out/max(out)
    out<-cbind.data.frame(x,dist=out)})
  
  out<-do.call(rbind.data.frame,out)
  return(out)
}))


### skl####

numCores <- parallel::detectCores()
cl<-parallel::makeCluster(numCores)
parallel::clusterExport(cl = cl,
                        varlist = c('rst',
                                    'geo_list_RU_annot',
                                    'skl'))

RU_dist<-do.call(rbind.data.frame,parallel::parLapply(uids,function(u){
  out<-exactextractr::exact_extract(x = rst[[u]],
                                    y= geo_list_RU_annot[geo_list_RU_annot$uid==u,],
                                    include_xy=T,
                                    include_cols='SOM')
  out<-lapply(out,function(x){
    out<-skl(x[,c('x','y')])
    out<-cbind.data.frame(x,dist=out)})
  out<-do.call(rbind.data.frame,out)
  return(out)
},cl=cl))

parallel::stopCluster(cl)

cl<-parallel::makeCluster(numCores)
parallel::clusterExport(cl = cl,
                        varlist = c('rst',
                                    'geo_list_BO_annot',
                                    'skl'))

BO_dist<-do.call(rbind.data.frame,parallel::parLapply(uids,function(u){
  out<-exactextractr::exact_extract(x = rst[[u]],
                                    y= geo_list_BO_annot[geo_list_BO_annot$uid==u,],
                                    include_xy=T,
                                    include_cols='SOM')
  out<-lapply(out,function(x){
    out<-skl(x[,c('x','y')])
    out<-cbind.data.frame(x,dist=out)})
  out<-do.call(rbind.data.frame,out)
  return(out)
},cl=cl))

parallel::stopCluster(cl)

### skl double norm####

numCores <- parallel::detectCores()
cl<-parallel::makeCluster(numCores)
parallel::clusterExport(cl = cl,
                        varlist = c('rst',
                                    'geo_list_RU_annot',
                                    'skl',
                                    'colToUse'))

RU_dist<-do.call(rbind.data.frame,parallel::parLapply(uids,function(u,dt=geo_list_RU_annot,cutoff=0.01){
  dtTable<-exactextractr::exact_extract(x = rst[[u]],
                                        y= dt[dt$uid==u,],
                                        include_xy=T,
                                        include_cols='SOM')
  dtTable<-lapply(1:length(dtTable),function(x){
    cbind.data.frame(dtTable[[x]],poly=x,dist=NA)
  })
  dtTable<-do.call(rbind.data.frame,dtTable)
  dtTable[,colToUse]<-scales::modulus_trans(0)$transform(dtTable[,colToUse])
  browser()
  clean<-apply(dtTable[,colToUse],2,function(x){
    out<-quantile(x,c(cutoff,1-cutoff))
    out<-x>=out[1] & x<=out[2]
    x[!out]<-NA
    return(x)
  })
  dtTable[,colToUse]<-clean
  out_SD<-apply(dtTable[,colToUse],2,sd,na.rm=T)
  out_MEANS<-aggregate(dtTable[,colToUse],list(poly=dtTable[,'poly']),mean,na.rm=T)
  for (i in out_MEANS$poly){
    dtTable[dtTable$poly==i,colToUse]<-scale(dtTable[dtTable$poly==i,colToUse],
                                             center = unlist(out_MEANS[out_MEANS$poly==i,colToUse,drop=T]),
                                             scale = out_SD)
    dist<-skl(dtTable[dtTable$poly==i,c('x','y')])
    # dist<-(dist-min(dist))/(max(dist)-min(dist))
    dtTable[dtTable$poly==i,'dist']<-dist
  }
  return(dtTable)
},cl = cl))

parallel::stopCluster(cl)

cl<-parallel::makeCluster(numCores)
parallel::clusterExport(cl = cl,
                        varlist = c('rst',
                                    'geo_list_BO_annot',
                                    'skl',
                                    'colToUse'))

BO_dist<-do.call(rbind.data.frame,parallel::parLapply(uids,function(u,dt=geo_list_BO_annot,cutoff=0.01){
  dtTable<-exactextractr::exact_extract(x = rst[[u]],
                                        y= dt[dt$uid==u,],
                                        include_xy=T,
                                        include_cols='SOM')
  dtTable<-lapply(1:length(dtTable),function(x){
    cbind.data.frame(dtTable[[x]],poly=x,dist=NA)
  })
  dtTable<-do.call(rbind.data.frame,dtTable)
  dtTable[,colToUse]<-scales::modulus_trans(0)$transform(dtTable[,colToUse])
  clean<-apply(dtTable[,colToUse],2,function(x){
    out<-quantile(x,c(cutoff,1-cutoff))
    out<-x>=out[1] & x<=out[2]
    x[!out]<-NA
    return(x)
  })
  dtTable[,colToUse]<-clean
  out_SD<-apply(dtTable[,colToUse],2,sd,na.rm=T)
  out_MEANS<-aggregate(dtTable[,colToUse],list(poly=dtTable[,'poly']),mean,na.rm=T)
  for (i in out_MEANS$poly){
    dtTable[dtTable$poly==i,colToUse]<-scale(dtTable[dtTable$poly==i,colToUse],
                                             center = unlist(out_MEANS[out_MEANS$poly==i,colToUse,drop=T]),
                                             scale = out_SD)
    dist<-skl(dtTable[dtTable$poly==i,c('x','y')])
    # dist<-(dist-min(dist))/(max(dist)-min(dist))
    dtTable[dtTable$poly==i,'dist']<-dist
  }
  return(dtTable)
},cl = cl))

parallel::stopCluster(cl)







lm_out<-lapply(setNames(levels(SOM_Label),levels(SOM_Label)),function(sl){
  RU_dist_sub<-RU_dist[RU_dist$SOM==sl,]
  RU_dist_sub[,colToUse]<-trns$transform(RU_dist_sub[,colToUse])
  RU_dist_sub[,colToUse]<-scale(RU_dist_sub[,colToUse])
  ex<-expression(paste0('dist~',paste(colToUse,collapse='+')))
  TEST<-lm(eval(ex),RU_dist_sub)
  RU_out1<-summary(TEST)
  RU_out2<-car::Anova(TEST)
  RU_out3<-anova(TEST)
  
  
  BO_dist_sub<-BO_dist[BO_dist$SOM==sl,]
  colnames(BO_dist)
  BO_dist_sub[,colToUse]<-trns$transform(BO_dist_sub[,colToUse])
  BO_dist_sub[,colToUse]<-scale(BO_dist_sub[,colToUse])
  ex<-expression(paste0('dist~',paste(colToUse,collapse='+')))
  TEST<-lm(eval(ex),BO_dist_sub)
  BO_out1<-summary(TEST)
  BO_out2<-car::Anova(TEST)
  BO_out3<-anova(TEST)
  
  out<-list(BO=list(LM=BO_out1,
                    AN1=BO_out2,
                    AN2=BO_out3),
            RU=list(LM=RU_out1,
                    AN1=RU_out2,
                    AN2=RU_out3))
  return(out)
})

TEMP<-lapply(lm_out,function(x){
  matrix(c(gtools::stars.pval(x$BO$AN1$`Pr(>F)`),
           gtools::stars.pval(x$RU$AN1$`Pr(>F)`)),
         ncol = 2, byrow = F)
})

TEMP<-lapply(lm_out,function(x){
  matrix(c(x$BO$LM$coefficients[,1]),
         ncol = 1, byrow = F)
})

TEMP<-lapply(lm_out,function(x){
  matrix(c(x$BO$LM$adj.r.squared,
           x$RU$LM$adj.r.squared),
         ncol = 2, byrow = T)
})

boot.function<-function(dt,id,ctu,trns){
  dt_sub<-dt[id,]
  # dt_sub[,ctu]<-trns$transform(dt_sub[,ctu])
  # dt_sub[,ctu]<-scale(dt_sub[,ctu],center = T,scale = T)
  ex<-expression(paste0('dist~',paste(ctu,collapse='+')))
  TEST.lm<-lm(eval(ex),dt_sub)
  out<-c(coef(TEST.lm),summary(TEST.lm)$adj.r.squared)
  return(out)
}


lm_out<-lapply(setNames(levels(SOM_Label),levels(SOM_Label)),function(sl){
  BO_out<-boot::boot(
    data = BO_dist[BO_dist$SOM==sl & BO_dist$dist!=0,],
    statistic = boot.function,
    R = 100,
    parallel ='snow',
    ctu=colToUse,
    trns = scales::modulus_trans(0),
    ncpus = 24)
  
  RU_out<-boot::boot(
    data = RU_dist[RU_dist$SOM==sl & RU_dist$dist!=0,],
    statistic = boot.function,
    R = 100,
    parallel = 'snow',
    ctu=colToUse,
    trns = scales::modulus_trans(0),
    ncpus = 24)
  
  return(list(BO=BO_out,
              RU=RU_out))})


lapply(lm_out,function(x){c(x$BO$t0[32],x$RU$t0[32])})
compare<-lapply(levels(SOM_Label),function(sl){
  sapply(1:ncol(lm_out$`01`$BO$t),function(x) t.test(lm_out[[sl]]$BO$t[,x],
                                                     lm_out[[sl]]$RU$t[,x])$p.value)
})

compare<-do.call(cbind,compare)

compare<-matrix(p.adjust(compare,'bonferroni'),ncol=ncol(compare),nrow=nrow(compare))


compare<-lapply(levels(SOM_Label),function(sl){
  bo<-as.vector(lm_out[[sl]]$BO$t)
  bo<-data.frame(val = bo,
                 param = unlist(lapply(c('q',colToUse,'rsq'),function(x)rep(x,nrow(lm_out[[sl]]$BO$t)))),
                 meth = rep('bo',length(bo)))
  
  ru<-as.vector(lm_out[[sl]]$RU$t)
  ru<-data.frame(val = ru,
                 param = unlist(lapply(c('q',colToUse,'rsq'),function(x)rep(x,nrow(lm_out[[sl]]$RU$t)))),
                 meth = rep('ru',length(ru)))
  tot<-rbind.data.frame(bo,ru)
  l<-lm(val~param*meth,tot)
  a<-aov(l)
  tk<-TukeyHSD(a,conf.level = 0.999)
})

compare<-do.call(rbind.data.frame,lapply(levels(SOM_Label),function(sl){
  bo<-as.vector(lm_out[[sl]]$BO$t[,32])
  bo<-data.frame(val = bo,
                 som = sl,
                 meth = rep('bo',length(bo)))
  
  ru<-as.vector(lm_out[[sl]]$RU$t[,32])
  ru<-data.frame(val = ru,
                 som = sl,
                 meth = rep('ru',length(ru)))
  tot<-rbind.data.frame(bo,ru)}))

l<-lm(val~som*meth,compare)
a<-aov(l)
tk<-TukeyHSD(a,conf.level = 0.999999)

BO_dist_sub<-BO_dist[BO_dist$SOM=='01',]
mh<-aggregate(x193ir.dna2.ir193di.~dist,BO_dist_sub,mean)
ms<-aggregate(x193ir.dna2.ir193di.~dist,BO_dist_sub,sd)
ns<-table(BO_dist_sub$dist)
ms$x193ir.dna2.ir193di.<-ms$x193ir.dna2.ir193di./sqrt(ns)
plot(mh)
for (i in 1:nrow(ms)){
  lines(matrix(c(ms[i,1],ms[i,1],mh[i,2]+ms[i,2],mh[i,2]-ms[i,2]),
               ncol = 2,
               byrow = F))
}

RU_dist_sub<-RU_dist[RU_dist$SOM=='01',]
ll<-lm(x191ir.dna1.ir191di.~dist,RU_dist_sub)
mh<-aggregate(x191ir.dna1.ir191di.~dist,RU_dist_sub,mean)
ms<-aggregate(x191ir.dna1.ir191di.~dist,RU_dist_sub,sd)
ns<-table(RU_dist_sub$dist)
# ms$x191ir.dna1.ir191di.<-ms$x191ir.dna1.ir191di./sqrt(ns)
plot(mh,xlim=c(0,1),ylim=c(min(RU_dist_sub$x191ir.dna1.ir191di.,na.rm = T),
                               max(RU_dist_sub$x191ir.dna1.ir191di.,na.rm=T)))
for (i in 1:nrow(ms)){
  lines(matrix(c(ms[i,1],ms[i,1],mh[i,2]+ms[i,2],mh[i,2]-ms[i,2]),
               ncol = 2,
               byrow = F))
}
abline(ll)

BO_dist_sub<-BO_dist[BO_dist$SOM=='01',]
ll<-lm(x193ir.dna2.ir193di.~dist,BO_dist_sub)
mh<-aggregate(x193ir.dna2.ir193di.~dist,BO_dist_sub,mean)
ms<-aggregate(x193ir.dna2.ir193di.~dist,BO_dist_sub,sd)
ns<-table(BO_dist_sub$dist)
ms$x193ir.dna2.ir193di.<-ms$x193ir.dna2.ir193di./sqrt(ns)
plot(mh,xlim=c(0,1),ylim=c(0,1),asp=1)
for (i in 1:nrow(ms)){
  lines(matrix(c(ms[i,1],ms[i,1],mh[i,2]+ms[i,2],mh[i,2]-ms[i,2]),
               ncol = 2,
               byrow = F))
}
abline(ll)




TEMP_mean<-aggregate(x191ir.dna1.ir191di.~dist,BO_dist[BO_dist$SOM=='01',],mean)
TEMP_sd<-aggregate(x191ir.dna1.ir191di.~dist,BO_dist[BO_dist$SOM=='01',],sd)
plot(TEMP_mean$dist,TEMP_mean$x191ir.dna1.ir191di.)

lm_out<-lapply(setNames(levels(SOM_Label),levels(SOM_Label)),function(sl){
  TEST1<-boot::boot(RU_dist[RU_dist$SOM==sl,],function(data,indices){
    RU_dist_sub<-data[indices,]
    RU_dist_sub[,colToUse]<-trns$transform(RU_dist_sub[,colToUse])
    RU_dist_sub[,colToUse]<-scale(RU_dist_sub[,colToUse])
    ex<-expression(paste0('dist~',paste(colToUse,collapse='+')))
    TEST<-lm(eval(ex),RU_dist_sub)
    RU_out1<-summary(TEST)
    return(RU_out1$adj.r.squared)
  },R = 10)
  TEST2<-boot::boot(BO_dist[BO_dist$SOM==sl,],function(data,indices){
    RU_dist_sub<-data[indices,]
    RU_dist_sub[,colToUse]<-trns$transform(RU_dist_sub[,colToUse])
    RU_dist_sub[,colToUse]<-scale(RU_dist_sub[,colToUse])
    ex<-expression(paste0('dist~',paste(colToUse,collapse='+')))
    TEST<-lm(eval(ex),RU_dist_sub)
    RU_out1<-summary(TEST)
    return(RU_out1$adj.r.squared)
  },R = 10)
  
  p.val<-t.test(TEST1$t,TEST2$t)['p.value']
  return(list(BO=TEST2,RU=TEST1,p.val=p.val))
})

p.adjust(unlist(lapply(lm_out,function(x) x$p.val)),'bonferroni')
order(unlist(lapply(lm_out,function(x){
  mean(x$BO$t)-mean(x$RU$t)
})))

unlist(lapply(lm_out,function(x){
  mean(x$BO$t)-mean(x$RU$t)
}))[order(unlist(lapply(lm_out,function(x){
  mean(x$BO$t)-mean(x$RU$t)
})))]

p.adjust(unlist(lapply(lm_out,function(x) x$p.val)),'bonferroni')[order(unlist(lapply(lm_out,function(x){
  mean(x$BO$t)-mean(x$RU$t)
})))]
