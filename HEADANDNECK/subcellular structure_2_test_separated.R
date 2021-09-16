targetFolder<-"C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/SUBCELLULAR_2"
source('C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/Scripts/skinlayer.R')

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

# geo_list_total<-as.matrix(rbind(geo_list_BO_trans,geo_list_RU_trans))

geo_list_total<-as.matrix(geo_list_RU_trans[geo_list_RU_trans$x158gd.ecad.gd158di.<0.5,])

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
                        # distf = 1
                        importance=importanceX)
fsom<-FlowSOM::BuildMST(fsom)

# saveRDS(fsom,
#         file = file.path(targetFolder,'som.R'  ))

fsomMFI<-FlowSOM::GetClusterMFIs(fsom)

fsomMFI<-scale(fsomMFI)
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

# postscript(file=file.path(targetFolder,'SOM_heatmap.eps'),
#            onefile = F,
#            horizontal = F,
#            paper = 'special',
#            width = 8,
#            height = 8,
#            bg = 'white',
#            pointsize = 18
# )
pheatmap::pheatmap(t(fsomMFI[,c('CD45ra','CD27','CD8','CD4','CD45ro','CD20','CD68')]),
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


numCores <- parallel::detectCores()
cl<-parallel::makeCluster(numCores)
parallel::clusterExport(cl = cl,
                        varlist = c('rst',
                                    'geo_list_RU_annot',
                                    'skl',
                                    'colToUse'))

RU_dist<-do.call(rbind.data.frame,parallel::parLapply(uids,function(u,dt=geo_list_RU_annot){
  dtTable<-exactextractr::exact_extract(x = rst[[u]],
                                        y= dt[dt$uid==u,],
                                        include_xy=T,
                                        include_cols='SOM')
  dtTable<-lapply(1:length(dtTable),function(x){
    cbind.data.frame(dtTable[[x]],poly=x,dist=NA)
  })
  dtTable<-do.call(rbind.data.frame,dtTable)
  dtTable[,colToUse]<-scales::modulus_trans(0)$transform(dtTable[,colToUse])
  
  out_SD<-apply(dtTable[,colToUse],2,sd,na.rm=T)
  out_MEANS<-aggregate(dtTable[,colToUse],list(poly=dtTable[,'poly']),mean,na.rm=T)
  
  dtTable<-do.call(rbind.data.frame,lapply(out_MEANS$poly,function(i){
    newDtTable<-dtTable[dtTable$poly==i,]
    newDtTable[,colToUse]<-scale(newDtTable[,colToUse],
                                 center = unlist(out_MEANS[out_MEANS$poly==i,colToUse,drop=T]),
                                 scale = out_SD)
    dist<-skl(newDtTable[,c('x','y')])
    newDtTable[,'dist']<-dist
    return(newDtTable)
  }))
  
  return(dtTable)
},cl=cl))

parallel::stopCluster(cl)

cl<-parallel::makeCluster(numCores)
parallel::clusterExport(cl = cl,
                        varlist = c('rst',
                                    'geo_list_BO_annot',
                                    'skl',
                                    'colToUse'))

BO_dist<-do.call(rbind.data.frame,parallel::parLapply(uids,function(u,dt=geo_list_BO_annot){
  dtTable<-exactextractr::exact_extract(x = rst[[u]],
                                        y= dt[dt$uid==u,],
                                        include_xy=T,
                                        include_cols='SOM')
  dtTable<-lapply(1:length(dtTable),function(x){
    cbind.data.frame(dtTable[[x]],poly=x,dist=NA)
  })
  dtTable<-do.call(rbind.data.frame,dtTable)
  dtTable[,colToUse]<-scales::modulus_trans(0)$transform(dtTable[,colToUse])
  
  out_SD<-apply(dtTable[,colToUse],2,sd,na.rm=T)
  out_MEANS<-aggregate(dtTable[,colToUse],list(poly=dtTable[,'poly']),mean,na.rm=T)
  
  dtTable<-do.call(rbind.data.frame,lapply(out_MEANS$poly,function(i){
    newDtTable<-dtTable[dtTable$poly==i,]
    newDtTable[,colToUse]<-scale(newDtTable[,colToUse],
                                 center = unlist(out_MEANS[out_MEANS$poly==i,colToUse,drop=T]),
                                 scale = out_SD)
    dist<-skl(newDtTable[,c('x','y')])
    newDtTable[,'dist']<-dist
    return(newDtTable)
  }))
  
  return(dtTable)
},cl=cl))

parallel::stopCluster(cl)


postscript(file=file.path(targetFolder,'hist_skin.eps'),
           onefile = F,
           horizontal = F,
           paper = 'special',
           width = 16,
           height = 8,
           bg = 'white',
           pointsize = 18
)

par(mfrow=c(1,2))

hist(BO_dist$dist,main = 'Reference',xlab = 'layer index',breaks = sort(unique(BO_dist$dist)))
hist(RU_dist$dist,main = 'Alternative',xlab = 'layer index',breaks = sort(unique(RU_dist$dist)))

dev.off()




xxlim<-0:7
size<-100
boo<-1:100

cl<-parallel::makeCluster(numCores)
parallel::clusterExport(cl = cl,
                        varlist = c('SOM_Label',
                                    'RU_dist',
                                    'xxlim',
                                    'size',
                                    'colToUse'))

lm_RU<-parallel::parLapply(boo,function(bo){
  lm_RU<-lapply(setNames(levels(SOM_Label),levels(SOM_Label)),function(sl,
                                                                       dtTable=RU_dist,
                                                                       xxxlim=xxlim,
                                                                       xxxsize=size,
                                                                       xxxcol=colToUse){
    
    dtTable<-dtTable[dtTable$SOM==sl,]
    actualX<-sort(unique(dtTable$dist))
    xxx<-xxxlim[xxxlim %in% actualX]
    sbst<-unlist(lapply(xxx,function(xl){
      
      out<-which(dtTable$dist==xl)
      out<-out[sample(1:length(out),xxxsize,replace = T)]
      return(out)
    }))
    dtTable<-dtTable[sbst,]
    ex<-expression(paste0('dist~',paste(xxxcol,collapse='+')))
    out<-lm(eval(ex),dtTable)
    out1<-caret::varImp(out)
    return(list(lm=out,imp=out1))
  })
  return(lm_RU)
},cl=cl)

parallel::stopCluster(cl)

# RU_R<-cbind(
#   apply(do.call(cbind,lapply(lm_RU,function(x)unlist(lapply(x,function(y)summary(y$lm)$adj.r.squared)))),1,mean),
#   apply(do.call(cbind,lapply(lm_RU,function(x)unlist(lapply(x,function(y)summary(y$lm)$adj.r.squared)))),1,sd)
# )

cl<-parallel::makeCluster(numCores)
parallel::clusterExport(cl = cl,
                        varlist = c('SOM_Label',
                                    'BO_dist',
                                    'xxlim',
                                    'size',
                                    'colToUse'))

lm_BO<-parallel::parLapply(boo,function(bo){
  lm_par<-lapply(setNames(levels(SOM_Label),levels(SOM_Label)),function(sl,
                                                                        dtTable=BO_dist,
                                                                        xxxlim=xxlim,
                                                                        xxxsize=size,
                                                                        xxxcol=colToUse){
    
    dtTable<-dtTable[dtTable$SOM==sl,]
    actualX<-sort(unique(dtTable$dist))
    xxx<-xxxlim[xxxlim %in% actualX]
    sbst<-unlist(lapply(xxx,function(xl){
      
      out<-which(dtTable$dist==xl)
      out<-out[sample(1:length(out),xxxsize,replace = T)]
      return(out)
    }))
    dtTable<-dtTable[sbst,]
    ex<-expression(paste0('dist~',paste(xxxcol,collapse='+')))
    out<-lm(eval(ex),dtTable)
    out1<-caret::varImp(out)
    return(list(lm=out,imp=out1))
  })
  return(lm_par)
},cl=cl)

parallel::stopCluster(cl)

# BO_R<-
# cbind(
#   apply(do.call(cbind,lapply(lm_BO,function(x)unlist(lapply(x,function(y)summary(y)$adj.r.squared)))),1,mean),
#   apply(do.call(cbind,lapply(lm_BO,function(x)unlist(lapply(x,function(y)summary(y)$adj.r.squared)))),1,sd)
# )
# 
# RU_R<-cbind(
#   apply(do.call(cbind,lapply(lm_RU,function(x)unlist(lapply(x,function(y)summary(y)$adj.r.squared)))),1,mean),
#   apply(do.call(cbind,lapply(lm_RU,function(x)unlist(lapply(x,function(y)summary(y)$adj.r.squared)))),1,sd)
# )

RU_R<-do.call(cbind,lapply(lm_RU,function(x)unlist(lapply(x,function(y)summary(y$lm)$adj.r.squared),recursive = F)))
BO_R<-do.call(cbind,lapply(lm_BO,function(x)unlist(lapply(x,function(y)summary(y$lm)$adj.r.squared),recursive = F)))

postscript(file=file.path(targetFolder,'R-squared.eps'),
           onefile = F,
           horizontal = F,
           paper = 'special',
           width = 8,
           height = 8,
           bg = 'white',
           pointsize = 18
)
par (mfrow=c(dimx,dimy),mar=c(0,2,1,0))
for (i in levels(SOM_Label)){
  tt<-t.test(BO_R[i,],RU_R[i,])$p.value
  tt<-p.adjust(tt,'bonferroni',n = dimx*dimy)
  tt<-gtools::stars.pval(tt)
  plot(NA,xlim=c(0.5,2.5),ylim=c(0,0.5),xaxt='n',main=paste(i,tt,sep=' '))
  points(x=1,mean(BO_R[i,]),pch=16)
  lines(matrix(c(1,mean(BO_R[i,])-sd(BO_R[i,]),
                 1,mean(BO_R[i,])+sd(BO_R[i,])),
               ncol = 2,byrow = T))
  points(x=2,mean(RU_R[i,]),pch=17)
  lines(matrix(c(2,mean(RU_R[i,])-sd(RU_R[i,]),
                 2,mean(RU_R[i,])+sd(RU_R[i,])),
               ncol = 2,byrow = T))
}

dev.off()


postscript(file=file.path(targetFolder,'R-squared_alternative.eps'),
           onefile = F,
           horizontal = F,
           paper = 'special',
           width = 15,
           height = 4,
           bg = 'white',
           pointsize = 18
)
par (mfrow=c(3,9),mar=c(1,2,1,1))
for (i in levels(SOM_Label)){
  tt<-t.test(BO_R[i,],RU_R[i,])$p.value
  tt<-p.adjust(tt,'bonferroni',n = dimx*dimy)
  tt<-gtools::stars.pval(tt)
  plot(NA,xlim=c(0.5,2.5),ylim=c(0,0.35),xaxt='n',main=paste(i,tt,sep=' '))
  points(x=1,mean(BO_R[i,]),pch=1)
  lines(matrix(c(0.75,mean(BO_R[i,])-sd(BO_R[i,]),
                 0.75,mean(BO_R[i,])+sd(BO_R[i,])),
               ncol = 2,byrow = T))
  lines(matrix(c(0.75-0.05,mean(BO_R[i,])+sd(BO_R[i,]),
                 0.75+0.05,mean(BO_R[i,])+sd(BO_R[i,])),
               ncol = 2,byrow = T))
  lines(matrix(c(0.75-0.05,mean(BO_R[i,])-sd(BO_R[i,]),
                 0.75+0.05,mean(BO_R[i,])-sd(BO_R[i,])),
               ncol = 2,byrow = T))
  
  points(x=2,mean(RU_R[i,]),pch=2)
  lines(matrix(c(1.75,mean(RU_R[i,])-sd(RU_R[i,]),
                 1.75,mean(RU_R[i,])+sd(RU_R[i,])),
               ncol = 2,byrow = T))
  lines(matrix(c(1.75-0.05,mean(RU_R[i,])-sd(RU_R[i,]),
                 1.75+0.05,mean(RU_R[i,])-sd(RU_R[i,])),
               ncol = 2,byrow = T))
  lines(matrix(c(1.75-0.05,mean(RU_R[i,])+sd(RU_R[i,]),
                 1.75+0.05,mean(RU_R[i,])+sd(RU_R[i,])),
               ncol = 2,byrow = T))
}

dev.off()









RU_P<-lapply(lm_RU,function(x)lapply(x,function(y)summary(y$lm)$coefficients[,'Estimate']))
BO_P<-lapply(lm_BO,function(x)lapply(x,function(y)summary(y$lm)$coefficients[,'Estimate']))

RU_P<-lapply(setNames(levels(SOM_Label),levels(SOM_Label)),function(sl){
  do.call(cbind,
          lapply(RU_P,function(x){
            x[[sl]]
          }))
})

BO_P<-lapply(setNames(levels(SOM_Label),levels(SOM_Label)),function(sl){
  do.call(cbind,
          lapply(BO_P,function(x){
            x[[sl]]
          }))
})

postscript(file=file.path(targetFolder,'parameters.eps'),
           onefile = F,
           horizontal = F,
           paper = 'special',
           width = 16,
           height = 8,
           bg = 'white',
           pointsize = 18
)

par (mfrow=c(dimx*dimy+1,length(colToUse)+1),mar=c(0,0,0,0))
plot(NA,xlim=c(0.5,2.5),ylim=c(-1,1),xaxt='n',yaxt='n',bty='n')
for (mrkr in colToUse){
  plot(NA,xlim=c(0.5,2.5),ylim=c(-1,1),xaxt='n',yaxt='n')
  text(1.5,0,strsplit(mrkr,'.',fixed=T)[[1]][2])
}

for (i in levels(SOM_Label)){
  plot(NA,xlim=c(0.5,2.5),ylim=c(-1,1),xaxt='n',yaxt='n',)
  text(1.5,0,i)
  for (mrkr in colToUse){
    
    tt_BO<-t.test(BO_P[[i]][mrkr,],mu=0)$p.value
    tt_BO<-p.adjust(tt_BO,'bonferroni',n = dimx*dimy*length(colToUse))
    tt_RU<-t.test(RU_P[[i]][mrkr,],mu=0)$p.value
    tt_RU<-p.adjust(tt_RU,'bonferroni',n = dimx*dimy*length(colToUse))
    mean_BO<-mean(BO_P[[i]][mrkr,])
    upper_BO<-mean(BO_P[[i]][mrkr,])+sd(BO_P[[i]][mrkr,])
    lower_BO<-mean(BO_P[[i]][mrkr,])-sd(BO_P[[i]][mrkr,])
    mean_RU<-mean(RU_P[[i]][mrkr,])
    upper_RU<-mean(RU_P[[i]][mrkr,])+sd(RU_P[[i]][mrkr,])
    lower_RU<-mean(RU_P[[i]][mrkr,])-sd(RU_P[[i]][mrkr,])
    upper_plot<-max(c(upper_BO,upper_RU))*1.2
    lower_plot<-min(c(lower_BO,lower_RU))*1.2
    
    plot(NA,xlim=c(0.5,2.5),ylim=c(lower_plot,upper_plot),xaxt='n',yaxt='n')
    if (tt_BO<0.05){
      if (mean(BO_P[[i]][mrkr,])<0){
      polygon(x=c(0.5,0.5,1.5),
              y = c(lower_plot,upper_plot,lower_plot),
              border = NA,
              col='red')
      } else {
        polygon(x=c(0.5,0.5,1.5),
                y = c(upper_plot,lower_plot,upper_plot),
                border = NA,
                col='red')
      }
      
    }
    if (tt_RU<0.05){
      if (mean(RU_P[[i]][mrkr,])<0){
        polygon(x=c(2.5,2.5,1.5),
                y = c(lower_plot,upper_plot,lower_plot),
                border = NA,
                col='blue')
      } else {
        polygon(x=c(2.5,2.5,1.5),
                y = c(upper_plot,lower_plot,upper_plot),
                border = NA,
                col='blue')
      }
      # rect(0.5,-1,1.5,1,col='red')
    }
    # if (tt_RU<0.05){
    #   rect(1.5,-1,2.5,1,col='blue')
    # }
    abline(h=0,col='gray')
    points(x=1,mean_BO,pch=16)
    lines(matrix(c(1,upper_BO,
                   1,lower_BO),
                 ncol = 2,byrow = T))
    points(x=2,mean_RU,pch=17)
    lines(matrix(c(2,upper_RU,
                   2,lower_RU),
                 ncol = 2,byrow = T))
  }
}

dev.off()

RU_IMP<-lapply(lm_RU,function(x)lapply(x,function(y)y$imp))
BO_IMP<-lapply(lm_BO,function(x)lapply(x,function(y)y$imp))

RU_IMP<-lapply(setNames(levels(SOM_Label),levels(SOM_Label)),function(sl){
  do.call(cbind,
          lapply(RU_IMP,function(x){
            x[[sl]]
          }))
})

BO_IMP<-lapply(setNames(levels(SOM_Label),levels(SOM_Label)),function(sl){
  do.call(cbind,
          lapply(BO_IMP,function(x){
            x[[sl]]
          }))
})

TOT_IMP<-lapply(setNames(colToUse,colToUse),function(mrkr){
  lapply(setNames(levels(SOM_Label),levels(SOM_Label)),function(sl){
    mean_BO<-mean(unlist(BO_IMP[[sl]][mrkr,]))
    mean_RU<-mean(unlist(RU_IMP[[sl]][mrkr,]))
    t.mean<-t.test(BO_IMP[[sl]][mrkr,],RU_IMP[[sl]][mrkr,])
    return(list(BO=mean_BO,RU=mean_RU,t=t.mean$p.value))
})
})

TOT_MAX<-ceiling(max(unlist(TOT_IMP)))
TOT_Color<-colorRampPalette(c('blue','green','yellow','red'))(TOT_MAX*2+1)

postscript(file=file.path(targetFolder,'key_importance.eps'),
           onefile = F,
           horizontal = F,
           paper = 'special',
           width = 5,
           height = 5,
           bg = 'white',
           pointsize = 18
)

plot(NA,xlim=c(0,TOT_MAX*2+2),ylim=c(0,TOT_MAX*2+2),asp=1)
for (i in seq(1,TOT_MAX*2+1,1)){
  rect(1,i,2,i+1,col=TOT_Color[i])
  text(4,i+0.5,format((i-1)/2,nsmall = 1),cex=0.5)
}

dev.off()

postscript(file=file.path(targetFolder,'importance.eps'),
           onefile = F,
           horizontal = F,
           paper = 'special',
           width = 16,
           height = 8,
           bg = 'white',
           pointsize = 18
)

par (mfrow=c(dimx*dimy+1,length(colToUse)+1),mar=c(0,0,0,0))
plot(NA,xlim=c(0.5,2.5),ylim=c(-1,1),xaxt='n',yaxt='n')
for (mrkr in colToUse){
  plot(NA,xlim=c(0.5,2.5),ylim=c(-1,1),xaxt='n',yaxt='n')
  text(1.5,0,strsplit(mrkr,'.',fixed=T)[[1]][2])
}

for (i in levels(SOM_Label)){
  plot(NA,xlim=c(0.5,2.5),ylim=c(-1,1),xaxt='n',yaxt='n',)
  text(1.5,0,i)
  for (mrkr in colToUse){
    
    plot(NA,xlim=c(0,2),ylim=c(0,2),xaxt='n',yaxt='n')
    
    BO_TEMP<-TOT_IMP[[mrkr]][[i]]$BO
    RU_TEMP<-TOT_IMP[[mrkr]][[i]]$RU
    rect(0.2,0.2,1,1.5,col = TOT_Color[ceiling(BO_TEMP*2)+1])
    rect(1,0.2,1.8,1.5,col = TOT_Color[ceiling(RU_TEMP*2)+1])
    text(1,1.8,gtools::stars.pval(TOT_IMP[[mrkr]][[i]]$t),cex=1.5)
  }
}

dev.off()

#### flagplot######
#### 
#### 



postscript(file=file.path(targetFolder,'flagplot.eps'),
           onefile = F,
           horizontal = F,
           paper = 'special',
           width = 16,
           height = 8,
           bg = 'white',
           pointsize = 18
)

par (mfrow=c(dimx*dimy+1,length(colToUse)+1),mar=c(0,0,0,0),bty='n')
plot(NA,xlim=c(0.5,2.5),ylim=c(-1,1),xaxt='n',yaxt='n')
for (mrkr in colToUse){
  plot(NA,xlim=c(0.5,2.5),ylim=c(-1,1),xaxt='n',yaxt='n')
  text(1.5,0,strsplit(mrkr,'.',fixed=T)[[1]][2])
}

for (i in levels(SOM_Label)){
  plot(NA,xlim=c(0.5,2.5),ylim=c(-1,1),xaxt='n',yaxt='n',)
  text(1.5,0,i)
  for (mrkr in colToUse){
    
    tt_BO<-t.test(BO_P[[i]][mrkr,],mu=0)$p.value
    tt_BO<-p.adjust(tt_BO,'bonferroni',n = dimx*dimy*length(colToUse))
    tt_RU<-t.test(RU_P[[i]][mrkr,],mu=0)$p.value
    tt_RU<-p.adjust(tt_RU,'bonferroni',n = dimx*dimy*length(colToUse))
    mean_BO<-mean(BO_P[[i]][mrkr,])
    upper_BO<-mean(BO_P[[i]][mrkr,])+sd(BO_P[[i]][mrkr,])
    lower_BO<-mean(BO_P[[i]][mrkr,])-sd(BO_P[[i]][mrkr,])
    mean_RU<-mean(RU_P[[i]][mrkr,])
    upper_RU<-mean(RU_P[[i]][mrkr,])+sd(RU_P[[i]][mrkr,])
    lower_RU<-mean(RU_P[[i]][mrkr,])-sd(RU_P[[i]][mrkr,])
    upper_plot<-max(c(upper_BO,upper_RU))*1.2
    lower_plot<-min(c(lower_BO,lower_RU))*1.2
    BO_TEMP<-TOT_IMP[[mrkr]][[i]]$BO
    RU_TEMP<-TOT_IMP[[mrkr]][[i]]$RU
    
    
    
    plot(NA,xlim=c(0,10),ylim=c(0,10),xaxt='n',yaxt='n',asp=1,bty='o')
    
    if (tt_BO<0.05){
      if (mean(BO_P[[i]][mrkr,])<0){
        rect(1,1,5,3,col = TOT_Color[ceiling(BO_TEMP*2)+1],border='black')
        polygon(x=c(1,1,5),
                y = c(3,9,3),
                border = 'black',
                col='red')
      } else {
        rect(1,9,5,7,col = TOT_Color[ceiling(BO_TEMP*2)+1],border='black')
        polygon(x=c(1,1,5),
                y = c(7,1,7),
                border = 'black',
                col='red')
      }
      
    }
    if (tt_RU<0.05){
      if (mean(RU_P[[i]][mrkr,])<0){
        rect(9,1,5,3,col = TOT_Color[ceiling(RU_TEMP*2)+1],border='black')
        polygon(x=c(9,9,5),
                y = c(3,9,3),
                border = 'black',
                col='blue')
      } else {
        rect(9,9,5,7,col = TOT_Color[ceiling(RU_TEMP*2)+1],border='black')
        polygon(x=c(9,9,5),
                y = c(7,1,7),
                border = 'black',
                col='blue')
      }
    }
  }
}

dev.off()

postscript(file=file.path(targetFolder,'flagplot_horiz.eps'),
           onefile = F,
           horizontal = F,
           paper = 'special',
           width = 16,
           height = 8,
           bg = 'white',
           pointsize = 18
)

par (mfrow=c(length(colToUse)+1,dimx*dimy+1),mar=c(0,0,0,0),bty='n')
plot(NA,xlim=c(0.5,2.5),ylim=c(-1,1),xaxt='n',yaxt='n')
for (i in levels(SOM_Label)){
  plot(NA,xlim=c(0.5,2.5),ylim=c(-1,1),xaxt='n',yaxt='n')
  text(1.5,0,i)
  
}

for (mrkr in colToUse){
  plot(NA,xlim=c(0.5,2.5),ylim=c(-1,1),xaxt='n',yaxt='n',)
  text(1.5,0,strsplit(mrkr,'.',fixed=T)[[1]][2])
  for (i in levels(SOM_Label)){
    tt_BO<-t.test(BO_P[[i]][mrkr,],mu=0)$p.value
    tt_BO<-p.adjust(tt_BO,'bonferroni',n = dimx*dimy*length(colToUse))
    tt_RU<-t.test(RU_P[[i]][mrkr,],mu=0)$p.value
    tt_RU<-p.adjust(tt_RU,'bonferroni',n = dimx*dimy*length(colToUse))
    mean_BO<-mean(BO_P[[i]][mrkr,])
    upper_BO<-mean(BO_P[[i]][mrkr,])+sd(BO_P[[i]][mrkr,])
    lower_BO<-mean(BO_P[[i]][mrkr,])-sd(BO_P[[i]][mrkr,])
    mean_RU<-mean(RU_P[[i]][mrkr,])
    upper_RU<-mean(RU_P[[i]][mrkr,])+sd(RU_P[[i]][mrkr,])
    lower_RU<-mean(RU_P[[i]][mrkr,])-sd(RU_P[[i]][mrkr,])
    upper_plot<-max(c(upper_BO,upper_RU))*1.2
    lower_plot<-min(c(lower_BO,lower_RU))*1.2
    BO_TEMP<-TOT_IMP[[mrkr]][[i]]$BO
    RU_TEMP<-TOT_IMP[[mrkr]][[i]]$RU
    
    
    
    plot(NA,xlim=c(0,10),ylim=c(0,10),xaxt='n',yaxt='n',asp=1,bty='o')
    
    if (tt_BO<0.05){
      if (mean(BO_P[[i]][mrkr,])<0){
        rect(1,1,5,3,col = TOT_Color[ceiling(BO_TEMP*2)+1],border='black')
        polygon(x=c(1,1,5),
                y = c(3,9,3),
                border = 'black',
                col='red')
      } else {
        rect(1,9,5,7,col = TOT_Color[ceiling(BO_TEMP*2)+1],border='black')
        polygon(x=c(1,1,5),
                y = c(7,1,7),
                border = 'black',
                col='red')
      }
      
    }
    if (tt_RU<0.05){
      if (mean(RU_P[[i]][mrkr,])<0){
        rect(9,1,5,3,col = TOT_Color[ceiling(RU_TEMP*2)+1],border='black')
        polygon(x=c(9,9,5),
                y = c(3,9,3),
                border = 'black',
                col='blue')
      } else {
        rect(9,9,5,7,col = TOT_Color[ceiling(RU_TEMP*2)+1],border='black')
        polygon(x=c(9,9,5),
                y = c(7,1,7),
                border = 'black',
                col='blue')
      }
    }
  }
}

dev.off()
