library(barbieHistologist)

rootFolder<-"C:/Users/luigi/Documents/BH_TEST/BH_TEST2.1/"
targetFolder<-'RESULTS_FAMILY'
brakeName<-'complex_'
rpttn<-"REP_1"
trainingSampl<-c(3)
brks<-c(1:5)

db<-factor(c('C1','C2','C3','C4','C5'),levels=c('C1','C2','C3','C4','C5'))

geo_list_BO<-sf::st_read(
  file.path(rootFolder,rpttn,"TEST_ROUND1/analysis/BODEN1/ex.sqlite")
)
if (any(!sf::st_is_valid(geo_list_BO))){
  geo_list_BO<-sf::st_make_valid(geo_list_BO)
}

geo_list_BO<-geo_list_BO[,c('uid','id','GEOMETRY')]

geo_list_RU<-sf::st_read(
  file.path(rootFolder,rpttn,"TEST_ROUND1/analysis/RUNIMC1/topLayers.sqlite")
)
if (any(!sf::st_is_valid(geo_list_RU))){
  geo_list_RU<-sf::st_make_valid(geo_list_RU)
}

geo_list_RU<-geo_list_RU[,c('uid','splitp_id','GEOMETRY')]
colnames(geo_list_RU)[2]<-'id'

uids<-unique(geo_list_RU$uid)
entrp<-list()
IOU<-list()
geo_list_GR_sub_total<-list()

for (densityBrake in brks){
  
  familyRaster<-bh_loadFamily(file.path(rootFolder,rpttn,paste0(brakeName,densityBrake,"/FAM_",densityBrake)))
  
  geo_list_GR_sub_total[[densityBrake]]<-RUNIMC::lazyCatMap(familyRaster$ID,fn_indexToExclude = 0,fn_uid = uids[densityBrake])
  colnames(geo_list_GR_sub_total[[densityBrake]])<-c('uid','id','GEOMETRY')
  sf::st_geometry(geo_list_GR_sub_total[[densityBrake]])<-'GEOMETRY'
  geo_list_GR_sub<-geo_list_GR_sub_total[[densityBrake]]
  
  
  #### BODEN ####
  geo_list_BO_sub<-geo_list_BO[geo_list_BO$uid==uids[densityBrake],]
  
  bigTable_BO<-bh_extractFamily(family = familyRaster,
                                sfc = geo_list_BO_sub,
                                sfcPrimaryKey = 'id')
  
  whichID<-unique(bigTable_BO$id)
  whichID<-whichID[whichID!=0]
  IOU[['BO']][[densityBrake]]<-sapply(whichID,function(idX){
    
    trueID<-bigTable_BO$ID[bigTable_BO$id==idX]
    trueID.table<-table(trueID)
    trueID.table<-trueID.table[names(trueID.table)!='0']
    trueID.winner<-as.numeric(names(trueID.table)[which.max(trueID.table)])
    
    predicted.geo<-geo_list_BO_sub$GEOMETRY[geo_list_BO_sub$id==idX]
    tru.geo<-geo_list_GR_sub$GEOMETRY[geo_list_GR_sub$id==trueID.winner]
    
    area.intersect<-sf::st_area(sf::st_intersection(predicted.geo,tru.geo))
    if (length(area.intersect)==0) area.intersect<-0
    area.union<-sf::st_area(sf::st_union(predicted.geo,tru.geo))
    out<-area.intersect/area.union
    return(out)
  },simplify = T,USE.NAMES = F)
  
  whichID<-unique(bigTable_BO$id)
  whichID<-whichID[whichID!=0]
  entrp[['BO']][[densityBrake]]<-sapply(whichID,function(idX){
    trueID<-bigTable_BO$ID[bigTable_BO$id==idX]
    trueID.table<-table(trueID)
    
    trueID.table<-trueID.table/sum(trueID.table)
    nNumber<-length(trueID.table)
    
    if (nNumber==1) out<-0 else {
      out<-(-sum(trueID.table*log(trueID.table)/log(nNumber)))}
    return(out)
  },simplify = T,USE.NAMES = F)
  
  
  #### RUNIMC ####
  geo_list_RU_sub<-geo_list_RU[geo_list_RU$uid==uids[densityBrake],]
  
  bigTable_RU<-bh_extractFamily(family = familyRaster,
                                sfc = geo_list_RU_sub,
                                sfcPrimaryKey = 'id')
  
  whichID<-unique(bigTable_RU$id)
  whichID<-whichID[whichID!=0]
  IOU[['RU']][[densityBrake]]<-sapply(whichID,function(idX){
    trueID<-bigTable_RU$ID[bigTable_RU$id==idX]
    trueID.table<-table(trueID)
    trueID.table<-trueID.table[names(trueID.table)!='0']
    trueID.winner<-as.numeric(names(trueID.table)[which.max(trueID.table)])
    
    predicted.geo<-geo_list_RU_sub$GEOMETRY[geo_list_RU_sub$id==idX]
    tru.geo<-geo_list_GR_sub$GEOMETRY[geo_list_GR_sub$id==trueID.winner]
    
    area.intersect<-sf::st_area(sf::st_intersection(predicted.geo,tru.geo))
    if (length(area.intersect)==0) area.intersect<-0
    area.union<-sf::st_area(sf::st_union(predicted.geo,tru.geo))
    out<-area.intersect/area.union
    return(out)
  },simplify = T,USE.NAMES = F)
  
  
  whichID<-unique(bigTable_RU$id)
  whichID<-whichID[whichID!=0]
  entrp[['RU']][[densityBrake]]<-sapply(whichID,function(idX){
    trueID<-bigTable_RU$ID[bigTable_RU$id==idX]
    trueID.table<-table(trueID)
    
    trueID.table<-trueID.table/sum(trueID.table)
    nNumber<-length(trueID.table)
    
    if (nNumber==1) out<-0 else {
      out<-(-sum(trueID.table*log(trueID.table)/log(nNumber)))}
    return(out)
  },simplify = T,USE.NAMES = F)
  
  
  
}



entrp_long<-lapply(names(entrp),function(i){
  out<-lapply(brks,function(ii){
    
    cbind.data.frame(entrp[[i]][[ii]],data.frame(Method=i,density=db[ii]))
    
  })
  do.call(rbind.data.frame,out)
})

entrp_long<-do.call(rbind.data.frame,entrp_long)

colnames(entrp_long)[1]<-'entropy'

entrp_long$Method<-factor(entrp_long$Method,levels=c('RU','BO'))
dirRes<-file.path(rootFolder,rpttn,targetFolder)
dir.create(file.path(rootFolder,rpttn,targetFolder))

write.csv(entrp_long,
          file = file.path(dirRes,'entropy.csv'))

ggp<-ggplot2::ggplot()+ggplot2::geom_boxplot(mapping =ggplot2::aes(y=entropy,
                                                                   x=density,
                                                                   fill = Method),
                                             position = 'dodge',
                                             data = entrp_long)+
  ggplot2::scale_fill_manual(values = c('gray50','white','gray25'))+
  ggplot2::ylim(c(0,1))+
  ggplot2::theme_classic()

ggplot2::ggsave(plot = ggp,
                filename = file.path(dirRes,"entropy.eps"),
                device = 'eps',
                width = 5,
                height = 5,
                dpi = 600)

IOU_long<-lapply(names(IOU),function(i){
  out<-lapply(brks,function(ii){
    
    cbind.data.frame(IOU[[i]][[ii]],data.frame(Method=i,density=db[ii]))
    
  })
  do.call(rbind.data.frame,out)
})

IOU_long<-do.call(rbind.data.frame,IOU_long)

colnames(IOU_long)[1]<-'IOU'
IOU_long$Method<-factor(IOU_long$Method,levels=c('RU','BO'))

write.csv(IOU_long,
          file = file.path(dirRes,'IOU.csv'))

ggp<-ggplot2::ggplot()+ggplot2::geom_boxplot(mapping =ggplot2::aes(y=IOU,
                                                                   x=density,
                                                                   fill = Method),
                                             position = 'dodge',
                                             data = IOU_long)+
  ggplot2::scale_fill_manual(values = c('gray50','white'))+
  ggplot2::ylim(c(0,1))+
  ggplot2::theme_classic()

ggplot2::ggsave(plot = ggp,
                filename = file.path(dirRes,"IOU.eps"),
                device = 'eps',
                width = 5,
                height = 5,
                dpi = 600)


stepPace=30


for (i in brks){
  
  area_GR<-data.frame(area=sf::st_area(geo_list_GR_sub_total[[i]]),Method='GRT')
  area_RU<-data.frame(area=sf::st_area(geo_list_RU[geo_list_RU$uid==unique(geo_list_RU$uid)[i],]),Method='RU',stringsAsFactors = F)
  area_BO<-data.frame(area=sf::st_area(geo_list_BO[geo_list_BO$uid==unique(geo_list_BO$uid)[i],]),Method='BO',stringsAsFactors = F)
  
  minArea<-0
  maxArea<-max(c(area_GR$area,area_BO$area,area_RU$area))
  
  stepCuts<-seq(minArea,maxArea,by=stepPace)[-1]
  
  area_GR<-area_GR[order(area_GR$area),]
  step_GR<-vapply(stepCuts,function(x){
    test<-area_GR$area<x
    if (!any(test)) return(1)
    max(which(test))
  },c(out=0))
  area_GR$area<-sum(area_GR$area)-cumsum(area_GR$area)
  step_GR<-matrix(c(step_GR,area_GR$area[step_GR]),ncol=2,byrow = F)
  
  
  area_BO<-area_BO[order(area_BO$area),]
  step_BO<-vapply(stepCuts,function(x){
    test<-area_BO$area<x
    if (!any(test)) return(1)
    max(which(test))
  },c(out=0))
  area_BO$area<-sum(area_BO$area)-cumsum(area_BO$area)
  step_BO<-matrix(c(step_BO,area_BO$area[step_BO]),ncol=2,byrow = F)
  
  
  area_RU<-area_RU[order(area_RU$area),]
  step_RU<-vapply(stepCuts,function(x){
    test<-area_RU$area<x
    if (!any(test)) return(1)
    max(which(test))
  },c(out=0))
  area_RU$area<-sum(area_RU$area)-cumsum(area_RU$area)
  step_RU<-matrix(c(step_RU,area_RU$area[step_RU]),ncol=2,byrow = F)
  
  xlm<-max(length(area_GR$area),length(area_BO$area),length(area_RU$area))
  ylm<-max(c(area_GR$area,area_BO$area,area_RU$area))
  
  step_total<-list()
  for (ii in 1:nrow(step_RU)){
    step_total[[ii]]<-rbind(step_RU[ii,],step_BO[ii,],step_GR[ii,])
    # interpol<-apply(step_total[[ii]],1,function(x){
    #   y1<-stats::approx(area_RU$area,
    #                     1:nrow(area_RU),
    #                     xout = x[2],
    #                     yleft = +Inf,
    #                     yright = -Inf)[[2]]
    #   y2<-stats::approx(area_BO$area,
    #                     1:nrow(area_BO),
    #                     xout = x[2],
    #                     yleft = +Inf,
    #                     yright = -Inf)[[2]]
    #   y3<-stats::approx(area_GR$area,
    #                     1:nrow(area_GR),
    #                     xout = x[2],
    #                     yleft = +Inf,
    #                     yright = -Inf)[[2]]
    #   out<-c(y1,y2,y3)
    # })
    # interpol<-vapply(1:ncol(interpol),function(iii){
    #   colVal<-interpol[,iii]
    #   colSort<-sort(colVal)
    #   out<-which(colSort==colVal[iii])
    #   if (length(out)!=1) out<-min(out)
    #   return(out)
    # },c(0))
    # step_total[[ii]]<-step_total[[ii]][interpol,]
    step_total[[ii]]<-step_total[[ii]][order(step_total[[ii]][,1]),]
    step_total[[ii]]<-rbind(
      # matrix(c(0,step_total[[ii]][1,2]),ncol=2),
      step_total[[ii]],
      matrix(c(xlm,step_total[[ii]][nrow(step_total[[ii]]),2]),ncol=2))
  }
  
  
  
  postscript(file = file.path(dirRes,paste0("area_",i,".eps")),
             onefile = T,
             bg = 'white',
             horizontal = F,
             height = 10,
             width = 5)
  par(mar=c(3,3,3,5),xpd=T)
  plot(NA,xlim=c(0,xlm),ylim=c(0,ylm),xaxt='n',yaxt='n',las=1,xlab=NA,ylab=NA)
  xBase<-10^floor(log(xlm,10))
  yBase<-10^floor(log(ylm,10))
  xSpan<-seq(0,xlm,by=xBase)
  ySpan<-seq(0,ylm,by=yBase)
  xTick<-xSpan/xBase
  yTick<-ySpan/yBase
  axis(side = 1,at=xSpan,labels = xTick,las=1)
  axis(side = 2,at=ySpan,labels = yTick,las=1)
  mtext(text = paste0('n of cells (X ',xBase,')'),
        side = 1,
        line = 2)
  mtext(text = paste0('area (X ',yBase,')'),
        side = 2,
        line = 2)
  colLin<-gray(seq(0.25,0.75,length.out=length(step_total)))
  for (ii in 1:length(step_total)){
    lines(step_total[[ii]],col=colLin[ii],lwd=8)
  }
  lines(area_GR$area)
  lines(area_BO$area)
  lines(area_RU$area)
  nbr<-floor(log10(length(area_GR$area)))-1
  points(step_GR,pch=0,cex=2)
  points(step_BO,pch=1,cex=2)
  points(step_RU,pch=2,cex=2)
  axisCuts<-vapply(step_total, function(x) {x[nrow(x),2]}, FUN.VALUE = c(y=0))
  whichCuts<-rle(axisCuts)
  positions<-cumsum(whichCuts$lengths)
  labelsPositions<-whichCuts$values
  labelLabels<-vapply(seq_along(labelsPositions),
                      function(lp){
                        if (whichCuts$lengths[lp]==1) return(as.character(stepCuts[lp]))
                        paste0(stepCuts[positions[lp]-whichCuts$lengths[lp]+1],
                               ' % ',
                               stepCuts[positions[lp]])
                      },
                      FUN.VALUE = c(''))
  axis(side = 4, at = labelsPositions,labels = labelLabels,las=1,xlab='area per cell')
  mtext(text = 'area X cell',
        side = 4,
        line = 3)
  
  dev.off()
  
}

postscript(file = file.path(dirRes,paste0("area_legend.eps")),
           onefile = T,
           bg = 'white',
           horizontal = F,
           height = 10,
           width = 5)

plot.new()
legend("center", legend = c("GrounTruth", "Reference", "alternative"), pch = c(0, 1, 2,NA))
dev.off()
