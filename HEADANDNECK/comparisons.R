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



dimx=8
dimy=8
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
        file = "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001//som"  )

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

pdf(file="C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/heatmap.pdf",
    width = 8,
    height = 5
)
pheatmap::pheatmap(t(fsomMFI),scale = 'none')

dev.off()

SOM_Label<-FlowSOM::GetClusters(fsom)

SOM_Label<-factor(as.character(formatC(SOM_Label,flag='0',format = 'd',digits=1)),levels = rownames(fsomMFI))


geo_list_BO_annot<-dplyr::bind_cols(geo_list_BO,data.frame(SOM=SOM_Label[1:nrow(geo_list_BO)]))


geo_list_RU_annot<-dplyr::bind_cols(geo_list_RU,data.frame(SOM=factor(as.character(SOM_Label[(1+nrow(geo_list_BO)):(nrow(geo_list_BO)+nrow(geo_list_RU))]))))

uids<-unique(geo_list_BO$uid)

for (u in uids){
postscript(file = paste0("C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/SOM_MAPS/",u,"_BO_map.eps"),
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
postscript(file = paste0("C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/SOM_MAPS/",u,"_RU_map.eps"),
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
rst<-RUNIMC:::IMCstackOpen(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/rasterStacks/Sample1_1_ROI003.txt.stk"
)

nms<-names(rst)

RGBmat<-matrix(c(1,0,0,
                 0,1,0,
                 0,0,1),nrow=3,byrow = F,dimnames = list(c('R','G','B'),nms[c(3,16,29)]))

rst2<-lapply(names(rst)[c(3,16,29)],function(x){
  rst<-RUNIMC::quantNorm(rst[[x]],0.95)
  RGBrst<-lapply(RGBmat[,x],function(y){
    rst*y
  })
})

rst3<-lapply(c('R','G','B'),function(clx){
  newrst<-lapply(rst2,'[[',clx)
  newrst<-raster::calc(raster::stack(newrst),sum)
  newrst[newrst>1]<-1
  return(newrst)
})

tiff(
  file = file.path("C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/GR_map.tiff"),
  units = 'in',
  res = 600,
  width = 25,
  height = 20,
  compression = 'lzw'
)
par(mar=c(0,0,0,0))
raster::plotRGB(raster::stack(rst3),scale=1)
dev.off()


prePostTable<-table(sf::st_drop_geometry(geo_list_RU_annot)[,c('primary_id','SOM')])

prePostTable<-apply(prePostTable,2,function(x)x/sum(x))

pdf(file="C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/heatmap_prePost.pdf",
    width = 10,
    height = 3.5
)


pp<-pheatmap::pheatmap(prePostTable,
                   scale = 'none',
                   angle_col = 0)

dev.off()

abb_BO<-table(geo_list_BO_annot$SOM)/sum(table(geo_list_BO_annot$SOM))*100
abb_RU<-table(geo_list_RU_annot$SOM)/sum(table(geo_list_RU_annot$SOM))*100
abb_BO_code<-cut(abb_BO,
                 breaks = c(-Inf,0.01,0.05,0.1,0.5,1,5,+Inf),
                 labels = c('<0.01%','0.05%','0.1%','0.5%','1%','5%','>5%'),
                 include.lowest = F)
abb_RU_code<-cut(abb_RU,
                 breaks = c(-Inf,0.01,0.05,0.1,0.5,1,5,+Inf),
                 labels = c('<0.01%','0.05%','0.1%','0.5%','1%','5%','>5%'),
                 include.lowest = F)
abb<-data.frame(BO=abb_BO_code[pp$tree_col$order],RU=abb_RU_code[pp$tree_col$order])
colnames(abb)<-c('BO','RU')
rownames(abb)<-levels(SOM_Label)[pp$tree_col$order]

newcol<-setNames(colorRampPalette(c('red','green','blue'))(length(levels(abb$BO))),levels(abb$BO))
pheatmap::pheatmap(t(fsomMFI)[,pp$tree_col$order],scale = 'none',cluster_cols = F,
                   annotation_col = abb,
                   annotation_colors = list(BO=newcol,RU=newcol),
                   drop_levels = F)
                   



# uids<-unique(geo_list_BO_annot$uid)
# TEMP_BO<-geo_list_BO_annot[geo_list_BO_annot$uid==uids[2],]
# TEMP_RU<-geo_list_RU_annot[geo_list_RU_annot$uid==uids[2],]
# TEMP_intersection<-sf::st_intersects(TEMP_BO,TEMP_RU)
# 
# TEMP_eval<-lapply(1:length(TEMP_intersection),function(x){
#   if (length(TEMP_intersection[[x]])!=0){
#     INT<-sf::st_area(sf::st_intersection(TEMP_BO[x,'GEOMETRY'],TEMP_RU[TEMP_intersection[[x]],'GEOMETRY']))
#     UNI<-sf::st_area(sf::st_union(TEMP_BO[x,'GEOMETRY'],TEMP_RU[TEMP_intersection[[x]],'GEOMETRY'],by_feature = T))
#     IOU<-INT/UNI
#   }
# })
# 
# TEMP_eval<-unlist(TEMP_eval)

# hist(TEMP_eval)

BO_area<-sf::st_area(geo_list_BO_annot)
BO_area<-cbind.data.frame(area=BO_area,sf::st_drop_geometry(geo_list_BO_annot)[,c('uid','SOM')])
RU_area<-sf::st_area(geo_list_RU_annot)
RU_area<-cbind.data.frame(area=RU_area,sf::st_drop_geometry(geo_list_RU_annot)[,c('uid','SOM')])

TEMP_rst<-RUNIMC:::retrieve.RsCollection(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/rasterStacks"
)
TEMP_bb<-do.call(rbind,lapply(TEMP_rst,function(x){
  uid<-x@uid
  totExt<-x@ncols*x@nrows
  return(data.frame(uid=uid,area=totExt))
}))
BO_tot<-aggregate(area~uid,BO_area,sum)
RU_tot<-aggregate(area~uid,RU_area,sum)
BO_perc<-BO_tot
RU_perc<-RU_tot

for (uids in TEMP_bb$uid){
  BO_perc$area[BO_perc$uid==uids]<-BO_perc$area[BO_perc$uid==uids]/TEMP_bb$area[TEMP_bb$uid==uids]*100
  RU_perc$area[RU_perc$uid==uids]<-RU_perc$area[RU_perc$uid==uids]/TEMP_bb$area[TEMP_bb$uid==uids]*100
}

postscript(file = "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/SOM_MAPS/Total_area.eps",
           onefile = F,
           paper = 'special',
           bg = 'white',
           horizontal = F,
           height = 5,
           width = 5)

plot(NA,xlim=c(0.5,2.5),ylim=c(0,100),xaxt='n',ylab='Area percentage')
axis(1,c(1,2),labels = c('Reference','Alternative'))

for (uids in TEMP_bb$uid){
 points(x=1,
        y=BO_perc$area[BO_perc$uid==uids],
        pch=16)
  points(x=2,
         y=RU_perc$area[RU_perc$uid==uids],
         pch=17)
  lines(matrix(c(1,
                 BO_perc$area[BO_perc$uid==uids],
                 2,
                 RU_perc$area[RU_perc$uid==uids]),
               ncol=2,
               byrow = T))
}

dev.off()


t.test(BO_perc$area,RU_perc$area,paired = T)
wilcox.test(BO_perc$area,RU_perc$area,paired = T)

BO_CLASS<-aggregate(area~SOM+uid,BO_area,sum)
RU_CLASS<-aggregate(area~SOM+uid,RU_area,sum)
BO_CLASS_perc<-BO_CLASS
RU_CLASS_perc<-RU_CLASS

for (uids in TEMP_bb$uid){
  BO_CLASS_perc$area[BO_perc$uid==uids]<-BO_CLASS_perc$area[BO_perc$uid==uids]/sum(BO_CLASS_perc$area[BO_perc$uid==uids])*100
  RU_CLASS_perc$area[RU_perc$uid==uids]<-RU_CLASS_perc$area[RU_perc$uid==uids]/sum(RU_CLASS_perc$area[RU_perc$uid==uids])*100
}

lls<-levels(SOM_Label)
for (uids in TEMP_bb$uid){
  BO_SOM<-lls[!(lls %in% BO_CLASS_perc$SOM[BO_CLASS_perc$uid==uids])]
  BO_CLASS_perc<-rbind.data.frame(BO_CLASS_perc,data.frame(SOM=BO_SOM,
                                                           uid=rep(uids,length(BO_SOM)),
                                                           area=rep(0,length(BO_SOM))))
  RU_SOM<-lls[!(lls %in% RU_CLASS_perc$SOM[RU_CLASS_perc$uid==uids])]
  RU_CLASS_perc<-rbind.data.frame(RU_CLASS_perc,data.frame(SOM=RU_SOM,
                                                           uid=rep(uids,length(RU_SOM)),
                                                                   area=rep(0,length(RU_SOM))))
}

BO_CLASS_perc<-BO_CLASS_perc[order(BO_CLASS_perc$uid,BO_CLASS_perc$SOM),]
RU_CLASS_perc<-RU_CLASS_perc[order(RU_CLASS_perc$uid,RU_CLASS_perc$SOM),]

tot_CLASS_perc<-rbind.data.frame(cbind.data.frame(BO_CLASS_perc,method='BO'),
                                 cbind.data.frame(RU_CLASS_perc,method='RU'))

tot_CLASS_perc$SOM<-as.factor(tot_CLASS_perc$SOM)
tot_CLASS_perc$uid<-as.factor(tot_CLASS_perc$uid)
tot_CLASS_perc$method<-as.factor(tot_CLASS_perc$method)
library(rstatix)
TEMP<-group_by(tot_CLASS_perc,SOM)
TEMP<-t_test(TEMP,area~method,paired = T,detailed = F,p.adjust.method = 'none')
p.adjust(TEMP$p,method = 'fdr')

postscript(file = "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/SOM_MAPS/Class_area.eps",
           onefile = F,
           paper = 'special',
           bg = 'white',
           horizontal = F,
           height = 10,
           width = 10)

lls<-levels(SOM_Label)

plot.new()
par(mfrow=c(5,5),mar=c(0,2,1,0))



for (cl in 1:length(lls)){
  plot(NA,
       xlim=c(0.5,2.5),
       ylim=c(0,max(c(BO_CLASS_perc$area[BO_CLASS_perc$SOM==lls[cl]],
                      RU_CLASS_perc$area[BO_CLASS_perc$SOM==lls[cl]]))),
       xaxt='n',
       main=paste0('SOM: ',lls[cl]))
  
  for (uids in TEMP_bb$uid){
  BO_y<-BO_CLASS_perc$area[BO_CLASS_perc$uid==uids & BO_CLASS_perc$SOM==lls[cl]]
  RU_y<-RU_CLASS_perc$area[RU_CLASS_perc$uid==uids & RU_CLASS_perc$SOM==lls[cl]]
  if  (length(BO_y)==0) BO_y<-0
  if  (length(RU_y)==0) RU_y<-0
  points(x=1,
         y=BO_y,
         pch=16)
  points(x=2,
         y=RU_y,
         pch=17)
  lines(matrix(c(1,
                 BO_y,
                 2,
                 RU_y),
               ncol=2,
               byrow = T))
}
}
dev.off()





for(uids in unique(BO_area$uid)){
  tiff(file.path(
    "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/analysis/RUNIMC1/test/Density_area",
    paste0(uids,'.tiff')),
    width = 10,
    height = 25,
    units = 'in',
    res = 300,
    compression='lzw')
  par(mfrow=c(25,1),mar=c(0,0,0,0))
  i=1
  for (sm in unique(BO_area$SOM)){
    if (length(BO_area[BO_area$uid==uids & BO_area$SOM==sm,'BO_area'])>3){
      BO_D<-density(BO_area[BO_area$uid==uids & BO_area$SOM==sm,'BO_area'])
    } else {(BO_D<-list(x=c(0,300),y=c(1,1)))}
    if (length(RU_area[RU_area$uid==uids & RU_area$SOM==sm,'RU_area'])>3){
      RU_D<-density(RU_area[RU_area$uid==uids & RU_area$SOM==sm,'RU_area'])
    } else {(RU_D<-list(x=c(0,300),y=c(1,1)))}
    plot(NA,xlim=c(0,300),ylim=c(0,1),ann=F)
    lines(BO_D$x,BO_D$y/max(BO_D$y),col='blue')
    lines(RU_D$x,RU_D$y/max(RU_D$y),col='red')
    
  }
  dev.off()
}



total_area<-dplyr::full_join(BO_area_aggr,
                             RU_area_aggr,
                             by=c('uid','SOM'))


rst<-RUNIMC:::IMCstackOpen(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/rasterStacks/Sample2_1_ROI002.txt.stk"
)

nms<-names(rst)

RGBmat<-matrix(c(1,0,0,
                 0,1,0,
                 0,0,1),nrow=3,byrow = F,dimnames = list(c('R','G','B'),nms[c(29,29,29)]))

rst2<-lapply(names(rst)[c(29,29,29)],function(x){
  rst<-RUNIMC::quantNorm(rst[[x]],0.95)
  RGBrst<-lapply(RGBmat[,x],function(y){
    rst*y
  })
})

rst3<-lapply(c('R','G','B'),function(clx){
  newrst<-lapply(rst2,'[[',clx)
  newrst<-raster::calc(raster::stack(newrst),sum)
  newrst[newrst>1]<-1
  return(newrst)
})

par(mar=c(3,3,3,3))

raster::plotRGB(raster::stack(rst3),scale=1)

plot(geo_list_BO_annot[
  geo_list_BO_annot$uid=='uid.65a13c3a5e3370b8e7ce28894957e0c0' ,'GEOMETRY'],col=NA,border='cyan',add=T)

plot(geo_list_RU_annot[
  geo_list_RU_annot$uid=='uid.65a13c3a5e3370b8e7ce28894957e0c0' ,'GEOMETRY'],col=NA,border='yellow',add=T)

plot(geo_list_RU_annot[
  geo_list_RU_annot$uid=='uid.65a13c3a5e3370b8e7ce28894957e0c0'&
    geo_list_RU_annot$SOM%in% c('05','10','20','25'),]['GEOMETRY'],col=NA,border='red')

plot(geo_list_BO_annot[
  geo_list_BO_annot$uid=='uid.65a13c3a5e3370b8e7ce28894957e0c0'&
    geo_list_BO_annot$SOM%in% c('05','10','20','25'),]['GEOMETRY'],col=NA,border='blue',add=T)


####### chi squared ######
count_BO_total<-(aggregate(SOM~uid,geo_list_BO_annot,table))

count_RU_total<-aggregate(SOM~uid,geo_list_RU_annot,table)

uids<-unique(count_BO_total$uid)
names(uids)<-uids
counts_total<-lapply(uids,function(x){
  out<-(rbind(count_BO_total$SOM[count_BO_total$uid==x],
              count_RU_total$SOM[count_RU_total$uid==x]))
  rownames(out)<-c('BO','RU')
  colnames(out)<-colnames(count_BO_total[1,'SOM'])
  return(out)
})

counts_array<-array(unlist(counts_total),
                    dim=c(2,25,14),
                    dimnames = list(c('BO','RU'),
                                    colnames(count_BO_total[1,'SOM']),
                                    
                                    names(counts_total)))
mantelhaen.test(counts_array,exact = T)


ratios<-log(colMeans(count_BO_total[,-1])/colMeans(count_RU_total[,-1],2))

chisq.test(as.vector(count_BO_total[4,-1]),as.vector(count_RU_total[4,-1]))

dev.off()

getLevel <- function(x,y,prob=c(0.75,0.5,0.25)) {
  kk <- MASS::kde2d(x,y)
  dx <- diff(kk$x[1:2])
  dy <- diff(kk$y[1:2])
  sz <- sort(kk$z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}


pop<-'24'
TEMP<-getLevel(geo_list_BO_trans$x160gd.cd68.gd160di.[geo_list_BO_annot$SOM==pop],
               geo_list_BO_trans$x191ir.dna1.ir191di.[geo_list_BO_annot$SOM==pop])

TEMP1<-MASS::kde2d(geo_list_BO_trans$x160gd.cd68.gd160di.[geo_list_BO_annot$SOM==pop],
                   geo_list_BO_trans$x191ir.dna1.ir191di.[geo_list_BO_annot$SOM==pop])

contour(TEMP1,levels = TEMP,drawlabels = F)

TEMP<-getLevel(geo_list_RU_trans$x160gd.cd68.gd160di.[geo_list_RU_annot$SOM==pop],
               geo_list_RU_trans$x191ir.dna1.ir191di.[geo_list_RU_annot$SOM==pop])

TEMP1<-MASS::kde2d(geo_list_RU_trans$x160gd.cd68.gd160di.[geo_list_RU_annot$SOM==pop],
                   geo_list_RU_trans$x191ir.dna1.ir191di.[geo_list_RU_annot$SOM==pop])

contour(TEMP1,levels = TEMP,drawlabels = F,col = 'red',add = T)


contour(TEMP,)
contour(TEMP1,add=T)

plot(geo_list_BO_trans$x156gd.cd4.gd156di.[geo_list_BO_annot$SOM=='11'],
     geo_list_BO_trans$x162dy.cd8a.dy162di.[geo_list_BO_annot$SOM=='11'])

plot(geo_list_RU_trans$x156gd.cd4.gd156di.[geo_list_RU_annot$SOM=='11'],
     geo_list_RU_trans$x162dy.cd8a.dy162di.[geo_list_RU_annot$SOM=='11'])

plot(geo_list_BO_trans$x156gd.cd4.gd156di.[geo_list_BO_annot$SOM=='06'],
     geo_list_BO_trans$x162dy.cd8a.dy162di.[geo_list_BO_annot$SOM=='06'])

plot(geo_list_RU_trans$x156gd.cd4.gd156di.[geo_list_RU_annot$SOM=='06'],
     geo_list_RU_trans$x162dy.cd8a.dy162di.[geo_list_RU_annot$SOM=='06'])





TEMP<-geo_list_BO_annot[geo_list_BO_annot$SOM %in% c('10','05','25','20','23','22'),]

ggp<-ggplot2::ggplot(TEMP)+
  ggplot2::geom_point(ggplot2::aes(x=x156gd.cd4.gd156di.,
                                   y=x162dy.cd8a.dy162di.),size=1)+
  ggplot2::geom_density2d(ggplot2::aes(x=x156gd.cd4.gd156di.,
                                       y=x162dy.cd8a.dy162di.,
                                       col=SOM),size=1)+
  ggplot2::scale_x_continuous(trans = trns,limits = c(0,10))+
  ggplot2::scale_y_continuous(trans = trns,limits = c(0,20))+
  ggplot2::scale_colour_brewer(palette = "Set1")+
  ggplot2::theme_classic()


ggplot2::ggsave(filename = "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/BO_CD4_CD8.eps",
                plot = ggp,
                device = 'eps',
                width = 5,
                height = 5)

ggp<-ggplot2::ggplot(TEMP)+
  ggplot2::geom_point(ggplot2::aes(x=x170er.cd3.er170di.,
                                   y=x161dy.cd20.dy161di.),size=1)+
  ggplot2::geom_density2d(ggplot2::aes(x=x170er.cd3.er170di.,
                                       y=x161dy.cd20.dy161di.,
                                       col=SOM),size=1)+
  ggplot2::scale_x_continuous(trans = trns,limits = c(0,10))+
  ggplot2::scale_y_continuous(trans = trns,limits = c(0,10))+
  ggplot2::scale_colour_brewer(palette = "Set1")+
  ggplot2::theme_classic()

ggplot2::ggsave(filename = "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/BO_CD3_CD20.eps",
                plot = ggp,
                device = 'eps',
                width = 5,
                height = 5)

TEMP<-geo_list_RU_annot[geo_list_RU_annot$SOM %in% c(c('10','05','25','20','23','22')),]

ggp<-ggplot2::ggplot(TEMP)+
  ggplot2::geom_point(ggplot2::aes(x=x156gd.cd4.gd156di.,
                                   y=x162dy.cd8a.dy162di.),size=1)+
  ggplot2::geom_density2d(ggplot2::aes(x=x156gd.cd4.gd156di.,
                                       y=x162dy.cd8a.dy162di.,
                                       col=SOM),size=1)+
  ggplot2::scale_x_continuous(trans = trns,limits = c(0,10))+
  ggplot2::scale_y_continuous(trans = trns,limits = c(0,20))+
  ggplot2::scale_colour_brewer(palette = "Set1")+
  ggplot2::theme_classic()

ggplot2::ggsave(filename = "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/RU_CD4_CD8.eps",
                plot = ggp,
                device = 'eps',
                width = 5,
                height = 5)

ggp<-ggplot2::ggplot(TEMP)+
  ggplot2::geom_point(ggplot2::aes(x=x170er.cd3.er170di.,
                                   y=x161dy.cd20.dy161di.),size=1)+
  ggplot2::geom_density2d(ggplot2::aes(x=x170er.cd3.er170di.,
                                       y=x161dy.cd20.dy161di.,
                                       col=SOM),size=1)+
  ggplot2::scale_x_continuous(trans = trns,limits = c(0,10))+
  ggplot2::scale_y_continuous(trans = trns,limits = c(0,10))+
  ggplot2::scale_colour_brewer(palette = "Set1")+
  ggplot2::theme_classic()

ggplot2::ggsave(filename = "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/RU_CD3_CD20.eps",
                plot = ggp,
                device = 'eps',
                width = 5,
                height = 5)




TEMP<-geo_list_BO_annot[geo_list_BO_annot$SOM %in% c('11','12','16','17','21'),]

ggp<-ggplot2::ggplot(TEMP)+
  ggplot2::geom_point(ggplot2::aes(x=x155gd.foxp3.gd155di.,
                                   y=x175lu.cd25.lu175di.),size=1)+
  ggplot2::geom_density2d(ggplot2::aes(x=x155gd.foxp3.gd155di.,
                                       y=x175lu.cd25.lu175di.,
                                       col=SOM),size=1)+
  ggplot2::scale_x_continuous(trans = trns,limits = c(0,10))+
  ggplot2::scale_y_continuous(trans = trns,limits = c(0,20))+
  ggplot2::scale_colour_brewer(palette = "Set1")+
  ggplot2::theme_classic()


ggplot2::ggsave(filename = "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/BO_CD4_CD8.eps",
                plot = ggp,
                device = 'eps',
                width = 5,
                height = 5)


TEMP<-geo_list_RU_annot[geo_list_RU_annot$SOM %in% c('11','12','16','17','21'),]

ggp<-ggplot2::ggplot(TEMP)+
  ggplot2::geom_point(ggplot2::aes(x=x155gd.foxp3.gd155di.,
                                   y=x175lu.cd25.lu175di.),size=1)+
  ggplot2::geom_density2d(ggplot2::aes(x=x155gd.foxp3.gd155di.,
                                       y=x175lu.cd25.lu175di.,
                                       col=SOM),size=1)+
  ggplot2::scale_x_continuous(trans = trns,limits = c(0,10))+
  ggplot2::scale_y_continuous(trans = trns,limits = c(0,20))+
  ggplot2::scale_colour_brewer(palette = "Set1")+
  ggplot2::theme_classic()


ggplot2::ggsave(filename = "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/BO_CD4_CD8.eps",
                plot = ggp,
                device = 'eps',
                width = 5,
                height = 5)




TEMP<-geo_list_RU_annot[geo_list_RU_annot$SOM %in% c('20','10','15','04','05','09'),]

ggp<-ggplot2::ggplot(TEMP)+
  ggplot2::geom_density(ggplot2::aes(x=area,color=SOM,y=..density..),size=1)+
  ggplot2::scale_x_continuous(limits = c(0,300))+
  ggplot2::scale_colour_brewer(palette = "Set1")+
  ggplot2::theme_classic()

ggplot2::ggsave(filename = "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/area_RU.eps",
                plot = ggp,
                device = 'eps',
                width = 5,
                height = 5)

TEMP<-geo_list_BO_annot[geo_list_BO_annot$SOM %in% c('20','10','15','04','05','09','21','16'),]

ggp<-ggplot2::ggplot(TEMP)+
  ggplot2::geom_density(ggplot2::aes(x=area,color=SOM,y=..density..),size=1)+
  ggplot2::scale_x_continuous(limits = c(0,300))+
  ggplot2::scale_colour_brewer(palette = "Set1")+
  ggplot2::theme_classic()

ggplot2::ggsave(filename = "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/area_BO.eps",
                plot = ggp,
                device = 'eps',
                width = 5,
                height = 5)

TEMP1<-sf::st_read(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/RUN_001/analysis/RUNIMC1/training/polygons/Sample1_1_ROI001.txt.sqlite"  
)

TEMP2<-sf::st_read(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/RUN_001/analysis/RUNIMC1/training/polygons/Sample1_1_ROI002.txt.sqlite"  
)

TEMP3<-sf::st_read(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/RUN_001/analysis/RUNIMC1/training/polygons/Sample1_1_ROI003.txt.sqlite"  
)

geo_list_gr<-dplyr::bind_rows(TEMP1,TEMP2,TEMP3)

rm(TEMP1)
rm(TEMP2)
rm(TEMP3)

GR_area<-sf::st_area(geo_list_gr)

geo_list_gr<-dplyr::bind_cols(geo_list_gr,data.frame(area=GR_area))

ggp<-ggplot2::ggplot(geo_list_gr)+
  ggplot2::geom_density(ggplot2::aes(x=area,color=label,y=..density..),size=1)+
  ggplot2::scale_x_continuous(limits = c(0,300))+
  ggplot2::scale_colour_brewer(palette = "Set1")+
  ggplot2::theme_classic()

ggplot2::ggsave(filename = "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/area_GR.eps",
                plot = ggp,
                device = 'eps',
                width = 5,
                height = 5)


BO<-aggregate(area~SOM,geo_Flist_BO_annot,sum)
RU<-aggregate(area~SOM,geo_list_RU_annot,sum)
BO<-cbind.data.frame(BO,Method='BO')
RU<-cbind.data.frame(RU,Method='RU')
TOT<-rbind.data.frame(BO,RU)
RATIO<-data.frame(SOM=BO$SOM,RATIO=log2(BO$area/RU$area))
ggp<-ggplot2::ggplot(RATIO)+ggplot2::geom_col(ggplot2::aes(x=SOM,y=RATIO))+
  ggplot2::theme_classic()

ggplot2::ggsave(filename = "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/RATIO.eps",
                plot = ggp,
                device = 'eps',
                width = 5,
                height = 2)


TEMP<-geo_list_BO_annot[geo_list_BO_annot$SOM %in% c('13'),]

ggp<-ggplot2::ggplot(TEMP)+
  ggplot2::geom_point(ggplot2::aes(x=x160gd.cd68.gd160di.,
                                   y=x191ir.dna1.ir191di.),size=1)+
  ggplot2::geom_density2d(ggplot2::aes(x=x160gd.cd68.gd160di.,
                                       y=x191ir.dna1.ir191di.,
                                       col=SOM),size=1)+
  ggplot2::scale_x_continuous(trans = trns,limits = c(0,40))+
  ggplot2::scale_y_continuous(trans = trns,limits = c(0,40))+
  ggplot2::scale_colour_brewer(palette = "Set1")+
  ggplot2::theme_classic()

ggplot2::ggsave(filename = "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/BO_CD68_DNA.eps",
                plot = ggp,
                device = 'eps',
                width = 5,
                height = 5)

TEMP<-geo_list_RU_annot[geo_list_RU_annot$SOM %in% c('13'),]

ggp<-ggplot2::ggplot(TEMP)+
  ggplot2::geom_point(ggplot2::aes(x=x160gd.cd68.gd160di.,
                                   y=x191ir.dna1.ir191di.),size=1)+
  ggplot2::geom_density2d(ggplot2::aes(x=x160gd.cd68.gd160di.,
                                       y=x191ir.dna1.ir191di.,
                                       col=SOM),size=1)+
  ggplot2::scale_x_continuous(trans = trns,limits = c(0,40))+
  ggplot2::scale_y_continuous(trans = trns,limits = c(0,40))+
  ggplot2::scale_colour_brewer(palette = "Set1")+
  ggplot2::theme_classic()

ggplot2::ggsave(filename = "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/RU_CD68_DNA.eps",
                plot = ggp,
                device = 'eps',
                width = 5,
                height = 5)



TEMP<-geo_list_RU_annot[geo_list_RU_annot$SOM %in% c('20','10','15','04','05','09','25','13','24','21','16','06','11'),]

NEWSOM<-vector('character',length = nrow(TEMP))

NEWSOM[TEMP$SOM %in% c('10','15','20')]<-'CD4'
NEWSOM[TEMP$SOM %in% c('05','09')]<-'CD8'
NEWSOM[TEMP$SOM %in% c('04')]<-'CD20'
NEWSOM[TEMP$SOM %in% c('25','13','24')]<-'Myeloid'
NEWSOM[TEMP$SOM %in% c('21','16','06','11')]<-'Epithelium'
NEWSOM<-factor(NEWSOM)
TEMP$SOM<-NEWSOM
ggp<-ggplot2::ggplot(TEMP)+
  ggplot2::geom_density(ggplot2::aes(x=area,color=SOM,y=..density..),size=1)+
  ggplot2::scale_x_continuous(limits = c(0,200))+
  ggplot2::scale_colour_brewer(palette = "Set1")+
  ggplot2::theme_classic()

ggplot2::ggsave(filename = "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/area_RU_coa.eps",
                plot = ggp,
                device = 'eps',
                width = 7,
                height = 5)

TEMP<-geo_list_BO_annot[geo_list_BO_annot$SOM %in% c('20','10','15','04','05','09','25','13','24','21','16','06','11'),]

NEWSOM<-vector('character',length = nrow(TEMP))

NEWSOM[TEMP$SOM %in% c('10','15','20')]<-'CD4'
NEWSOM[TEMP$SOM %in% c('05','09')]<-'CD8'
NEWSOM[TEMP$SOM %in% c('04')]<-'CD20'
NEWSOM[TEMP$SOM %in% c('25','13','24')]<-'Myeloid'
NEWSOM[TEMP$SOM %in% c('21','16','06','11')]<-'Epithelium'
NEWSOM<-factor(NEWSOM)
TEMP$SOM<-NEWSOM
ggp<-ggplot2::ggplot(TEMP)+
  ggplot2::geom_density(ggplot2::aes(x=area,color=SOM,y=..density..),size=1)+
  ggplot2::scale_x_continuous(limits = c(0,200))+
  ggplot2::scale_colour_brewer(palette = "Set1")+
  ggplot2::theme_classic()

ggplot2::ggsave(filename = "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/area_BO_coa.eps",
                plot = ggp,
                device = 'eps',
                width = 7,
                height = 5)


TEMP<-geo_list_gr[geo_list_gr$label %in% c('CD4','CD8','CD20','MLYD','TMR'),]
NEWLABEL<-vector('character',nrow((TEMP)))
NEWLABEL[TEMP$label=='CD4']<-'CD4'
NEWLABEL[TEMP$label=='CD8']<-'CD8'
NEWLABEL[TEMP$label=='CD20']<-'CD20'
NEWLABEL[TEMP$label=='MLYD']<-'Myeloid'
NEWLABEL[TEMP$label=='TMR']<-'Epithelium'
NEWLABEL<-factor(NEWLABEL)
TEMP$label<-NEWLABEL
ggp<-ggplot2::ggplot(TEMP)+
  ggplot2::geom_density(ggplot2::aes(x=area,color=label,y=..density..),size=1)+
  ggplot2::scale_x_continuous(limits = c(0,200))+
  ggplot2::scale_colour_brewer(palette = "Set1")+
  ggplot2::theme_classic()

ggplot2::ggsave(filename = "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/area_GR_coa.eps",
                plot = ggp,
                device = 'eps',
                width = 7,
                height = 5)
