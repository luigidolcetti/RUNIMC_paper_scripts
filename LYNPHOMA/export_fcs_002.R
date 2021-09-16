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

colnames(geo_list_BO_trans)<-unlist(lapply(strsplit(colnames(geo_list_BO_trans),'.',fixed = T),'[',2))

geo_list_BO_trans<-cbind(uid=geo_list_BO$uid[drop=T],geo_list_BO_trans,area=sf::st_area(geo_list_BO))

fcs_BO<-lapply(unique(geo_list_BO_trans$uid),function(x,y=geo_list_BO_trans){
  TEMP<-as.matrix(y[y$uid==x,-1])
  out<-flowCore::flowFrame(TEMP)})

colnames(geo_list_RU_trans)<-unlist(lapply(strsplit(colnames(geo_list_RU_trans),'.',fixed = T),'[',2))

geo_list_RU_trans<-cbind(uid=geo_list_RU$uid[drop=T],geo_list_RU_trans,area=sf::st_area(geo_list_RU))

fcs_RU<-lapply(unique(geo_list_RU_trans$uid),function(x,y=geo_list_RU_trans){
  TEMP<-as.matrix(y[y$uid==x,-1])
  out<-flowCore::flowFrame(TEMP)})

for (i in 1:3){
  flowCore::write.FCS(fcs_BO[[i]],
                      paste0("C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/FCS/BO_",i))
  flowCore::write.FCS(fcs_RU[[i]],
                      paste0("C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/FCS/RU_",i))
}


TEMP<-geo_list_RU_trans[geo_list_RU_trans$uid==unique(geo_list_RU_trans$uid)[1],c('cd68','cd20')]
TEMP<-geo_list_BO_trans[geo_list_BO_trans$uid==unique(geo_list_BO_trans$uid)[1],c('cd68','cd20')]

cln<-colnames(geo_list_BO_trans)[2:17]

cln<-c('cd3','cd8','cd4','cd20','cd21','cd68','cd31','asma','cd16')
trsh<-c(1.25, 1.75, 1.5, 0.75, 1.5, 2.75, 1.25, 1.5, 0.8)
names(trsh)<-cln

for (i in 1:3){
  tiff(
    filename = file.path(
      "C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/scatters",
      paste0('RU_alternative',i,'.tiff')
    ),
    width = 5,
    height = 5,
    units = 'in',
    res = 300,
    compression = 'lzw'
  )
  par(mfrow=c(length(cln),length(cln)))
  par(mai=c(0.02,0.02,0.02,0.02))
  cln_name<-c('CD3','CD8','CD4','CD20','CD21','CD68','CD31','ACTA2','CD16')
  cln_name<-setNames(cln_name,cln)
  for (cc in cln){
    for (rr in cln){
      if (cc!=rr){
        
        TEMP<-geo_list_RU_trans[geo_list_RU_trans$uid==unique(geo_list_RU_trans$uid)[i],c(rr,cc)]
        # TEMP<-TEMP[sample(1:nrow(TEMP),1000),]
        
        TEMPD<-ks::kde(TEMP)
        
        newZ<-predict(object = TEMPD,x=TEMP)
        
        newZ<-cut(newZ,
                  include.lowest = T,
                  breaks = c(-Inf,TEMPD$cont[c(T,F,F,F,F,F,F,F,F,F,F)],+Inf),
                  labels=F)
        colZ<-setNames(1:max(newZ),colorRampPalette(c('blue','green','red'))(max(newZ)))
        newZ<-names(colZ)[newZ]
        plot(TEMP,col=newZ,cex=0.2,xaxt='n',yaxt='n')
        abline(v=trsh[rr],h = trsh[cc],lty=3,col='gray10')
      } else {
        plot(NA,xaxt='n',yaxt='n',xlim=c(0,1),ylim=c(0,1))
        text(labels=cln_name[cc],x=0.5,y=0.5,cex=1.5)
      }
      
    }
  }
  dev.off()
}

for (i in 1:3){
  tiff(
    filename = file.path(
      "C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/scatters",
      paste0('BO_alternative',i,'.tiff')
    ),
    width = 5,
    height = 5,
    units = 'in',
    res = 300,
    compression = 'lzw'
  )
  par(mfrow=c(length(cln),length(cln)))
  par(mai=c(0.02,0.02,0.02,0.02))
  
  cln_name<-c('CD3','CD8','CD4','CD20','CD21','CD68','CD31','ACTA2','CD16')
  cln_name<-setNames(cln_name,cln)
  
  for (cc in cln){
    for (rr in cln){
      if (cc!=rr){
        
        TEMP<-geo_list_BO_trans[geo_list_BO_trans$uid==unique(geo_list_BO_trans$uid)[i],c(rr,cc)]
        # TEMP<-TEMP[sample(1:nrow(TEMP),1000),]
        
        TEMPD<-ks::kde(TEMP)
        
        newZ<-predict(object = TEMPD,x=TEMP)
        
        newZ<-cut(newZ,
                  include.lowest = T,
                  breaks = c(-Inf,TEMPD$cont[c(T,F,F,F,F,F,F,F,F,F,F)],+Inf),
                  labels=F)
        colZ<-setNames(1:max(newZ),colorRampPalette(c('blue','green','red'))(max(newZ)))
        newZ<-names(colZ)[newZ]
        plot(TEMP,col=newZ,cex=0.2,xaxt='n',yaxt='n')
        abline(v=trsh[rr],h = trsh[cc],lty=3,col='gray10')
      } else {
        plot(NA,xaxt='n',yaxt='n',xlim=c(0,1),ylim=c(0,1))
        text(labels=cln_name[cc],x=0.5,y=0.5,cex=1.5)
      }
      
    }
  }
  dev.off()
}


uids<-unique(geo_list_BO_trans$uid)
geo_list_BO_Class<-dplyr::bind_cols(geo_list_BO[,'uid'],CLASS='UNDF')
geo_list_RU_Class<-dplyr::bind_cols(geo_list_RU[,'uid'],CLASS='UNDF')

stmnt<-expression(cd3>trsh['cd3'] & cd4>trsh['cd4'] & cd8<trsh['cd8'])

geo_list_BO_Class$CLASS[with(geo_list_BO_trans,eval(stmnt)) & geo_list_BO_Class$CLASS=='UNDF']<-'1_CD4-Lymph'
geo_list_RU_Class$CLASS[with(geo_list_RU_trans,eval(stmnt)) & geo_list_RU_Class$CLASS=='UNDF']<-'1_CD4-Lymph'

stmnt<-expression(cd3>trsh['cd3'] & cd4<trsh['cd4'] & cd8>trsh['cd8'])

geo_list_BO_Class$CLASS[with(geo_list_BO_trans,eval(stmnt)) & geo_list_BO_Class$CLASS=='UNDF']<-'2_CD8-Lymph'
geo_list_RU_Class$CLASS[with(geo_list_RU_trans,eval(stmnt)) & geo_list_RU_Class$CLASS=='UNDF']<-'2_CD8-Lymph'

stmnt<-expression(cd3<trsh['cd3'] & cd68>trsh['cd68'])

geo_list_BO_Class$CLASS[with(geo_list_BO_trans,eval(stmnt)) & geo_list_BO_Class$CLASS=='UNDF']<-'7_mph'
geo_list_RU_Class$CLASS[with(geo_list_RU_trans,eval(stmnt)) & geo_list_RU_Class$CLASS=='UNDF']<-'7_mph'

stmnt<-expression(cd3<trsh['cd3'] & cd68<trsh['cd68'] & asma>trsh['asma'])

geo_list_BO_Class$CLASS[with(geo_list_BO_trans,eval(stmnt)) & geo_list_BO_Class$CLASS=='UNDF']<-'8_fbr'
geo_list_RU_Class$CLASS[with(geo_list_RU_trans,eval(stmnt)) & geo_list_RU_Class$CLASS=='UNDF']<-'8_fbr'

stmnt<-expression(cd3<trsh['cd3'] & cd68<trsh['cd68'] & cd31>trsh['cd31'])

geo_list_BO_Class$CLASS[with(geo_list_BO_trans,eval(stmnt)) & geo_list_BO_Class$CLASS=='UNDF']<-'9_endt'
geo_list_RU_Class$CLASS[with(geo_list_RU_trans,eval(stmnt)) & geo_list_RU_Class$CLASS=='UNDF']<-'9_endt'

stmnt<-expression(cd20<trsh['cd20'] & cd21>trsh['cd21'])
geo_list_BO_Class$CLASS[with(geo_list_BO_trans,eval(stmnt)) & geo_list_BO_Class$CLASS=='UNDF']<-'5_CD21+'
geo_list_RU_Class$CLASS[with(geo_list_RU_trans,eval(stmnt)) & geo_list_RU_Class$CLASS=='UNDF']<-'5_CD21+'

stmnt<-expression(cd20>trsh['cd20'] & cd21>trsh['cd21'])
geo_list_BO_Class$CLASS[with(geo_list_BO_trans,eval(stmnt)) & geo_list_BO_Class$CLASS=='UNDF']<-'4_B-cellCD21+'
geo_list_RU_Class$CLASS[with(geo_list_RU_trans,eval(stmnt)) & geo_list_RU_Class$CLASS=='UNDF']<-'4_B-cellCD21+'

stmnt<-expression(cd16>trsh['cd16'] )
geo_list_BO_Class$CLASS[with(geo_list_BO_trans,eval(stmnt))& geo_list_BO_Class$CLASS=='UNDF']<-'6_nk'
geo_list_RU_Class$CLASS[with(geo_list_RU_trans,eval(stmnt))& geo_list_RU_Class$CLASS=='UNDF']<-'6_nk'

stmnt<-expression(cd20>trsh['cd20'])
geo_list_BO_Class$CLASS[with(geo_list_BO_trans,eval(stmnt)) & geo_list_BO_Class$CLASS=='UNDF']<-'3_B-Cells'
geo_list_RU_Class$CLASS[with(geo_list_RU_trans,eval(stmnt)) & geo_list_RU_Class$CLASS=='UNDF']<-'3_B-Cells'

ppal<-c('#fce96a','#8eff78','#85edeb','#85c3ed','#8faaf7','#cd84f0','#ed8cad','#fc3232','#ffb04a','#787773')

for (i in 1:3){
  postscript(file = paste0("C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/scatters/BO_map_",i,".eps"),
             height = 30,
             width = 30,
             onefile = F,
             horizontal = F,
             bg = 'white',
             paper = 'special',
             pointsize = 16)
  par(mar=c(1,0,0,0))
  plot(geo_list_BO_Class[geo_list_BO_Class$uid==uids[i],'CLASS'],key.pos = 1,key.length = 1.1,lwd=0.1,pal=ppal)
  dev.off()
  postscript(file = paste0("C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/scatters/RU_map_",i,".eps"),
             height = 30,
             width = 30,
             onefile = F,
             horizontal = F,
             bg = 'white',
             paper = 'special',
             pointsize = 16)
  plot(geo_list_RU_Class[geo_list_RU_Class$uid==uids[i],'CLASS'],key.pos = 1,key.length = 1.1,lwd=0.1,pal=ppal)
  dev.off()
}

geo_list_RU_Class<-dplyr::bind_cols(geo_list_RU_Class,area=sf::st_area(geo_list_RU_Class))
geo_list_BO_Class<-dplyr::bind_cols(geo_list_BO_Class,area=sf::st_area(geo_list_BO_Class))

BO_C_tot<-aggregate(area~uid,geo_list_BO_Class,sum)
RU_C_tot<-aggregate(area~uid,geo_list_RU_Class,sum)


BO_C_CLASS<-aggregate(area~uid+CLASS,geo_list_BO_Class,sum)
RU_C_CLASS<-aggregate(area~uid+CLASS,geo_list_RU_Class,sum)

for (i in uids){
  BO_C_CLASS$area[BO_C_CLASS$uid==i]<-BO_C_CLASS$area[BO_C_CLASS$uid==i]/BO_C_tot$area[BO_C_tot$uid==i]*100
  RU_C_CLASS$area[RU_C_CLASS$uid==i]<-RU_C_CLASS$area[RU_C_CLASS$uid==i]/RU_C_tot$area[RU_C_tot$uid==i]*100
  
  
}

TOT_CLASS<-rbind.data.frame(cbind.data.frame(RU_C_CLASS,method = 'RU'),
                            cbind.data.frame(BO_C_CLASS,method = 'BO'))


anv<-aov(area~as.factor(CLASS)*as.factor(method),data = TOT_CLASS)
summary(anv)
TukeyHSD(anv)

tt<-unlist(lapply(setNames(unique(TOT_CLASS$CLASS),unique(TOT_CLASS$CLASS)),function(x){
  t.test(TOT_CLASS$area[TOT_CLASS$CLASS==x & TOT_CLASS$method=='BO'],
         TOT_CLASS$area[TOT_CLASS$CLASS==x & TOT_CLASS$method=='RU'],paired = T)$p.value
}))

tt<-p.adjust(tt,method = 'bonferroni')
tt
# plot(y=rep(NA,length(TOT_CLASS$CLASS[TOT_CLASS$method=='BO'])),
#      x=as.factor(TOT_CLASS$CLASS[TOT_CLASS$method=='BO']),ylim=c(0,50))
# points(y=TOT_CLASS$area[TOT_CLASS$method=='BO'],
#      x=as.factor(TOT_CLASS$CLASS[TOT_CLASS$method=='BO']),col='red')
# 
# points(y=TOT_CLASS$area[TOT_CLASS$method=='RU'],
#        x=as.factor(TOT_CLASS$CLASS[TOT_CLASS$method=='RU']),col='blue')
classes<-setNames(unique(TOT_CLASS$CLASS),unique(TOT_CLASS$CLASS))
met<-setNames(c('BO','RU'),c('BO','RU'))
postscript(file = paste0("C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/scatters/classes.eps"),
           height = 7,
           width = 12,
           onefile = F,
           horizontal = F,
           bg = 'white',
           paper = 'special',
           pointsize = 12)

plot(NA,xlim=c(0.5,10.5),ylim=c(0,60),xaxt='n',xlab=NA,ylab='Area percentage')
axis(side = 1,
     at = 1:length(classes),
     labels = classes,
     las=1)
ccl<-setNames(c('red','green','blue'),uids)
for (cl in 1:length(classes)){
  for (u in uids){
    points(
      x=cl-0.2,
      y=TOT_CLASS$area[TOT_CLASS$CLASS==classes[cl] &
                         TOT_CLASS$uid==u &
                         TOT_CLASS$method==met[1]],
      pch=16,
      col=ccl[u])
    points(
      x=cl+0.2,
        y=TOT_CLASS$area[TOT_CLASS$CLASS==classes[cl] &
                         TOT_CLASS$uid==u &
                         TOT_CLASS$method==met[2]],
      pch=17,
      col=ccl[u])
    lines(matrix(c(cl-0.2,
                   TOT_CLASS$area[TOT_CLASS$CLASS==classes[cl] &
                                    TOT_CLASS$uid==u &
                                    TOT_CLASS$method==met[1]],
                   cl+0.2,
                   TOT_CLASS$area[TOT_CLASS$CLASS==classes[cl] &
                                    TOT_CLASS$uid==u &
                                    TOT_CLASS$method==met[2]]),
                 ncol=2,
                 byrow = T),
          col=ccl[u])
  }
}
dev.off()




###### area comparisons ####
###### 
###### 


BO_area<-sf::st_area(geo_list_BO)
BO_area<-cbind.data.frame(area=BO_area,sf::st_drop_geometry(geo_list_BO)[,c('uid'),drop=F])
RU_area<-sf::st_area(geo_list_RU)
RU_area<-cbind.data.frame(area=RU_area,sf::st_drop_geometry(geo_list_RU)[,c('uid'),drop=F])

TEMP_rst<-RUNIMC:::retrieve.RsCollection(
  "C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/rasterStacks"
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

postscript(file = "C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/scatters/Total_area.eps",
           onefile = F,
           paper = 'special',
           bg = 'white',
           horizontal = F,
           height = 5,
           width = 5)

plot(NA,xlim=c(0.5,2.5),ylim=c(70,100),xaxt='n',ylab='Area percentage')
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
