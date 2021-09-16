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

colnames(geo_list_RU_trans)<-unlist(lapply(strsplit(colnames(geo_list_RU_trans),'.',fixed = T),'[',2))

geo_list_RU_trans<-cbind(uid=geo_list_RU$uid[drop=T],geo_list_RU_trans,area=sf::st_area(geo_list_RU))

cln<-colnames(geo_list_BO_trans)[2:17]

cln<-c('cd3','cd8','cd4','cd20','cd21','cd68','cd31','asma','cd16')
trsh<-c(1.25, 1.75, 1.5, 0.75, 1.5, 2.75, 1.25, 1.5, 0.8)
names(trsh)<-cln

for (i in 1:3){
  tiff(
    filename = file.path(
      "C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/scatters_DNA",
      paste0('RU_alternative',i,'.tiff')
    ),
    width = 2,
    height = 18,
    units = 'in',
    res = 300,
    compression = 'lzw'
  )
  par(mfrow=c(length(cln),1))
  par(mar=c(4,4,1,1))
  cln_name<-c('CD3','CD8','CD4','CD20','CD21','CD68','CD31','ACTA2','CD16')
  cln_name<-setNames(cln_name,cln)
  for (cc in 'dna'){
    for (rr in cln){
        TEMP<-geo_list_RU_trans[geo_list_RU_trans$uid==unique(geo_list_RU_trans$uid)[i],c(rr,cc)]
      
          TEMPD<-ks::kde(TEMP)
        
        newZ<-predict(object = TEMPD,x=TEMP)
        
        newZ<-cut(newZ,
                  include.lowest = T,
                  breaks = c(-Inf,TEMPD$cont[c(T,F,F,F,F,F,F,F,F,F,F)],+Inf),
                  labels=F)
        colZ<-setNames(1:max(newZ),colorRampPalette(c('blue','green','red'))(max(newZ)))
        newZ<-names(colZ)[newZ]
        
        plot(TEMP,col=newZ,cex=0.1,cex.lab=0.8,xlab=cln_name[rr],ylab='DNA',ylim=c(0,2.5))
        abline(v=trsh[rr],h=1,lty=3,col='gray10')
      
    }
  }
  dev.off()
}

for (i in 1:3){
  tiff(
    filename = file.path(
      "C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/scatters_DNA",
      paste0('BO_alternative',i,'.tiff')
    ),
    width = 2,
    height = 18,
    units = 'in',
    res = 300,
    compression = 'lzw'
  )
  par(mfrow=c(length(cln),1))
  par(mar=c(4,4,1,1))
  cln_name<-c('CD3','CD8','CD4','CD20','CD21','CD68','CD31','ACTA2','CD16')
  cln_name<-setNames(cln_name,cln)
  for (cc in 'dna'){
    for (rr in cln){
      TEMP<-geo_list_BO_trans[geo_list_BO_trans$uid==unique(geo_list_BO_trans$uid)[i],c(rr,cc)]
      
      TEMPD<-ks::kde(TEMP)
      
      newZ<-predict(object = TEMPD,x=TEMP)
      
      newZ<-cut(newZ,
                include.lowest = T,
                breaks = c(-Inf,TEMPD$cont[c(T,F,F,F,F,F,F,F,F,F,F)],+Inf),
                labels=F)
      colZ<-setNames(1:max(newZ),colorRampPalette(c('blue','green','red'))(max(newZ)))
      newZ<-names(colZ)[newZ]
      
      plot(TEMP,col=newZ,cex=0.1,cex.lab=0.8,xlab=cln_name[rr],ylab='DNA',ylim=c(0,2.5))
      abline(v=trsh[rr],h=1,lty=3,col='gray10')
      
    }
  }
  dev.off()
}
