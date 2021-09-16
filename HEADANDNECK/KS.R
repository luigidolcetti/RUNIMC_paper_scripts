BO_TEMP<-geo_list_BO_annot
colnames(BO_TEMP)
BO_TEMP<-BO_TEMP[,c('uid','id','SOM')]
sf::write_sf(BO_TEMP,
             "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/analysis/BODEN/BP_SOM.sqlite")

RU_TEMP<-geo_list_RU_annot
colnames(RU_TEMP)
RU_TEMP<-RU_TEMP[,c(1:7,39)]
sf::write_sf(RU_TEMP,
             "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/analysis/RUNIMC1/RU_SOM.sqlite")
# 

BO_TEMP<-sf::read_sf(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/analysis/BODEN/BP_SOM.sqlite"
)

RU_TEMP<-sf::read_sf(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/analysis/RUNIMC1/RU_SOM.sqlite"
)

mystudy<-RUNIMC::retrieve(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/study.xml"
)

uids<-RUNIMC::st_uids(mystudy)

BO_SD<-lapply(uids,function(u){
  TEMP<-BO_TEMP[BO_TEMP$uid==u,]
  TEMPRST<-mystudy$raster[[u]]
  out<-exactextractr::exact_extract(TEMPRST,
                                    TEMP,
                                    fun='stdev')
})

BO_SD<-do.call(rbind.data.frame,BO_SD)
BO_SD<-cbind.data.frame(BO_TEMP,BO_SD)

RU_SD<-lapply(uids,function(u){
  TEMP<-RU_TEMP[RU_TEMP$uid==u,]
  TEMPRST<-mystudy$raster[[u]]
  out<-exactextractr::exact_extract(TEMPRST,
                                    TEMP,
                                    fun='stdev')
})

RU_SD<-do.call(rbind.data.frame,RU_SD)
RU_SD<-cbind.data.frame(RU_TEMP,RU_SD)


BO_density<-density(BO_SD$stdev.x161dy.cd20.dy161di.[BO_SD$som=='05'])
RU_density<-density(RU_SD$stdev.x161dy.cd20.dy161di.[RU_SD$som=='05'])
plot(BO_density)
lines(RU_density,col='red')

t.test(BO_SD$stdev.x161dy.cd20.dy161di.[BO_SD$som=='05'],RU_SD$stdev.x161dy.cd20.dy161di.[RU_SD$som=='05'])

ks.test(RU_SD$stdev.x160gd.cd68.gd160di.[RU_SD$som=='23'],BO_SD$stdev.x160gd.cd68.gd160di.[BO_SD$som=='23'],alternative = 'g')
ks.test(RU_SD$stdev.x161dy.cd20.dy161di.[RU_SD$som=='05'],BO_SD$stdev.x161dy.cd20.dy161di.[BO_SD$som=='05'],alternative = 'g')

TEMP<-ks.test(RU_SD$stdev.x162dy.cd8a.dy162di.[RU_SD$som=='10'],BO_SD$stdev.x162dy.cd8a.dy162di.[BO_SD$som=='10'],alternative = 'g')

SOM_I<-sort(unique(RU_SD$som))
names(SOM_I)<-SOM_I
CHL_I<-colnames(RU_SD)[10:39]
names(CHL_I)<-CHL_I
ks<-lapply(SOM_I,function(sI){
  vapply(CHL_I,function(cI){
    if (length(RU_SD[RU_SD$som==sI,cI])>3 &
        length(BO_SD[BO_SD$som==sI,cI])>3){
    out<-ks.test(RU_SD[RU_SD$som==sI,cI],BO_SD[BO_SD$som==sI,cI],alternative = 'g')
    return(out$p.value)} else {return(NULL)}
  },c(0))
})

ks<-do.call(rbind,ks)

RU_ks.adj<-matrix(p.adjust(ks,
             method = 'bonferroni'),
             ncol=ncol(ks),
             nrow=nrow(ks),
             byrow = F,
             dimnames = list(rownames(ks),colnames(ks)))



SOM_I<-sort(unique(RU_SD$som))
names(SOM_I)<-SOM_I
CHL_I<-colnames(RU_SD)[10:39]
names(CHL_I)<-CHL_I
ks<-lapply(SOM_I,function(sI){
  vapply(CHL_I,function(cI){
    if (length(RU_SD[RU_SD$som==sI,cI])>3 &
        length(BO_SD[BO_SD$som==sI,cI])>3){
      out<-ks.test(BO_SD[BO_SD$som==sI,cI],RU_SD[RU_SD$som==sI,cI],alternative = 'g')
      return(out$p.value)} else {return(NULL)}
  },c(0))
})

ks<-do.call(rbind,ks)

BO_ks.adj<-matrix(p.adjust(ks,
                           method = 'bonferroni'),
                  ncol=ncol(ks),
                  nrow=nrow(ks),
                  byrow = F,
                  dimnames = list(rownames(ks),colnames(ks)))

BO_val<-sf::st_drop_geometry(geo_list_BO_annot)
chnls<-colnames(BO_val)[3:32]
names(chnls)<-chnls
BO_val<-lapply(chnls,function(x){
  out<-aggregate(BO_val[,c(x)],list(SOM=BO_val[,'SOM']),mean)
  colnames(out)<-c('SOM',x)
  return(out[,2,drop=F])
})
BO_val<-do.call(cbind,BO_val)

BO_cor<-lapply(chnls,function(x){cor(BO_val[,x],BO_ks.adj[,paste0('stdev.',x)])})


RU_val<-sf::st_drop_geometry(geo_list_RU_annot)
chnls<-colnames(RU_val)[9:38]
names(chnls)<-chnls
RU_val<-lapply(chnls,function(x){
  out<-aggregate(RU_val[,c(x)],list(SOM=RU_val[,'SOM']),mean)
  colnames(out)<-c('SOM',x)
  return(out[,2,drop=F])
})
RU_val<-do.call(cbind,RU_val)

RU_cor<-lapply(chnls,function(x){cor(RU_val[,x],RU_ks.adj[,paste0('stdev.',x)])})

plot(RU_val[,chnls[5]],RU_ks.adj[,paste0('stdev.',chnls[5])])

plot(BO_val[,chnls[5]],BO_ks.adj[,paste0('stdev.',chnls[5])])
