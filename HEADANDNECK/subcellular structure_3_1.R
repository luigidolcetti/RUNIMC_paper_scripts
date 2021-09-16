targetFolder<-"C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/SUBCELLULAR_3"
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

geo_list_total<-as.matrix(rbind(geo_list_BO_trans,geo_list_RU_trans))

dimx=5
dimy=5
meta=10

fsom<-readRDS(file = file.path(targetFolder,'som.R'  ))

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

geo_list_BO_annot<-dplyr::bind_cols(geo_list_BO,data.frame(SOM=SOM_Label[1:nrow(geo_list_BO)]))

geo_list_RU_annot<-dplyr::bind_cols(geo_list_RU,data.frame(SOM=factor(as.character(SOM_Label[(1+nrow(geo_list_BO)):(nrow(geo_list_BO)+nrow(geo_list_RU))]))))


geo_list_BO_annot<-dplyr::bind_cols(geo_list_BO_annot,data.frame(area=sf::st_area(geo_list_BO_annot)))
geo_list_RU_annot$area<-sf::st_area(geo_list_RU_annot)
                                    
Tot_BO_area<-aggregate(area~uid,sf::st_drop_geometry(geo_list_BO_annot),sum)
Tot_RU_area<-aggregate(area~uid,sf::st_drop_geometry(geo_list_RU_annot),sum)

par_BO_area<-aggregate(area~SOM+uid,sf::st_drop_geometry(geo_list_BO_annot),sum)
par_RU_area<-aggregate(area~SOM+uid,sf::st_drop_geometry(geo_list_RU_annot),sum)

for (u in unique(par_BO_area$uid)){
  par_BO_area$area[par_BO_area$uid==u]<-par_BO_area$area[par_BO_area$uid==u]/Tot_BO_area$area[Tot_BO_area$uid==u]*100
  par_RU_area$area[par_RU_area$uid==u]<-par_RU_area$area[par_RU_area$uid==u]/Tot_RU_area$area[Tot_RU_area$uid==u]*100
}


BO_long<-tidyr::pivot_wider(data = par_BO_area,
                             names_from = SOM,
                            values_from = area)
BO_long[is.na(BO_long)]<-0
RU_long<-tidyr::pivot_wider(data = par_RU_area,
                            names_from = SOM,
                            values_from = area)
RU_long[is.na(RU_long)]<-0
tt<-unlist(lapply(levels(SOM_Label),function(x) t.test(BO_long[,x,drop=T],RU_long[,x,drop=T],paired=T)$p.value))
tt<-p.adjust(tt,'bonferroni')
names(tt)<-levels(SOM_Label)
ss<-gtools::stars.pval(tt)

SL<-levels(SOM_Label)

postscript(file=file.path(targetFolder,'SOMxROI.eps'),
           onefile = F,
           horizontal = F,
           paper = 'special',
           width = 16,
           height = 8,
           bg = 'white',
           pointsize = 18
)
plot(NA,xlim=c(1,25),ylim=c(0,max(rbind.data.frame(BO_long[,-1],RU_long[,-1]))+3),bty='l',
     xaxt='n',xlab='FlowSom label',ylab='Area percentage')
axis(1,at=1:25,labels = SL,)


for (cc in 1:25){
  text(cc,max(rbind.data.frame(BO_long[,-1],RU_long[,-1])+3),labels = ss[cc])
  for (rr in 1:14){
    points(cc-0.2,BO_long[rr,SL[cc]],pch=16)
    points(cc+0.2,RU_long[rr,SL[cc]],pch=17)
    lines(matrix(c(cc-0.2,BO_long[rr,SL[cc]],cc+0.2,RU_long[rr,SL[cc]]),
                 ncol=2,
                 byrow = T))
  }
}

dev.off()
