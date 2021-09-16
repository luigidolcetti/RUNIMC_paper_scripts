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


rst<-RUNIMC:::retrieve.RsCollection(
  "C:/Users/luigi/Documents/BH_TEST/HEADANDNECK/RUN_001/rasterStacks"
)
uids<-unique(geo_list_BO$uid)


TEMP<-(RUNIMC::quantNorm(rst[[uids[8]]]$x160gd.cd68.gd160di.,0.99))-
             (RUNIMC::quantNorm(rst[[uids[8]]]$x162dy.cd8a.dy162di.,0.99))

tiff(filename=file.path(targetFolder,'CD68vsCD8.tiff'),
     width = 10,
     height = 10,
     units = 'in',
     pointsize = 12,
     compression = 'lzw',
     res = 600,
     bg = 'white')

raster::plot(TEMP,col=colorRampPalette(c('green','black','red'))(255),
             xlim=c(200,400),ylim=c(0,200),
             xaxt='n',
             yaxt='n',
             legend=F,
             bty='n')
dev.off()

postscript(file=file.path(targetFolder,'CD68vsCD8_BO.eps'),
           onefile = F,
           horizontal = F,
           paper = 'special',
           width = 16,
           height = 16,
           bg = 'white',
           pointsize = 18
)
raster::plot(TEMP,col=colorRampPalette(c('green','black','red'))(255),
             xlim=c(200,400),ylim=c(0,200),
             xaxt='n',
             yaxt='n',
             legend=F,
             bty='n')

plot(geo_list_BO_annot[geo_list_BO_annot$uid==uids[8] &
                       geo_list_BO_annot$SOM %in% c('20'),'GEOMETRY'],add=T,border='cyan',lwd=2)

dev.off()

postscript(file=file.path(targetFolder,'CD68vsCD8_RU.eps'),
           onefile = F,
           horizontal = F,
           paper = 'special',
           width = 16,
           height = 16,
           bg = 'white',
           pointsize = 18
)

raster::plot(TEMP,col=colorRampPalette(c('green','black','red'))(255),
             xlim=c(200,400),ylim=c(0,200),
             xaxt='n',
             yaxt='n',
             legend=F,
             bty='n')

plot(geo_list_RU_annot[geo_list_RU_annot$uid==uids[8] &
                         geo_list_RU_annot$SOM %in% c('20'),'GEOMETRY'],add=T,border='cyan',lwd=2)

dev.off()

TEMP<-(RUNIMC::quantNorm(rst[[uids[10]]]$x158gd.ecad.gd158di.,0.80))-
  (RUNIMC::quantNorm(rst[[uids[10]]]$x191ir.dna1.ir191di.,0.99))



tiff(filename=file.path(targetFolder,'ecadvsnuc.tiff'),
     width = 10,
     height = 10,
     units = 'in',
     pointsize = 12,
     compression = 'lzw',
     res = 600,
     bg = 'white')


raster::plot(TEMP,col=colorRampPalette(c('green','black','red'))(255),
             xlim=c(300,450),ylim=c(250,400),
             xaxt='n',
             yaxt='n',
             legend=F,
             bty='n')

dev.off()

postscript(file=file.path(targetFolder,'ecadvsnuc_BO.eps'),
           onefile = F,
           horizontal = F,
           paper = 'special',
           width = 16,
           height = 16,
           bg = 'white',
           pointsize = 18
)
raster::plot(TEMP,col=colorRampPalette(c('green','black','red'))(255),
             xlim=c(300,450),ylim=c(250,400),
             xaxt='n',
             yaxt='n',
             legend=F,
             bty='n')

plot(geo_list_BO_annot[geo_list_BO_annot$uid==uids[10] &
                         geo_list_BO_annot$SOM %in% c('01','02','03','04','05','07','08'),'GEOMETRY'],add=T,border='cyan',lwd=2)

dev.off()

postscript(file=file.path(targetFolder,'ecadvsnuc_RU.eps'),
           onefile = F,
           horizontal = F,
           paper = 'special',
           width = 16,
           height = 16,
           bg = 'white',
           pointsize = 18
)
raster::plot(TEMP,col=colorRampPalette(c('green','black','red'))(255),
             xlim=c(300,450),ylim=c(250,400),
             xaxt='n',
             yaxt='n',
             legend=F,
             bty='n')

plot(geo_list_RU_annot[geo_list_RU_annot$uid==uids[10] &
                         geo_list_RU_annot$SOM %in% c('01','02','03','04','05','07','08'),'GEOMETRY'],add=T,border='cyan',lwd=2)

dev.off()


aggregate(area~uid+SOM,geo_list_RU_annot,length)
