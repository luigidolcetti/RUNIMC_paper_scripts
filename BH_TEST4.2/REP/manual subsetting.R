library(RUNIMC)
library(barbieHistologist)

rootFolder<-"C:/Users/luigi/Documents/BH_TEST/BH_TEST4.2"
targetFolder<-'RESULTS_SUBSETTING'
rpttn<-"REP_1"
trainingSampl<-0
brks<-c(1:4)
db<-factor(c('R1','R2','R3','R4'),levels=c('R1','R2','R3','R4'))


dir.create(file.path(rootFolder,rpttn,targetFolder))

mystudy<-retrieve(file.path(rootFolder,rpttn,'TEST_ROUND1/study.xml'))

TEMP_transform<-scales::modulus_trans(0)

geo_list_BO<-sf::st_read(
  file.path(rootFolder,rpttn,"TEST_ROUND1/analysis/BODEN1/ex.sqlite")
)
if (any(!sf::st_is_valid(geo_list_BO))){
  geo_list_BO<-sf::st_make_valid(geo_list_BO)
}

geo_list_RU<-sf::st_read(
  file.path(rootFolder,rpttn,"TEST_ROUND1/analysis/RUNIMC1/topLayers.sqlite")
)
if (any(!sf::st_is_valid(geo_list_RU))){
  geo_list_RU<-sf::st_make_valid(geo_list_RU)
}


geo_list_BO<-dplyr::bind_cols(geo_list_BO,data.frame(cell=rep('UNDF',nrow(geo_list_BO)),
                                                     newCell=rep('UNDF',nrow(geo_list_BO)),
                                                     Method=rep('BO',nrow(geo_list_BO))))
geo_list_RU<-dplyr::bind_cols(geo_list_RU,data.frame(cell=rep('UNDF',nrow(geo_list_RU)),
                                                     newCell=rep('UNDF',nrow(geo_list_RU)),
                                                     Method=rep('RU',nrow(geo_list_RU))))

geo_list_RU<-geo_list_RU[,-c(2:8)]

geo_list_GR_sub<-list()

for (densityBrake in brks){
  
  familyRaster<-bh_loadFamily(file.path(rootFolder,rpttn,paste0("sample_",densityBrake,"/FAM_",densityBrake)))
  familyVector<-sf::st_read(file.path(rootFolder,rpttn,paste0('sample_',densityBrake,'/GroundTruth_',densityBrake,'.sqlite')))
  familyVector<-familyVector[familyVector$compartment=='cytoplasm',]
  familyVector<-sf::st_drop_geometry(familyVector)
  familyVector<-familyVector[c(1,4)]
  colnames(familyVector)<-c('cell','id')
  
  geo_list_GR_sub[[densityBrake]]<-RUNIMC::lazyCatMap(familyRaster$ID,fn_indexToExclude = 0,fn_uid =st_uids(mystudy)[densityBrake])
  geo_list_GR_sub[[densityBrake]]<-extractMeanPixel(fn_raster = mystudy$raster,fn_polygons = geo_list_GR_sub[[densityBrake]])
  colnames(geo_list_GR_sub[[densityBrake]])[2]<-'id'
  colnames(geo_list_GR_sub[[densityBrake]])[15]<-'GEOMETRY'
  sf::st_geometry(geo_list_GR_sub[[densityBrake]])<-'GEOMETRY'
  geo_list_GR_sub[[densityBrake]]<-dplyr::inner_join(geo_list_GR_sub[[densityBrake]],familyVector,by='id')
  geo_list_GR_sub[[densityBrake]]<-dplyr::bind_cols(geo_list_GR_sub[[densityBrake]],data.frame(newCell=rep('UNDF',nrow(geo_list_GR_sub[[densityBrake]])),
                                                                                               Method=rep('GR',nrow(geo_list_GR_sub[[densityBrake]]))))
}

geo_list_GR<-do.call(dplyr::bind_rows,geo_list_GR_sub)
  
  
  
  geo_list_total<-dplyr::bind_rows(geo_list_GR,
                                   geo_list_BO,
                                   geo_list_RU)
  
  cell_types<-sort(unique(sf::st_drop_geometry(geo_list_GR)[,'cell']))
  channels<-colnames(sf::st_drop_geometry(geo_list_GR))[3:14]
  
  lower_bound<-matrix(-Inf,nrow=length(cell_types),ncol = length(channels),dimnames = list(cell_types,channels))
  upper_bound<-matrix(+Inf,nrow=length(cell_types),ncol = length(channels),dimnames = list(cell_types,channels))
  
  TEMP<-sf::st_drop_geometry(geo_list_GR)
  TEMP[3:14]<-TEMP_transform$transform(TEMP[3:14])
  
  lower_bound[1:12,channels[1]]<-2
  upper_bound[-(1:12),channels[1]]<-2
  
  lower_bound[2:9,channels[2]]<-2.1
  upper_bound[-(2:9),channels[2]]<-2.1
  
  lower_bound[2:5,channels[3]]<-0.9
  upper_bound[-(2:5),channels[3]]<-0.9
  
  lower_bound[6:9,channels[4]]<-0.35
  upper_bound[-(6:9),channels[4]]<-0.35
  
  lower_bound[c(4,5,8,9),channels[5]]<-2
  upper_bound[-c(4,5,8,9),channels[5]]<-2
  
  lower_bound[c(1,3,5,7,9),channels[6]]<-5.2
  upper_bound[-c(1,3,5,7,9),channels[6]]<-5.2
  
  lower_bound[c(1),channels[7]]<-2
  upper_bound[-c(1),channels[7]]<-2
  
  lower_bound[c(10,12),channels[8]]<-2.5
  upper_bound[-c(10,12),channels[8]]<-2.5
  
  lower_bound[c(11,12),channels[9]]<-0.15
  upper_bound[-c(11,12),channels[9]]<-0.15
  
  lower_bound[c(10,11,12),channels[10]]<-1.8
  upper_bound[-c(10,11,12),channels[10]]<-1.8
  
  lower_bound[13,channels[11]]<-0.5
  
  
  GLT<-sf::st_drop_geometry(geo_list_total)
  GLT[3:14]<-TEMP_transform$transform(GLT[3:14])
  for (i in cell_types){
    uB<-paste(paste0('GLT[,"',channels,'"]','<',upper_bound[i,]),collapse = ' & ')
    lB<-paste(paste0('GLT[,"',channels,'"]','>',lower_bound[i,]),collapse = ' & ')
    logiEXP<-paste0('(',uB,') & (',lB,')')
    geo_list_total[eval(parse(text=logiEXP)),'newCell']<-i
    GLT[eval(parse(text=logiEXP)),'newCell']<-i
  }
 
  
  
  
  area<-sf::st_area(geo_list_total)
  area<-cbind.data.frame(area,sf::st_drop_geometry(geo_list_total)[,c('uid','cell','newCell','Method')])
  areaSUB<-aggregate(area~Method+uid+newCell,area,sum)
  areaSUB$area<-areaSUB$area/1000
  areaSUB$Method<-factor(areaSUB$Method,levels = c('BO','RU','GR'))
  
  
  
  ggp<-ggplot2::ggplot(areaSUB)+
    ggplot2::geom_boxplot(ggplot2::aes(x=Method,y=area,fill=Method),
                          position = ggplot2::position_dodge2(preserve = "single"),notch=F)+
    ggplot2::facet_wrap(ggplot2::vars(newCell),scales = 'free',nrow = 2)+
    ggplot2::scale_fill_manual(values = c('gray100','gray50','gray75'))+
    ggplot2::scale_size_manual(values=c(1,1,1))+
    ggplot2::labs(y = 'Area (x1000 pixels)')+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank())
  
  ggplot2::ggsave(plot = ggp,
                  filename =  file.path(rootFolder,rpttn,targetFolder,"Tot_area.eps"),
                  device = 'eps',
                  width = 14,
                  height = 3,
                  dpi = 600)
  
  write.csv(areaSUB,
            file =  file.path(rootFolder,rpttn,targetFolder,"Tot_area.csv"),
            row.names = F)
  
  
  counts<-sf::st_drop_geometry(geo_list_total)[,c('uid','cell','newCell','Method')]
  counts$newCell<-factor(counts$newCell)
  counts<-aggregate(newCell~uid+Method,counts,table)
  TEMP_norm<-as.data.frame(counts$newCell)
 
  counts<-cbind(counts[,c(1,2)],as.data.frame(TEMP_norm))
  counts<-tidyr::pivot_longer(counts,
                             names_to = 'newCell',
                             values_to='Counts',
                             cols=colnames(counts)[3:16])
  
  counts$Method<-factor(counts$Method,levels = c('BO','RU','GR'))
  
  
  
  ggp<-ggplot2::ggplot(counts)+
    ggplot2::geom_boxplot(ggplot2::aes(x=Method,y=Counts,fill=Method),
                          position = ggplot2::position_dodge2(preserve = "single"),notch=F)+
    ggplot2::facet_wrap(ggplot2::vars(newCell),scales = 'free',nrow = 2)+
    ggplot2::scale_fill_manual(values = c('gray100','gray50','gray75'))+
    ggplot2::scale_size_manual(values=c(1,1,1))+
    ggplot2::labs(y = 'Counts')+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank())
  
  ggplot2::ggsave(plot = ggp,
                  filename =  file.path(rootFolder,rpttn,targetFolder,"Tot_counts.eps"),
                  device = 'eps',
                  width = 14,
                  height = 3,
                  dpi = 600)
  
  write.csv(counts,
            file =  file.path(rootFolder,rpttn,targetFolder,"Tot_counts.csv"),
            row.names = F)
  
  
  
  newColors<-c('cyan','green1','green2','green3','green4','orange1','orange2','orange3','orange4','purple1','purple2','purple3','gray80','red')
  names(newColors)<-sort(unique(geo_list_total$newCell))
  postscript(file = file.path(rootFolder,rpttn,targetFolder,"GR.eps"),
             width = 40,
             height = 40,
             bg = 'white',
             horizontal = F,
             pointsize = 10)
  par(oma=c(0,0,0,0))
  plot(geo_list_total[geo_list_total$uid==st_uids(mystudy)[1] & geo_list_total$Method=='GR',]['newCell'],main='GR',key.pos = 2,key.length = lcm(8),
       col=newColors[sf::st_drop_geometry(geo_list_total)[geo_list_total$uid==st_uids(mystudy)[1] & geo_list_total$Method=='GR','newCell']])
  dev.off()
  
  postscript(file = file.path(rootFolder,rpttn,targetFolder,"BO.eps"),
             width = 40,
             height = 40,
             bg = 'white',
             horizontal = F,
             pointsize = 10)
  par(oma=c(0,0,0,0))
  plot(geo_list_total[geo_list_total$uid==st_uids(mystudy)[1] & geo_list_total$Method=='BO',]['newCell'],main='BO',key.pos = 2,key.length = lcm(8),
       col=newColors[sf::st_drop_geometry(geo_list_total)[geo_list_total$uid==st_uids(mystudy)[1] & geo_list_total$Method=='BO','newCell']])
  dev.off()
  
  postscript(file = file.path(rootFolder,rpttn,targetFolder,"RU.eps"),
             width = 40,
             height = 40,
             bg = 'white',
             horizontal = F,
             pointsize = 10)
  par(oma=c(0,0,0,0))
  plot(geo_list_total[geo_list_total$uid==st_uids(mystudy)[1] & geo_list_total$Method=='RU',]['newCell'],main='RU',key.pos = 2,key.length = lcm(8),
       col=newColors[sf::st_drop_geometry(geo_list_total)[geo_list_total$uid==st_uids(mystudy)[1] & geo_list_total$Method=='RU','newCell']])
  dev.off()
  
  
  postscript(file =  file.path(rootFolder,rpttn,targetFolder,"colors.eps"),
             width = 40,
             height = 40,
             bg = 'white',
             horizontal = F,
             pointsize = 10)
  plot(rep(1,length(newColors)),pch=15,col=newColors,cex=4)
  dev.off()
  
  
  