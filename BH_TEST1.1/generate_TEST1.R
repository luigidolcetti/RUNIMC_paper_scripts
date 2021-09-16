library(barbieHistologist)

CD45<-bh_defineMarker(markerName = 'cd45',
                      channelName = 'ch1',
                      patternFunction = '.pattern.random',
                      patternModifier = list(mean = 5, sd=5),
                      compartment = 'cytoplasm')
CD4<-bh_defineMarker(markerName = 'cd4',
                     channelName = 'ch2',
                     patternFunction = '.pattern.skin2',
                     patternModifier = list(mean = 2, sd=3,k=5,X0=0.2),
                     compartment = 'cytoplasm')
CD8<-bh_defineMarker(markerName = 'cd8',
                     channelName = 'ch3',
                     patternFunction = '.pattern.skin2',
                     patternModifier = list(mean = 4, sd=4,k=1,X0=0.5 ),
                     compartment = 'cytoplasm')
CD11c<-bh_defineMarker(markerName = 'cd11c',
                       channelName = 'ch4',
                       patternFunction = '.pattern.random',
                       patternModifier = list(mean = 4, sd=100),
                       compartment = 'cytoplasm')
CD1Low<-bh_defineMarker(markerName = 'cd1000',
                        channelName = 'ch5',
                        patternFunction = '.pattern.random',
                        patternModifier = list(mean = 2, sd=3),
                        compartment = 'cytoplasm')
CD1high<-bh_defineMarker(markerName = 'cd1000',
                         channelName = 'ch5',
                         patternFunction = '.pattern.random',
                         patternModifier = list(mean = 8, sd=4),
                         compartment = 'cytoplasm')
DNA<-bh_defineMarker(markerName = 'dna',
                     channelName = 'ch6',
                     patternFunction = '.pattern.random',
                     patternModifier = list(mean = 10, sd=5),
                     compartment = 'nucleus')


cytop1<-bh_defineShape(majorAxis = c(5,1),
                       minorAxis = c(5,1),
                       roundness = c(1,0),
                       nArms = 30,
                       fixedArms = T,
                       orientation =  NULL,
                       armExt = c(0.3,0.1),
                       armElbow = 2,
                       armSwing = 0,
                       armTrig = c(-2,0.1))

cytop2<-bh_defineShape(majorAxis = c(5,1),
                       minorAxis = c(5,1),
                       roundness = c(0.5,0),
                       nArms = 12,
                       fixedArms = F,
                       orientation = NULL,
                       armExt = c(0.3,0.1),
                       armElbow = 2,
                       armSwing = 0,
                       armTrig = c(-2,0.1))

nuc1<-bh_defineShape(majorAxis = c(3,0.5),
                     minorAxis = c(3,0.5),
                     roundness = c(1,0),
                     nArms = 40,
                     fixedArms = T,
                     orientation = NULL,
                     armExt = c(0.5,0),
                     armElbow = 2,
                     armSwing = 0,
                     armTrig = c(-1,0.1))

cell1<-bh_defineCell(name = 'CD4_lo',
                     cytoplasm = cytop1,
                     nucleus = nuc1,
                     organelle = NULL,
                     markers = list(CD45,CD4,DNA,CD1Low))
cell2<-bh_defineCell(name = 'CD4_hi',
                     cytoplasm = cytop1,
                     nucleus = nuc1,
                     organelle = NULL,
                     markers = list(CD45,CD4,DNA,CD1high))

cell3<-bh_defineCell(name = 'CD8_lo',
                     cytoplasm = cytop1,
                     nucleus = nuc1,
                     organelle = NULL,
                     markers = list(CD45,CD8,DNA,CD1Low))
cell4<-bh_defineCell(name = 'CD8_hi',
                     cytoplasm = cytop1,
                     nucleus = nuc1,
                     organelle = NULL,
                     markers = list(CD45,CD8,DNA,CD1high))

cell5<-bh_defineCell(name = 'dendr',
                     cytoplasm = cytop2,
                     nucleus = nuc1,
                     organelle = NULL,
                     markers = list(CD45,DNA,CD11c))


##### just a test ####

# maxdim = 100
# tissue1<-bh_defineTissue(coords = c(0,maxdim,0,maxdim),
#                          resolution = 1,
#                          bg = 0,
#                          markers = list(CD45,CD4,CD8,CD11c,CD1high,DNA))
# 
# TEMP_population<-bh_populate_byInteract(cellPrototype = list(cell1,cell2,cell3,cell4,cell5),
#                              proportion = c(0.2,0.2,0.2,0.2,0.2),
#                              tissue = tissue1,
#                              require_cytoplasm = T,
#                              require_nucleus = T,
#                              require_organelle = F,
#                              minDist_interact = 0,
#                              minDist_enclose = 0.5,
#                              areaTresh = 0.99,
#                              cropToMesure = T)
# 
# GEOM_list<-bh_asSFC(cells = TEMP_population)
# 
# plot(GEOM_list[1],col=NA,border=c('blue','red','green'))
# 
# TEMP<-sf::st_area(GEOM_list[GEOM_list$compartment=='cytoplasm',])
# hist(TEMP,breaks = 20)
# summary(TEMP)
# 
# TEMP_pic<-bh_engrave(tissue = tissue1,
#                      cells = TEMP_population$CD4_lo)
# TEMP_pic<-bh_engrave(tissue = TEMP_pic,
#                      cells = TEMP_population$CD4_hi)
# TEMP_pic<-bh_engrave(tissue = TEMP_pic,
#                      cells = TEMP_population$CD8_lo)
# TEMP_pic<-bh_engrave(tissue = TEMP_pic,
#                      cells = TEMP_population$CD8_hi)
# TEMP_pic<-bh_engrave(tissue = TEMP_pic,
#                      cells = TEMP_population$dendr)
# 
# TEMP_mod<-TEMP_pic
# 
# for (i in names(TEMP_mod)){
#   V<-raster::values(TEMP_mod[[i]])
#   M<-mean(V[V!=0])
#   S<-sd(V[V!=0])
#   TEMP_mod[[i]]<-bh_field_perturbation(raster = TEMP_mod[[i]],
#                                        fun = .perturbation.constant,
#                                        fun.param = list(field = c(0,maxdim,0,maxdim),
#                                                         amount = 0.1,
#                                                         mean = M,
#                                                         sd = S))
# }
#######

dirList<-list.dirs(
  "C:/Users/k1343421/Documents/BH2/BH_TEST1.2",recursive = F
)

for (ii in 1:3){
densityList<-c(0.05,0.1,0.25,0.5,0.5,0.75,0.99)

fileRoot= dirList[ii]
maxdim = 500
dir.create(fileRoot)

for (dd in seq_along(densityList)){
  
  folderDensity<-file.path(fileRoot,paste0('densityBrake_',dd))
  dir.create(folderDensity)
  
  tissue1<-bh_defineTissue(coords = c(0,maxdim,0,maxdim),
                           resolution = 1,
                           bg = 0,
                           markers = list(CD45,CD4,CD8,CD11c,CD1high,DNA))
  
  TEMP_population<-bh_populate_byInteract(cellPrototype = list(cell1,cell2,cell3,cell4,cell5),
                                          proportion = c(0.2,0.2,0.2,0.2,0.2),
                                          tissue = tissue1,
                                          require_cytoplasm = T,
                                          require_nucleus = T,
                                          require_organelle = F,
                                          areaTresh = densityList[dd],
                                          cropToMesure = T)
  
  bh_savePopulation(TEMP_population,
                    file=file.path(folderDensity,paste0('population_',dd,'.R')))
  
  TEMP_pic<-bh_engrave(tissue = tissue1,
                       cells = TEMP_population$CD4_lo)
  TEMP_pic<-bh_engrave(tissue = TEMP_pic,
                       cells = TEMP_population$CD4_hi)
  TEMP_pic<-bh_engrave(tissue = TEMP_pic,
                       cells = TEMP_population$CD8_lo)
  TEMP_pic<-bh_engrave(tissue = TEMP_pic,
                       cells = TEMP_population$CD8_hi)
  TEMP_pic<-bh_engrave(tissue = TEMP_pic,
                       cells = TEMP_population$dendr)
  
  TEMP_mod<-TEMP_pic
  
  for (i in names(TEMP_mod)){
    V<-raster::values(TEMP_mod[[i]])
    M<-mean(V[V!=0])
    S<-sd(V[V!=0])
    TEMP_mod[[i]]<-bh_field_perturbation(raster = TEMP_mod[[i]],
                                         fun = .perturbation.constant,
                                         fun.param = list(field = c(0,maxdim,0,maxdim),
                                                          amount = 0.05,
                                                          mean = M,
                                                          sd = S))
  }
  
  folderTiff<-file.path(folderDensity,paste0('TIFF_',dd))
  dir.create(folderTiff)
  
  for (i in names(TEMP_mod)){
    bh_saveTiff(TEMP_mod[[i]],
                filePath = file.path(folderTiff,
                                     paste0(i,'.tif')))
  }
  
  XYZ_table<-bh_asXYZ(tissue = TEMP_mod)
  
  colnames(XYZ_table)<-c('X','Y','ch1-CD45(ch1)','ch2-CD4(ch2)','ch3-CD8(ch3)','ch4-CD11c(ch4)','ch5-CDx(ch5)','ch6-DNA(ch6)')
  
  write.table(XYZ_table,
              file.path(folderDensity,paste0('XYZtable_',dd,'.txt')),
              row.names = F,
              quote = F,
              sep='\t')
  
  GEOM_list<-bh_asSFC(cells = TEMP_population)
  
  bh_saveSFC(GEOM_list,
             file.path(folderDensity,paste0('GroundTruth_',dd,'.sqlite')))
  
  FAM<-bh_familyPicture(tissue = tissue1,
                        sfc = GEOM_list,
                        compartment = 'cytoplasm')
  
  bh_saveFamily(familyPicture = FAM,
                familyName = paste0('FAM_',dd),
                filePath = folderDensity)
  
}

}
