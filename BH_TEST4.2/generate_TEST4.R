library(barbieHistologist)

m_CD45<-bh_defineMarker(markerName = 'cd45',
                        channelName = 'ch01',
                        patternFunction = '.pattern.random',
                        patternModifier = list(mean = 3, sd=0.5),
                        compartment = 'cytoplasm')

m_CD3<-bh_defineMarker(markerName = 'cd3',
                       channelName = 'ch02',
                       patternFunction = '.pattern.random',
                       patternModifier = list(mean = 3, sd=0.7),
                       compartment = 'cytoplasm')

m_CD4<-bh_defineMarker(markerName = 'cd4',
                       channelName = 'ch03',
                       patternFunction = '.pattern.skin2',
                       patternModifier = list(mean = 4, sd=0.3,k=0.5,X0=0.5),
                       compartment = 'cytoplasm')

m_CD8<-bh_defineMarker(markerName = 'cd8',
                       channelName = 'ch04',
                       patternFunction = '.pattern.skin2',
                       patternModifier = list(mean = 4, sd=0.3,k=1,X0=0.3 ),
                       compartment = 'cytoplasm')

m_CCR7_hi<-bh_defineMarker(markerName = 'ccr7_hi',
                           channelName = 'ch05',
                           patternFunction = '.pattern.random',
                           patternModifier = list(mean = 3, sd=0.3),
                           compartment = 'cytoplasm')

m_CCR7_low<-bh_defineMarker(markerName = 'ccr7_low',
                            channelName = 'ch05',
                            patternFunction = '.pattern.random',
                            patternModifier = list(mean = 1, sd=0.3),
                            compartment = 'cytoplasm')

m_CD45RO_hi<-bh_defineMarker(markerName = 'cd45ro_hi',
                             channelName = 'ch06',
                             patternFunction = '.pattern.random',
                             patternModifier = list(mean = 6, sd=0.3),
                             compartment = 'cytoplasm')

m_CD45RO_low<-bh_defineMarker(markerName = 'cd45ro_low',
                              channelName = 'ch06',
                              patternFunction = '.pattern.random',
                              patternModifier = list(mean = 3, sd=0.3),
                              compartment = 'cytoplasm')

m_CD45RO_sup<-bh_defineMarker(markerName = 'cd45ro_sup',
                              channelName = 'ch06',
                              patternFunction = '.pattern.random',
                              patternModifier = list(mean = 7, sd=0.1),
                              compartment = 'cytoplasm')

m_CD20<-bh_defineMarker(markerName = 'cd20',
                        channelName = 'ch07',
                        patternFunction = '.pattern.skin2',
                        patternModifier = list(mean = 10, sd=0.1,k=0.5,X0=0.2 ),
                        compartment = 'cytoplasm')

m_CD68<-bh_defineMarker(markerName = 'cd68',
                        channelName = 'ch08',
                        patternFunction = '.pattern.random',
                        patternModifier = list(mean = 4, sd=0.3),
                        compartment = 'cytoplasm')

m_CD33<-bh_defineMarker(markerName = 'cd33',
                        channelName = 'ch09',
                        patternFunction = '.pattern.skin2',
                        patternModifier = list(mean = 1, sd=0.3,k=0.1,X0=0.3 ),
                        compartment = 'cytoplasm')

m_CD74<-bh_defineMarker(markerName = 'cd74',
                        channelName = 'ch10',
                        patternFunction = '.pattern.random',
                        patternModifier = list(mean = 3, sd=0.3),
                        compartment = 'cytoplasm')

m_ECAD<-bh_defineMarker(markerName = 'ecad',
                        channelName = 'ch11',
                        patternFunction = '.pattern.skin2',
                        patternModifier = list(mean = 5, sd=0.1,k=1.5,X0=0.3 ),
                        compartment = 'cytoplasm')

DNA1<-bh_defineMarker(markerName = 'dna1',
                      channelName = 'ch12',
                      patternFunction = '.pattern.random',
                      patternModifier = list(mean = 5, sd=0.3),
                      compartment = 'nucleus')

DNA2<-bh_defineMarker(markerName = 'dna2',
                      channelName = 'ch12',
                      patternFunction = '.pattern.random',
                      patternModifier = list(mean = 6, sd=0.3),
                      compartment = 'nucleus')


cytop_imm<-bh_defineShape(majorAxis = c(2,0.1),
                          minorAxis = c(2,0.1),
                          roundness = c(1,0),
                          nArms = 6,
                          fixedArms = T,
                          orientation = NULL,
                          armExt = c(1.1,0.1),
                          armElbow = 3,
                          armSwing = 0,
                          armTrig = c(-1,0.1))

cytop_tum<-bh_defineShape(majorAxis = c(5,0.5),
                          minorAxis = c(5,0.5),
                          roundness = c(0.4,0),
                          nArms = 5,
                          fixedArms = T,
                          orientation = 90,
                          armExt = c(1.1,0.1),
                          armElbow = 3,
                          armSwing = 0,
                          armTrig = c(-2,0.1))

nuc_imm<-bh_defineShape(majorAxis = c(2,0.05),
                        minorAxis = c(1,0.05),
                        roundness = c(1,0),
                        nArms = 6,
                        fixedArms = T,
                        orientation = NULL,
                        armExt = c(1,0),
                        armElbow = 2,
                        armSwing = 0,
                        armTrig = c(-1,0.1))

nuc_tum<-bh_defineShape(majorAxis = c(3,0.01),
                        minorAxis = c(2,0.01),
                        roundness = c(1,0),
                        nArms = 6,
                        fixedArms = T,
                        orientation = NULL,
                        armExt = c(1,0),
                        armElbow = 2,
                        armSwing = 0,
                        armTrig = c(-1,0.1))




cell1<-bh_defineCell(name = 'CD4_UL',
                     cytoplasm = cytop_imm,
                     nucleus = nuc_imm,
                     organelle = NULL,
                     markers = list(m_CD45,
                                    m_CD3,
                                    m_CD4,
                                    m_CD45RO_low,
                                    m_CCR7_hi,
                                    DNA1))
cell2<-bh_defineCell(name = 'CD4_UR',
                     cytoplasm = cytop_imm,
                     nucleus = nuc_imm,
                     organelle = NULL,
                     markers = list(m_CD45,
                                    m_CD3,
                                    m_CD4,
                                    m_CD45RO_hi,
                                    m_CCR7_hi,
                                    DNA1))
cell3<-bh_defineCell(name = 'CD4_LR',
                     cytoplasm = cytop_imm,
                     nucleus = nuc_imm,
                     organelle = NULL,
                     markers = list(m_CD45,
                                    m_CD3,
                                    m_CD4,
                                    m_CD45RO_hi,
                                    m_CCR7_low,
                                    DNA1))
cell4<-bh_defineCell(name = 'CD4_LL',
                     cytoplasm = cytop_imm,
                     nucleus = nuc_imm,
                     organelle = NULL,
                     markers = list(m_CD45,
                                    m_CD3,
                                    m_CD4,
                                    m_CD45RO_low,
                                    m_CCR7_low,
                                    DNA1))

########
cell5<-bh_defineCell(name = 'CD8_UL',
                     cytoplasm = cytop_imm,
                     nucleus = nuc_imm,
                     organelle = NULL,
                     markers = list(m_CD45,
                                    m_CD3,
                                    m_CD8,
                                    m_CD45RO_low,
                                    m_CCR7_hi,
                                    DNA1))
cell6<-bh_defineCell(name = 'CD8_UR',
                     cytoplasm = cytop_imm,
                     nucleus = nuc_imm,
                     organelle = NULL,
                     markers = list(m_CD45,
                                    m_CD3,
                                    m_CD8,
                                    m_CD45RO_hi,
                                    m_CCR7_hi,
                                    DNA1))
cell7<-bh_defineCell(name = 'CD8_LR',
                     cytoplasm = cytop_imm,
                     nucleus = nuc_imm,
                     organelle = NULL,
                     markers = list(m_CD45,
                                    m_CD3,
                                    m_CD8,
                                    m_CD45RO_hi,
                                    m_CCR7_low,
                                    DNA1))
cell8<-bh_defineCell(name = 'CD8_LL',
                     cytoplasm = cytop_imm,
                     nucleus = nuc_imm,
                     organelle = NULL,
                     markers = list(m_CD45,
                                    m_CD3,
                                    m_CD8,
                                    m_CD45RO_low,
                                    m_CCR7_low,
                                    DNA1))

####################

cell9<-bh_defineCell(name = 'CD20',
                     cytoplasm = cytop_imm,
                     nucleus = nuc_imm,
                     organelle = NULL,
                     markers = list(m_CD45,
                                    m_CD20,
                                    m_CD45RO_sup,
                                    DNA1))

cell10<-bh_defineCell(name = 'MAC',
                      cytoplasm = cytop_imm,
                      nucleus = nuc_imm,
                      organelle = NULL,
                      markers = list(m_CD45,
                                     m_CD68,
                                     m_CD74,
                                     DNA1))

cell11<-bh_defineCell(name = 'MONO',
                      cytoplasm = cytop_imm,
                      nucleus = nuc_imm,
                      organelle = NULL,
                      markers = list(m_CD45,
                                     m_CD33,
                                     m_CD74,
                                     DNA1))

cell12<-bh_defineCell(name = 'MONOMAC',
                      cytoplasm = cytop_imm,
                      nucleus = nuc_imm,
                      organelle = NULL,
                      markers = list(m_CD45,
                                     m_CD68,
                                     m_CD33,
                                     m_CD74,
                                     DNA1))

cell13<-bh_defineCell(name = 'TUM',
                      cytoplasm = cytop_tum,
                      nucleus = nuc_tum,
                      organelle = NULL,
                      markers = list(m_ECAD,
                                     DNA2))

dirList<-list.dirs(
  "C:/Users/luigi/Documents/BH_TEST/BH_TEST4.2",recursive = F
)

for (ii in 1:3){
  
  densityList<-c(0.98,0.98,0.98,0.98)
  fileRoot= dirList[ii]
  maxdim = 500
  dir.create(fileRoot)
  
  for (dd in seq_along(densityList)){
    
    folderDensity<-file.path(fileRoot,paste0('sample_',dd))
    dir.create(folderDensity)
    
    tissue1<-bh_defineTissue(coords = c(0,maxdim,0,maxdim),
                             resolution = 1,
                             bg = 0,
                             markers = list(m_CD45,
                                            m_CD3,
                                            m_CD4,
                                            m_CD8,
                                            m_CCR7_hi,
                                            m_CD45RO_hi,
                                            m_CD20,
                                            m_CD68,
                                            m_CD33,
                                            m_CD74,
                                            m_ECAD,
                                            DNA1))
    
    TEMP_population<-bh_populate_byInteract(cellPrototype = list(cell1,
                                                                 cell2,
                                                                 cell3,
                                                                 cell4,
                                                                 cell5,
                                                                 cell6,
                                                                 cell7,
                                                                 cell8,
                                                                 cell9,
                                                                 cell10,
                                                                 cell11,
                                                                 cell12,
                                                                 cell13),
                                            proportion = c(0.01,0.05,0.01,0.03,
                                                           0.003,0.003,0.003,0.001,
                                                           0.1,
                                                           0.03,0.03,0.04,
                                                           0.69),
                                            tissue = tissue1,
                                            require_cytoplasm = T,
                                            require_nucleus = T,
                                            require_organelle = F,
                                            cropToMesure = F,
                                            areaTresh = densityList[dd])
    
    bh_savePopulation(TEMP_population,
                      file=file.path(folderDensity,paste0('population_',dd,'.R')))
    
    TEMP_pic<-bh_engrave(tissue = tissue1,
                         cells = TEMP_population$CD4_UL)
    TEMP_pic<-bh_engrave(tissue = TEMP_pic,
                         cells = TEMP_population$CD4_UR)
    TEMP_pic<-bh_engrave(tissue = TEMP_pic,
                         cells = TEMP_population$CD4_LR)
    TEMP_pic<-bh_engrave(tissue = TEMP_pic,
                         cells = TEMP_population$CD4_LL)
    
    TEMP_pic<-bh_engrave(tissue = TEMP_pic,
                         cells = TEMP_population$CD8_UL)
    TEMP_pic<-bh_engrave(tissue = TEMP_pic,
                         cells = TEMP_population$CD8_UR)
    TEMP_pic<-bh_engrave(tissue = TEMP_pic,
                         cells = TEMP_population$CD8_LR)
    TEMP_pic<-bh_engrave(tissue = TEMP_pic,
                         cells = TEMP_population$CD8_LL)
    
    TEMP_pic<-bh_engrave(tissue = TEMP_pic,
                         cells = TEMP_population$CD20)
    
    TEMP_pic<-bh_engrave(tissue = TEMP_pic,
                         cells = TEMP_population$MONO)
    TEMP_pic<-bh_engrave(tissue = TEMP_pic,
                         cells = TEMP_population$MAC)
    TEMP_pic<-bh_engrave(tissue = TEMP_pic,
                         cells = TEMP_population$MONOMAC)
    
    TEMP_pic<-bh_engrave(tissue = TEMP_pic,
                         cells = TEMP_population$TUM)
    
    TEMP_mod<-TEMP_pic
    # MM<-c(3,3,4,4,2,6,10,4,1,3,5,5.5)
    # SS<-c(0.5,0.7,0.5,0.3,0.3,0.3,0.1,0.3,0.3,0.3,0.3,0.3)
    for (i in seq_along(names(TEMP_mod))){
      V<-raster::values(TEMP_mod[[i]])
      M<-mean(V)
      S<-sd(V)
      TEMP_mod[[i]]<-bh_field_perturbation(raster = TEMP_mod[[i]],
                                           fun = .perturbation.constant,
                                           fun.param = list(field = c(0,maxdim,0,maxdim),
                                                            amount = 0.05,
                                                            mean = M,
                                                            sd = S))
    }
    
    TEMP_TRANS<-scales::modulus_trans(0)
    for (i in names(TEMP_mod)){
      TEMP_mod[[i]]<-TEMP_TRANS$inverse(TEMP_mod[[i]])
    }
    
    folderTiff<-file.path(folderDensity,paste0('TIFF_',dd))
    dir.create(folderTiff)
    
    for (i in names(TEMP_mod)){
      bh_saveTiff(TEMP_mod[[i]],
                  filePath = file.path(folderTiff,
                                       paste0(i,'.tif')))
    }
    
    XYZ_table<-bh_asXYZ(tissue = TEMP_mod)
    
    colnames(XYZ_table)<-c('X','Y',
                           'ch01-cd45(ch01)',
                           'ch02-cd3(ch02)',
                           'ch03-cd4(ch03)',
                           'ch04-cd8(ch04)',
                           'ch05-ccr7(ch05)',
                           'ch06-cd45ro(ch06)',
                           'ch07-cd20(ch07)',
                           'ch08-cd68(ch08)',
                           'ch09-cd33(ch09)',
                           'ch10-cd74(ch10)',
                           'ch11-ecad(ch11)',
                           'ch12-dna(ch12)')
    
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
