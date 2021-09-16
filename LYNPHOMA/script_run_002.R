library(RUNIMC)

mystudy<-initStudy(fn_studyName = 'RUN_002',
                   fn_rootFolder = "C:/Users/luigi/Documents/BH_TEST/MYELOMA"  ,
                   fn_rawDataFolder = "C:/Users/luigi/Documents/BH_TEST/MYELOMA/RAWDATA"  )


exportTiffImage(mystudy,
                what = 'raw',
                fileFolder = "C:/Users/luigi/Documents/BH_TEST/MYELOMA/TIFF_002")

newAnalysis(x = mystudy,
            analysisName = 'RUNIMC1')

runIMC_V2(mystudy)

addFilter(x = mystudy,
          filter = 'vanvliet',
          parameters = list(list(sigma=1,order=0,axis='x'),
                            list(sigma=3,order=0,axis='x'),
                            list(sigma=1,order=1,axis='x'),
                            list(sigma=3,order=1,axis='x'),
                            list(sigma=1,order=2,axis='x'),
                            list(sigma=3,order=2,axis='x'),
                            list(sigma=1,order=0,axis='y'),
                            list(sigma=3,order=0,axis='y'),
                            list(sigma=1,order=1,axis='y'),
                            list(sigma=3,order=1,axis='y'),
                            list(sigma=1,order=2,axis='y'),
                            list(sigma=3,order=2,axis='y')),
          channels = ch_Rnames(mystudy),
          append = F)

addFilter(x = mystudy,
          filter = 'deriche',
          parameters = list(list(sigma=1,order=0,axis='x'),
                            list(sigma=3,order=0,axis='x'),
                            list(sigma=1,order=1,axis='x'),
                            list(sigma=3,order=1,axis='x'),
                            list(sigma=1,order=2,axis='x'),
                            list(sigma=3,order=2,axis='x'),
                            list(sigma=1,order=0,axis='y'),
                            list(sigma=3,order=0,axis='y'),
                            list(sigma=1,order=1,axis='y'),
                            list(sigma=3,order=1,axis='y'),
                            list(sigma=1,order=2,axis='y'),
                            list(sigma=3,order=2,axis='y')),
          channels = ch_Rnames(mystudy),
          append = T)

addFilter(mystudy,'blur_anisotropic',list(list(amplitude=3),list(amplitude=1)),channels = ch_Rnames(mystudy),append=T)

deployFilters(mystudy)

addExtractionDirectives(mystudy,c(1,0.25),'core',append=F)

extractTrainingFeatures(mystudy)


# mystudy<-retrieve("C:/Users/k1343421/Documents/BH/RUN/TEST_ROUND/study.xml"  )
availableTF<-tf_featureList(mystudy)

addClassificationDirectives(x = mystudy,
                            method = 'randomForest',
                            methodParameters = list(responseVariable = 'label',
                                                    predictiveFeatures = availableTF,
                                                    PvalueTreshold = 0.2,
                                                    ntree=5000,
                                                    mtry=21,
                                                    importance=F,
                                                    nodesize=1,
                                                    do.trace=1000))

makeClassificationModel(mystudy,
                        method = 'randomForest')

classify(mystudy,saveToDisk = T,method = 'randomForest')

localCorrection(mystudy,
                matrixExtent = 3) ### runned but using just label

archive(mystudy)

##### first failed
mystudy<-retrieve(
  "C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/study.xml"
)

addClassificationDirectives(x = mystudy,
                            method = 'randomOnions',
                            methodParameters = list(responseVariable = 'SLI',
                                                    classificationLyr = 'label',
                                                    predictiveFeatures = availableTF,
                                                    prefix = 'sko_',
                                                    labels = tf_labelList(mystudy),
                                                    importance = F,
                                                    mtry = 130,
                                                    ntree = 3000,
                                                    do.trace=1000))

makeClassificationModel(x=mystudy,method = 'randomOnions')

classify(mystudy,saveToDisk = T,method = 'randomOnions')

addSegmentationDirectives(mystudy,method = 'pandaMap')

mystudy$currentAnalysis$segmentationDirectives

mystudy$currentAnalysis$segmentationDirectives@methodParameters$lowerQuantile = 0
mystudy$currentAnalysis$segmentationDirectives@methodParameters$upperQuantile = 1
mystudy$currentAnalysis$segmentationDirectives@methodParameters$lowerAreaLimit = 5
mystudy$currentAnalysis$segmentationDirectives@methodParameters$numberOfCores = 24
mystudy$currentAnalysis$segmentationDirectives@methodParameters$ClampDetectionDirection = '4'
mystudy$currentAnalysis$segmentationDirectives@methodParameters$overlapExtent = 40:60
mystudy$currentAnalysis$segmentationDirectives@methodParameters$movingWindowDimension = 100:200
mystudy$currentAnalysis$segmentationDirectives@methodParameters$nOfCutBrakes = 10
mystudy$currentAnalysis$segmentationDirectives@methodParameters$numberOfCores = 24

archive(mystudy)

segment(mystudy,
        labelLayer = names(mystudy$currentAnalysis$classification[[1]])[c(11:18)])

archive(mystudy)

topLayers<-lapply(st_uids(mystudy),function(uids){
  out<-RUNIMC:::.topLayer(x = mystudy$currentAnalysis$exprs,
                          uid = uids)})

topLayers<-do.call(dplyr::bind_rows,topLayers)


sf::write_sf(topLayers,
             "C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/analysis/RUNIMC1/topLayers.sqlite",
             append = F,
             delete_layer = T)


####################### BODEN ############


newAnalysis(x = mystudy,
            analysisName = 'BODEN')

importTiffMask(x = mystudy,
               fileFolder = "C:/Users/luigi/Documents/BH_TEST/MYELOMA/TIFF_002/CP/cpToImg"  ,
               layerName = 'bo',
)

addSegmentationDirectives(mystudy,method = 'lazyCatMap')

mystudy$currentAnalysis$segmentationDirectives@methodParameters$indexToExclude<-0

segment(mystudy,labelLayer = names(mystudy$currentAnalysis$classification[[1]]))


archive(mystudy)

