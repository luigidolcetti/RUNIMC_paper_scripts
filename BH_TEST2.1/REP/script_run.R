library(RUNIMC)

mystudy<-initStudy(fn_studyName = 'TEST_ROUND1',
                   fn_rootFolder = "C:/Users/luigi/Documents/BH_TEST/BH_TEST2.1/REP_1",
                   fn_rawDataFolder = "C:/Users/luigi/Documents/BH_TEST/BH_TEST2.1/REP_1/XYZfiles"  )

newAnalysis(x = mystudy,
            analysisName = 'RUNIMC1')
# runIMC()

addFilter(x = mystudy,
          filter = 'vanvliet',
          parameters = list(list(sigma=1,order=0,axis='x'),
                            list(sigma=3,order=0,axis='x'),
                            list(sigma=1,order=1,axis='x'),
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
                                                    PvalueTreshold = 0.25,
                                                    ntree=5000,
                                                    mtry=12,
                                                    importance=F,
                                                    nodesize=1,
                                                    do.trace=500))

makeClassificationModel(mystudy,
                        method = 'randomForest')
classify(mystudy,saveToDisk = T,method = 'randomForest')

localCorrection(mystudy,
                matrixExtent = 3)

addClassificationDirectives(x = mystudy,
                            method = 'randomOnions',
                            methodParameters = list(responseVariable = 'SLI',
                                                    classificationLyr = 'label_clean',
                                                    predictiveFeatures = availableTF,
                                                    prefix = 'sko_',
                                                    labels = tf_labelList(mystudy)[1:3],
                                                    importance = F,
                                                    mtry = 50,
                                                    ntree = 5000,
                                                    do.trace=500))



makeClassificationModel(x=mystudy,method = 'randomOnions')


archive(mystudy)

classify(mystudy,saveToDisk = T,method = 'randomOnions')

addSegmentationDirectives(mystudy,method = 'pandaMap')

mystudy$currentAnalysis$segmentationDirectives

mystudy$currentAnalysis$segmentationDirectives@methodParameters$lowerQuantile = 0
mystudy$currentAnalysis$segmentationDirectives@methodParameters$upperQuantile = 1
mystudy$currentAnalysis$segmentationDirectives@methodParameters$lowerAreaLimit = 4
mystudy$currentAnalysis$segmentationDirectives@methodParameters$numberOfCores = 24
mystudy$currentAnalysis$segmentationDirectives@methodParameters$ClampDetectionDirection = '4'
mystudy$currentAnalysis$segmentationDirectives@methodParameters$overlapExtent = 40:60
mystudy$currentAnalysis$segmentationDirectives@methodParameters$movingWindowDimension = 100:200
mystudy$currentAnalysis$segmentationDirectives@methodParameters$nOfCutBrakes = 7

archive(mystudy)

segment(mystudy,
        labelLayer = names(mystudy$currentAnalysis$classification[[1]])[c(7:9)])

archive(mystudy)

topLayers<-lapply(st_uids(mystudy),function(uids){
  out<-RUNIMC:::.topLayer(x = mystudy$currentAnalysis$exprs,
                          uid = uids)})

topLayers<-do.call(dplyr::bind_rows,topLayers)


sf::write_sf(topLayers,
             "C:/Users/luigi/Documents/BH_TEST/BH_TEST2.1/REP_1/TEST_ROUND1/analysis/RUNIMC1/topLayers.sqlite")


####################### BODEN_noShrink ############

library(RUNIMC)
memory.limit(32000)

# mystudy<-retrieve("C:/Users/luigi/Documents/BH_TEST/BH_TEST2.1/REP_1/TEST_ROUND1/study.xml" )

newAnalysis(x = mystudy,
            analysisName = 'BODEN1')

importTiffMask(x = mystudy,
               fileFolder = "C:/Users/luigi/Documents/BH_TEST/BH_TEST2.1/REP_1/cellprofiler/cpTotImg"  ,
               layerName = 'bo_',
)

addSegmentationDirectives(mystudy,method = 'lazyCatMap')

mystudy$currentAnalysis$segmentationDirectives@methodParameters$indexToExclude<-0

segment(mystudy,labelLayer = 'bo__BDN')



archive(mystudy)



