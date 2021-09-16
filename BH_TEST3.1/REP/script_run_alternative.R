library(RUNIMC)

exFile_study<-"C:/Users/luigi/Documents/BH_TEST/BH_TEST3.1/REP_1/TEST_ROUND1/study.xml"
exFile_analysis<-"C:/Users/luigi/Documents/BH_TEST/BH_TEST3.1/REP_1/TEST_ROUND1/analysis/RUNIMC1/analysis.xml"

mystudy<-retrieve(exFile_study)
mystudy$currentAnalysis<-RUNIMC:::retrieve.xml(exFile_analysis)

addSegmentationDirectives(mystudy,method = 'pandaMap')

mystudy$currentAnalysis$segmentationDirectives

mystudy$currentAnalysis$segmentationDirectives@methodParameters$lowerQuantile = 0
mystudy$currentAnalysis$segmentationDirectives@methodParameters$upperQuantile = 1
mystudy$currentAnalysis$segmentationDirectives@methodParameters$lowerAreaLimit = 16
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
             "C:/Users/luigi/Documents/BH_TEST/BH_TEST3.1/REP_1/TEST_ROUND1/analysis/RUNIMC1/topLayers.sqlite")


