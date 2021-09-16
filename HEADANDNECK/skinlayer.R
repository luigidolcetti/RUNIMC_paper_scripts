skl<-function(x){
  workMatrix<-x
  neigbours<-matrix(c(0,1,0,-1,1,0,-1,0),ncol = 2,byrow = T)
  skinLayer<-vector(mode = 'numeric',length = nrow(workMatrix))
  skinLayerIndex<-0
  flag<-nrow(workMatrix)
  while(flag>0){
    testPixel<-apply(workMatrix,1,function(rr){
      scorePixel<-apply(neigbours,1,function(ngbrs){
        newNgbrs<-rr+ngbrs
        NgbrsCompare<-(newNgbrs[1]==workMatrix[,1]) & (newNgbrs[2]==workMatrix[,2])
        if (any(NgbrsCompare)) return(1) else return(0)
      })
      
      scorePixel<-sum(scorePixel)
      if (scorePixel<4) return(T) else return(F)
    })
    
    skinLayer[testPixel]<-skinLayerIndex
    workMatrix[testPixel,]<-c(Inf,Inf)
    flag<-flag-sum(testPixel)
    skinLayerIndex<-skinLayerIndex-1
  }
  skinLayer<-skinLayer-(min(skinLayer))
  return(skinLayer)
}
