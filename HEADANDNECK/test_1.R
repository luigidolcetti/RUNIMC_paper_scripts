RU_TEMP<-table(geo_list_RU_annot$SOM[geo_list_RU_annot$uid==uids[5],drop=T])
BO_TEMP<-table(geo_list_BO_annot$SOM[geo_list_BO_annot$uid==uids[5],drop=T])
RU_TEMP<-RU_TEMP/sum(RU_TEMP)
BO_TEMP<-BO_TEMP/sum(BO_TEMP)
-sum(RU_TEMP*log(RU_TEMP))
-sum(BO_TEMP*log(BO_TEMP))

RU<-vapply(uids,function(x){
  RU_TEMP<-table(geo_list_RU_annot$SOM[geo_list_RU_annot$uid==x,drop=T])
  RU_TEMP<-RU_TEMP/sum(RU_TEMP)
  -sum(RU_TEMP*log(RU_TEMP))
},c(1))


BO<-vapply(uids,function(x){
  BO_TEMP<-table(geo_list_BO_annot$SOM[geo_list_BO_annot$uid==x,drop=T])
  BO_TEMP<-BO_TEMP/sum(BO_TEMP)
  -sum(BO_TEMP*log(BO_TEMP))
},c(1))

mean(RU)
mean(BO)
t.test(RU,BO,paired = T)
