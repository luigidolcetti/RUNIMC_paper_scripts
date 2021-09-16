RU_TEMP<-table(geo_list_RU_annot$SOM[geo_list_RU_annot$uid==uids[5],drop=T])
BO_TEMP<-table(geo_list_BO_annot$SOM[geo_list_BO_annot$uid==uids[5],drop=T])
RU_TEMP<-RU_TEMP/sum(RU_TEMP)
BO_TEMP<-BO_TEMP/sum(BO_TEMP)
-sum(RU_TEMP*log(RU_TEMP))
-sum(BO_TEMP*log(BO_TEMP))
a=0
RU<-vapply(uids,function(x){
  RU_TEMP<-table(geo_list_RU_annot$SOM[geo_list_RU_annot$uid==x,drop=T])
  RU_TEMP<-RU_TEMP/sum(RU_TEMP)
  if (a==1) { out<-exp(-sum(RU_TEMP*log(RU_TEMP)))} else {
    out<-sum(RU_TEMP^a)^(1/(1-a))
  }
},c(1))


BO<-vapply(uids,function(x){
  BO_TEMP<-table(geo_list_BO_annot$SOM[geo_list_BO_annot$uid==x,drop=T])
  BO_TEMP<-BO_TEMP/sum(BO_TEMP)
  if (a==1) { out<-exp(-sum(BO_TEMP*log(BO_TEMP)))} else {
    out<-sum(BO_TEMP^a)^(1/(1-a))
  }
},c(1))

mean(RU)
mean(BO)
  wilcox.test(RU,BO,paired = T)
  t.test(RU,BO,paired = T)
boxplot(cbind(RU,BO))
shapiro.test(BO)


cbind(
table(geo_list_RU_annot$SOM[geo_list_RU_annot$uid=='uid.031837e3a256dfc52113ad3845cf19bc',drop=T]),
table(geo_list_BO_annot$SOM[geo_list_BO_annot$uid=='uid.031837e3a256dfc52113ad3845cf19bc',drop=T])
)

