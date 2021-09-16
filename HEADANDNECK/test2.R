library(SpadeR)
RU_TEMP<-aggregate(SOM~uid,geo_list_RU_annot,table)
rownames(RU_TEMP)<-RU_TEMP[,1]
RU_TEMP<-t(RU_TEMP[,-1])
# Diversity(data = RU_TEMP,
#           datatype = 'abundance')
RU_TEMP<-apply(RU_TEMP,2,function(x){
  TEMP<-x[x!=0]
  TEMP<-TEMP/sum(TEMP)
  -sum(TEMP*log(TEMP))
})

BO_TEMP<-aggregate(SOM~uid,geo_list_BO_annot,table)
rownames(BO_TEMP)<-BO_TEMP[,1]
BO_TEMP<-t(BO_TEMP[,-1])
# Diversity(data = BO_TEMP,
#           datatype = 'abundance')
BO_TEMP<-apply(BO_TEMP,2,function(x){
  TEMP<-x[x!=0]
  TEMP<-TEMP/sum(TEMP)
  -sum(TEMP*log(TEMP))
})

t.test(RU_TEMP,BO_TEMP,paired = T)

plot(NA,xlim=c(0.5,2.5),ylim=c(0.6,1))
points(rep(1,length(BO_TEMP)),BO_TEMP,pch=16)
points(rep(2,length(RU_TEMP)),RU_TEMP,pch=17)
lines(matrix(c(0.8,mean(BO_TEMP),1.2,mean(BO_TEMP)),
             ncol=2,
             byrow = T))
lines(matrix(c(1.8,mean(RU_TEMP),2.2,mean(RU_TEMP)),
             ncol=2,
             byrow = T))

for (i in 1:length(BO_TEMP)){
  lines(matrix(c(1,BO_TEMP[i],2,RU_TEMP[i]),
               ncol=2,
               byrow = T))
}
