rootFolder<-"C:/Users/luigi/Documents/BH_TEST/BH_TEST5.1"
targetFolder<-'RESULTS_SUBSETTING'
rpttn<-c("REP_01","REP_02","REP_03")
trainingSampl<-0
brks<-c(1:4)
db<-factor(c('R1','R2','R3','R4'),levels=c('R1','R2','R3','R4'))

pathToRepeat<-file.path(rootFolder,rpttn,targetFolder)
area_list<-vector('list',length(rpttn))
count_list<-vector('list',length(rpttn))

for (ii in seq_along(rpttn)){
  area_list[[ii]]<-read.csv(file.path(pathToRepeat[[ii]],'Tot_area.csv'))
  area_list[[ii]]<-cbind.data.frame(area_list[[ii]],Rep=ii)
  count_list[[ii]]<-read.csv(file.path(pathToRepeat[[ii]],'Tot_counts.csv'))
  count_list[[ii]]<-cbind.data.frame(count_list[[ii]],Rep=ii)
}

area_tot<-do.call(rbind.data.frame,area_list)
counts_tot<-do.call(rbind.data.frame,count_list)

##### area ####
areaXrep<-aggregate(area~Rep+newCell+uid,area_tot,function(x){
  BO<-(x[1]-x[2])^2
  RU<-(x[3]-x[2])^2
  out<-c(BO,RU)
  names(out)<-c('BO','RU')
  out
})

areaXrep<-aggregate(area~uid+Rep,areaXrep,function(x)sqrt(sum(x)/length(x)))


areaT<-t.test(areaXrep$BO,areaXrep$RU,paired = T)
areaP<-gtools::stars.pval(areaT$p.value)

ymn<-min(areaXrep$BO,areaXrep$RU)
ymx<-max(areaXrep$BO,areaXrep$RU)
ypad<-(ymx-ymn)*0.25

postscript(file = file.path(rootFolder,'area_RSEM.eps'),
           onefile = T,
           bg = 'white',
           horizontal = F,
           height = 5,
           width = 5)
plot(x=c(rnorm(nrow(areaXrep),1,0.1),
         rnorm(nrow(areaXrep),2,0.1)),
     y=c(areaXrep$BO,
         areaXrep$RU),
     xlim = c(0,3),
     ylim = c(ymn-ypad,ymx+ypad),
     xlab = NA,
     xaxt = 'n',
     las = 1,
     bty = 'l',
     ylab = 'Area RSEM')

axis(1,c(1,2),c('Reference','Alternative'))
axis(3,c(1,2),c(NA,NA))
mtext(areaP,3,1)

points(x=c(1,2),colMeans(areaXrep[,c('BO','RU')]),pch=c(16,17))
lines(x=c(1,1),c(mean(areaXrep$BO)+sd(areaXrep$BO),mean(areaXrep$BO)-sd(areaXrep$BO)))
lines(x=c(0.8,1.2),c(mean(areaXrep$BO)+sd(areaXrep$BO),mean(areaXrep$BO)+sd(areaXrep$BO)))
lines(x=c(0.8,1.2),c(mean(areaXrep$BO)-sd(areaXrep$BO),mean(areaXrep$BO)-sd(areaXrep$BO)))

lines(x=c(2,2),c(mean(areaXrep$RU)+sd(areaXrep$RU),mean(areaXrep$RU)-sd(areaXrep$RU)))
lines(x=c(1.8,2.2),c(mean(areaXrep$RU)+sd(areaXrep$RU),mean(areaXrep$RU)+sd(areaXrep$RU)))
lines(x=c(1.8,2.2),c(mean(areaXrep$RU)-sd(areaXrep$RU),mean(areaXrep$RU)-sd(areaXrep$RU)))

dev.off()

##### counts ####
countsXrep<-aggregate(Counts~Rep+newCell+uid,counts_tot,function(x){
  BO<-(x[1]-x[2])^2
  RU<-(x[3]-x[2])^2
  out<-c(BO,RU)
  names(out)<-c('BO','RU')
  out
})

countsXrep<-aggregate(Counts~uid+Rep,countsXrep,function(x)sqrt(sum(x)/length(x)))

countsT<-t.test(countsXrep$BO,countsXrep$RU,paired = T)
countsP<-gtools::stars.pval(countsT$p.value)

ymn<-min(countsXrep$BO,countsXrep$RU)
ymx<-max(countsXrep$BO,countsXrep$RU)
ypad<-(ymx-ymn)*0.25

postscript(file = file.path(rootFolder,'counts_RSEM.eps'),
           onefile = T,
           bg = 'white',
           horizontal = F,
           height = 5,
           width = 5)
plot(x=c(rnorm(nrow(countsXrep),1,0.1),
         rnorm(nrow(countsXrep),2,0.1)),
     y=c(countsXrep$BO,
         countsXrep$RU),
     xlim = c(0,3),
     ylim = c(ymn-ypad,ymx+ypad),
     xlab = NA,
     xaxt = 'n',
     las = 1,
     bty = 'l',
     ylab = 'Counts RSEM')

axis(1,c(1,2),c('Reference','Alternative'),)
axis(3,c(1,2),c(NA,NA))
mtext(countsP,3,1)

points(x=c(1,2),colMeans(countsXrep[,c('BO','RU')]),pch=c(16,17))
lines(x=c(1,1),c(mean(countsXrep$BO)+sd(countsXrep$BO),mean(countsXrep$BO)-sd(countsXrep$BO)))
lines(x=c(0.8,1.2),c(mean(countsXrep$BO)+sd(countsXrep$BO),mean(countsXrep$BO)+sd(countsXrep$BO)))
lines(x=c(0.8,1.2),c(mean(countsXrep$BO)-sd(countsXrep$BO),mean(countsXrep$BO)-sd(countsXrep$BO)))

lines(x=c(2,2),c(mean(countsXrep$RU)+sd(countsXrep$RU),mean(countsXrep$RU)-sd(countsXrep$RU)))
lines(x=c(1.8,2.2),c(mean(countsXrep$RU)+sd(countsXrep$RU),mean(countsXrep$RU)+sd(countsXrep$RU)))
lines(x=c(1.8,2.2),c(mean(countsXrep$RU)-sd(countsXrep$RU),mean(countsXrep$RU)-sd(countsXrep$RU)))

dev.off()
