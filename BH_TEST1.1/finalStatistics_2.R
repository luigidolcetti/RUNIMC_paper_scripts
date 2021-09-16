library(barbieHistologist)

rootFolder<-"C:/Users/luigi/Documents/BH_TEST/BH_TEST1.1/"
targetFolder<-'RESULTS_FAMILY'
brakeName<-'densityBrake_'
dirRes<-"C:/Users/luigi/Documents/BH_TEST/BH_TEST1.1/Final_stat"
rpttn<-c("REP_1","REP_2","REP_3")
trainingSampl<-4
brks<-c(1:3,5:7)

db<-factor(c('D1','D2','D3','x','D4','D5','D6'),levels=c('D1','D2','D3','x','D4','D5','D6'))

JD<-lapply(rpttn,function(rep){

geo_list_BO<-sf::st_read(
  file.path(rootFolder,rep,"TEST_ROUND1/analysis/BODEN1/ex.sqlite")
)
if (any(!sf::st_is_valid(geo_list_BO))){
  geo_list_BO<-sf::st_make_valid(geo_list_BO)
}

geo_list_BO<-geo_list_BO[,c('uid','id','GEOMETRY')]

geo_list_RU<-sf::st_read(
  file.path(rootFolder,rep,"TEST_ROUND1/analysis/RUNIMC1/topLayers.sqlite")
)
if (any(!sf::st_is_valid(geo_list_RU))){
  geo_list_RU<-sf::st_make_valid(geo_list_RU)
}

geo_list_RU<-geo_list_RU[,c('uid','splitp_id','GEOMETRY')]
colnames(geo_list_RU)[2]<-'id'

uids<-unique(geo_list_RU$uid)
geo_list_BO<-geo_list_BO[geo_list_BO$uid %in% uids[brks],]
geo_list_RU<-geo_list_RU[geo_list_RU$uid %in% uids[brks],]

geo_list_GR<-lapply(brks,function(densityBrake){
  
  familyRaster<-bh_loadFamily(file.path(rootFolder,rep,paste0(brakeName,densityBrake,"/FAM_",densityBrake)))
  
  geo_list_GR_sub<-RUNIMC::lazyCatMap(familyRaster$ID,fn_indexToExclude = 0,fn_uid = uids[densityBrake])
  colnames(geo_list_GR_sub)<-c('uid','id','GEOMETRY')
  sf::st_geometry(geo_list_GR_sub)<-'GEOMETRY'
  geo_list_GR_sub<-dplyr::bind_cols(geo_list_GR_sub,smp=densityBrake)
})
geo_list_GR<-do.call(dplyr::bind_rows,geo_list_GR)

area.GR<-cbind.data.frame(uid=sf::st_drop_geometry(geo_list_GR)[,'uid'],area=sf::st_area(geo_list_GR),Method='GR')
area.RU<-cbind.data.frame(uid=sf::st_drop_geometry(geo_list_RU)[,'uid'],area=sf::st_area(geo_list_RU),Method='RU')
area.BO<-cbind.data.frame(uid=sf::st_drop_geometry(geo_list_BO)[,'uid'],area=sf::st_area(geo_list_BO),Method='BO')

area.tot<-rbind.data.frame(area.GR,area.BO,area.RU)

jd<-lapply(unique(area.tot$uid),function(uids){
  subMat<-area.tot[area.tot$uid==uids,]
  TEMP_MAX<-max(subMat$area)+10
  BO<-hist(subMat$area[subMat$Method=='BO'],breaks = seq(0,TEMP_MAX,10),plot=F)$counts
  GR<-hist(subMat$area[subMat$Method=='GR'],breaks = seq(0,TEMP_MAX,10),plot=F)$counts
  RU<-hist(subMat$area[subMat$Method=='RU'],breaks = seq(0,TEMP_MAX,10),plot=F)$counts
  BO.jd<-philentropy::JSD(rbind(BO,GR),est.prob = 'empirical')
  RU.jd<-philentropy::JSD(rbind(RU,GR),est.prob = 'empirical')
  out<-data.frame(rep=rep,x=(which(unique(area.tot$uid)==uids)),BO=BO.jd,RU=RU.jd)
  return(out)
})

jd<-do.call(rbind.data.frame,jd)

})


JD<-do.call(rbind.data.frame,JD)

BO.mean<-aggregate(BO~x,JD,mean)
BO.sd<-aggregate(BO~x,JD,sd)
RU.mean<-aggregate(RU~x,JD,mean)
RU.sd<-aggregate(RU~x,JD,sd)


ymn<-min(c(JD$BO,JD$RU))
ymx<-max(c(JD$BO,JD$RU))
ypad<-(ymx-ymn)*0.25

postscript(file = file.path(dirRes,'JSD.eps'),
           onefile = T,
           bg = 'white',
           horizontal = F,
           height = 5,
           width = 5)
plot(NA,
     xlim = c(0,7),
     ylim = c(ymn-ypad,ymx+ypad),
     xlab = 'Cell density',
     xaxt = 'n',
     las = 1,
     bty = 'l',
     ylab = 'Jensen-Shannon divergence')

axis(1,1:length(brks),db[brks])
# axis(3,c(1,2),c(NA,NA))
# mtext(p.entropy,3,1)

points(BO.mean$x-0.1,BO.mean$BO,pch=16)
points(RU.mean$x+0.1,RU.mean$RU,pch=17)
lines(BO.mean$x-0.1,BO.mean$BO)
lines(RU.mean$x+0.1,RU.mean$RU)

for(i in 1:nrow(BO.mean)){
lines(x=c(BO.mean$x[i]-0.1,BO.mean$x[i]-0.1),c(BO.mean$BO[i]+BO.sd$BO[i],BO.mean$BO[i]-BO.sd$BO[i]))
}
for(i in 1:nrow(BO.mean)){
  lines(x=c(BO.mean$x[i]-0.1-0.3,BO.mean$x[i]-0.1+0.3),c(BO.mean$BO[i]+BO.sd$BO[i],BO.mean$BO[i]+BO.sd$BO[i]))
  lines(x=c(BO.mean$x[i]-0.1-0.3,BO.mean$x[i]-0.1+0.3),c(BO.mean$BO[i]-BO.sd$BO[i],BO.mean$BO[i]-BO.sd$BO[i]))
}

for(i in 1:nrow(RU.mean)){
  lines(x=c(RU.mean$x[i]+0.1,RU.mean$x[i]+0.1),c(RU.mean$RU[i]+RU.sd$RU[i],RU.mean$RU[i]-RU.sd$RU[i]))
}
for(i in 1:nrow(RU.mean)){
  lines(x=c(RU.mean$x[i]+0.1-0.3,RU.mean$x[i]+0.1+0.3),c(RU.mean$RU[i]+RU.sd$RU[i],RU.mean$RU[i]+RU.sd$RU[i]))
  lines(x=c(RU.mean$x[i]+0.1-0.3,RU.mean$x[i]+0.1+0.3),c(RU.mean$RU[i]-RU.sd$RU[i],RU.mean$RU[i]-RU.sd$RU[i]))
}

dev.off()

cumulative.BO<-aggregate(BO~rep,JD,sum)
cumulative.RU<-aggregate(RU~rep,JD,sum)

cumMean.BO<-mean(cumulative.BO$BO)
cumSd.BO<-sd(cumulative.BO$BO)
cumMean.RU<-mean(cumulative.RU$RU)
cumSd.RU<-sd(cumulative.RU$RU)

t.CUM<-t.test(cumulative.BO$BO,cumulative.RU$RU,paired = T)
p.CUM<-gtools::stars.pval(t.CUM$p.value)

ymn<-min(cumulative.BO$BO,cumulative.RU$RU)
ymx<-max(cumulative.BO$BO,cumulative.RU$RU)
ypad<-(ymx-ymn)*0.25

postscript(file = file.path(dirRes,'cumulative JSD.eps'),
           onefile = T,
           bg = 'white',
           horizontal = F,
           height = 5,
           width = 5)
plot(x=c(rnorm(nrow(cumulative.BO),1,0.1),
         rnorm(nrow(cumulative.RU),2,0.1)),
     y=c(cumulative.BO$BO,
         cumulative.RU$RU),
     xlim = c(0,3),
     ylim = c(ymn-ypad,ymx+ypad),
     xlab = NA,
     xaxt = 'n',
     las = 1,
     bty = 'l',
     ylab = 'Cumulative JSD')

axis(1,c(1,2),c('Reference','Alternative'))
axis(3,c(1,2),c(NA,NA))
mtext(p.CUM,3,1)

points(x=c(1,2),c(cumMean.BO,cumMean.RU),pch=c(16,17))
lines(x=c(1,1),c(cumMean.BO+cumSd.BO,cumMean.BO-cumSd.BO))
lines(x=c(0.8,1.2),c(cumMean.BO+cumSd.BO,cumMean.BO+cumSd.BO))
lines(x=c(0.8,1.2),c(cumMean.BO-cumSd.BO,cumMean.BO-cumSd.BO))

lines(x=c(2,2),c(cumMean.RU+cumSd.RU,cumMean.RU-cumSd.RU))
lines(x=c(1.8,2.2),c(cumMean.RU+cumSd.RU,cumMean.RU+cumSd.RU))
lines(x=c(1.8,2.2),c(cumMean.RU-cumSd.RU,cumMean.RU-cumSd.RU))

dev.off()
