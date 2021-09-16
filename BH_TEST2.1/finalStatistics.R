rep_1<-read.csv(
  "C:/Users/luigi/Documents/BH_TEST/BH_TEST2.1/REP_1/RESULTS_FAMILY/entropy.csv"  
)
rep_2<-read.csv(
  "C:/Users/luigi/Documents/BH_TEST/BH_TEST2.1/REP_2/RESULTS_FAMILY/entropy.csv"  
)
rep_3<-read.csv(
  "C:/Users/luigi/Documents/BH_TEST/BH_TEST2.1/REP_3/RESULTS_FAMILY/entropy.csv"  
)

dirRes<-"C:/Users/luigi/Documents/BH_TEST/BH_TEST2.1/Final_stat"
dir.create(dirRes)

rep_total<-rbind(cbind(rep_1,rep=1),
                 cbind(rep_2,rep=2),
                 cbind(rep_3,rep=3))


aggr_entropy<-aggregate(entropy~Method+density+rep,rep_total,median)
aggr_entropy<-aggregate(entropy~Method+rep,aggr_entropy,sum)

aggr_entropy<-data.frame(rep=unique(aggr_entropy$rep),BO=aggr_entropy$entropy[aggr_entropy$Method=='BO'],RU=aggr_entropy$entropy[aggr_entropy$Method=='RU'])

t.entropy<-t.test(aggr_entropy$BO,aggr_entropy$RU,paired = T)
p.entropy<-gtools::stars.pval(t.entropy$p.value)

ymn<-min(aggr_entropy$BO,aggr_entropy$RU)
ymx<-max(aggr_entropy$BO,aggr_entropy$RU)
ypad<-(ymx-ymn)*0.25

postscript(file = file.path(dirRes,'cumulative entropy.eps'),
           onefile = T,
           bg = 'white',
           horizontal = F,
           height = 5,
           width = 5)
plot(x=c(rnorm(nrow(aggr_entropy),1,0.1),
         rnorm(nrow(aggr_entropy),2,0.1)),
     y=c(aggr_entropy$BO,
         aggr_entropy$RU),
     xlim = c(0,3),
     ylim = c(ymn-ypad,ymx+ypad),
     xlab = NA,
     xaxt = 'n',
     las = 1,
     bty = 'l',
     ylab = 'Cumulative entropy')

axis(1,c(1,2),c('Reference','Alternative'))
axis(3,c(1,2),c(NA,NA))
mtext(p.entropy,3,1)

points(x=c(1,2),colMeans(aggr_entropy[,c('BO','RU')]),pch=c(16,17))
lines(x=c(1,1),c(mean(aggr_entropy$BO)+sd(aggr_entropy$BO),mean(aggr_entropy$BO)-sd(aggr_entropy$BO)))
lines(x=c(0.8,1.2),c(mean(aggr_entropy$BO)+sd(aggr_entropy$BO),mean(aggr_entropy$BO)+sd(aggr_entropy$BO)))
lines(x=c(0.8,1.2),c(mean(aggr_entropy$BO)-sd(aggr_entropy$BO),mean(aggr_entropy$BO)-sd(aggr_entropy$BO)))

lines(x=c(2,2),c(mean(aggr_entropy$RU)+sd(aggr_entropy$RU),mean(aggr_entropy$RU)-sd(aggr_entropy$RU)))
lines(x=c(1.8,2.2),c(mean(aggr_entropy$RU)+sd(aggr_entropy$RU),mean(aggr_entropy$RU)+sd(aggr_entropy$RU)))
lines(x=c(1.8,2.2),c(mean(aggr_entropy$RU)-sd(aggr_entropy$RU),mean(aggr_entropy$RU)-sd(aggr_entropy$RU)))

dev.off()



rep_1<-read.csv(
  "C:/Users/luigi/Documents/BH_TEST/BH_TEST2.1/REP_1/RESULTS_FAMILY/IOU.csv"  
)
rep_2<-read.csv(
  "C:/Users/luigi/Documents/BH_TEST/BH_TEST2.1/REP_2/RESULTS_FAMILY/IOU.csv"  
)
rep_3<-read.csv(
  "C:/Users/luigi/Documents/BH_TEST/BH_TEST2.1/REP_3/RESULTS_FAMILY/IOU.csv"  
)

dirRes<-"C:/Users/luigi/Documents/BH_TEST/BH_TEST2.1/Final_stat"


rep_total<-rbind(cbind(rep_1,rep=1),
                 cbind(rep_2,rep=2),
                 cbind(rep_3,rep=3))


aggr_IOU<-aggregate(IOU~Method+density+rep,rep_total,median)
aggr_IOU<-aggregate(IOU~Method+rep,aggr_IOU,sum)

aggr_IOU<-data.frame(rep=unique(aggr_IOU$rep),BO=aggr_IOU$IOU[aggr_IOU$Method=='BO'],RU=aggr_IOU$IOU[aggr_IOU$Method=='RU'])

t.IOU<-t.test(aggr_IOU$BO,aggr_IOU$RU,paired = T)
p.IOU<-gtools::stars.pval(t.IOU$p.value)

ymn<-min(aggr_IOU$BO,aggr_IOU$RU)
ymx<-max(aggr_IOU$BO,aggr_IOU$RU)
ypad<-(ymx-ymn)*0.25

postscript(file = file.path(dirRes,'cumulative IOU.eps'),
           onefile = T,
           bg = 'white',
           horizontal = F,
           height = 5,
           width = 5)
plot(x=c(rnorm(nrow(aggr_IOU),1,0.1),
         rnorm(nrow(aggr_IOU),2,0.1)),
     y=c(aggr_IOU$BO,
         aggr_IOU$RU),
     xlim = c(0,3),
     ylim = c(ymn-ypad,ymx+ypad),
     xlab = NA,
     xaxt = 'n',
     las = 1,
     bty = 'l',
     ylab = 'Cumulative IOU')

axis(1,c(1,2),c('Reference','Alternative'))
axis(3,c(1,2),c(NA,NA))
mtext(p.IOU,3,1)

points(x=c(1,2),colMeans(aggr_IOU[,c('BO','RU')]),pch=c(16,17))
lines(x=c(1,1),c(mean(aggr_IOU$BO)+sd(aggr_IOU$BO),mean(aggr_IOU$BO)-sd(aggr_IOU$BO)))
lines(x=c(0.8,1.2),c(mean(aggr_IOU$BO)+sd(aggr_IOU$BO),mean(aggr_IOU$BO)+sd(aggr_IOU$BO)))
lines(x=c(0.8,1.2),c(mean(aggr_IOU$BO)-sd(aggr_IOU$BO),mean(aggr_IOU$BO)-sd(aggr_IOU$BO)))

lines(x=c(2,2),c(mean(aggr_IOU$RU)+sd(aggr_IOU$RU),mean(aggr_IOU$RU)-sd(aggr_IOU$RU)))
lines(x=c(1.8,2.2),c(mean(aggr_IOU$RU)+sd(aggr_IOU$RU),mean(aggr_IOU$RU)+sd(aggr_IOU$RU)))
lines(x=c(1.8,2.2),c(mean(aggr_IOU$RU)-sd(aggr_IOU$RU),mean(aggr_IOU$RU)-sd(aggr_IOU$RU)))

dev.off()









