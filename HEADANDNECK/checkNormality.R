library(RUNIMC)

mystudy<-retrieve(
  "C:/Users/k1343421/Documents/RUNIMC_PAPER/RUNIMC/RUN_001/study.xml"  
)

raster::plot(mystudy$raster$uid.acac583e59b91dd2b55f8bc72c751e7d$x191ir.dna1.ir191di.)

TEMP<-raster::getValues(mystudy$raster$uid.acac583e59b91dd2b55f8bc72c751e7d$x143nd.vimentin.nd143di.)

hist(TEMP)

censoredZero<-TEMP

trnsf<-scales::modulus_trans(0)

# plot(censoredZero,trnsf$transform(censoredZero))

qqnorm(trnsf$transform(censoredZero[sample(1:length(censoredZero),100)]))
qqnorm((censoredZero[sample(1:length(censoredZero),100)]))
hist(trnsf$transform(censoredZero))
