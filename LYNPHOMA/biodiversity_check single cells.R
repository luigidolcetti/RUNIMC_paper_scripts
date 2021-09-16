geo_list_BO<-sf::st_read(
  "C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/analysis/BODEN/ex.sqlite"
)

geo_list_RU<-sf::st_read(
  "C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/analysis/RUNIMC1/topLayers.sqlite"
)


trns<-scales::modulus_trans(0)

colToUse<-colnames(geo_list_BO)[3:17]

geo_list_BO_trans<-sf::st_drop_geometry(geo_list_BO)[,colToUse]
geo_list_BO_trans<-trns$transform(geo_list_BO_trans)

geo_list_RU_trans<-sf::st_drop_geometry(geo_list_RU)[,colToUse]
geo_list_RU_trans<-trns$transform(geo_list_RU_trans)

colnames(geo_list_BO_trans)<-unlist(lapply(strsplit(colnames(geo_list_BO_trans),'.',fixed = T),'[',2))

geo_list_BO_trans<-cbind(uid=geo_list_BO$uid[drop=T],geo_list_BO_trans,area=sf::st_area(geo_list_BO))

fcs_BO<-lapply(unique(geo_list_BO_trans$uid),function(x,y=geo_list_BO_trans){
  TEMP<-as.matrix(y[y$uid==x,-1])
  out<-flowCore::flowFrame(TEMP)})

colnames(geo_list_RU_trans)<-unlist(lapply(strsplit(colnames(geo_list_RU_trans),'.',fixed = T),'[',2))

geo_list_RU_trans<-cbind(uid=geo_list_RU$uid[drop=T],geo_list_RU_trans,area=sf::st_area(geo_list_RU))


cln<-colnames(geo_list_BO_trans)[2:17]

cln<-c('cd3','cd8','cd4','cd20','cd21','cd68','cd31','asma','cd16')
trsh<-c(1.25, 1.75, 1.5, 0.75, 1.5, 2.75, 1.25, 1.5, 0.8)
names(trsh)<-cln



uids<-unique(geo_list_BO_trans$uid)
geo_list_BO_Class<-dplyr::bind_cols(geo_list_BO[,'uid'],CLASS='UNDF')
geo_list_RU_Class<-dplyr::bind_cols(geo_list_RU[,'uid'],CLASS='UNDF')

stmnt<-expression(cd3>trsh['cd3'] & cd4>trsh['cd4'] & cd8<trsh['cd8'])

geo_list_BO_Class$CLASS[with(geo_list_BO_trans,eval(stmnt)) & geo_list_BO_Class$CLASS=='UNDF']<-'1_CD4-Lymph'
geo_list_RU_Class$CLASS[with(geo_list_RU_trans,eval(stmnt)) & geo_list_RU_Class$CLASS=='UNDF']<-'1_CD4-Lymph'

stmnt<-expression(cd3>trsh['cd3'] & cd4<trsh['cd4'] & cd8>trsh['cd8'])

geo_list_BO_Class$CLASS[with(geo_list_BO_trans,eval(stmnt)) & geo_list_BO_Class$CLASS=='UNDF']<-'2_CD8-Lymph'
geo_list_RU_Class$CLASS[with(geo_list_RU_trans,eval(stmnt)) & geo_list_RU_Class$CLASS=='UNDF']<-'2_CD8-Lymph'

stmnt<-expression(cd3<trsh['cd3'] & cd68>trsh['cd68'])

geo_list_BO_Class$CLASS[with(geo_list_BO_trans,eval(stmnt)) & geo_list_BO_Class$CLASS=='UNDF']<-'7_mph'
geo_list_RU_Class$CLASS[with(geo_list_RU_trans,eval(stmnt)) & geo_list_RU_Class$CLASS=='UNDF']<-'7_mph'

stmnt<-expression(cd3<trsh['cd3'] & cd68<trsh['cd68'] & asma>trsh['asma'])

geo_list_BO_Class$CLASS[with(geo_list_BO_trans,eval(stmnt)) & geo_list_BO_Class$CLASS=='UNDF']<-'8_fbr'
geo_list_RU_Class$CLASS[with(geo_list_RU_trans,eval(stmnt)) & geo_list_RU_Class$CLASS=='UNDF']<-'8_fbr'

stmnt<-expression(cd3<trsh['cd3'] & cd68<trsh['cd68'] & cd31>trsh['cd31'])

geo_list_BO_Class$CLASS[with(geo_list_BO_trans,eval(stmnt)) & geo_list_BO_Class$CLASS=='UNDF']<-'9_endt'
geo_list_RU_Class$CLASS[with(geo_list_RU_trans,eval(stmnt)) & geo_list_RU_Class$CLASS=='UNDF']<-'9_endt'

stmnt<-expression(cd20<trsh['cd20'] & cd21>trsh['cd21'])
geo_list_BO_Class$CLASS[with(geo_list_BO_trans,eval(stmnt)) & geo_list_BO_Class$CLASS=='UNDF']<-'5_CD21+'
geo_list_RU_Class$CLASS[with(geo_list_RU_trans,eval(stmnt)) & geo_list_RU_Class$CLASS=='UNDF']<-'5_CD21+'

stmnt<-expression(cd20>trsh['cd20'] & cd21>trsh['cd21'])
geo_list_BO_Class$CLASS[with(geo_list_BO_trans,eval(stmnt)) & geo_list_BO_Class$CLASS=='UNDF']<-'4_B-cellCD21+'
geo_list_RU_Class$CLASS[with(geo_list_RU_trans,eval(stmnt)) & geo_list_RU_Class$CLASS=='UNDF']<-'4_B-cellCD21+'

stmnt<-expression(cd16>trsh['cd16'] )
geo_list_BO_Class$CLASS[with(geo_list_BO_trans,eval(stmnt))& geo_list_BO_Class$CLASS=='UNDF']<-'6_nk'
geo_list_RU_Class$CLASS[with(geo_list_RU_trans,eval(stmnt))& geo_list_RU_Class$CLASS=='UNDF']<-'6_nk'

stmnt<-expression(cd20>trsh['cd20'])
geo_list_BO_Class$CLASS[with(geo_list_BO_trans,eval(stmnt)) & geo_list_BO_Class$CLASS=='UNDF']<-'3_B-Cells'
geo_list_RU_Class$CLASS[with(geo_list_RU_trans,eval(stmnt)) & geo_list_RU_Class$CLASS=='UNDF']<-'3_B-Cells'

rts<-RUNIMC:::retrieve.RsCollection(
  "C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/rasterStacks"
)

bbbox<-lapply(rts,raster::extent)

biodiversity_BO<-lapply(unique(geo_list_BO_Class$uid),function(u){
  newpoly<-geo_list_BO_Class[geo_list_BO_Class$uid==u,]
  grd<-sf::st_make_grid(x = newpoly,
                        cellsize = 50,
                        square = F)
  
  bd<-lapply(grd,function(x){
    insidePoly<-sf::st_covers(x,newpoly)[[1]]
    out<-table(newpoly[insidePoly,'CLASS',drop=T])
    # out<-out/sum(out)
    # out<--sum(out*log(out))
    # out<-sf::st_sf(nPoly=length(insidePoly),
    #                meanArea=mean(sf::st_area(newpoly[insidePoly,])),
    #                entropy=out,
    #                perplexity=exp(out),
    #                GEOMETRY=sf::st_sfc(x),
    #                sf_column_name = 'GEOMETRY')

  })
  # out<-do.call(dplyr::bind_rows,bd)
  # return(out)
})

# plot(biodiversity_BO[[1]]['nPoly'],breaks = c(0,10,20,30,40,50,60))


biodiversity_RU<-lapply(unique(geo_list_RU_Class$uid),function(u){
  newpoly<-geo_list_RU_Class[geo_list_RU_Class$uid==u,]
  grd<-sf::st_make_grid(x = newpoly,
                        cellsize = 50,
                        square = F)
  
  bd<-lapply(grd,function(x){
    insidePoly<-sf::st_covers(x,newpoly)[[1]]
    out<-table(newpoly[insidePoly,'CLASS',drop=T])
    # out<-out/sum(out)
    # out<--sum(out*log(out))
    # out<-sf::st_sf(nPoly=length(insidePoly),
    #                meanArea=mean(sf::st_area(newpoly[insidePoly,])),
    #                entropy=out,
    #                perplexity=exp(out),
    #                GEOMETRY=sf::st_sfc(x),
    #                sf_column_name = 'GEOMETRY')
    # 
  })
  # out<-do.call(dplyr::bind_rows,bd)
  # return(out)
})


biodifference<-lapply(1:3,function(u){
  tile<-length(biodiversity_RU[[u]])
  out<-lapply(1:tile,function(tl){
    # browser()
    cellType_BO<-names(biodiversity_BO[[u]][[tl]])
    cellType_RU<-names(biodiversity_RU[[u]][[tl]])
    in_BO<-cellType_BO[!(cellType_BO %in% cellType_RU)]
    in_RU<-cellType_RU[!(cellType_RU %in% cellType_BO)]
    if (length(in_BO)>0) out_BO<-data.frame(direction=rep(-1,length(in_BO)),cellType=in_BO) else out_BO<-data.frame(direction=numeric(0),cellType=character(0))
    if (length(in_RU)>0) out_RU<-data.frame(direction=rep(1,length(in_RU)),cellType=in_RU) else out_RU<-data.frame(direction=numeric(0),cellType=character(0))
  out<-rbind.data.frame(out_RU,out_BO)
    })
  out<-do.call(rbind.data.frame,out)
  })

for (i in 1:3){
TEMP<-table(biodifference[[i]])
postscript(file = paste0("C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/Ecology/exclusive_cells_",i,".eps"),
           height = 10,
           width = 20,
           onefile = F,
           horizontal = F,
           bg = 'white',
           paper = 'special',
           pointsize = 12)
plot(NA,xlim=c(0,ncol(TEMP)+1),ylim=c(0,300),xaxt='na',ylab='number of tiles',cex=0.8,las=1)
axis(1,at=1:(length(colnames(TEMP))),labels = colnames(TEMP),cex=0.5)
for (u in 1:(ncol(TEMP))){
  rect(u-0.4,0,u,TEMP[1,u],col = 'white',border = 'black')
  rect(u,0,u+0.4,TEMP[2,u],col='gray50',border='black')
  xaxt='na'
}
dev.off()
}

