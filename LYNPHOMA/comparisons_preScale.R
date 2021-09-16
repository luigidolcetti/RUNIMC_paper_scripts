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


for (i in unique(geo_list_BO$uid)){
geo_list_BO_trans[geo_list_BO$uid==i,]<-scale(geo_list_BO_trans[geo_list_BO$uid==i,])
geo_list_RU_trans[geo_list_RU$uid==i,]<-scale(geo_list_RU_trans[geo_list_RU$uid==i,])
}

geo_list_total<-as.matrix(rbind(geo_list_BO_trans,geo_list_RU_trans))


# importanceX<-c(0.1,0.1,0.3,
#                0.8,0.8,0.8,
#                0.1,0.8,0.8,
#                0.5,0.8,0.5,
#                0.8,0.8,1,
#                0.8,0.8,1,
#                1,0.1,0.8,
#                0.8,0.8,0.3,
#                1,0.8,0.8,
#                0.8,0.1,0.1)
# 


dimx=5
dimy=5
meta=10


fsom<-FlowSOM::ReadInput(geo_list_total,
                         compensate = F,
                         transform = F,
                         scale = F,
                         silent = F)
fsom<-FlowSOM::BuildSOM(fsom = fsom,
                        colsToUse = colToUse,
                        silent = F,
                        xdim=dimx,
                        ydim=dimy)
# rlen=10,
# mst = 30,
# distf = 1,
# importance=importanceX)
fsom<-FlowSOM::BuildMST(fsom)

saveRDS(fsom,
        file = "C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/COMPARISONS_scaled/5x5_som"  )

fsomMFI<-FlowSOM::GetMFIs(fsom)

fsomMFI<-apply(fsomMFI,2,function(x){(x-min(x))/(max(x)-min(x))})

colnames(fsomMFI)<-unlist(lapply(strsplit(colnames(fsomMFI),'.',fixed = T),'[',2))

rownames(fsomMFI)<-as.character(formatC(1:nrow(fsomMFI),flag='0',digits = 1,format = 'd'))

pdf(file="C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/COMPARISONS_scaled/5x5_heatmap.pdf",
    width = 5,
    height = 5
)
pheatmap::pheatmap(t(fsomMFI),scale = 'none')

dev.off()

SOM_Label<-FlowSOM::GetClusters(fsom)

SOM_Label<-factor(as.character(formatC(SOM_Label,flag='0',format = 'd',digits=1)),levels = rownames(fsomMFI))


geo_list_BO_annot<-dplyr::bind_cols(geo_list_BO,data.frame(SOM=SOM_Label[1:nrow(geo_list_BO)]))


geo_list_RU_annot<-dplyr::bind_cols(geo_list_RU,data.frame(SOM=factor(as.character(SOM_Label[(1+nrow(geo_list_BO)):(nrow(geo_list_BO)+nrow(geo_list_RU))]))))

uids<-unique(geo_list_BO$uid)

for (i in 1:3){
  postscript(file = paste0("C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/COMPARISONS_scaled/5x5_BO_map_",i,".eps"),
             height = 10,
             width = 50,
             horizontal = F,
             bg = 'white',
             pointsize = 10)
  plot(geo_list_BO_annot[geo_list_BO_annot$uid==uids[i],]['SOM'],
       key.pos = 1,
       key.length = 1)
  dev.off()
  postscript(file = paste0("C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/COMPARISONS_scaled/5x5_RU_map_",i,".eps"),
             height = 10,
             width = 50,
             horizontal = F,
             bg = 'white',
             pointsize = 10)
  plot(geo_list_RU_annot[geo_list_RU_annot$uid==uids[i],]['SOM'],
       key.pos = 1,
       key.length = 1)
  dev.off()
}
prePostTable<-table(sf::st_drop_geometry(geo_list_RU_annot)[,c('primary_id','SOM')])

prePostTable<-apply(prePostTable,2,function(x)x/sum(x))

pdf(file="C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/COMPARISONS_scaled/5X5_heatmap_prePost.pdf",
    width = 10,
    height = 3.5
)


pheatmap::pheatmap(prePostTable,
                   scale = 'none',
                   angle_col = 0)

dev.off()


tableSOM<-pbapply::pblapply(1:3,function(i){
  BO<-table(geo_list_BO_annot[geo_list_BO_annot$uid==uids[i],]$SOM)
  RU<-table(geo_list_RU_annot[geo_list_RU_annot$uid==uids[i],]$SOM)
  BO_area<-unlist(lapply(levels(SOM_Label),function(sl){
    sum(sf::st_area(geo_list_BO_annot[geo_list_BO_annot$uid==uids[i] & geo_list_BO_annot$SOM==sl,]))
  }))
  RU_area<-unlist(lapply(levels(SOM_Label),function(sl){
    sum(sf::st_area(geo_list_RU_annot[geo_list_RU_annot$uid==uids[i] & geo_list_RU_annot$SOM==sl,]))
  }))
  cbind(BO,BO_area,RU,RU_area)
})


for (ii in 1:3){
  postscript(file = paste0("C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/COMPARISONS_scaled/5x5_BO_map_split_",ii,".eps"),
             height = 50,
             width = 50,
             horizontal = F,
             bg = 'white',
             pointsize = 10)
 
  par(mar=c(0,0,1,0),mfrow=c(5,5))
  for (i in levels(SOM_Label)){
    TEMP<-geo_list_BO_annot$GEOMETRY[geo_list_BO_annot$uid==uids[ii] & geo_list_BO_annot$SOM ==i]
    if (length(TEMP)==0){
      plot.new()
    } else {
      plot(TEMP,
           xlim = sf::st_bbox(geo_list_BO_annot)[c(1, 3)],
           ylim = sf::st_bbox(geo_list_BO_annot)[c(2, 4)],
           main = i)
    }
    
  }
  dev.off()
  postscript(file = paste0("C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/COMPARISONS_scaled/5x5_RU_map_split_",ii,".eps"),
             height = 50,
             width = 50,
             horizontal = F,
             bg = 'white',
             pointsize = 10)
  
  par(mar=c(0,0,1,0),mfrow=c(5,5))
  for (i in levels(SOM_Label)){
    TEMP<-geo_list_RU_annot$GEOMETRY[geo_list_RU_annot$uid==uids[ii] & geo_list_RU_annot$SOM ==i]
    if (length(TEMP)==0){
      plot.new()
    } else {
      plot(TEMP,
           xlim = sf::st_bbox(geo_list_RU_annot)[c(1, 3)],
           ylim = sf::st_bbox(geo_list_RU_annot)[c(2, 4)],
           main = i)
    }
    
  }
  dev.off()
}

lbl<-'11'
hist(sf::st_area(geo_list_RU_annot[geo_list_RU_annot$uid==uids[1] & geo_list_RU_annot$SOM %in% lbl,]['GEOMETRY']),breaks = 100)
hist(sf::st_area(geo_list_BO_annot[geo_list_BO_annot$uid==uids[1] & geo_list_BO_annot$SOM %in% lbl,]['GEOMETRY']),breaks = 100)

hist(geo_list_RU_trans[geo_list_RU_annot$uid==uids[2] & geo_list_RU_annot$SOM %in% lbl,]$x161dy.cd20.dy161di.)
hist(geo_list_BO_trans[geo_list_BO_annot$uid==uids[2] & geo_list_BO_annot$SOM %in% lbl,]$x161dy.cd20.dy161di.)


TEMP_class<-RUNIMC:::retrieve.RsCollection(
  "C:/Users/luigi/Documents/BH_TEST/MYELOMA/RUN_002/rasterStacks"
)
raster::plot(TEMP_class[[1]]$x161dy.cd20.dy161di.,xlim=c(0,100),ylim=c(0,100))
plot(geo_list_RU[geo_list_RU$uid==unique(geo_list_RU$uid)[1],'GEOMETRY'],col=NA,add=T)
