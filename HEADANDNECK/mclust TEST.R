library(mclust)
library(factoextra)
BO_TEMP<-scale(geo_list_BO_trans)
RU_TEMP<-scale(geo_list_RU_trans)
BO_BIC<-mclust::mclustBIC(BO_TEMP,G = 1:25,modelNames = 'VVV')
RU_BIC<-mclust::mclustBIC(RU_TEMP,G = 1:25,modelNames = 'VVV')
plot(BO_BIC)
plot(RU_BIC)

BO_MC<-mclust::Mclust(data = BO_TEMP,
                      x = BO_BIC  )

RU_MC<-mclust::Mclust(data = RU_TEMP,
                      x = RU_BIC  )

fviz_mclust(BO_MC, what = "classification", pointsize = 0.1)
fviz_mclust(RU_MC, "classification", geom = "point", pointsize = 0.1)

BO_SOM_Label<-BO_MC$classification

BO_SOM_Label<-factor(as.character(formatC(BO_SOM_Label,flag='0',format = 'd',digits=1)))

geo_list_BO_annot<-dplyr::bind_cols(geo_list_BO,data.frame(SOM=BO_SOM_Label[1:nrow(geo_list_BO)]))


RU_SOM_Label<-RU_MC$classification

RU_SOM_Label<-factor(as.character(formatC(RU_SOM_Label,flag='0',format = 'd',digits=1)))

geo_list_RU_annot<-dplyr::bind_cols(geo_list_RU,data.frame(SOM=RU_SOM_Label[1:nrow(geo_list_RU)]))

uids<-unique(geo_list_BO$uid)

plot(geo_list_BO_annot[geo_list_BO_annot$uid==uids[3],'SOM'])
plot(geo_list_RU_annot[geo_list_RU_annot$uid==uids[3],'SOM'])

BO_MFI<-geo_list_BO_annot[,3:33,drop=T]
BO_MFI<-aggregate(.~SOM,BO_MFI,median)
BO_MFI<-BO_MFI[,-1]
BO_MFI<-apply(BO_MFI,2,function(x){(x-min(x))/(max(x)-min(x))})
colnames(BO_MFI)<-c('CD38','SMA','VIM',
                        'CD14','CD33','CD16',
                        'PANK','CD11b','CD15',
                        'IgD','CD45','CD24',
                        'CD11c','FoxP3','CD4',
                        'E-CAD','CD68','CD20',
                        'CD8','Arg1','CD45ra',
                        'CD74','CD127','COLL',
                        'CD3','CD27','CD45ro',
                        'CD25','DNA1','DNA2')
pheatmap::pheatmap(t(BO_MFI))


RU_MFI<-geo_list_RU_annot[,9:39,drop=T]
RU_MFI<-aggregate(.~SOM,RU_MFI,median)
RU_MFI<-RU_MFI[,-1]
RU_MFI<-apply(RU_MFI,2,function(x){(x-min(x))/(max(x)-min(x))})
colnames(RU_MFI)<-c('CD38','SMA','VIM',
                    'CD14','CD33','CD16',
                    'PANK','CD11b','CD15',
                    'IgD','CD45','CD24',
                    'CD11c','FoxP3','CD4',
                    'E-CAD','CD68','CD20',
                    'CD8','Arg1','CD45ra',
                    'CD74','CD127','COLL',
                    'CD3','CD27','CD45ro',
                    'CD25','DNA1','DNA2')
pheatmap::pheatmap(t(RU_MFI))
