ent<-read.csv(
  "C:/Users/luigi/Documents/BH_TEST/BH_TEST2.1/REP_1/RESULTS_FAMILY/entropy.csv"  
)

ent$Method<-factor(ent$Method,levels=c('BO','RU'))

ggp<-ggplot2::ggplot()+ggplot2::geom_boxplot(mapping =ggplot2::aes(y=entropy,
                                                                   x=density,
                                                                   fill = Method),
                                             position = 'dodge',
                                             data = ent)+
  ggplot2::scale_fill_manual(values = c('white','gray50'))+
  ggplot2::ylim(c(0,1))+
  ggplot2::theme_classic()

ggplot2::ggsave(plot = ggp,
                filename = "C:/Users/luigi/Documents/BH_TEST/BH_TEST2.1/REP_1/RESULTS_FAMILY/entropy_BS.eps",
                device = 'eps',
                width = 5,
                height = 5,
                dpi = 600)


IOU<-read.csv(
  "C:/Users/luigi/Documents/BH_TEST/BH_TEST2.1/REP_1/RESULTS_FAMILY/IOU.csv"  
)

IOU$Method<-factor(IOU$Method,levels=c('BO','RU'))

ggp<-ggplot2::ggplot()+ggplot2::geom_boxplot(mapping =ggplot2::aes(y=IOU,
                                                                   x=density,
                                                                   fill = Method),
                                             position = 'dodge',
                                             data = IOU)+
  ggplot2::scale_fill_manual(values = c('white','gray50'))+
  ggplot2::ylim(c(0,1))+
  ggplot2::theme_classic()

ggplot2::ggsave(plot = ggp,
                filename = "C:/Users/luigi/Documents/BH_TEST/BH_TEST2.1/REP_1/RESULTS_FAMILY/IOU_BS.eps",
                device = 'eps',
                width = 5,
                height = 5,
                dpi = 600)
