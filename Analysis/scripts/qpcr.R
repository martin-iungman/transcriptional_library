library(tidyverse)
qpcr=readxl::read_xlsx("Experimental_data/qPCR_Gates.xlsx")
qpcr=qpcr%>%filter(Task=="Unknown", Sample!="NTC")%>%mutate(Cq=as.numeric(Cq), Gate=str_remove(Sample,"_rep_.")%>%str_remove("Gate_")%>%as.numeric()%>%as_factor(), Rep=str_remove(Sample, "Gate_._"))
qpcr$tech=rep(c(1,2,3), times=7*3*2)
qpcr_wide=qpcr%>%select(Cq,Target,Sample, Rep,tech,Gate)%>%
  pivot_wider(values_from=Cq, names_from=Target)%>%mutate(delta=EGFP-GAPDH)

df=qpcr_wide%>%filter(Rep!="rep_1", Gate!=1)
df=df%>%group_by(Rep,Gate)%>%mutate(mean_GAPDH=mean(GAPDH), delta=2^-(EGFP-mean_GAPDH))
mean_g2=df_poster%>%filter(Gate==2)%>%select(Rep,delta)%>%group_by(Rep)%>%summarise(mean_g2=mean(delta))  
df_poster=df_poster%>%left_join(mean_g2)%>%mutate(FC=delta/mean_g2)%>%mutate(Gate=as.factor(as.numeric(Gate)-1))
p2=df_poster%>%filter(Rep=="rep_2")%>%
  ggplot(aes(Gate, FC, fill=Gate, col=Gate))+geom_line(linewidth=2)+geom_point(size=5)+geom_point(shape = 1,size=5,colour = "black")+ggtitle("Rep 1")+scale_y_log10()+geom_point()+ylab("")+rcartocolor::scale_color_carto_d(palette="Emrld")+theme(legend.position = "none", text=element_text(size=15))
p3=df_poster%>%filter(Rep=="rep_3")%>%
  ggplot(aes(Gate, FC, fill=Gate, col=Gate))+geom_line(linewidth=2)+geom_point(size=5)+geom_point(shape = 1,size=5,colour = "black")+ggtitle("Rep 2")+scale_y_log10()+geom_point()+ylab("")+rcartocolor::scale_color_carto_d(palette="Emrld")+theme(legend.position = "none", text=element_text(size=15))
cowplot::plot_grid(p2,p3) #el posta


