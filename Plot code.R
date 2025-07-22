library("ggplot2")
library("caret")
library("scales")
library("pheatmap")
library("gridExtra")
library("cowplot")
library("PNWColors")
library("MetBrewer")
library("MASS")
library("reshape2")
library("dplyr")
library("corrplot")
library("psych")
library("splines")
library("qcc")
library("mgcv")
library("car")
library("dlnm")
library(wk)
library(s2)
library(deldir)
library(sf)
library(sp)
library(spData)
library(spdep)
library(openxlsx)
library(tsModel)
library("metafor")
library("forestplot")
library("mvmeta")
library(ggspatial)
library(ggsn)
library(rgdal)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(grid)
library(mapproj)
library(openxlsx)
library(cowplot)
library(showtext)
library(paletteer)
library(stringr)
library(rgeos)
library(ggnewscale)
library(gtable)
library(ggtext)
library(ggh4x)
library("trend")
getwd()

setwd("D:/~")

#######################Spatial auto-correlation##################################################
#Produce the age-specific, and PLAD-specific incidence
setwd("D:/~")
Provincename<-c("Beijing","Tianjin","Hebei","Shanxi","Inner Mongolia",
                "Liaoning","Jilin","Heilongjiang","Shanghai","Jiangsu","Zhejiang",
                "Anhui","Fujian","Jiangxi","Shandong","Henan","Hubei","Hunan",
                "Guangdong","Guangxi","Hainan","Chongqing","Sichuan","Guizhou",
                "Yunnan","Tibet","Shaanxi","Gansu","Qinghai","Ningxia","Xinjiang")
#23 age groups
for (i in 1:23)
{Name<-paste0("Agegroup",i,".xlsx")
Data<-read.xlsx(Name)
Incidence<-aggregate(Case~Year+Order,data=Data,sum)
Incidence$Order<-rep(Provincename,each=15)
Incidence$Population<-aggregate(Population.book~Year+Order,data=Data,mean)[,3]
Incidence$Incidence<-Incidence$Case/Incidence$Population*100000
assign(paste("Incidence.",i,sep=""),Incidence)
}

#Merge the incidence in the age groups 0–5
Incidence.target<-Incidence.1
Incidence.target$Case<-Incidence.1$Case+Incidence.2$Case+
  Incidence.3$Case+Incidence.4$Case+Incidence.5$Case
Incidence.target$Population<-Incidence.1$Population+Incidence.2$Population+
  Incidence.3$Population+Incidence.4$Population+Incidence.5$Population
Incidence.target$Incidence<-Incidence.target$Case/Incidence.target$Population*100000

#Merge the incidence in the age groups 6–14
Incidence.study<-Incidence.6
Incidence.study$Case<-Incidence.6$Case+
  Incidence.7$Case+Incidence.8$Case+
  Incidence.9$Case+Incidence.10$Case+
  Incidence.11$Case+Incidence.12$Case+
  Incidence.13$Case+Incidence.14$Case
Incidence.study$Population<-Incidence.6$Population+
  Incidence.7$Population+Incidence.8$Population+
  Incidence.9$Population+Incidence.10$Population+
  Incidence.11$Population+Incidence.12$Population+
  Incidence.13$Population+Incidence.14$Population
Incidence.study$Incidence<-Incidence.study$Case/Incidence.study$Population*100000

#Merge the incidence in the age groups 15–22
Incidence.nontarget<-Incidence.15
Incidence.nontarget$Case<-Incidence.15$Case+Incidence.16$Case+
  Incidence.17$Case+Incidence.18$Case+
  Incidence.19$Case+Incidence.20$Case+
  Incidence.21$Case+Incidence.22$Case
Incidence.nontarget$Population<-Incidence.15$Population+Incidence.16$Population+
  Incidence.17$Population+Incidence.18$Population+
  Incidence.19$Population+Incidence.20$Population+
  Incidence.21$Population+Incidence.22$Population
Incidence.nontarget$Incidence<-Incidence.nontarget$Case/Incidence.nontarget$Population*100000

#The whole population
library(reshape2)
setwd("D:/~")
Incidence.0<-read.xlsx("各省份全年病毒性肝炎发病率.xlsx")
Case<-Incidence.0[,c(1,(2+1):(2+31))]
Incidence<-Incidence.0[,c(1,(35+1):(35+31))]
Incidence.whole<-melt(Case,id.vars="year",variable.name="Order",value.name="Case")
Incidence.whole$Incidence<-melt(Incidence,id.vars="year",variable.name="Order",value.name="Incidence")[,3]
colnames(Incidence.whole)<-c("Year","Order","Case","Incidence")

#Tests for spatial auto-correlation
setwd("D:/~")
Weight<-read.xlsx("空间权重矩阵.xlsx")
Weight<-as.matrix(Weight[,-1])
#The incidence for the whole population (Incidence.whole)
Incidence.whole<-read.xlsx("Agegroup1-22.xlsx")
#Global spatial auto-correlation
Spatial.three.populations<-as.data.frame(matrix(0,ncol=8,nrow=15))
w_contact<-mat2listw(Weight, style="W")
for (i in 2005:2019){
  inf<-Incidence.whole[which(Incidence.whole$Year==i),4]
  moran_test_contact<-moran.test(inf,w_contact,alternative="two.sided")
  Spatial.three.populations[(i-2004),1]<-moran_test_contact$estimate[1]
  Spatial.three.populations[(i-2004),2]<-moran_test_contact$p.value
}

#The incidence in the age groups 0–4
Incidence.target<-read.xlsx("Agegroup1-4.xlsx")
for (i in 2005:2019){
  inf<-Incidence.target[which(Incidence.target$Year==i),5]
  moran_test_contact<-moran.test(inf,w_contact,alternative="two.sided")
  Spatial.three.populations[(i-2004),3]<-moran_test_contact$estimate[1]
  Spatial.three.populations[(i-2004),4]<-moran_test_contact$p.value
}

#The incidence in the age groups 5–14
Incidence.study<-read.xlsx("Agegroup5-14.xlsx")
for (i in 2005:2019){
  inf<-Incidence.study[which(Incidence.study$Year==i),5]
  moran_test_contact<-moran.test(inf,w_contact,alternative="two.sided")
  Spatial.three.populations[(i-2004),5]<-moran_test_contact$estimate[1]
  Spatial.three.populations[(i-2004),6]<-moran_test_contact$p.value
}

#The incidence in the age groups 15–22
Incidence.nontarget<-read.xlsx("Agegroup15-22.xlsx")
for (i in 2005:2019){
  inf<-Incidence.nontarget[which(Incidence.nontarget$Year==i),5]
  moran_test_contact<-moran.test(inf,w_contact,alternative="two.sided")
  Spatial.three.populations[(i-2004),7]<-moran_test_contact$estimate[1]
  Spatial.three.populations[(i-2004),8]<-moran_test_contact$p.value
}
colnames(Spatial.three.populations)<-c("I1","P1",
                                       "I2","P2",
                                       "I3","P3",
                                       "I4","P4")

#######################Temporal trend################################################
setwd("D:/~")
#Temporal trends of the incidence of four population
Mannkendall.three.populations<-as.data.frame(matrix(0,ncol=8+1,nrow=32))
Code<-c("1-22","1-5 update","6-14 update","15-22")
for(j in 2:4){
  Name<-paste0("Agegroup",Code[j],".xlsx")
  Incidence<-read.xlsx(Name)
  Incidence$Number<-rep(c(1:31),each=15)
  for(i in 1:31)
  {result<-mk.test(Incidence$Incidence[which(Incidence$Number==i)])
  Mannkendall.three.populations[i,j*2]<-result$statistic
  Mannkendall.three.populations[i,j*2+1]<-result$p.value
  }
  Nation<-aggregate(Case~Year,Incidence,sum)
  Nation$Population<-aggregate(Population~Year,Incidence,sum)[,2]
  Nation$Incidence<-Nation$Case/Nation$Population*100000
  result<-mk.test(Nation$Incidence)
  Mannkendall.three.populations[32,j*2]<-result$statistic
  Mannkendall.three.populations[32,j*2+1]<-result$p.value
}
j=1
Name<-paste0("Agegroup",Code[j],".xlsx")
Incidence<-read.xlsx(Name)
Incidence$Number<-rep(c(1:31),each=15)
for(i in 1:31)
{result<-mk.test(Incidence$Incidence[which(Incidence$Number==i)])
Mannkendall.three.populations[i,j*2]<-result$statistic
Mannkendall.three.populations[i,j*2+1]<-result$p.value
}
Nation<-read.xlsx("D:/研究生/乙肝干预/全国乙肝数据/整理后的数据/乙肝/各省份全年病毒性肝炎发病率.xlsx")
result<-mk.test(Nation$China_Incidence)
Mannkendall.three.populations[32,j*2]<-result$statistic
Mannkendall.three.populations[32,j*2+1]<-result$p.value
Mannkendall.three.populations$V1<-c(Provincename,"Nation")
colnames(Mannkendall.three.populations)<-c("Province",
                                           "MK1","P1",
                                           "MK2","P2",
                                           "MK3","P3",
                                           "MK4","P4")

#######################Monthly incidence of hepatitis B in different age groups (Figure S3)################################################
setwd("D:/~")
Incidence.month<-read.xlsx("Monthly incidence.xlsx")
Incidence.month<-Incidence.month[which(Incidence.month$Group=="Whole population"|
                                         Incidence.month$Group=="Non-intervention population"),]
Incidence.month$Time<-as.Date(Incidence.month$Time,origin="1899-12-30")
Incidence.month$Category[which(Incidence.month$Category=="Whole population")]<-c("Three different populations")
Incidence.other<-read.xlsx("Figure S3-S6 update.xlsx")
Incidence.other$Time<-as.Date(paste0(Incidence.other$Year,'-',
                                     ifelse(Incidence.other$Month<10,paste0("0",Incidence.other$Month),Incidence.other$Month),
                                     '-01'))
Incidence.study<-aggregate(Incidence.other$Case,by=list(Incidence.other$Year,
                                                        Incidence.other$Month),sum)
colnames(Incidence.study)<-c("Year","Month","Case")
Incidence.study$Population<-aggregate(Incidence.other$Population,by=list(Incidence.other$Year,
                                                                         Incidence.other$Month),sum)$x
Incidence.study$Incidence<-Incidence.study$Case/Incidence.study$Population*100000
Incidence.study$Group<-c("Study population")
Incidence.study$Category<-c("Three different populations")
Incidence.study$Time<-as.Date(paste0(Incidence.study$Year,'-',
                                     ifelse(Incidence.study$Month<10,paste0("0",Incidence.study$Month),Incidence.study$Month),
                                     '-01'))
Incidence.other<-Incidence.other[,c("Time","Incidence","Group","Category")]
Incidence.study<-Incidence.study[,c("Time","Incidence","Group","Category")]
Incidence.other<-rbind(Incidence.other,Incidence.study)
Incidence.month<-rbind(Incidence.month,Incidence.other)
Incidence.month$Group<-factor(Incidence.month$Group,
                              levels=c("Whole population",
                                       "Study population","Non-intervention population",
                                       "5-year old children","6-year old children",
                                       "7-year old children","8-year old children",
                                       "9-year old children","Children aged 10–14",
                                       "Population aged 15–19","Population aged 20–24",
                                       "Population aged 25–29"))
Incidence.month$Category<-factor(Incidence.month$Category,
                                 levels=c("Three different populations",
                                          "Study populations aged 5",
                                          "Study populations aged 6–7",
                                          "Study populations aged 8–9",
                                          "Study populations aged 10–19",
                                          "Study populations aged 20–29"))
plot2<-ggplot(Incidence.month,aes(y=Incidence,x=Time,colour=Group))+
  geom_line(aes(group=Group),linewidth=0.8)+
  geom_point(aes(group=Group),size=1.2)+
  scale_color_manual(values=c(met.brewer("Demuth",10)[c(7,10)],met.brewer("Redon",12)))+
  scale_x_date(date_breaks="1 year",date_labels="%Y")+
  scale_y_continuous(position="left")+#limits=c(0,15.5),breaks=c(0,5,10,15),
  xlab("Year")+ylab("Incidence (cases per 100,000 individuals)")+
  theme_bw()+
  labs(fill="Group")+
  facet_wrap(Category~.,nrow=3,ncol=2,scales = "free_y")+
  theme(legend.position="top",  
        legend.title=element_blank(),
        legend.key.height=unit(10, "mm"),
        legend.key.width=unit(10, "mm"),
        legend.background=element_blank(),
        legend.text=element_text(size=15,family='serif'),
        plot.title=element_text(size=15,color = 'black',family ='serif'),
        title=element_text(size=15,family ='serif'),
        axis.title=element_text(size=15,family ='serif'),
        axis.text=element_text(size=14,face="plain",color="black",family ='serif'),
        strip.background=element_blank(),
        strip.text=element_text(face="bold",size=15,family="serif"))+
  guides(color=guide_legend(ncol=5))+
  facetted_pos_scales(
    y = list(Category=="Whole population"~scale_y_continuous(limits=c(0,9),breaks=c(0,3,6,9)),
             Category=="Three different populations"~scale_y_continuous(limits=c(0,12),breaks=c(0,4,8,12)),
             #Category=="Study population aged 4–5"~scale_y_continuous(limits=c(0,3),breaks=c(0,1,2,3)),
             Category=="Study population aged 6–7"~scale_y_continuous(limits=c(0,3),breaks=c(0,1,2,3)),
             Category=="Study population aged 8–9"~scale_y_continuous(limits=c(0,3),breaks=c(0,1,2,3)),
             Category=="Study population aged 10–29"~scale_y_continuous(limits=c(0,15),breaks=c(0,5,10,15))))
setwd("D:/New results")
tiff(file="Figure S3.tiff",width=420,height=270,units="mm",res=300,compression="lzw")
print(plot2)
dev.off()

png(file="Figure S3.png",width=420,height=270,units="mm",res=300)
print(plot2)
dev.off()

#######################The seasonality of hepatitis B cases in different age groups (Figure S4)############################
#Whole population
col<-met.brewer("Troy",8)
setwd("D:/~")
Season1<-read.xlsx("Agegroup.xlsx")
#Set thousand separator
Season1$Label<-format(Season1$std.Cases,big.mark=",")
#Bar plot
Figure1<-ggplot()+
  geom_col(data=Season1,aes(x=Order,y=std.Cases),color='black', width=0.7, fill=col[2],alpha=0.6)+
  theme_bw()+
  scale_x_continuous(limits=c(0.5,12.5), breaks=seq(1,12,1),labels=c("Jan","Feb","Mar","Apr",
                                                                     "May","Jun","Jul","Aug",
                                                                     "Sep","Oct","Nov","Dec"))+
  geom_text(data=Season1,aes(x=Order,y=std.Cases,label=Label),vjust=-0.5,size=6.0,fontface='bold',family="serif")+
  scale_y_continuous(expand = c(0,0),limits=c(0,120000),breaks=c(0,40000,80000,120000),labels=c("0","40,000","80,000","120,000"),
                     sec.axis = sec_axis(~./400,name = 'Standardized incidence of hepatitis B (1/100,000)',
                                         breaks=c(0,100,200,300)))+
  geom_line(data=Season1,aes(x=Order,y=std.Incidence*400),linewidth=1,color=col[1])+
  geom_point(data=Season1,aes(x=Order,y=std.Incidence*400),size=2,color=col[1])+
  scale_linetype_manual(values=c("solid"))+
  xlab("Month")+ylab("Standardized cases of hepatitis B")+
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text = element_text(colour = "black",size = 20,family="serif"),
        axis.title = element_text(size = 20,color = 'black',family="serif"),
        plot.title=element_text(size=22,color = 'black',family ='serif'),
        plot.margin = margin(t = 1,
                             r = 0.5,
                             b = 0.5,
                             l = 0.5,
                             unit = "cm"))+
  geom_rect(aes(xmin=9.5,xmax=10.5,ymin=0.9*120000,ymax=0.95*120000),
            fill=col[2],color='black',alpha=0.6)+
  annotate(geom='text',x=11.55,y=0.925*120000,label='No. of cases',size=8,family="serif"
  )+
  annotate('segment',x=9.5,xend=10.5,y=0.85*120000,yend=0.85*120000,
           linetype=1,cex=1,color=col[1])+
  annotate('text',x=10,y=0.851*120000,label='•',
           size=10,color=col[1])+
  annotate('text',x=11.5,y=0.85*120000,label='Incidence',family="serif",
           size=8)+
  ggtitle("A")

Year1<-read.xlsx("Year Agegroup.xlsx")
Year1$Year<-factor(Year1$Year)
Cumulative.cases1<-aggregate(Year1$Cases,by=list(Year1$Year),sum)
colnames(Cumulative.cases1)<-c("Year","Cumulative.cases")
Year1<-left_join(Year1,Cumulative.cases1,by="Year")
Year1$Proportion<-Year1$Cases/Year1$Cumulative.cases*100
min(Year1$Proportion)
max(Year1$Proportion)
#Heat plot
Figure2<-ggplot(Year1,aes(x=Season,y=Year))+
  geom_tile(aes(fill=Proportion))+
  scale_fill_gradientn(colors=rev(col),
                       limit=c(6,10.3),
                       breaks=c(6,8,10),labels=c("6.0","8.0","10.0"))+
  scale_y_discrete(limits=rev(c("2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019")),position = "right")+
  scale_x_discrete(limits=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),position = "bottom")+
  labs(x = '',y="",fill="")+
  ggtitle("B")+
  theme_minimal()+
  labs(fill="Proportion\n(%)")+
  theme(axis.ticks.x = element_blank(),axis.ticks.y=element_blank(),
        axis.text = element_text(colour = "black",size = 20,family ='serif'),
        axis.text.x=element_text(family ='serif'),
        legend.text = element_text(colour="black",size=20,family ='serif'),
        legend.title=element_text(size=20,color="black",family ='serif'),
        plot.title = element_text(size = 22,color = 'black',family ='serif'))

#Study population
col<-met.brewer("Troy",8)
setwd("D:/~")
Incidence<-read.xlsx("D:/Figure S3-S6 update.xlsx")
Season2<-aggregate(Incidence$Case,by=list(Incidence$Month),sum)
colnames(Season2)<-c("Order","Cases")
Season2$Season<-Season1$Season
Season2$std.Cases<-round(Season2$Cases/c(31,28.25,31,30,31,30,31,31,30,31,30,31)*30.5/15)
Population<-aggregate(Incidence$Population,by=list(Incidence$Year,Incidence$Month),sum)
Population<-aggregate(Population$x,by=list(Population$Group.2),mean)$x/c(31,28.25,31,30,31,30,31,31,30,31,30,31)*30.5/15
Season2$Population<-Population
Season2$std.Incidence<-Season2$std.Cases/Season2$Population*100000
Season2$Label<-format(Season2$std.Cases,big.mark=",")
Figure3<-ggplot()+
  geom_col(data=Season2,aes(x=Order,y=std.Cases),color='black', width=0.7, fill=col[5],alpha=0.6)+
  theme_bw()+
  scale_x_continuous(limits=c(0.5,12.5), breaks=seq(1,12,1),labels=c("Jan","Feb","Mar","Apr",
                                                                     "May","Jun","Jul","Aug",
                                                                     "Sep","Oct","Nov","Dec"))+
  geom_text(data=Season2,aes(x=Order,y=std.Cases,label=Label),vjust=-0.5,size=6.0,fontface='bold',family="serif")+
  scale_y_continuous(expand = c(0,0),limits=c(0,30000),breaks=c(0,10000,20000,30000),labels=c("0","10,000","20,000","30,000"),
                     sec.axis = sec_axis(~./200,name = 'Standardized incidence of hepatitis B (1/100,000)',
                                         breaks=c(0,75,150)))+
  geom_line(data=Season2,aes(x=Order,y=std.Incidence*200),linewidth=1,color=col[6])+
  geom_point(data=Season2,aes(x=Order,y=std.Incidence*200),size=2,color=col[6])+
  scale_linetype_manual(values=c("solid"))+
  xlab("Month")+ylab("Standardized cases of hepatitis B")+
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text = element_text(colour = "black",size = 20,family="serif"),
        axis.title = element_text(size = 20,color = 'black',family="serif"),
        plot.title=element_text(size=22,color = 'black',family ='serif'),
        plot.margin = margin(t = 1,  
                             r = 0.5,  
                             b = 0.5,  
                             l = 0.5,  
                             unit = "cm"))+
  geom_rect(aes(xmin=9.5,xmax=10.5,ymin=0.9*30000,ymax=0.95*30000),
            fill=col[5],color='black',alpha=0.6)+
  annotate(geom='text',x=11.5,y=0.925*30000,label='No. of cases',size=8,family="serif"
  )+
  annotate('segment',x=9.5,xend=10.5,y=0.85*30000,yend=0.85*30000,
           linetype=1,cex=1,color=col[6])+
  annotate('text',x=10,y=0.851*30000,label='•',
           size=10,color=col[6])+
  annotate('text',x=11.5,y=0.85*30000,label='Incidence',family="serif",
           size=8)+
  ggtitle("C")

Year2<-aggregate(Incidence$Case,by=list(Incidence$Year,Incidence$Month),sum)
colnames(Year2)<-c("Year","Month","Cases")
Year2$Season<-rep(c("Jan","Feb","Mar","Apr","May","Jun",
                    "Jul","Aug","Sep","Oct","Nov","Dec"),each=15)
Year2$Year<-factor(Year2$Year)
Cumulative.cases2<-aggregate(Year2$Cases,by=list(Year2$Year),sum)
colnames(Cumulative.cases2)<-c("Year","Cumulative.cases")
Year2<-left_join(Year2,Cumulative.cases2,by="Year")
Year2$Proportion<-Year2$Cases/Year2$Cumulative.cases*100
min(Year2$Proportion)
max(Year2$Proportion)
Figure4<-ggplot(Year2,aes(x=Season,y=Year))+
  geom_tile(aes(fill=Proportion))+
  scale_fill_gradientn(colors=rev(col),#col1
                       limit=c(6,10.5),breaks=c(6,8,10),labels=c("6.0","8.0","10.0"))+
  scale_y_discrete(limits=rev(c("2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019")),position = "right")+
  scale_x_discrete(limits=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),position = "bottom")+
  labs(x = '',y="",fill="")+
  ggtitle("D")+
  theme_minimal()+
  labs(fill="Proportion\n(%)")+
  theme(axis.ticks.x = element_blank(),axis.ticks.y=element_blank(),
        axis.text = element_text(colour = "black",size = 20,family ='serif'),
        axis.text.x=element_text(family ='serif'),
        legend.text = element_text(colour="black",size=20,family ='serif'),
        legend.title=element_text(size=20,color="black",family ='serif'),
        plot.title = element_text(size = 22,color = 'black',family ='serif'))

#Non-intervention population
col<-met.brewer("Troy",8)
setwd("D:/~")
Season3<-read.xlsx("Agegroup15-22.xlsx")
Season3$Label<-format(Season3$std.Cases,big.mark=",")
Figure5<-ggplot()+
  geom_col(data=Season3,aes(x=Order,y=std.Cases),color='black', width=0.7, fill=col[4],alpha=0.6)+
  theme_bw()+
  scale_x_continuous(limits=c(0.5,12.5), breaks=seq(1,12,1),labels=c("Jan","Feb","Mar","Apr",
                                                                     "May","Jun","Jul","Aug",
                                                                     "Sep","Oct","Nov","Dec"))+
  geom_text(data=Season3,aes(x=Order,y=std.Cases,label=Label),vjust=-0.5,size=6.0,fontface='bold',family="serif")+
  scale_y_continuous(expand = c(0,0),limits=c(0,80000),breaks=c(0,40000,80000),labels=c("0","40,000","80,000"),
                     sec.axis = sec_axis(~./200,name = 'Standardized incidence of hepatitis B (1/100,000)',
                                         breaks=c(0,200,400)))+
  geom_line(data=Season3,aes(x=Order,y=std.Incidence*225),linewidth=1,color=col[3])+
  geom_point(data=Season3,aes(x=Order,y=std.Incidence*225),size=2,color=col[3])+
  scale_linetype_manual(values=c("solid"))+
  xlab("Month")+ylab("Standardized cases of hepatitis B")+
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text = element_text(colour = "black",size = 20,family="serif"),
        axis.title = element_text(size = 20,color = 'black',family="serif"),
        plot.title=element_text(size=22,color = 'black',family ='serif'),
        plot.margin = margin(t = 1,  
                             r = 0.5,  
                             b = 0.5,  
                             l = 0.5,  
                             unit = "cm"))+
  geom_rect(aes(xmin=9.5,xmax=10.5,ymin=0.9*80000,ymax=0.95*80000),
            fill=col[4],color='black',alpha=0.6)+
  annotate(geom='text',x=11.5,y=0.925*80000,label='No. of cases',size=8,family="serif"
  )+
  annotate('segment',x=9.5,xend=10.5,y=0.85*80000,yend=0.85*80000,
           linetype=1,cex=1,color=col[3])+
  annotate('text',x=10,y=0.851*80000,label='•',
           size=10,color=col[3])+
  annotate('text',x=11.5,y=0.85*80000,label='Incidence',family="serif",
           size=8)+
  ggtitle("E")

Year3<-read.xlsx("Year Agegroup15-22.xlsx")
Year3$Year<-factor(Year3$Year)
Cumulative.cases3<-aggregate(Year3$Cases,by=list(Year3$Year),sum)
colnames(Cumulative.cases3)<-c("Year","Cumulative.cases")
Year3<-left_join(Year3,Cumulative.cases3,by="Year")
Year3$Proportion<-Year3$Cases/Year3$Cumulative.cases*100
min(Year3$Proportion)
max(Year3$Proportion)
Figure6<-ggplot(Year3,aes(x=Season,y=Year))+
  geom_tile(aes(fill=Proportion))+
  scale_fill_gradientn(colors=rev(col),#col1
                       limit=c(6,10.5),breaks=c(6,8,10),labels=c("6.0","8.0","10.0"))+
  scale_y_discrete(limits=rev(c("2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019")),position = "right")+
  scale_x_discrete(limits=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),position = "bottom")+
  labs(x = '',y="",fill="")+
  ggtitle("F")+
  theme_minimal()+
  labs(fill="Proportion\n(%)")+
  theme(axis.ticks.x = element_blank(),axis.ticks.y=element_blank(),
        axis.text = element_text(colour = "black",size = 20,family ='serif'),
        axis.text.x=element_text(family ='serif'),
        legend.text = element_text(colour="black",size=20,family ='serif'),
        legend.title=element_text(size=20,color="black",family ='serif'),
        plot.title = element_text(size = 22,color = 'black',family ='serif'))

plot1<-plot_grid(Figure1,Figure2,
                 Figure3,Figure4,
                 Figure5,Figure6,
                 ncol=2)#,align = "vh"
setwd("D:/~")
tiff(file="Figure S4.tiff",width=660,height=500,units="mm",res=300,compression="lzw")
print(plot1)
dev.off()

png(file="Figure S4.png",width=660,height=500,units="mm",res=300)
print(plot1)
dev.off()

#######################Seasonality of the incidence with wavelet analyses (Figure S5, S6)########################################
setwd("D:/~")
Incidence.study.population<-Incidence.month[which(Incidence.month$Group=="Study population"),]
Incidence.study.population<-arrange(Incidence.study.population,Time)
Incidence.study.population<-Incidence.study.population[,c(1,2)]
#Wavelet analysis
index.ticks<-seq(1,180,by=24)
index.labels<-c("2005-01","2007-01","2009-01","2011-01","2013-01",
                "2015-01","2017-01","2019-01")
wavelet_result<-analyze.wavelet(Incidence.study.population,my.series=2,loess.span=0)
png(filename="Local study population.png",width=175,height=140,units="mm",res = 300)
par(family="serif")
wt.image(wavelet_result,main="Study population",#color.palette="heat.colors",
         legend.params=list(lab="Power levels",label.digits=3),
         periodlab="Period (months)",
         timelab="Time elapsed (months)",
         spec.time.axis=list(at=index.ticks,labels=index.labels),
         siglvl=0.05,plot.ridge=FALSE)

png(filename="Global study population.png",width=175,height=140,units="mm",res = 300)
par(family="serif")
wavelet.plot<-wt.avg(wavelet_result,main="Study population",#color.palette="heat.colors",
                     periodlab="Period (months)",
                     averagelab="Average wavelet power")
dev.off()

#Composite figure
setwd("D:/~")
Fig1<-image_read("Local 1.png")
Fig2<-image_read("Local study population.png")
Fig3<-image_read("Local 4.png")
Fig4<-image_read("Local 6.png")
Fig5<-image_read("Local 7.png")
Fig6<-image_read("Local 8.png")
Fig7<-image_read("Local 9.png")
Fig8<-image_read("Local 10.png")
Fig9<-image_read("Local 11.png")
Fig10<-image_read("Local 12.png")
Fig11<-image_read("Local 13.png")
Fig12<-image_read("Local 14.png")
First.row<-image_append(c(Fig1,Fig2,Fig3,Fig4),stack=F)
Second.row<-image_append(c(Fig5,Fig6,Fig7,Fig8),stack=F)
Third.row<-image_append(c(Fig9,Fig10,Fig11,Fig12),stack=F)
Final.image<-image_append(c(First.row,Second.row,Third.row),stack=T)
image_write(Final.image,"Fig S5.png")

setwd("D:/~")
Fig1<-image_read("Global 1.png")
Fig2<-image_read("Global study population.png")
Fig3<-image_read("Global 4.png")
Fig4<-image_read("Global 6.png")
Fig5<-image_read("Global 7.png")
Fig6<-image_read("Global 8.png")
Fig7<-image_read("Global 9.png")
Fig8<-image_read("Global 10.png")
Fig9<-image_read("Global 11.png")
Fig10<-image_read("Global 12.png")
Fig11<-image_read("Global 13.png")
Fig12<-image_read("Global 14.png")
First.row<-image_append(c(Fig1,Fig2,Fig3,Fig4),stack=F)
Second.row<-image_append(c(Fig5,Fig6,Fig7,Fig8),stack=F)
Third.row<-image_append(c(Fig9,Fig10,Fig11,Fig12),stack=F)
Final.image<-image_append(c(First.row,Second.row,Third.row),stack=T)
image_write(Final.image,"Fig S6.png")

###########################Map of the incidence of study population (Figure 2)###############################
Agegroup<-c("5-year old children","6-year old children",
            "7-year old children","8-year old children",
            "9-year old children","Children aged 10–14",
            "Population aged 15–19","Population aged 20–24",
            "Population aged 25–29")
i<-6
Agegroup.file<-paste0("agegroup",i,".xlsx")
Agegroup.sub<-read.xlsx(Agegroup.file)
Incidence.year.province<-aggregate(Agegroup.sub$Case,by=list(Agegroup.sub$Year,Agegroup.sub$Province),sum)
colnames(Incidence.year.province)<-c("Year","Province","Case")
Incidence.year.province$Population<-aggregate(Agegroup.sub$Population.book,by=list(Agegroup.sub$Year,Agegroup.sub$Province),mean)$x
Incidence.year.province$Incidence<-Incidence.year.province$Case/Incidence.year.province$Population*100000
Incidence.year.province$Group<-Agegroup[i-5]
Incidence.year.province<-arrange(Incidence.year.province,Year,Province)
for(i in 7:14)
{Agegroup.file<-paste0("agegroup",i,".xlsx")
Agegroup.sub<-read.xlsx(Agegroup.file)
Incidence<-aggregate(Agegroup.sub$Case,by=list(Agegroup.sub$Year,Agegroup.sub$Province),sum)
colnames(Incidence)<-c("Year","Province","Case")
Incidence$Population<-aggregate(Agegroup.sub$Population.book,by=list(Agegroup.sub$Year,Agegroup.sub$Province),mean)$x
Incidence$Incidence<-Incidence$Case/Incidence$Population*100000
Incidence$Group<-Agegroup[i-5]
Incidence<-arrange(Incidence,Year,Province)
Incidence.year.province<-rbind(Incidence.year.province,Incidence)
}
Incidence.year.province<-arrange(Incidence.year.province,Group,Province,Year)

setwd("D:/~")
#Shape file of China
China<-read_sf("province_2004.shp")
China<-China[-which(is.na(China$NAME2004)),]
ggplot()+geom_sf(data=China)
#Boundary
Boundary<-read_sf("国界线.shp")
Boundary<-st_transform(Boundary,CRS("+proj=longlat +datum=WGS84 "))
ggplot()+geom_sf(data=China)+geom_sf(data=Boundary)

#2005
Incidence<-Incidence.year.province%>%
  filter(Year==2005)
Incidence.2005<-aggregate(Incidence$Case,by=list(Incidence$Year,Incidence$Province),sum)
colnames(Incidence.2005)<-c("Year","Province","Case")
Incidence.2005$Population<-aggregate(Incidence$Population,by=list(Incidence$Year,Incidence$Province),sum)$x
Incidence.2005$Incidence<-Incidence.2005$Case/Incidence.2005$Population*100000
Incidence.2005$Province.name<-c("安徽","北京","重庆","福建","甘肃",
                                "广东","广西","贵州","海南","河北",
                                "黑龙江","河南","湖北","湖南","内蒙古",
                                "江苏","江西","吉林","辽宁","宁夏",
                                "青海","陕西","山东","上海","山西",
                                "四川","天津","西藏","新疆","云南",
                                "浙江")
China.2005<-left_join(China,Incidence.2005,by=c("NAME2004"="Province.name"))

summary(Incidence.2005$Incidence)
#Range of the shape file
China.range<-st_bbox(China)
#Representative fraction
rangelongchina<-China.range[3]-China.range[1]
lengthFor600km<-600*rangelongchina/5200
startlong<-85 
startlat<-15 
heightlat<-0.8 
scalebarsub1long<-c(0,0,lengthFor600km,lengthFor600km)+startlong
scalebarsub1lat<-c(0,heightlat,heightlat,0)+startlat
scalebar<-data.frame(long=scalebarsub1long,lat=scalebarsub1lat)
Scalebar.black<-st_polygon(list(matrix(c(scalebar$long[1],scalebar$lat[1],
                                         scalebar$long[2],scalebar$lat[2],
                                         scalebar$long[3],scalebar$lat[3],
                                         scalebar$long[4],scalebar$lat[4],
                                         scalebar$long[1],scalebar$lat[1]),
                                       ncol=2,byrow=TRUE)))
Scalebar.black<-st_sfc(Scalebar.black,crs="+proj=longlat +datum=WGS84 ")
Scalebar.white<-st_polygon(list(matrix(c(scalebar$long[1]+lengthFor600km,scalebar$lat[1],
                                         scalebar$long[2]+lengthFor600km,scalebar$lat[2],
                                         scalebar$long[3]+lengthFor600km,scalebar$lat[3],
                                         scalebar$long[4]+lengthFor600km,scalebar$lat[4],
                                         scalebar$long[1]+lengthFor600km,scalebar$lat[1]),
                                       ncol=2,byrow=TRUE)))
Scalebar.white<-st_sfc(Scalebar.white,crs="+proj=longlat +datum=WGS84 ")
#Compass
scalenorth<-0.7 
startlongn<-85
startlatn<-53
scalenorthpart1<-data.frame(long=c(startlongn,startlongn,startlongn-2*scalenorth),
                            lat=c(startlatn,startlatn+2*scalenorth,startlatn-1*scalenorth))
scalenorthpart2<-data.frame(long=c(startlongn,startlongn,startlongn+2*scalenorth),
                            lat=c(startlatn,startlatn+2*scalenorth,startlatn-1*scalenorth))
Scalenorth.black<-st_polygon(list(matrix(c(scalenorthpart1$long[1],scalenorthpart1$lat[1],
                                           scalenorthpart1$long[2],scalenorthpart1$lat[2],
                                           scalenorthpart1$long[3],scalenorthpart1$lat[3],
                                           scalenorthpart1$long[1],scalenorthpart1$lat[1]),
                                         ncol=2,byrow=TRUE)))
Scalenorth.black<-st_sfc(Scalenorth.black,crs="+proj=longlat +datum=WGS84 ")
Scalenorth.white<-st_polygon(list(matrix(c(scalenorthpart2$long[1],scalenorthpart2$lat[1],
                                           scalenorthpart2$long[2],scalenorthpart2$lat[2],
                                           scalenorthpart2$long[3],scalenorthpart2$lat[3],
                                           scalenorthpart2$long[1],scalenorthpart2$lat[1]),
                                         ncol=2,byrow=TRUE)))
Scalenorth.white<-st_sfc(Scalenorth.white,crs="+proj=longlat +datum=WGS84 ")


plot.2005<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2005,aes(fill=Incidence),colour="black", size=.5)+
  theme_bw()+
  scale_fill_gradientn(colours=met.brewer("Hiroshige",30)[16:30],
                       limit=c(0,242),breaks=c(0,60,120,180,240),na.value = "white")+
  theme_bw()+
  xlab("Longitude")+
  ylab("Latitude")+
  scale_x_continuous(breaks=c(80,90,100,110,120,130),labels=c("80°E","90°E","100°E","110°E","120°E","130°E"))+
  scale_y_continuous(limit=c(13,55),breaks=c(20,30,40,50),labels=c("20°N","30°N","40°N","50°N"))+
  labs(fill="Incidence\n(1/100,000)")+
  ggtitle(" \nA\n                                       2005")+
  geom_sf(data=Scalebar.black,fill="black",colour="black",size=0.3)+
  geom_sf(data=Scalebar.white,fill="white",colour="black",size=0.3)+
  annotate("text",x=c(startlong,startlong+lengthFor600km,startlong+2*lengthFor600km), 
           y=c(startlat, startlat, startlat)-1,label=c("0","600","1200 km"),family="serif")+
  geom_sf(data=Scalenorth.black,fill="black",color="black",size=0.3)+
  geom_sf(data=Scalenorth.white,fill="white",color="black",size=0.3)+
  annotate("text",x=c(startlongn),y=c(startlatn-1.5),
           label=c("N"),size=4,family='serif')+
  theme(axis.title=element_text(size=15,family ='serif'), #family="RMN",
        axis.text=element_text(size=12,face="plain",color="black",family ='serif'), #family="RMN",
        title=element_text(size=16,family ='serif'), #family="RMN",
        legend.position=c(0.14,0.2),
        legend.text=element_text(size=12,family ='serif'), #family="RMN",
        legend.title=element_text(size=13,color="black",family ='serif'), #family="RMN",
        legend.key.height=unit(4, "mm"),
        legend.key.width=unit(4, "mm"),
        legend.background=element_blank())


nineline.2005<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2005,aes(fill=Incidence),colour="black", size=.5)+
  theme_bw()+
  scale_fill_gradientn(colours=met.brewer("Hiroshige",30)[16:30],
                       limit=c(0,242),breaks=c(0,60,120,180,240),na.value = "white")+        theme_bw()+
  ylim(5,22)+
  xlim(107,122)+
  guides(fill="none")+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"))

Fig.2005<-ggdraw()+
  draw_plot(plot.2005)+
  draw_plot(nineline.2005,x=0.81,y=0.05,width=0.13,height=0.39)

#2012
Incidence<-Incidence.year.province%>%
  filter(Year==2012)
Incidence.2012<-aggregate(Incidence$Case,by=list(Incidence$Year,Incidence$Province),sum)
colnames(Incidence.2012)<-c("Year","Province","Case")
Incidence.2012$Population<-aggregate(Incidence$Population,by=list(Incidence$Year,Incidence$Province),sum)$x
Incidence.2012$Incidence<-Incidence.2012$Case/Incidence.2012$Population*100000
Incidence.2012$Province.name<-c("安徽","北京","重庆","福建","甘肃",
                                "广东","广西","贵州","海南","河北",
                                "黑龙江","河南","湖北","湖南","内蒙古",
                                "江苏","江西","吉林","辽宁","宁夏",
                                "青海","陕西","山东","上海","山西",
                                "四川","天津","西藏","新疆","云南",
                                "浙江")
China.2012<-left_join(China,Incidence.2012,by=c("NAME2004"="Province.name"))
summary(Incidence.2012$Incidence)

plot.2012<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2012,aes(fill=Incidence),colour="black", size=.5)+
  theme_bw()+
  scale_fill_gradientn(colours=met.brewer("Hiroshige",30)[16:30],
                       limit=c(0,242),breaks=c(0,60,120,180,240),na.value = "white")+        theme_bw()+
  xlab("Longitude")+
  ylab("Latitude")+
  scale_x_continuous(breaks=c(80,90,100,110,120,130),labels=c("80°E","90°E","100°E","110°E","120°E","130°E"))+
  scale_y_continuous(limit=c(13,55),breaks=c(20,30,40,50),labels=c("20°N","30°N","40°N","50°N"))+
  labs(fill="Incidence\n(1/100,000)")+
  ggtitle(" \nB\n                                       2012")+
  geom_sf(data=Scalebar.black,fill="black",colour="black",size=0.3)+
  geom_sf(data=Scalebar.white,fill="white",colour="black",size=0.3)+
  annotate("text",x=c(startlong,startlong+lengthFor600km,startlong+2*lengthFor600km), 
           y=c(startlat, startlat, startlat)-1,label=c("0","600","1200 km"),family="serif")+
  geom_sf(data=Scalenorth.black,fill="black",color="black",size=0.3)+
  geom_sf(data=Scalenorth.white,fill="white",color="black",size=0.3)+
  annotate("text",x=c(startlongn),y=c(startlatn-1.5),
           label=c("N"),size=4,family='serif')+
  theme(axis.title=element_text(size=15,family ='serif'), #family="RMN",
        axis.text=element_text(size=12,face="plain",color="black",family ='serif'), #family="RMN",
        title=element_text(size=16,family ='serif'), #family="RMN",
        legend.position=c(0.14,0.2),
        legend.text=element_text(size=12,family ='serif'), #family="RMN",
        legend.title=element_text(size=13,color="black",family ='serif'), #family="RMN",
        legend.key.height=unit(4, "mm"),
        legend.key.width=unit(4, "mm"),
        legend.background=element_blank())


nineline.2012<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2012,aes(fill=Incidence),colour="black", size=.5)+
  theme_bw()+
  scale_fill_gradientn(colours=met.brewer("Hiroshige",30)[16:30],
                       limit=c(0,242),breaks=c(0,60,120,180,240),na.value = "white")+        theme_bw()+
  ylim(5,22)+
  xlim(107,122)+
  guides(fill="none")+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"))

Fig.2012<-ggdraw()+
  draw_plot(plot.2012)+
  draw_plot(nineline.2012,x=0.81,y=0.05,width=0.13,height=0.39)


#2019
Incidence<-Incidence.year.province%>%
  filter(Year==2019)
Incidence.2019<-aggregate(Incidence$Case,by=list(Incidence$Year,Incidence$Province),sum)
colnames(Incidence.2019)<-c("Year","Province","Case")
Incidence.2019$Population<-aggregate(Incidence$Population,by=list(Incidence$Year,Incidence$Province),sum)$x
Incidence.2019$Incidence<-Incidence.2019$Case/Incidence.2019$Population*100000
Incidence.2019$Province.name<-c("安徽","北京","重庆","福建","甘肃",
                                "广东","广西","贵州","海南","河北",
                                "黑龙江","河南","湖北","湖南","内蒙古",
                                "江苏","江西","吉林","辽宁","宁夏",
                                "青海","陕西","山东","上海","山西",
                                "四川","天津","西藏","新疆","云南",
                                "浙江")
China.2019<-left_join(China,Incidence.2019,by=c("NAME2004"="Province.name"))
summary(Incidence.2019$Incidence)

plot.2019<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2019,aes(fill=Incidence),colour="black", size=.5)+
  theme_bw()+
  scale_fill_gradientn(colours=met.brewer("Hiroshige",30)[16:30],
                       limit=c(0,242),breaks=c(0,60,120,180,240),na.value = "white")+        theme_bw()+
  xlab("Longitude")+
  ylab("Latitude")+
  scale_x_continuous(breaks=c(80,90,100,110,120,130),labels=c("80°E","90°E","100°E","110°E","120°E","130°E"))+
  scale_y_continuous(limit=c(13,55),breaks=c(20,30,40,50),labels=c("20°N","30°N","40°N","50°N"))+
  labs(fill="Incidence\n(1/100,000)")+
  ggtitle(" \nC\n                                       2019")+
  geom_sf(data=Scalebar.black,fill="black",colour="black",size=0.3)+
  geom_sf(data=Scalebar.white,fill="white",colour="black",size=0.3)+
  annotate("text",x=c(startlong,startlong+lengthFor600km,startlong+2*lengthFor600km), 
           y=c(startlat, startlat, startlat)-1,label=c("0","600","1200 km"),family="serif")+
  geom_sf(data=Scalenorth.black,fill="black",color="black",size=0.3)+
  geom_sf(data=Scalenorth.white,fill="white",color="black",size=0.3)+
  annotate("text",x=c(startlongn),y=c(startlatn-1.5),
           label=c("N"),size=4,family='serif')+
  theme(axis.title=element_text(size=15,family ='serif'), #family="RMN",
        axis.text=element_text(size=12,face="plain",color="black",family ='serif'), #family="RMN",
        title=element_text(size=16,family ='serif'), #family="RMN",
        legend.position=c(0.14,0.2),
        legend.text=element_text(size=12,family ='serif'), #family="RMN",
        legend.title=element_text(size=13,color="black",family ='serif'), #family="RMN",
        legend.key.height=unit(4, "mm"),
        legend.key.width=unit(4, "mm"),
        legend.background=element_blank())


nineline.2019<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2019,aes(fill=Incidence),colour="black", size=.5)+
  theme_bw()+
  scale_fill_gradientn(colours=met.brewer("Hiroshige",30)[16:30],
                       limit=c(0,242),breaks=c(0,60,120,180,240),na.value = "white")+        theme_bw()+
  ylim(5,22)+
  xlim(107,122)+
  guides(fill="none")+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"))

Fig.2019<-ggdraw()+
  draw_plot(plot.2019)+
  draw_plot(nineline.2019,x=0.81,y=0.05,width=0.13,height=0.39)

#The average incidence
Incidence<-Incidence.year.province
Case.preparation<-aggregate(Incidence$Case,by=list(Incidence$Year,Incidence$Province),sum)
colnames(Case.preparation)<-c("Year","Province","Case")
Incidence.average<-aggregate(Case.preparation$Case,by=list(Case.preparation$Province),mean)
colnames(Incidence.average)<-c("Province","Case")
Population.preparation<-aggregate(Incidence$Population,by=list(Incidence$Year,Incidence$Province),sum)
colnames(Population.preparation)<-c("Year","Province","Population")
Incidence.average$Population<-aggregate(Population.preparation$Population,by=list(Population.preparation$Province),mean)$x
Incidence.average$Incidence<-Incidence.average$Case/Incidence.average$Population*100000
Incidence.average$Province.name<-c("安徽","北京","重庆","福建","甘肃",
                                   "广东","广西","贵州","海南","河北",
                                   "黑龙江","河南","湖北","湖南","内蒙古",
                                   "江苏","江西","吉林","辽宁","宁夏",
                                   "青海","陕西","山东","上海","山西",
                                   "四川","天津","西藏","新疆","云南",
                                   "浙江")
China.average<-left_join(China,Incidence.average,by=c("NAME2004"="Province.name"))
summary(Incidence.average$Incidence)

plot.average<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.average,aes(fill=Incidence),colour="black", size=.5)+
  theme_bw()+
  scale_fill_gradientn(colours=met.brewer("Hiroshige",30)[16:30],
                       limit=c(0,242),breaks=c(0,60,120,180,240),na.value = "white")+        theme_bw()+
  xlab("Longitude")+
  ylab("Latitude")+
  scale_x_continuous(breaks=c(80,90,100,110,120,130),labels=c("80°E","90°E","100°E","110°E","120°E","130°E"))+
  scale_y_continuous(limit=c(13,55),breaks=c(20,30,40,50),labels=c("20°N","30°N","40°N","50°N"))+
  labs(fill="Incidence\n(1/100,000)")+
  ggtitle(" \nD\n                                   Average")+
  geom_sf(data=Scalebar.black,fill="black",colour="black",size=0.3)+
  geom_sf(data=Scalebar.white,fill="white",colour="black",size=0.3)+
  annotate("text",x=c(startlong,startlong+lengthFor600km,startlong+2*lengthFor600km), 
           y=c(startlat, startlat, startlat)-1,label=c("0","600","1200 km"),family="serif")+
  geom_sf(data=Scalenorth.black,fill="black",color="black",size=0.3)+
  geom_sf(data=Scalenorth.white,fill="white",color="black",size=0.3)+
  annotate("text",x=c(startlongn),y=c(startlatn-1.5),
           label=c("N"),size=4,family='serif')+
  theme(axis.title=element_text(size=15,family ='serif'), #family="RMN",
        axis.text=element_text(size=12,face="plain",color="black",family ='serif'), #family="RMN",
        title=element_text(size=16,family ='serif'), #family="RMN",
        legend.position=c(0.14,0.2),
        legend.text=element_text(size=12,family ='serif'), #family="RMN",
        legend.title=element_text(size=13,color="black",family ='serif'), #family="RMN",
        legend.key.height=unit(4, "mm"),
        legend.key.width=unit(4, "mm"),
        legend.background=element_blank())

nineline.average<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.average,aes(fill=Incidence),colour="black", size=.5)+
  theme_bw()+
  scale_fill_gradientn(colours=met.brewer("Hiroshige",30)[16:30],
                       limit=c(0,242),breaks=c(0,60,120,180,240),na.value = "white")+        theme_bw()+
  ylim(5,22)+
  xlim(107,122)+
  guides(fill="none")+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"))

Fig.average<-ggdraw()+
  draw_plot(plot.average)+
  draw_plot(nineline.average,x=0.81,y=0.05,width=0.13,height=0.39)

Fig.study<-plot_grid(Fig.2005,Fig.2012,Fig.2019,Fig.average,
                     ncol=2,align = "vh")
tiff(filename="Figure 1.tiff", width = 330, height = 330,units = "mm", res = 300, compression = "lzw")
print(Fig.study)
dev.off()
png(filename="Figure 1.png", width = 330, height = 330,units = "mm", res = 300)
print(Fig.study)
dev.off()
#######################Four measures of intervention-related effects (Figure 3)######################################
#Slope changes
setwd("D:/~")
#Shape file of China
China<-read_sf("province_2004.shp")
China<-China[-which(is.na(China$NAME2004)),]
ggplot()+geom_sf(data=China)
#Boundary of China
Boundary<-read_sf("国界线.shp")
Boundary<-st_transform(Boundary,CRS("+proj=longlat +datum=WGS84 "))
ggplot()+geom_sf(data=China)+geom_sf(data=Boundary)

#2002
setwd("D:/~")
Map.2002<-read.xlsx("2002年效应值地图数据 update.xlsx")
China.2002<-left_join(China,Map.2002,by=c("NAME2004"="province"))

#The range of the shape file
China.range<-st_bbox(China)
#Representative Fraction
rangelongchina<-China.range[3]-China.range[1]
lengthFor600km<-600*rangelongchina/5200
startlong<-85 
startlat<-15 
heightlat<-0.8 
scalebarsub1long<-c(0,0,lengthFor600km,lengthFor600km)+startlong
scalebarsub1lat<-c(0,heightlat,heightlat,0)+startlat
scalebar<-data.frame(long=scalebarsub1long,lat=scalebarsub1lat)
Scalebar.black<-st_polygon(list(matrix(c(scalebar$long[1],scalebar$lat[1],
                                         scalebar$long[2],scalebar$lat[2],
                                         scalebar$long[3],scalebar$lat[3],
                                         scalebar$long[4],scalebar$lat[4],
                                         scalebar$long[1],scalebar$lat[1]),
                                       ncol=2,byrow=TRUE)))
Scalebar.black<-st_sfc(Scalebar.black,crs="+proj=longlat +datum=WGS84 ")
Scalebar.white<-st_polygon(list(matrix(c(scalebar$long[1]+lengthFor600km,scalebar$lat[1],
                                         scalebar$long[2]+lengthFor600km,scalebar$lat[2],
                                         scalebar$long[3]+lengthFor600km,scalebar$lat[3],
                                         scalebar$long[4]+lengthFor600km,scalebar$lat[4],
                                         scalebar$long[1]+lengthFor600km,scalebar$lat[1]),
                                       ncol=2,byrow=TRUE)))
Scalebar.white<-st_sfc(Scalebar.white,crs="+proj=longlat +datum=WGS84 ")
#Compass
scalenorth<-0.7 
startlongn<-85
startlatn<-53
scalenorthpart1<-data.frame(long=c(startlongn,startlongn,startlongn-2*scalenorth),
                            lat=c(startlatn,startlatn+2*scalenorth,startlatn-1*scalenorth))
scalenorthpart2<-data.frame(long=c(startlongn,startlongn,startlongn+2*scalenorth),
                            lat=c(startlatn,startlatn+2*scalenorth,startlatn-1*scalenorth))
Scalenorth.black<-st_polygon(list(matrix(c(scalenorthpart1$long[1],scalenorthpart1$lat[1],
                                           scalenorthpart1$long[2],scalenorthpart1$lat[2],
                                           scalenorthpart1$long[3],scalenorthpart1$lat[3],
                                           scalenorthpart1$long[1],scalenorthpart1$lat[1]),
                                         ncol=2,byrow=TRUE)))
Scalenorth.black<-st_sfc(Scalenorth.black,crs="+proj=longlat +datum=WGS84 ")
Scalenorth.white<-st_polygon(list(matrix(c(scalenorthpart2$long[1],scalenorthpart2$lat[1],
                                           scalenorthpart2$long[2],scalenorthpart2$lat[2],
                                           scalenorthpart2$long[3],scalenorthpart2$lat[3],
                                           scalenorthpart2$long[1],scalenorthpart2$lat[1]),
                                         ncol=2,byrow=TRUE)))
Scalenorth.white<-st_sfc(Scalenorth.white,crs="+proj=longlat +datum=WGS84 ")

#Calculating centroids
dataprovincename<-read.xlsx("省份名称标注经纬度.xlsx")
dataprovincename$long<-sapply(dataprovincename$经纬度, function(x){strsplit(x, split = ",")[[1]][1]}) %>% as.numeric()
dataprovincename$lat<-sapply(dataprovincename$经纬度, function(x){strsplit(x, split = ",")[[1]][2]}) %>% as.numeric()
dataprovincename<-dataprovincename[!is.na(dataprovincename$long),]
dataprovincename<-as.data.frame(dataprovincename)
dataprovincename$Order<-c(12,34,1,13,28,19,20,24,21,3,16,8,17,18,7,10,14,6,5,30,
                          29,15,4,27,9,23,32,2,26,33,31,25,11,22)
dataprovincename<-arrange(dataprovincename,Order)
dataprovincename<-dataprovincename[,c(2:4)]
colnames(dataprovincename)<-c("province","longitude","latitude")
China.2002<-left_join(China.2002,dataprovincename,by=c("NAME2004"="province"))

Change.2002<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2002,aes(fill=Change),colour="black", size=.5)+
  theme_bw()+
  scale_fill_manual(values=brewer.pal(11,"RdBu")[c(9,8,3,4)],
                    breaks=c("Decrease faster","Decrease slower","Increase faster","Increase slower"),
                    labels=c("Decrease faster","Decrease slower","Increase faster","Increase slower"),na.value = "white")+
  xlab("Longitude")+
  ylab("Latitude")+
  scale_x_continuous(breaks=c(80,90,100,110,120,130),labels=c("80°E","90°E","100°E","110°E","120°E","130°E"))+
  scale_y_continuous(limit=c(13,55),breaks=c(20,30,40,50),labels=c("20°N","30°N","40°N","50°N"))+
  labs(fill="Trend change")+
  ggtitle("  \nA\n                                NIP")+
  geom_sf(data=Scalebar.black,fill="black",colour="black",size=0.3)+
  geom_sf(data=Scalebar.white,fill="white",colour="black",size=0.3)+
  annotate("text",x=c(startlong,startlong+lengthFor600km,startlong+2*lengthFor600km), 
           y=c(startlat, startlat, startlat)-1,label=c("0","600","1200 km"),family="serif")+
  geom_sf(data=Scalenorth.black,fill="black",color="black",size=0.3)+
  geom_sf(data=Scalenorth.white,fill="white",color="black",size=0.3)+
  annotate("text",x=c(startlongn),y=c(startlatn-1.5),
           label=c("N"),size=4,family='serif')+
  geom_text(data=China.2002,
            aes(x=longitude,y=latitude,label=Symbol.change),
            family="serif",size=4)+
  theme(axis.title=element_text(size=17,family ='serif'), #family="RMN",
        axis.text=element_text(size=15,face="plain",color="black",family ='serif'), #family="RMN",
        title=element_text(size=19,family ='serif'), #family="RMN",
        legend.position=c(0.16,0.22),
        legend.text=element_text(size=14,family ='serif'), #family="RMN",
        legend.title=element_text(size=15,color="black",family ='serif'), #family="RMN",
        legend.key.height=unit(4, "mm"),
        legend.key.width=unit(4, "mm"),
        legend.background=element_blank())

nineline.2002<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2002,aes(fill=Change),colour="black", size=.5)+
  theme_bw()+
  scale_fill_manual(values=brewer.pal(11,"RdBu")[c(9,8,3,4)],
                    breaks=c("Decrease faster","Decrease slower","Increase faster","Increase slower"),
                    labels=c("Decrease faster","Decrease slower","Increase faster","Increase slower"),na.value = "white")+
  ylim(5,22)+
  xlim(107,122)+
  guides(fill="none")+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"))

Change.2002<-ggdraw()+
  draw_plot(Change.2002)+
  draw_plot(nineline.2002,x=0.8,y=0.1,width=0.13,height=0.39)

#2009
setwd("D:/~")
Map.2009<-read.xlsx("2009年效应值地图数据 update.xlsx")
China.2009<-left_join(China,Map.2009,by=c("NAME2004"="province"))
China.2009<-left_join(China.2009,dataprovincename,by=c("NAME2004"="province"))
Change.2009<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2009,aes(fill=Change),colour="black", size=.5)+
  theme_bw()+
  scale_fill_manual(values=brewer.pal(11,"RdBu")[c(9,8,3,4)],
                    breaks=c("Decrease faster","Decrease slower","Increase faster","Increase slower"),
                    labels=c("Decrease faster","Decrease slower","Increase faster","Increase slower"),na.value = "white")+
  xlab("Longitude")+
  ylab("Latitude")+
  scale_x_continuous(breaks=c(80,90,100,110,120,130),labels=c("80°E","90°E","100°E","110°E","120°E","130°E"))+
  scale_y_continuous(limit=c(13,55),breaks=c(20,30,40,50),labels=c("20°N","30°N","40°N","50°N"))+
  labs(fill="Trend change")+
  ggtitle(" \n  \n                    Catch-up vaccination")+
  geom_sf(data=Scalebar.black,fill="black",colour="black",size=0.3)+
  geom_sf(data=Scalebar.white,fill="white",colour="black",size=0.3)+
  annotate("text",x=c(startlong,startlong+lengthFor600km,startlong+2*lengthFor600km), 
           y=c(startlat, startlat, startlat)-1,label=c("0","600","1200 km"),family="serif")+
  geom_sf(data=Scalenorth.black,fill="black",color="black",size=0.3)+
  geom_sf(data=Scalenorth.white,fill="white",color="black",size=0.3)+
  annotate("text",x=c(startlongn),y=c(startlatn-1.5),
           label=c("N"),size=4,family='serif')+
  geom_text(data=China.2009,
            aes(x=longitude,y=latitude,label=Symbol.change),
            family="serif",size=4)+
  theme(axis.title=element_text(size=17,family ='serif'), #family="RMN",
        axis.text=element_text(size=15,face="plain",color="black",family ='serif'), #family="RMN",
        title=element_text(size=19,family ='serif'), #family="RMN",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        legend.position=c(0.16,0.22),
        legend.text=element_text(size=14,family ='serif'), #family="RMN",
        legend.title=element_text(size=15,color="black",family ='serif'), #family="RMN",
        legend.key.height=unit(4, "mm"),
        legend.key.width=unit(4, "mm"),
        legend.background=element_blank())

nineline.2009<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2009,aes(fill=Change),colour="black", size=.5)+
  theme_bw()+
  scale_fill_manual(values=brewer.pal(11,"RdBu")[c(9,8,3,4)],
                    breaks=c("Decrease faster","Decrease slower","Increase faster","Increase slower"),
                    labels=c("Decrease faster","Decrease slower","Increase faster","Increase slower"),na.value = "white")+
  ylim(5,22)+
  xlim(107,122)+
  guides(fill="none")+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"))

Change.2009<-ggdraw()+
  draw_plot(Change.2009)+
  draw_plot(nineline.2009,x=0.75,y=0.1,width=0.13,height=0.39)

#ER
China.2002$ER[which(China.2002$ER>100)]<-100
ER.2002<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2002,aes(fill=ER),colour="black", size=.5)+
  theme_bw()+
  scale_fill_gradientn(colours=rev(met.brewer("Hokusai2",31)),
                       limit=c(-100,100),breaks=c(-100,-50,0,50,100),
                       labels=c("-100","-50","0","50","100+"),na.value = "white")+
  xlab("Longitude")+
  ylab("Latitude")+
  scale_x_continuous(breaks=c(80,90,100,110,120,130),labels=c("80°E","90°E","100°E","110°E","120°E","130°E"))+
  scale_y_continuous(limit=c(13,55),breaks=c(20,30,40,50),labels=c("20°N","30°N","40°N","50°N"))+
  labs(fill="ER (%)")+
  ggtitle(" \nB\n                                NIP")+
  geom_sf(data=Scalebar.black,fill="black",colour="black",size=0.3)+
  geom_sf(data=Scalebar.white,fill="white",colour="black",size=0.3)+
  annotate("text",x=c(startlong,startlong+lengthFor600km,startlong+2*lengthFor600km), 
           y=c(startlat, startlat, startlat)-1,label=c("0","600","1200 km"),family="serif")+
  geom_sf(data=Scalenorth.black,fill="black",color="black",size=0.3)+
  geom_sf(data=Scalenorth.white,fill="white",color="black",size=0.3)+
  annotate("text",x=c(startlongn),y=c(startlatn-1.5),
           label=c("N"),size=4,family='serif')+
  geom_text(data=China.2002,
            aes(x=longitude,y=latitude,label=Symbol.ER),
            family="serif",size=4)+
  theme(axis.title=element_text(size=17,family ='serif'), #family="RMN",
        axis.text=element_text(size=15,face="plain",color="black",family ='serif'), #family="RMN",
        title=element_text(size=19,family ='serif'), #family="RMN",
        legend.position=c(0.14,0.2),
        legend.text=element_text(size=14,family ='serif'), #family="RMN",
        legend.title=element_text(size=15,color="black",family ='serif'), #family="RMN",
        legend.key.height=unit(4, "mm"),
        legend.key.width=unit(4, "mm"),
        legend.background=element_blank())

nineline.2002<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2002,aes(fill=ER),colour="black", size=.5)+
  theme_bw()+
  scale_fill_gradientn(colours=rev(met.brewer("Hokusai2",31)),
                       limit=c(-100,100),breaks=c(-100,-50,0,50,100),
                       labels=c("-100","-50","0","50","100+"),na.value = "white")+
  ylim(5,22)+
  xlim(107,122)+
  guides(fill="none")+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"))

ER.2002<-ggdraw()+
  draw_plot(ER.2002)+
  draw_plot(nineline.2002,x=0.8,y=0.1,width=0.13,height=0.39)

#2009
China.2009$ER[which(China.2009$ER>0)]<-0
ER.2009<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2009,aes(fill=ER),colour="black", size=.5)+
  theme_bw()+
  scale_fill_gradientn(colours=rev(met.brewer("Hokusai2",31)),
                       limit=c(-100,0),breaks=c(-100,-75,-50,-25,0),
                       labels=c("-100","-75","-50","-25","0+"),na.value = "white")+
  xlab("Longitude")+
  ylab("Latitude")+
  scale_x_continuous(breaks=c(80,90,100,110,120,130),labels=c("80°E","90°E","100°E","110°E","120°E","130°E"))+
  scale_y_continuous(limit=c(13,55),breaks=c(20,30,40,50),labels=c("20°N","30°N","40°N","50°N"))+
  labs(fill="ER (%)")+
  ggtitle(" \n  \n                    Catch-up vaccination")+
  geom_sf(data=Scalebar.black,fill="black",colour="black",size=0.3)+
  geom_sf(data=Scalebar.white,fill="white",colour="black",size=0.3)+
  annotate("text",x=c(startlong,startlong+lengthFor600km,startlong+2*lengthFor600km), 
           y=c(startlat, startlat, startlat)-1,label=c("0","600","1200 km"),family="serif")+
  geom_sf(data=Scalenorth.black,fill="black",color="black",size=0.3)+
  geom_sf(data=Scalenorth.white,fill="white",color="black",size=0.3)+
  annotate("text",x=c(startlongn),y=c(startlatn-1.5),
           label=c("N"),size=4,family='serif')+
  geom_text(data=China.2009,
            aes(x=longitude,y=latitude,label=Symbol.ER),
            family="serif",size=4)+
  theme(axis.title=element_text(size=17,family ='serif'), #family="RMN",
        axis.text=element_text(size=15,face="plain",color="black",family ='serif'), #family="RMN",
        title=element_text(size=19,family ='serif'), #family="RMN",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        legend.position=c(0.16, 0.2),
        legend.text=element_text(size=14,family ='serif'), #family="RMN",
        legend.title=element_text(size=15,color="black",family ='serif'), #family="RMN",
        legend.key.height=unit(4, "mm"),
        legend.key.width=unit(4, "mm"),
        legend.background=element_blank())

nineline.2009<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2009,aes(fill=ER),colour="black", size=.5)+
  theme_bw()+
  scale_fill_gradientn(colours=rev(met.brewer("Hokusai2",31)),
                       limit=c(-100,0),breaks=c(-100,-75,-50,-25,0),
                       labels=c("-100","-75","-50","-25","0+"),na.value = "white")+
  ylim(5,22)+
  xlim(107,122)+
  guides(fill="none")+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"))

ER.2009<-ggdraw()+
  draw_plot(ER.2009)+
  draw_plot(nineline.2009,x=0.75,y=0.1,width=0.13,height=0.39)

#EC
China.2002$EC[which(China.2002$EC>0)]<-0
EC.2002<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2002,aes(fill=EC),colour="black", size=.5)+
  theme_bw()+
  scale_fill_gradientn(colours=rev(met.brewer("Hiroshige",30)[16:30]),
                       limit=c(-80,0),breaks=c(-80,-60,-40,-20,0),
                       labels=c("-80","-60","-40","-20","0+"),na.value = "white")+
  xlab("Longitude")+
  ylab("Latitude")+
  scale_x_continuous(breaks=c(80,90,100,110,120,130),labels=c("80°E","90°E","100°E","110°E","120°E","130°E"))+
  scale_y_continuous(limit=c(13,55),breaks=c(20,30,40,50),labels=c("20°N","30°N","40°N","50°N"))+
  labs(fill="EC\n(1,000 cases)")+
  ggtitle(" \nC\n                                NIP")+
  geom_sf(data=Scalebar.black,fill="black",colour="black",size=0.3)+
  geom_sf(data=Scalebar.white,fill="white",colour="black",size=0.3)+
  annotate("text",x=c(startlong,startlong+lengthFor600km,startlong+2*lengthFor600km), 
           y=c(startlat, startlat, startlat)-1,label=c("0","600","1200 km"),family="serif")+
  geom_sf(data=Scalenorth.black,fill="black",color="black",size=0.3)+
  geom_sf(data=Scalenorth.white,fill="white",color="black",size=0.3)+
  annotate("text",x=c(startlongn),y=c(startlatn-1.5),
           label=c("N"),size=4,family='serif')+
  geom_text(data=China.2002,
            aes(x=longitude,y=latitude,label=Symbol.EC),
            family="serif",size=4)+
  theme(axis.title=element_text(size=17,family ='serif'), #family="RMN",
        axis.text=element_text(size=15,face="plain",color="black",family ='serif'), #family="RMN",
        title=element_text(size=19,family ='serif'), #family="RMN",
        legend.position=c(0.14,0.2),
        legend.text=element_text(size=14,family ='serif'), #family="RMN",
        legend.title=element_text(size=15,color="black",family ='serif'), #family="RMN",
        legend.key.height=unit(4, "mm"),
        legend.key.width=unit(4, "mm"),
        legend.background=element_blank())

nineline.2002<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2002,aes(fill=EC),colour="black", size=.5)+
  theme_bw()+
  scale_fill_gradientn(colours=rev(met.brewer("Hiroshige",30)[16:30]),
                       limit=c(-80,0),breaks=c(-80,-60,-40,-20,0),
                       labels=c("-80","-60","-40","-20","0"),na.value = "white")+
  ylim(5,22)+
  xlim(107,122)+
  guides(fill="none")+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"))

EC.2002<-ggdraw()+
  draw_plot(EC.2002)+
  draw_plot(nineline.2002,x=0.8,y=0.1,width=0.13,height=0.39)

#2009
China.2009$EC[which(China.2009$EC>0)]<-0
EC.2009<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2009,aes(fill=EC),colour="black", size=.5)+
  theme_bw()+
  scale_fill_gradientn(colours=rev(met.brewer("Hiroshige",30)[16:30]),
                       limit=c(-80,0),breaks=c(-80,-60,-40,-20,0),
                       labels=c("-80","-60","-40","-20","0+"),na.value = "white")+
  xlab("Longitude")+
  ylab("Latitude")+
  scale_x_continuous(breaks=c(80,90,100,110,120,130),labels=c("80°E","90°E","100°E","110°E","120°E","130°E"))+
  scale_y_continuous(limit=c(13,55),breaks=c(20,30,40,50),labels=c("20°N","30°N","40°N","50°N"))+
  labs(fill="EC\n(1,000 cases)")+
  ggtitle(" \n  \n                    Catch-up vaccination")+
  geom_sf(data=Scalebar.black,fill="black",colour="black",size=0.3)+
  geom_sf(data=Scalebar.white,fill="white",colour="black",size=0.3)+
  annotate("text",x=c(startlong,startlong+lengthFor600km,startlong+2*lengthFor600km), 
           y=c(startlat, startlat, startlat)-1,label=c("0","600","1200 km"),family="serif")+
  geom_sf(data=Scalenorth.black,fill="black",color="black",size=0.3)+
  geom_sf(data=Scalenorth.white,fill="white",color="black",size=0.3)+
  annotate("text",x=c(startlongn),y=c(startlatn-1.5),
           label=c("N"),size=4,family='serif')+
  geom_text(data=China.2009,
            aes(x=longitude,y=latitude,label=Symbol.EC),
            family="serif",size=4)+
  theme(axis.title=element_text(size=17,family ='serif'), #family="RMN",
        axis.text=element_text(size=15,face="plain",color="black",family ='serif'), #family="RMN",
        title=element_text(size=19,family ='serif'), #family="RMN",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        legend.position=c(0.16, 0.2),
        legend.text=element_text(size=14,family ='serif'), #family="RMN",
        legend.title=element_text(size=15,color="black",family ='serif'), #family="RMN",
        legend.key.height=unit(4, "mm"),
        legend.key.width=unit(4, "mm"),
        legend.background=element_blank())

nineline.2009<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2009,aes(fill=ER),colour="black", size=.5)+
  theme_bw()+
  scale_fill_gradientn(colours=rev(met.brewer("Hiroshige",30)[16:30]),
                       limit=c(-80,0),breaks=c(-80,-60,-40,-20,0),
                       labels=c("-80","-60","-40","-20","0+"),na.value = "white")+
  ylim(5,22)+
  xlim(107,122)+
  guides(fill="none")+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"))

EC.2009<-ggdraw()+
  draw_plot(EC.2009)+
  draw_plot(nineline.2009,x=0.75,y=0.1,width=0.13,height=0.39)
#EIR
China.2002$EIR[which(China.2002$EIR>0)]<-0
EIR.2002<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2002,aes(fill=EIR),colour="black", size=.5)+
  theme_bw()+
  scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(9,"Blues")[3:9])(30)),
                       limit=c(-320,0),breaks=c(-320,-240,-160,-80,0),
                       labels=c("-320","-240","-160","-80","0+"),na.value = "white")+
  xlab("Longitude")+
  ylab("Latitude")+
  scale_x_continuous(breaks=c(80,90,100,110,120,130),labels=c("80°E","90°E","100°E","110°E","120°E","130°E"))+
  scale_y_continuous(limit=c(13,55),breaks=c(20,30,40,50),labels=c("20°N","30°N","40°N","50°N"))+
  labs(fill="EIR\n(1/100,000)")+
  ggtitle(" \nD\n                                NIP")+
  geom_sf(data=Scalebar.black,fill="black",colour="black",size=0.3)+
  geom_sf(data=Scalebar.white,fill="white",colour="black",size=0.3)+
  annotate("text",x=c(startlong,startlong+lengthFor600km,startlong+2*lengthFor600km), 
           y=c(startlat, startlat, startlat)-1,label=c("0","600","1200 km"),family="serif")+
  geom_sf(data=Scalenorth.black,fill="black",color="black",size=0.3)+
  geom_sf(data=Scalenorth.white,fill="white",color="black",size=0.3)+
  annotate("text",x=c(startlongn),y=c(startlatn-1.5),
           label=c("N"),size=4,family='serif')+
  geom_text(data=China.2002,
            aes(x=longitude,y=latitude,label=Symbol.EIR),
            family="serif",size=4)+
  theme(axis.title=element_text(size=17,family ='serif'), #family="RMN",
        axis.text=element_text(size=15,face="plain",color="black",family ='serif'), #family="RMN",
        title=element_text(size=19,family ='serif'), #family="RMN",
        legend.position=c(0.14,0.2),
        legend.text=element_text(size=14,family ='serif'), #family="RMN",
        legend.title=element_text(size=15,color="black",family ='serif'), #family="RMN",
        legend.key.height=unit(4, "mm"),
        legend.key.width=unit(4, "mm"),
        legend.background=element_blank())

nineline.2002<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2002,aes(fill=EIR),colour="black", size=.5)+
  theme_bw()+
  scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(9,"Blues")[3:9])(30)),
                       limit=c(-320,0),breaks=c(-320,-240,-160,-80,0),
                       labels=c("-320","-240","-160","-80","0+"),na.value = "white")+
  ylim(5,22)+
  xlim(107,122)+
  guides(fill="none")+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"))

EIR.2002<-ggdraw()+
  draw_plot(EIR.2002)+
  draw_plot(nineline.2002,x=0.8,y=0.1,width=0.13,height=0.39)

#2009
China.2009$EIR[which(China.2009$EIR>0)]<-0
EIR.2009<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2009,aes(fill=EIR),colour="black", size=.5)+
  theme_bw()+
  scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(9,"Blues")[3:9])(30)),
                       limit=c(-160,0),breaks=c(-160,-120,-80,-40,0),
                       labels=c("-160","-120","-80","-40","0+"),na.value = "white")+
  xlab("Longitude")+
  ylab("Latitude")+
  scale_x_continuous(breaks=c(80,90,100,110,120,130),labels=c("80°E","90°E","100°E","110°E","120°E","130°E"))+
  scale_y_continuous(limit=c(13,55),breaks=c(20,30,40,50),labels=c("20°N","30°N","40°N","50°N"))+
  labs(fill="EIR\n(1/100,000)")+
  ggtitle(" \n  \n                    Catch-up vaccination")+
  geom_sf(data=Scalebar.black,fill="black",colour="black",size=0.3)+
  geom_sf(data=Scalebar.white,fill="white",colour="black",size=0.3)+
  annotate("text",x=c(startlong,startlong+lengthFor600km,startlong+2*lengthFor600km), 
           y=c(startlat, startlat, startlat)-1,label=c("0","600","1200 km"),family="serif")+
  geom_sf(data=Scalenorth.black,fill="black",color="black",size=0.3)+
  geom_sf(data=Scalenorth.white,fill="white",color="black",size=0.3)+
  annotate("text",x=c(startlongn),y=c(startlatn-1.5),
           label=c("N"),size=4,family='serif')+
  geom_text(data=China.2009,
            aes(x=longitude,y=latitude,label=Symbol.EIR),
            family="serif",size=4)+
  theme(axis.title=element_text(size=17,family ='serif'), #family="RMN",
        axis.text=element_text(size=15,face="plain",color="black",family ='serif'), #family="RMN",
        title=element_text(size=19,family ='serif'), #family="RMN",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        legend.position=c(0.16, 0.2),
        legend.text=element_text(size=14,family ='serif'), #family="RMN",
        legend.title=element_text(size=15,color="black",family ='serif'), #family="RMN",
        legend.key.height=unit(4, "mm"),
        legend.key.width=unit(4, "mm"),
        legend.background=element_blank())

nineline.2009<-ggplot()+
  geom_sf(data=Boundary,colour="grey10",size=0.1)+
  geom_sf(data=China.2009,aes(fill=ER),colour="black", size=.5)+
  theme_bw()+
  scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(9,"Blues")[3:9])(30)),
                       limit=c(-160,0),breaks=c(-160,-120,-80,-40,0),
                       labels=c("-160","-120","-80","-40","0+"),na.value = "white")+
  ylim(5,22)+
  xlim(107,122)+
  guides(fill="none")+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"))

EIR.2009<-ggdraw()+
  draw_plot(EIR.2009)+
  draw_plot(nineline.2009,x=0.75,y=0.1,width=0.13,height=0.39)

Fig.effect<-plot_grid(Change.2002,Change.2009,
                      ER.2002,ER.2009,
                      EC.2002,EC.2009,
                      EIR.2002,EIR.2009,
                      ncol=2,align = "vh")
tiff(filename = "Figure 3 update II.tiff", width = 360, height = 640,units = "mm", res = 300, compression = "lzw")
print(Fig.effect)
dev.off()
png(filename = "DFigure 3 update II.png", width = 360, height = 640,units = "mm", res = 300)
print(Fig.effect)
dev.off()


#######################Forest plot (Figure 4)#####################################################
setwd("D:/~")
EIR<-read.xlsx("EIR forest plot update II.xlsx")
colnames(EIR)[c(5,9,13,17)]<-c("NIP\nEC (95% eCI)",
                               "NIP\nEIR (95% eCI)",
                               "Catch-up vaccination\nEC (95% eCI)",
                               "Catch-up vaccination\nEIR (95% eCI)")
EIR$` ` <- paste(rep(" ", 20), collapse = " ")
EIR<-EIR[,c(1:5,18,6:17)]
EIR$`  ` <- paste(rep(" ", 20), collapse = " ")
EIR<-EIR[,c(1:10,19,11:18)]
EIR$`   ` <- paste(rep(" ", 20), collapse = " ")
EIR<-EIR[,c(1:15,20,16:19)]
EIR$`    ` <- paste(rep(" ", 20), collapse = " ")
EIR[c(2,10,11,12,19,23,26,29),c(5,10)]<-"  "
EIR[c(2,3,4,5,6,12,19,23,26,29),c(15,20)]<-"  "
EIR$Point1<-EIR$Point1/1000
EIR$Point3<-EIR$Point3/1000
EIR$Low1<-EIR$Low1/1000
EIR$Low3<-EIR$Low3/1000
EIR$High1<-EIR$High1/1000
EIR$High3<-EIR$High3/1000
EIR.forest<-EIR[,c(1,5,6,10,11,15,16,20,21)]

color.sub<-met.brewer("Troy",8)

tm<-forest_theme(base_size=12,
                 base_family="serif",
                 refline_lty="solid",
                 ci_pch=c(1),
                 ci_col=c(color.sub[8]),
                 ci_lw=3,
                 refline_lwd=1,
                 refline_col="grey20")
Fig.forest<-forest(EIR.forest,
                   est=list(EIR$Point1,
                            EIR$Point2,
                            EIR$Point3,
                            EIR$Point4),
                   lower=list(EIR$Low1,
                              EIR$Low2,
                              EIR$Low3,
                              EIR$Low4),
                   upper=list(EIR$High1,
                              EIR$High2,
                              EIR$High3,
                              EIR$High4),
                   #sizes=dt$se,
                   ci_column = c(3,5,7,9),
                   ref_line = 0,
                   xlim=list(c(-360,180),
                             c(-123,60),
                             c(-540,1),
                             c(-90,6)),
                   ticks_at=list(c(-360,-180,0,180),
                                 c(-120,-60,0,60),
                                 c(-540,-360,-180,0),
                                 c(-90,-60,-30,0)),
                   theme=tm)
Fig.forest<-edit_plot(Fig.forest,
                      row=c(1,2,12,19,23,26,29),
                      gp = gpar(fontface = "bold"))
Fig.forest<-add_border(Fig.forest,part="header")
Fig.forest<-edit_plot(Fig.forest,
                      row=c(1:31),
                      col=c(5),
                      which="ci",
                      gp = gpar(col=color.sub[7]))
Fig.forest<-edit_plot(Fig.forest,
                      row=c(1:31),
                      col=c(7),
                      which="ci",
                      gp = gpar(col=color.sub[6]))
Fig.forest<-edit_plot(Fig.forest,
                      row=c(1:31),
                      col=c(9),
                      which="ci",
                      gp = gpar(col=color.sub[5]))

tiff(filename="Figure 4 update II.tiff", width = 480, height = 210,units = "mm", res = 300, compression = "lzw")
print(Fig.forest)
dev.off()

png(filename="Figure 4 update II.png", width = 480, height = 210,units = "mm", res = 300)
print(Fig.forest)
dev.off()

