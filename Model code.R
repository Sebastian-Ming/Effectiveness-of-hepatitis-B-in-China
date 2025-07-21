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
#library(ggsn)
#library(rgdal)
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
#library(rgeos)
library(ggnewscale)
library(gtable)
library(ggtext)
library(ggh4x)
getwd()

setwd("D:/~")
#######################Functions################################################
n.cir=10000
#2002 intervention NIP for each PLAD
Name<-c("Beijing","Tianjin","Hebei","Shanxi","Inner Mongolia",
        "Liaoning","Jilin","Heilongjiang","Shanghai","Jiangsu","Zhejiang",
        "Anhui","Fujian","Jiangxi","Shandong","Henan","Hubei","Hunan",
        "Guangdong","Guangxi","Hainan","Chongqing","Sichuan","Guizhou",
        "Yunnan","Tibet","Shaanxi","Gansu","Qinghai","Ningxia","Xinjiang")
#"model" is the corresponding statistical model
#"dataframe" is the original data
ER.fun.2002<-function(model,dataframe){
  #Length of coefficients
  lengthfun<-length(summary(model)$p.coeff)
  #Extracting coefficients of the intervention 
  coef.number<-as.matrix(t(t(as.data.frame(summary(model)$p.coeff)[c(2,lengthfun),1])))
  #Duration of post-intervention period
  line.row<-which(dataframe$Intervention.2002==1)
  length.row<-length(line.row)
  time.number<-model.matrix(model)[line.row,c(2,lengthfun)]
  #Covariance matrix
  cov.number<-vcov(model)[c(2,lengthfun),c(2,lengthfun)]
  pre.ER<-c()
  ER<-c()
  se.ER<-c()
  ER.low<-c()
  ER.high<-c()
  #Calculating ERs
  for(k in 1:length.row)
  {pre.ER[k]<-time.number[k,]%*%coef.number
  ER[k]=(exp(pre.ER[k])-1)*100
  se.ER[k]<-sqrt(time.number[k,]%*%cov.number%*%t(t(time.number[k,])))
  ER.low[k]=(exp(pre.ER[k]-1.96*se.ER[k])-1)*100
  ER.high[k]=(exp(pre.ER[k]+1.96*se.ER[k])-1)*100
  }
  result=data.frame(ER=ER,ER.low=ER.low,ER.high=ER.high)
  result$Date<-rev(seq.Date(from=as.Date("2019/12/01",format="%Y/%m/%d"),
                            by="-1 month",length.out=length.row))
  return(result)
}

#2009 intervention catch-up vaccination for each PLAD
ER.fun.2009<-function(model,dataframe){
  lengthfun<-length(summary(model)$p.coeff)
  coef.number<-as.matrix(t(t(as.data.frame(summary(model)$p.coeff)[c(2,lengthfun),1])))
  line.row<-which(dataframe$Intervention.2009==1)
  length.row<-length(line.row)
  time.number<-model.matrix(model)[line.row,c(2,lengthfun)]
  cov.number<-vcov(model)[c(2,lengthfun),c(2,lengthfun)]
  pre.ER<-c()
  ER<-c()
  se.ER<-c()
  ER.low<-c()
  ER.high<-c()
  for(k in 1:length.row)
  {pre.ER[k]<-time.number[k,]%*%coef.number
  ER[k]=(exp(pre.ER[k])-1)*100
  se.ER[k]<-sqrt(time.number[k,]%*%cov.number%*%t(t(time.number[k,])))
  ER.low[k]=(exp(pre.ER[k]-1.96*se.ER[k])-1)*100
  ER.high[k]=(exp(pre.ER[k]+1.96*se.ER[k])-1)*100
  }
  result=data.frame(ER=ER,ER.low=ER.low,ER.high=ER.high)
  result$Date<-rev(seq.Date(from=as.Date("2019/12/01",format="%Y/%m/%d"),
                            by="-1 month",length.out=length.row))
  #list(result=result)
  return(result)
}

#2002 intervention NIP for pooled effect
#"model" is the meta analysis model for each age group
#"i" is the theoratically lasting years of the intervention effect 
ER.meta.2002<-function(model,i){
  coef.meta<-as.matrix(t(model$coefficients))
  time.meta<-as.matrix(cbind(rep(1,180-i*12),c((i*12):179)-i*12))
  cov.meta<-vcov(model)
  pre.ER.meta<-c()
  ER.meta<-c()
  se.ER.meta<-c()
  ER.low.meta<-c()
  ER.high.meta<-c()
  for(l in 1:(180-i*12))
  {pre.ER.meta[l]<-time.meta[l,]%*%coef.meta
  ER.meta[l]<-(exp(pre.ER.meta[l])-1)*100
  se.ER.meta[l]<-sqrt(time.meta[l,]%*%cov.meta%*%t(t(time.meta[l,]))) 
  ER.low.meta[l]=(exp(pre.ER.meta[l]-1.96*se.ER.meta[l])-1)*100
  ER.high.meta[l]=(exp(pre.ER.meta[l]+1.96*se.ER.meta[l])-1)*100
  }
  ER.whole<-ER.meta
  ER.whole.low<-ER.low.meta
  ER.whole.upper<-ER.high.meta
  month.meta<-c((i*12+1):180)
  ER<-data.frame(Time=month.meta,ER=ER.meta,ER.lower=ER.low.meta,ER.upper=ER.high.meta)
  ER$ER<-sprintf("%0.2f",ER$ER)
  ER$ER.lower<-sprintf("%0.2f",ER$ER.lower)
  ER$ER.upper<-sprintf("%0.2f",ER$ER.upper)
  ER$Date<-rev(seq.Date(from=as.Date("2019/12/01",format="%Y/%m/%d"),
                        by="-1 month",length.out=length(month.meta)))
  return(ER)
}

#2009 intervention catch-up vaccination for pooled effect
#"initiation" is the start point of the catch-up vaccination for each age group
#"year" indicates the calendar year 2009, etc.
ER.meta.2009<-function(model,i,initiation,year){
  coef.meta<-as.matrix(t(model$coefficients))
  #Plus 1 if starting in December
  time.meta<-as.matrix(cbind(rep(1,((i-1)*12)+2),c((initiation-1):(initiation-2+((i-1)*12)+2))))
  cov.meta<-vcov(model)
  pre.ER.meta<-c()
  ER.meta<-c()
  se.ER.meta<-c()
  ER.low.meta<-c()
  ER.high.meta<-c()
  for(l in 1:(((i-1)*12)+2))
  {pre.ER.meta[l]<-time.meta[l,]%*%coef.meta
  ER.meta[l]<-(exp(pre.ER.meta[l])-1)*100
  se.ER.meta[l]<-sqrt(time.meta[l,]%*%cov.meta%*%t(t(time.meta[l,]))) 
  ER.low.meta[l]=(exp(pre.ER.meta[l]-1.96*se.ER.meta[l])-1)*100
  ER.high.meta[l]=(exp(pre.ER.meta[l]+1.96*se.ER.meta[l])-1)*100
  }
  ER.whole<-ER.meta
  ER.whole.low<-ER.low.meta
  ER.whole.upper<-ER.high.meta
  month.meta<-c((initiation-1):(initiation-2+((i-1)*12)+2))
  ER<-data.frame(Time=month.meta,ER=ER.meta,ER.lower=ER.low.meta,ER.upper=ER.high.meta)
  ER$ER<-sprintf("%0.2f",ER$ER)
  ER$ER.lower<-sprintf("%0.2f",ER$ER.lower)
  ER$ER.upper<-sprintf("%0.2f",ER$ER.upper)
  ER$Date<-seq.Date(from=as.Date(paste0(year,"/11/01"),format="%Y/%m/%d"),
                    by="month",length.out=length(month.meta))
  return(ER)
}

#Estimation of EIR
#2009
#"Sequence" is a list of PLADs to be modelled
#"database" is the original data
#"model" is the corresponding statistical model
#"Location" is the location of intervention effect in the model
#"intervention" is the affecting years of intervention
#The unit of intervention needs to be month.
EIR.2009<-function(Sequence,database,model,Location,intervention,n.cir){
  set.seed(202403)
  EMR.reserve<-matrix(0,nrow=31,ncol=8)
  EMR.reserve.diff<-c()
  MA.reserve<-matrix(0,nrow=n.cir,ncol=31)
  for(i in Sequence)
  {data.number<-database[which(database$Order==i),]
  model.number<-eval(parse(text=paste(model,".",i,sep="")))
  mod.number<-as.matrix(data.frame(model.matrix(model.number),offset=log(data.number$Population.book)))
  mod.number.0<-mod.number
  mod.number.0[,Location]<-0
  line.row<-which(data.number$Intervention.2009==1)
  length.row<-length(line.row)
  freq.number<-data.number$Case[line.row]
  pop.number<-data.number$Population.book[line.row]
  #1 indicates the offset
  coef.number<-c(summary(model.number)$p.coeff,1)
  diff.number<-exp(mod.number%*% coef.number)-exp(mod.number.0%*%coef.number)
  diff.number<-diff.number[which(diff.number!=0)]
  EIR.number<-sum(diff.number)/intervention/(sum(pop.number)/length(pop.number))*100000
  EMR.reserve[i,1]<-sum(diff.number)
  EMR.reserve[i,2]<-intervention
  EMR.reserve[i,3]<-sum(pop.number)/length(pop.number)
  EMR.reserve[i,4]<-EIR.number
  #The calculation of the interval
  coef.ini<-coef.number[Location]
  length.coef<-length(coef.ini)
  cov.number<-vcov(model.number)[Location,Location]
  eigen.number<-eigen(cov.number)
  r.norm<-rnorm(length.coef*n.cir)
  r.norm<-matrix(r.norm,n.cir)
  coef.sim<-coef.ini+eigen.number$vectors%*%diag(sqrt(eigen.number$values),length.coef)%*%t(r.norm)
  for(k in 1:n.cir){
    coef.number[Location[1]]<-coef.sim[1,k]
    coef.number[Location[2]]<-coef.sim[2,k]
    diff.number<-exp(mod.number%*%coef.number)-exp(mod.number.0%*%coef.number)
    diff.number<-diff.number[which(diff.number!=0)]
    EMR.reserve.diff[k]<-sum(diff.number)
  }
  MA.reserve[,i]<-EMR.reserve.diff
  if(length(which(is.infinite(EMR.reserve.diff)))==0)
  {EMR.reserve[i,5]<-sort(EMR.reserve.diff)[n.cir*0.025]
  EMR.reserve[i,7]<-sort(EMR.reserve.diff)[n.cir*0.975]}else{
    EMR.reserve[i,5]<-sort(EMR.reserve.diff[-which(is.infinite(EMR.reserve.diff))])[length(EMR.reserve.diff[-which(is.infinite(EMR.reserve.diff))])*0.025]
    EMR.reserve[i,7]<-sort(EMR.reserve.diff[-which(is.infinite(EMR.reserve.diff))])[length(EMR.reserve.diff[-which(is.infinite(EMR.reserve.diff))])*0.975]
  }
  EMR.reserve[i,6]<-EMR.reserve[i,5]/EMR.reserve[i,2]/EMR.reserve[i,3]*100000
  EMR.reserve[i,8]<-EMR.reserve[i,7]/EMR.reserve[i,2]/EMR.reserve[i,3]*100000
  }
  #The calculation of country level results
  EMR.mainland.point<-sum(EMR.reserve[Sequence,1])/mean(EMR.reserve[Sequence,2])/sum(EMR.reserve[Sequence,3])*100000
  EMR.mainland.low<-sort(rowSums(MA.reserve[,Sequence]))[n.cir*0.025]/mean(EMR.reserve[Sequence,2])/sum(EMR.reserve[Sequence,3])*100000
  EMR.mainland.high<-sort(rowSums(MA.reserve[,Sequence]))[n.cir*0.975]/mean(EMR.reserve[Sequence,2])/sum(EMR.reserve[Sequence,3])*100000
  Mainland<-c(sum(EMR.reserve[Sequence,1]),mean(EMR.reserve[Sequence,2]),sum(EMR.reserve[Sequence,3]),
              EMR.mainland.point,sort(rowSums(MA.reserve[,Sequence]))[n.cir*0.025],EMR.mainland.low,
              sort(rowSums(MA.reserve[,Sequence]))[n.cir*0.975],EMR.mainland.high) 
  EMR.reserve<-rbind(EMR.reserve,Mainland)
  EMR.reserve<as.data.frame(EMR.reserve)
  row.names(EMR.reserve)<-c(Name,"Mainland")
  colnames(EMR.reserve)<-c("EC.point","Intervention.year","Population",
                           "EIR.point","EC.low","EIR.low","EC.high","EIR.high")
  return(EMR.reserve)
}

#2002
EIR.2002<-function(Sequence,database,model,Location,intervention,n.cir){
  set.seed(202403)
  EMR.reserve<-matrix(0,nrow=31,ncol=8)
  EMR.reserve.diff<-c()
  MA.reserve<-matrix(0,nrow=n.cir,ncol=31)
  for(i in Sequence)
  {data.number<-database[which(database$Order==i),]
  model.number<-eval(parse(text=paste(model,".",i,sep="")))
  mod.number<-as.matrix(data.frame(model.matrix(model.number),offset=log(data.number$Population.book)))
  mod.number.0<-mod.number
  mod.number.0[,Location]<-0
  line.row<-which(data.number$Intervention.2002==1)
  length.row<-length(line.row)
  freq.number<-data.number$Case[line.row]
  pop.number<-data.number$Population.book[line.row]
  coef.number<-c(summary(model.number)$p.coeff,1)
  diff.number<-exp(mod.number%*% coef.number)-exp(mod.number.0%*%coef.number)
  diff.number<-diff.number[which(diff.number!=0)]
  EIR.number<-sum(diff.number)/intervention/(sum(pop.number)/length(pop.number))*100000
  EMR.reserve[i,1]<-sum(diff.number)
  EMR.reserve[i,2]<-intervention
  EMR.reserve[i,3]<-sum(pop.number)/length(pop.number)
  EMR.reserve[i,4]<-EIR.number
  coef.ini<-coef.number[Location]
  length.coef<-length(coef.ini)
  cov.number<-vcov(model.number)[Location,Location]
  eigen.number<-eigen(cov.number)
  r.norm<-rnorm(length.coef*n.cir)
  r.norm<-matrix(r.norm,n.cir)
  coef.sim<-coef.ini+eigen.number$vectors%*%diag(sqrt(eigen.number$values),length.coef)%*%t(r.norm)
  for(k in 1:n.cir){
    coef.number[Location[1]]<-coef.sim[1,k]
    coef.number[Location[2]]<-coef.sim[2,k]
    diff.number<-exp(mod.number%*%coef.number)-exp(mod.number.0%*%coef.number)
    diff.number<-diff.number[which(diff.number!=0)]
    EMR.reserve.diff[k]<-sum(diff.number)
  }
  MA.reserve[,i]<-EMR.reserve.diff
  EMR.reserve[i,5]<-sort(EMR.reserve.diff)[n.cir*0.025]
  EMR.reserve[i,7]<-sort(EMR.reserve.diff)[n.cir*0.975]
  EMR.reserve[i,6]<-EMR.reserve[i,5]/EMR.reserve[i,2]/EMR.reserve[i,3]*100000
  EMR.reserve[i,8]<-EMR.reserve[i,7]/EMR.reserve[i,2]/EMR.reserve[i,3]*100000
  }
  EMR.reserve[,2]<-intervention
  EMR.mainland.point<-sum(EMR.reserve[Sequence,1])/mean(EMR.reserve[Sequence,2])/sum(EMR.reserve[Sequence,3])*100000
  EMR.mainland.low<-sort(rowSums(MA.reserve[,Sequence]))[n.cir*0.025]/mean(EMR.reserve[Sequence,2])/sum(EMR.reserve[Sequence,3])*100000
  EMR.mainland.high<-sort(rowSums(MA.reserve[,Sequence]))[n.cir*0.975]/mean(EMR.reserve[Sequence,2])/sum(EMR.reserve[Sequence,3])*100000
  Mainland<-c(sum(EMR.reserve[Sequence,1]),mean(EMR.reserve[Sequence,2]),sum(EMR.reserve[Sequence,3]),
              EMR.mainland.point,sort(rowSums(MA.reserve[,Sequence]))[n.cir*0.025],EMR.mainland.low,
              sort(rowSums(MA.reserve[,Sequence]))[n.cir*0.975],EMR.mainland.high) 
  EMR.reserve<-rbind(EMR.reserve,Mainland)
  EMR.reserve<as.data.frame(EMR.reserve)
  row.names(EMR.reserve)<-c(Name,"Mainland")
  colnames(EMR.reserve)<-c("EC.point","Intervention.year","Population",
                           "EIR.point","EC.low","EIR.low","EC.high","EIR.high")
  return(EMR.reserve)
}

#######################Age group 6, only 2002 intervention###############################
setwd("D:/~")
Hepatitis.B.agegroup6<-read.xlsx("agegroup6.xlsx")
Hepatitis.B.agegroup6$Month.factor<-factor(Hepatitis.B.agegroup6$Month)
Hepatitis.B.agegroup6$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup6$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup6$Intervention.2002<-factor(rep(c(rep(0,12*2),rep(1,13*12)),31))
Hepatitis.B.agegroup6$Interaction.2002<-rep(c(rep(0,2*12),c(0:(180-2*12-1))),31)
regions<-as.character(unique(Hepatitis.B.agegroup6$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup6[Hepatitis.B.agegroup6$Province==x, ])
seq(datalist)
yori2<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori2<-vector("list", length(datalist)); names(Sori2) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup6[which(Hepatitis.B.agegroup6$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+Time+
              ns(Temperature,df=3)+
              Holiday+Interaction.2002,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori2[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun),1]
Sori2[[i]] <- vcov(mfirst)[c(2,lengthfun),c(2,lengthfun)]
assign(paste("model2.",i,sep=""),mfirst)
}

mvall2<-mvmeta(yori2[-c(1,2,6,9,26),]~1, Sori2[-c(1,2,6,9,26)], method="reml")
summary(mvall2)

ER2<-ER.meta.2002(mvall2,2)
ER2$ER.lower<-as.numeric(ER2$ER.lower)
ER2$ER.upper<-as.numeric(ER2$ER.upper)
ER2$ER<-as.numeric(ER2$ER)

EIR2<-EIR.2002(c(3:5,7,8,10:25,27:31),Hepatitis.B.agegroup6,"model2",c(2,9),13,10000)
EIR2<-data.frame(EIR2)
EIR2$EC.point<-sprintf("%0.0f",EIR2$EC.point)
EIR2$Population<-sprintf("%0.0f",EIR2$Population)
EIR2$EIR.point<-sprintf("%0.2f",EIR2$EIR.point)
EIR2$EC.low<-sprintf("%0.0f",EIR2$EC.low)
EIR2$EIR.low<-sprintf("%0.2f",EIR2$EIR.low)
EIR2$EC.high<-sprintf("%0.0f",EIR2$EC.high)
EIR2$EIR.high<-sprintf("%0.2f",EIR2$EIR.high)

Coefficient2<-matrix(0,nrow=32,ncol=6)
Coefficient2<-as.data.frame(Coefficient2)
colnames(Coefficient2)<-c("Age","PLADs","RR1 (95%CI)","RR2 (95%CI)","RR3 (95%CI)","Trend")
Coefficient2$Age<-c("4")
yori<-matrix(0,length(datalist),1,dimnames=list(regions, paste("beta",1,sep="")))
Sori<-vector("list", length(datalist)); names(Sori2) <- regions
for (i in 1:31)
{model.number<-eval(parse(text=paste("model2.",i,sep="")))
Point1<-model.number$coefficients[2]
Se1<-summary(model.number)$se[2]
Point2<-model.number$coefficients[3]
Se2<-summary(model.number)$se[3]
length.coef<-length(model.number$coefficients)
Point3<-model.number$coefficients[length.coef]
Se3<-summary(model.number)$se[length.coef]
Coefficient2[i,2]<-Name[i]
Coefficient2[i,3]<-paste0(sprintf("%0.3f",Point1)," (",sprintf("%0.3f",Point1-1.96*Se1),
                          ", ",sprintf("%0.3f",Point1+1.96*Se1),")")
Coefficient2[i,4]<-paste0(sprintf("%0.3f",Point2)," (",sprintf("%0.3f",Point2-1.96*Se2),
                          ", ",sprintf("%0.3f",Point2+1.96*Se2),")") 
Coefficient2[i,5]<-paste0(sprintf("%0.3f",Point3)," (",sprintf("%0.3f",Point3-1.96*Se3),
                          ", ",sprintf("%0.3f",Point3+1.96*Se3),")")
Coefficient2[i,6]<-round(summary(model.number)$p.pv[length.coef],3)
yori[i,]<-as.data.frame(summary(model.number)$p.coeff)[c(3),1]
Sori[[i]]<-vcov(model.number)[c(3),c(3)]
}
mvall<-mvmeta(yori[-c(1,2,6,9,26),]~1, Sori[-c(1,2,6,9,26)], method="reml")
mvall<-summary(mvall)
mvall2<-summary(mvall2)
Coefficient2[32,2]<-c("Chinese mainland")
Coefficient2[32,4]<-paste0(sprintf("%0.3f",mvall$coefficients[1])," (",
                           sprintf("%0.3f",mvall$coefficients[1]-1.96*mvall$coefficients[2]),
                           ", ",sprintf("%0.3f",mvall$coefficients[1]+1.96*mvall$coefficients[2]),")")
Coefficient2[32,3]<-paste0(sprintf("%0.3f",mvall2$coefficients[1,1])," (",
                           sprintf("%0.3f",mvall2$coefficients[1,1]-1.96*mvall2$coefficients[1,2]),
                           ", ",sprintf("%0.3f",mvall2$coefficients[1,1]+1.96*mvall2$coefficients[1,2]),")")
Coefficient2[32,5]<-paste0(sprintf("%0.3f",mvall2$coefficients[2,1])," (",
                           sprintf("%0.3f",mvall2$coefficients[2,1]-1.96*mvall2$coefficients[2,2]),
                           ", ",sprintf("%0.3f",mvall2$coefficients[2,1]+1.96*mvall2$coefficients[2,2]),")")
Coefficient2[32,6]<-round(mvall2$coefficients[2,4],3)

mvall.2<-mvall
mvall2<-mvmeta(yori2[-c(1,2,6,9,26),]~1, Sori2[-c(1,2,6,9,26)], method="reml")
#######################Age group 7, only 2002 intervention###############################
setwd("D:/~")
Hepatitis.B.agegroup7<-read.xlsx("agegroup7.xlsx")
Hepatitis.B.agegroup7$Month.factor<-factor(Hepatitis.B.agegroup7$Month)
Hepatitis.B.agegroup7$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup7$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup7$Intervention.2002<-factor(rep(c(rep(0,12*3),rep(1,12*12)),31))
Hepatitis.B.agegroup7$Interaction.2002<-rep(c(rep(0,3*12),c(0:(180-3*12-1))),31)
regions<-as.character(unique(Hepatitis.B.agegroup7$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup7[Hepatitis.B.agegroup7$Province==x, ])
seq(datalist)
yori3<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori3<-vector("list", length(datalist)); names(Sori3) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup7[which(Hepatitis.B.agegroup7$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+Time+
              ns(Temperature,df=3)+
              Holiday+Interaction.2002,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori3[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun),1]
Sori3[[i]] <- vcov(mfirst)[c(2,lengthfun),c(2,lengthfun)]
assign(paste("model3.",i,sep=""),mfirst)
}

mvall3<-mvmeta(yori3[-c(1,2,9,14,26),]~1, Sori3[-c(1,2,9,14,26)], method="reml")
summary(mvall3)

ER3<-ER.meta.2002(mvall3,3)
ER3$ER.lower<-as.numeric(ER3$ER.lower)
ER3$ER.upper<-as.numeric(ER3$ER.upper)
ER3$ER<-as.numeric(ER3$ER)

EIR3<-EIR.2002(c(3:8,10:13,15:25,27:31),Hepatitis.B.agegroup7,"model3",c(2,9),12,10000)
EIR3<-data.frame(EIR3)
EIR3$EC.point<-sprintf("%0.0f",EIR3$EC.point)
EIR3$Population<-sprintf("%0.0f",EIR3$Population)
EIR3$EIR.point<-sprintf("%0.2f",EIR3$EIR.point)
EIR3$EC.low<-sprintf("%0.0f",EIR3$EC.low)
EIR3$EIR.low<-sprintf("%0.2f",EIR3$EIR.low)
EIR3$EC.high<-sprintf("%0.0f",EIR3$EC.high)
EIR3$EIR.high<-sprintf("%0.2f",EIR3$EIR.high)

Coefficient3<-matrix(0,nrow=32,ncol=6)
Coefficient3<-as.data.frame(Coefficient3)
colnames(Coefficient3)<-c("Age","PLADs","RR1 (95%CI)","RR2 (95%CI)","RR3 (95%CI)","Trend")
Coefficient3$Age<-c("5")
yori<-matrix(0,length(datalist),1,dimnames=list(regions, paste("beta",1,sep="")))
Sori<-vector("list", length(datalist)); names(Sori) <- regions
for (i in 1:31)
{model.number<-eval(parse(text=paste("model3.",i,sep="")))
Point1<-model.number$coefficients[2]
Se1<-summary(model.number)$se[2]
Point2<-model.number$coefficients[3]
Se2<-summary(model.number)$se[3]
length.coef<-length(model.number$coefficients)
Point3<-model.number$coefficients[length.coef]
Se3<-summary(model.number)$se[length.coef]
Coefficient3[i,2]<-Name[i]
Coefficient3[i,3]<-paste0(sprintf("%0.3f",Point1)," (",sprintf("%0.3f",Point1-1.96*Se1),
                          ", ",sprintf("%0.3f",Point1+1.96*Se1),")")
Coefficient3[i,4]<-paste0(sprintf("%0.3f",Point2)," (",sprintf("%0.3f",Point2-1.96*Se2),
                          ", ",sprintf("%0.3f",Point2+1.96*Se2),")") 
Coefficient3[i,5]<-paste0(sprintf("%0.3f",Point3)," (",sprintf("%0.3f",Point3-1.96*Se3),
                          ", ",sprintf("%0.3f",Point3+1.96*Se3),")")
Coefficient3[i,6]<-round(summary(model.number)$p.pv[length.coef],3)
yori[i,]<-as.data.frame(summary(model.number)$p.coeff)[c(3),1]
Sori[[i]]<-vcov(model.number)[c(3),c(3)]
}
mvall<-mvmeta(yori[-c(1,2,9,14,26),]~1, Sori[-c(1,2,9,14,26)], method="reml")
mvall<-summary(mvall)
mvall3<-summary(mvall3)
Coefficient3[32,2]<-c("Chinese mainland")
Coefficient3[32,4]<-paste0(sprintf("%0.3f",mvall$coefficients[1])," (",
                           sprintf("%0.3f",mvall$coefficients[1]-1.96*mvall$coefficients[2]),
                           ", ",sprintf("%0.3f",mvall$coefficients[1]+1.96*mvall$coefficients[2]),")")
Coefficient3[32,3]<-paste0(sprintf("%0.3f",mvall3$coefficients[1,1])," (",
                           sprintf("%0.3f",mvall3$coefficients[1,1]-1.96*mvall3$coefficients[1,2]),
                           ", ",sprintf("%0.3f",mvall3$coefficients[1,1]+1.96*mvall3$coefficients[1,2]),")")
Coefficient3[32,5]<-paste0(sprintf("%0.3f",mvall3$coefficients[2,1])," (",
                           sprintf("%0.3f",mvall3$coefficients[2,1]-1.96*mvall3$coefficients[2,2]),
                           ", ",sprintf("%0.3f",mvall3$coefficients[2,1]+1.96*mvall3$coefficients[2,2]),")")
Coefficient3[32,6]<-round(mvall3$coefficients[2,4],3)

mvall.3<-mvall
mvall3<-mvmeta(yori3[-c(1,2,9,14,26),]~1, Sori3[-c(1,2,9,14,26)], method="reml")
#######################Age group 8, only 2002 intervention###############################
setwd("D:/~")
Hepatitis.B.agegroup8<-read.xlsx("agegroup8.xlsx")
Hepatitis.B.agegroup8$Month.factor<-factor(Hepatitis.B.agegroup8$Month)
Hepatitis.B.agegroup8$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup8$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup8$Intervention.2002<-factor(rep(c(rep(0,12*4),rep(1,12*11)),31))
Hepatitis.B.agegroup8$Interaction.2002<-rep(c(rep(0,4*12),c(0:(180-4*12-1))),31)
regions<-as.character(unique(Hepatitis.B.agegroup8$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup8[Hepatitis.B.agegroup8$Province==x, ])
seq(datalist)
yori4<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori4<-vector("list", length(datalist)); names(Sori4) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup8[which(Hepatitis.B.agegroup8$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+Time+
              ns(Temperature,df=3)+
              Holiday+Interaction.2002,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori4[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun),1]
Sori4[[i]] <- vcov(mfirst)[c(2,lengthfun),c(2,lengthfun)]
assign(paste("model4.",i,sep=""),mfirst)
}

mvall4<-mvmeta(yori4[-c(1,2,9,14,26),]~1, Sori4[-c(1,2,9,14,26)], method="reml")
summary(mvall4)

ER4<-ER.meta.2002(mvall4,4)
ER4$ER.lower<-as.numeric(ER4$ER.lower)
ER4$ER.upper<-as.numeric(ER4$ER.upper)
ER4$ER<-as.numeric(ER4$ER)

EIR4<-EIR.2002(c(3:8,10:13,15:25,27:31),Hepatitis.B.agegroup8,"model4",c(2,9),11,10000)
EIR4<-data.frame(EIR4)
EIR4$EC.point<-sprintf("%0.0f",EIR4$EC.point)
EIR4$Population<-sprintf("%0.0f",EIR4$Population)
EIR4$EIR.point<-sprintf("%0.2f",EIR4$EIR.point)
EIR4$EC.low<-sprintf("%0.0f",EIR4$EC.low)
EIR4$EIR.low<-sprintf("%0.2f",EIR4$EIR.low)
EIR4$EC.high<-sprintf("%0.0f",EIR4$EC.high)
EIR4$EIR.high<-sprintf("%0.2f",EIR4$EIR.high)

Coefficient4<-matrix(0,nrow=32,ncol=6)
Coefficient4<-as.data.frame(Coefficient4)
colnames(Coefficient4)<-c("Age","PLADs","RR1 (95%CI)","RR2 (95%CI)","RR3 (95%CI)","Trend")
Coefficient4$Age<-c("7")
yori<-matrix(0,length(datalist),1,dimnames=list(regions, paste("beta",1,sep="")))
Sori<-vector("list", length(datalist)); names(Sori) <- regions
for (i in 1:31)
{model.number<-eval(parse(text=paste("model4.",i,sep="")))
Point1<-model.number$coefficients[2]
Se1<-summary(model.number)$se[2]
Point2<-model.number$coefficients[3]
Se2<-summary(model.number)$se[3]
length.coef<-length(model.number$coefficients)
Point3<-model.number$coefficients[length.coef]
Se3<-summary(model.number)$se[length.coef]
Coefficient4[i,2]<-Name[i]
Coefficient4[i,3]<-paste0(sprintf("%0.3f",Point1)," (",sprintf("%0.3f",Point1-1.96*Se1),
                          ", ",sprintf("%0.3f",Point1+1.96*Se1),")")
Coefficient4[i,4]<-paste0(sprintf("%0.3f",Point2)," (",sprintf("%0.3f",Point2-1.96*Se2),
                          ", ",sprintf("%0.3f",Point2+1.96*Se2),")") 
Coefficient4[i,5]<-paste0(sprintf("%0.3f",Point3)," (",sprintf("%0.3f",Point3-1.96*Se3),
                          ", ",sprintf("%0.3f",Point3+1.96*Se3),")")
Coefficient4[i,6]<-round(summary(model.number)$p.pv[length.coef],3)
yori[i,]<-as.data.frame(summary(model.number)$p.coeff)[c(3),1]
Sori[[i]]<-vcov(model.number)[c(3),c(3)]
}
mvall<-mvmeta(yori[-c(1,2,9,14,26),]~1, Sori[-c(1,2,9,14,26)], method="reml")
mvall<-summary(mvall)
mvall4<-summary(mvall4)
Coefficient4[32,2]<-c("Chinese mainland")
Coefficient4[32,4]<-paste0(sprintf("%0.3f",mvall$coefficients[1])," (",
                           sprintf("%0.3f",mvall$coefficients[1]-1.96*mvall$coefficients[2]),
                           ", ",sprintf("%0.3f",mvall$coefficients[1]+1.96*mvall$coefficients[2]),")")
Coefficient4[32,3]<-paste0(sprintf("%0.3f",mvall4$coefficients[1,1])," (",
                           sprintf("%0.3f",mvall4$coefficients[1,1]-1.96*mvall4$coefficients[1,2]),
                           ", ",sprintf("%0.3f",mvall4$coefficients[1,1]+1.96*mvall4$coefficients[1,2]),")")
Coefficient4[32,5]<-paste0(sprintf("%0.3f",mvall4$coefficients[2,1])," (",
                           sprintf("%0.3f",mvall4$coefficients[2,1]-1.96*mvall4$coefficients[2,2]),
                           ", ",sprintf("%0.3f",mvall4$coefficients[2,1]+1.96*mvall4$coefficients[2,2]),")")
Coefficient4[32,6]<-round(mvall4$coefficients[2,4],3)

mvall.4<-mvall
mvall4<-mvmeta(yori4[-c(1,2,9,14,26),]~1, Sori4[-c(1,2,9,14,26)], method="reml")
#######################Age group 9, both 2002 and 2009 intervention########################
setwd("D:/~")
Hepatitis.B.agegroup9<-read.xlsx("agegroup9.xlsx")
Hepatitis.B.agegroup9$Month.factor<-factor(Hepatitis.B.agegroup9$Month)
Hepatitis.B.agegroup9$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup9$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
#2002 intervention term
Hepatitis.B.agegroup9$Intervention.2002<-factor(rep(c(rep(0,12*5),rep(1,12*10)),31))
#2009 intervention term
Hepatitis.B.agegroup9$Intervention.2009<-factor(rep(c(rep(0,12*4),rep(0,10),rep(1,2),rep(0,12*10)),31))
#2002 interaction term
Hepatitis.B.agegroup9$Interaction.2002<-rep(c(rep(0,5*12),c(0:(180-5*12-1))),31)
#2009 interaction term
Hepatitis.B.agegroup9$Interaction.2009<-rep(c(rep(0,4*12),rep(0,10),c(0:(2-1)),rep(0,10*12)),31)
regions<-as.character(unique(Hepatitis.B.agegroup9$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup9[Hepatitis.B.agegroup9$Province==x, ])
seq(datalist)
yori5.1<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori5.1<-vector("list", length(datalist)); names(Sori5.1) <- regions
yori5.2<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori5.2<-vector("list", length(datalist)); names(Sori5.2) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup9[which(Hepatitis.B.agegroup9$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+Intervention.2009+
              Time+ns(Temperature,df=3)+
              Holiday+Interaction.2002+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori5.1[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun-1),1]
Sori5.1[[i]] <- vcov(mfirst)[c(2,lengthfun-1),c(2,lengthfun-1)]
yori5.2[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(3,lengthfun),1]
Sori5.2[[i]] <- vcov(mfirst)[c(3,lengthfun),c(3,lengthfun)]
assign(paste("model5.",i,sep=""),mfirst)
}

#ER of 2002 intervention
mvall5.1<-mvmeta(yori5.1[-c(1,2,9,14,26),]~1, Sori5.1[-c(1,2,9,14,26)], method="reml")
summary(mvall5.1)

ER5.1<-ER.meta.2002(mvall5.1,5)
ER5.1$ER.lower<-as.numeric(ER5.1$ER.lower)
ER5.1$ER.upper<-as.numeric(ER5.1$ER.upper)
ER5.1$ER<-as.numeric(ER5.1$ER)

EIR5.1<-EIR.2002(c(3:8,10:13,15:25,27:31),Hepatitis.B.agegroup9,"model5",c(2,10),10,10000)
EIR5.1<-data.frame(EIR5.1)
EIR5.1$EC.point<-sprintf("%0.0f",EIR5.1$EC.point)
EIR5.1$Population<-sprintf("%0.0f",EIR5.1$Population)
EIR5.1$EIR.point<-sprintf("%0.2f",EIR5.1$EIR.point)
EIR5.1$EC.low<-sprintf("%0.0f",EIR5.1$EC.low)
EIR5.1$EIR.low<-sprintf("%0.2f",EIR5.1$EIR.low)
EIR5.1$EC.high<-sprintf("%0.0f",EIR5.1$EC.high)
EIR5.1$EIR.high<-sprintf("%0.2f",EIR5.1$EIR.high)

#Without the further estimation of the 2009 intervention effectiveness among the childrend aged 8y

#Save the cofficients β of 2002
Coefficient5.1<-matrix(0,nrow=32,ncol=6)
Coefficient5.1<-as.data.frame(Coefficient5.1)
colnames(Coefficient5.1)<-c("Age","PLADs","RR1 (95%CI)","RR2 (95%CI)","RR3 (95%CI)","Trend")
Coefficient5.1$Age<-c("8")
yori<-matrix(0,length(datalist),1,dimnames=list(regions, paste("beta",1,sep="")))
Sori<-vector("list", length(datalist)); names(Sori) <- regions
for (i in 1:31)
{model.number<-eval(parse(text=paste("model5.",i,sep="")))
Point1<-model.number$coefficients[2]
Se1<-summary(model.number)$se[2]
Point2<-model.number$coefficients[4]
Se2<-summary(model.number)$se[4]
length.coef<-length(model.number$coefficients)
Point3<-model.number$coefficients[length.coef-1]
Se3<-summary(model.number)$se[length.coef-1]
Coefficient5.1[i,2]<-Name[i]
Coefficient5.1[i,3]<-paste0(sprintf("%0.3f",Point1)," (",sprintf("%0.3f",Point1-1.96*Se1),
                            ", ",sprintf("%0.3f",Point1+1.96*Se1),")")
Coefficient5.1[i,4]<-paste0(sprintf("%0.3f",Point2)," (",sprintf("%0.3f",Point2-1.96*Se2),
                            ", ",sprintf("%0.3f",Point2+1.96*Se2),")") 
Coefficient5.1[i,5]<-paste0(sprintf("%0.3f",Point3)," (",sprintf("%0.3f",Point3-1.96*Se3),
                            ", ",sprintf("%0.3f",Point3+1.96*Se3),")")
Coefficient5.1[i,6]<-round(summary(model.number)$p.pv[length.coef],3)
yori[i,]<-as.data.frame(summary(model.number)$p.coeff)[c(4),1]
Sori[[i]]<-vcov(model.number)[c(4),c(4)]
}
mvall<-mvmeta(yori[-c(1,2,9,14,26),]~1, Sori[-c(1,2,9,14,26)], method="reml")
mvall<-summary(mvall)
mvall5.1<-summary(mvall5.1)
Coefficient5.1[32,2]<-c("Chinese mainland")
Coefficient5.1[32,4]<-paste0(sprintf("%0.3f",mvall$coefficients[1])," (",
                             sprintf("%0.3f",mvall$coefficients[1]-1.96*mvall$coefficients[2]),
                             ", ",sprintf("%0.3f",mvall$coefficients[1]+1.96*mvall$coefficients[2]),")")
Coefficient5.1[32,3]<-paste0(sprintf("%0.3f",mvall5.1$coefficients[1,1])," (",
                             sprintf("%0.3f",mvall5.1$coefficients[1,1]-1.96*mvall5.1$coefficients[1,2]),
                             ", ",sprintf("%0.3f",mvall5.1$coefficients[1,1]+1.96*mvall5.1$coefficients[1,2]),")")
Coefficient5.1[32,5]<-paste0(sprintf("%0.3f",mvall5.1$coefficients[2,1])," (",
                             sprintf("%0.3f",mvall5.1$coefficients[2,1]-1.96*mvall5.1$coefficients[2,2]),
                             ", ",sprintf("%0.3f",mvall5.1$coefficients[2,1]+1.96*mvall5.1$coefficients[2,2]),")")
Coefficient5.1[32,6]<-round(mvall5.1$coefficients[2,4],3)

mvall.5<-mvall
mvall5.1<-mvmeta(yori5.1[-c(1,2,9,14,26),]~1, Sori5.1[-c(1,2,9,14,26)], method="reml")
#######################Age group 10, both 2002 and 2009 intervention#######################
setwd("D:/~")
Hepatitis.B.agegroup10<-read.xlsx("agegroup10.xlsx")
Hepatitis.B.agegroup10$Month.factor<-factor(Hepatitis.B.agegroup10$Month)
Hepatitis.B.agegroup10$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup10$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup10$Intervention.2002<-factor(rep(c(rep(0,12*6),rep(1,12*9)),31))
Hepatitis.B.agegroup10$Intervention.2009<-factor(rep(c(rep(0,12*4),rep(0,10),rep(1,2+12),rep(0,12*9)),31))
Hepatitis.B.agegroup10$Interaction.2002<-rep(c(rep(0,6*12),c(0:(180-6*12-1))),31)
Hepatitis.B.agegroup10$Interaction.2009<-rep(c(rep(0,4*12),rep(0,10),c(0:(2+12-1)),rep(0,9*12)),31)
regions<-as.character(unique(Hepatitis.B.agegroup10$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup10[Hepatitis.B.agegroup10$Province==x, ])
seq(datalist)
yori6.1<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori6.1<-vector("list", length(datalist)); names(Sori6.1) <- regions
yori6.2<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori6.2<-vector("list", length(datalist)); names(Sori6.2) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+Intervention.2009+
              Time+ns(Temperature,df=3)+
              Holiday+Interaction.2002+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori6.1[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun-1),1]
Sori6.1[[i]] <- vcov(mfirst)[c(2,lengthfun-1),c(2,lengthfun-1)]
yori6.2[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(3,lengthfun),1]
Sori6.2[[i]] <- vcov(mfirst)[c(3,lengthfun),c(3,lengthfun)]
assign(paste("model6.",i,sep=""),mfirst)
}

mvall6.1<-mvmeta(yori6.1[-c(1,2,9,14,26),]~1, Sori6.1[-c(1,2,9,14,26)], method="reml")
summary(mvall6.1)

ER6.1<-ER.meta.2002(mvall6.1,6)
ER6.1$ER.lower<-as.numeric(ER6.1$ER.lower)
ER6.1$ER.upper<-as.numeric(ER6.1$ER.upper)
ER6.1$ER<-as.numeric(ER6.1$ER)

EIR6.1<-EIR.2002(c(3:8,10:13,15:25,27:31),Hepatitis.B.agegroup10,"model6",c(2,10),9,10000)
EIR6.1<-data.frame(EIR6.1)
EIR6.1$EC.point<-sprintf("%0.0f",EIR6.1$EC.point)
EIR6.1$Population<-sprintf("%0.0f",EIR6.1$Population)
EIR6.1$EIR.point<-sprintf("%0.2f",EIR6.1$EIR.point)
EIR6.1$EC.low<-sprintf("%0.0f",EIR6.1$EC.low)
EIR6.1$EIR.low<-sprintf("%0.2f",EIR6.1$EIR.low)
EIR6.1$EC.high<-sprintf("%0.0f",EIR6.1$EC.high)
EIR6.1$EIR.high<-sprintf("%0.2f",EIR6.1$EIR.high)

mvall6.2<-mvmeta(yori6.2[-c(1,2,9,14,26),]~1, Sori6.2[-c(1,2,9,14,26)], method="reml")
summary(mvall6.2)

ER6.2<-ER.meta.2009(mvall6.2,2,59,2009)
ER6.2$ER.lower<-as.numeric(ER6.2$ER.lower)
ER6.2$ER.upper<-as.numeric(ER6.2$ER.upper)
ER6.2$ER<-as.numeric(ER6.2$ER)

EIR6.2<-EIR.2009(c(3:8,10:13,15:25,27:31),Hepatitis.B.agegroup10,"model6",c(3,11),(2+1*12)/12,10000)
EIR6.2<-data.frame(EIR6.2)
EIR6.2$EC.point<-sprintf("%0.0f",EIR6.2$EC.point)
EIR6.2$Population<-sprintf("%0.0f",EIR6.2$Population)
EIR6.2$EIR.point<-sprintf("%0.2f",EIR6.2$EIR.point)
EIR6.2$EC.low<-sprintf("%0.0f",EIR6.2$EC.low)
EIR6.2$EIR.low<-sprintf("%0.2f",EIR6.2$EIR.low)
EIR6.2$EC.high<-sprintf("%0.0f",EIR6.2$EC.high)
EIR6.2$EIR.high<-sprintf("%0.2f",EIR6.2$EIR.high)

Coefficient6.1<-matrix(0,nrow=32,ncol=6)
Coefficient6.1<-as.data.frame(Coefficient6.1)
colnames(Coefficient6.1)<-c("Age","PLADs","RR1 (95%CI)","RR2 (95%CI)","RR3 (95%CI)","Trend")
Coefficient6.1$Age<-c("9")
yori<-matrix(0,length(datalist),1,dimnames=list(regions, paste("beta",1,sep="")))
Sori<-vector("list", length(datalist)); names(Sori) <- regions
for (i in 1:31)
{model.number<-eval(parse(text=paste("model6.",i,sep="")))
Point1<-model.number$coefficients[2]
Se1<-summary(model.number)$se[2]
Point2<-model.number$coefficients[4]
Se2<-summary(model.number)$se[4]
length.coef<-length(model.number$coefficients)
Point3<-model.number$coefficients[length.coef-1]
Se3<-summary(model.number)$se[length.coef-1]
Coefficient6.1[i,2]<-Name[i]
Coefficient6.1[i,3]<-paste0(sprintf("%0.3f",Point1)," (",sprintf("%0.3f",Point1-1.96*Se1),
                            ", ",sprintf("%0.3f",Point1+1.96*Se1),")")
Coefficient6.1[i,4]<-paste0(sprintf("%0.3f",Point2)," (",sprintf("%0.3f",Point2-1.96*Se2),
                            ", ",sprintf("%0.3f",Point2+1.96*Se2),")") 
Coefficient6.1[i,5]<-paste0(sprintf("%0.3f",Point3)," (",sprintf("%0.3f",Point3-1.96*Se3),
                            ", ",sprintf("%0.3f",Point3+1.96*Se3),")")
Coefficient6.1[i,6]<-round(summary(model.number)$p.pv[length.coef],3)
yori[i,]<-as.data.frame(summary(model.number)$p.coeff)[c(4),1]
Sori[[i]]<-vcov(model.number)[c(4),c(4)]
}
mvall<-mvmeta(yori[-c(1,2,9,14,26),]~1, Sori[-c(1,2,9,14,26)], method="reml")
mvall<-summary(mvall)
mvall6.1<-summary(mvall6.1)
Coefficient6.1[32,2]<-c("Chinese mainland")
Coefficient6.1[32,4]<-paste0(sprintf("%0.3f",mvall$coefficients[1])," (",
                             sprintf("%0.3f",mvall$coefficients[1]-1.96*mvall$coefficients[2]),
                             ", ",sprintf("%0.3f",mvall$coefficients[1]+1.96*mvall$coefficients[2]),")")
Coefficient6.1[32,3]<-paste0(sprintf("%0.3f",mvall6.1$coefficients[1,1])," (",
                             sprintf("%0.3f",mvall6.1$coefficients[1,1]-1.96*mvall6.1$coefficients[1,2]),
                             ", ",sprintf("%0.3f",mvall6.1$coefficients[1,1]+1.96*mvall6.1$coefficients[1,2]),")")
Coefficient6.1[32,5]<-paste0(sprintf("%0.3f",mvall6.1$coefficients[2,1])," (",
                             sprintf("%0.3f",mvall6.1$coefficients[2,1]-1.96*mvall6.1$coefficients[2,2]),
                             ", ",sprintf("%0.3f",mvall6.1$coefficients[2,1]+1.96*mvall6.1$coefficients[2,2]),")")
Coefficient6.1[32,6]<-round(mvall6.1$coefficients[2,4],3)

Coefficient6.2<-matrix(0,nrow=32,ncol=6)
Coefficient6.2<-as.data.frame(Coefficient6.2)
colnames(Coefficient6.2)<-c("Age","PLADs","RR1 (95%CI)","RR2 (95%CI)","RR3 (95%CI)","Trend")
Coefficient6.2$Age<-c("9")
yori<-matrix(0,length(datalist),1,dimnames=list(regions, paste("beta",1,sep="")))
Sori<-vector("list", length(datalist)); names(Sori) <- regions
for (i in 1:31)
{model.number<-eval(parse(text=paste("model6.",i,sep="")))
Point1<-model.number$coefficients[3]
Se1<-summary(model.number)$se[3]
Point2<-model.number$coefficients[4]
Se2<-summary(model.number)$se[4]
length.coef<-length(model.number$coefficients)
Point3<-model.number$coefficients[length.coef]
Se3<-summary(model.number)$se[length.coef]
Coefficient6.2[i,2]<-Name[i]
Coefficient6.2[i,3]<-paste0(sprintf("%0.3f",Point1)," (",sprintf("%0.3f",Point1-1.96*Se1),
                            ", ",sprintf("%0.3f",Point1+1.96*Se1),")")
Coefficient6.2[i,4]<-paste0(sprintf("%0.3f",Point2)," (",sprintf("%0.3f",Point2-1.96*Se2),
                            ", ",sprintf("%0.3f",Point2+1.96*Se2),")") 
Coefficient6.2[i,5]<-paste0(sprintf("%0.3f",Point3)," (",sprintf("%0.3f",Point3-1.96*Se3),
                            ", ",sprintf("%0.3f",Point3+1.96*Se3),")")
Coefficient6.2[i,6]<-round(summary(model.number)$p.pv[length.coef],3)
yori[i,]<-as.data.frame(summary(model.number)$p.coeff)[c(4),1]
Sori[[i]]<-vcov(model.number)[c(4),c(4)]
}
mvall<-mvmeta(yori[-c(1,2,9,14,26),]~1, Sori[-c(1,2,9,14,26)], method="reml")
mvall<-summary(mvall)
mvall6.2<-summary(mvall6.2)
Coefficient6.2[32,2]<-c("Chinese mainland")
Coefficient6.2[32,4]<-paste0(sprintf("%0.3f",mvall$coefficients[1])," (",
                             sprintf("%0.3f",mvall$coefficients[1]-1.96*mvall$coefficients[2]),
                             ", ",sprintf("%0.3f",mvall$coefficients[1]+1.96*mvall$coefficients[2]),")")
Coefficient6.2[32,3]<-paste0(sprintf("%0.3f",mvall6.2$coefficients[1,1])," (",
                             sprintf("%0.3f",mvall6.2$coefficients[1,1]-1.96*mvall6.2$coefficients[1,2]),
                             ", ",sprintf("%0.3f",mvall6.2$coefficients[1,1]+1.96*mvall6.2$coefficients[1,2]),")")
Coefficient6.2[32,5]<-paste0(sprintf("%0.3f",mvall6.2$coefficients[2,1])," (",
                             sprintf("%0.3f",mvall6.2$coefficients[2,1]-1.96*mvall6.2$coefficients[2,2]),
                             ", ",sprintf("%0.3f",mvall6.2$coefficients[2,1]+1.96*mvall6.2$coefficients[2,2]),")")
Coefficient6.2[32,6]<-round(mvall6.2$coefficients[2,4],3)

mvall.6<-mvall
mvall6.1<-mvmeta(yori6.1[-c(1,2,9,14,26),]~1, Sori6.1[-c(1,2,9,14,26)], method="reml")
mvall6.2<-mvmeta(yori6.2[-c(1,2,9,14,26),]~1, Sori6.2[-c(1,2,9,14,26)], method="reml")
#######################Age group 11, both 2002 and 2009 intervention#######################
setwd("D:/~")
Hepatitis.B.agegroup11<-read.xlsx("agegroup11.xlsx")
Hepatitis.B.agegroup11$Month.factor<-factor(Hepatitis.B.agegroup11$Month)
Hepatitis.B.agegroup11$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup11$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup11$Intervention.2002<-factor(rep(c(rep(0,12*7),rep(1,12*8)),31))
Hepatitis.B.agegroup11$Intervention.2009<-factor(rep(c(rep(0,12*4),rep(0,10),rep(1,2+2*12),rep(0,12*8)),31))
Hepatitis.B.agegroup11$Interaction.2002<-rep(c(rep(0,7*12),c(0:(180-7*12-1))),31)
Hepatitis.B.agegroup11$Interaction.2009<-rep(c(rep(0,4*12),rep(0,10),c(0:(2+2*12-1)),rep(0,8*12)),31)
regions<-as.character(unique(Hepatitis.B.agegroup11$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup11[Hepatitis.B.agegroup11$Province==x, ])
seq(datalist)
yori7.1<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori7.1<-vector("list", length(datalist)); names(Sori7.1) <- regions
yori7.2<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori7.2<-vector("list", length(datalist)); names(Sori7.2) <- regions
Percent<-as.data.frame(matrix(0,nrow=31,ncol=6))
Percent[,1]<-regions
year.pre<-rep(4,31)*12
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+Intervention.2009+
              Time+ns(Temperature,df=3)+
              Holiday+Interaction.2002+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori7.1[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun-1),1]
Sori7.1[[i]] <- vcov(mfirst)[c(2,lengthfun-1),c(2,lengthfun-1)]
yori7.2[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(3,lengthfun),1]
Sori7.2[[i]] <- vcov(mfirst)[c(3,lengthfun),c(3,lengthfun)]
assign(paste("model7.",i,sep=""),mfirst)
}

mvall7.1<-mvmeta(yori7.1~1, Sori7.1, method="reml")
summary(mvall7.1)

ER7.1<-ER.meta.2002(mvall7.1,7)
ER7.1$ER.lower<-as.numeric(ER7.1$ER.lower)
ER7.1$ER.upper<-as.numeric(ER7.1$ER.upper)
ER7.1$ER<-as.numeric(ER7.1$ER)

EIR7.1<-EIR.2002(c(1:31),Hepatitis.B.agegroup11,"model7",c(2,10),8,10000)
EIR7.1<-data.frame(EIR7.1)
EIR7.1$EC.point<-sprintf("%0.0f",EIR7.1$EC.point)
EIR7.1$Population<-sprintf("%0.0f",EIR7.1$Population)
EIR7.1$EIR.point<-sprintf("%0.2f",EIR7.1$EIR.point)
EIR7.1$EC.low<-sprintf("%0.0f",EIR7.1$EC.low)
EIR7.1$EIR.low<-sprintf("%0.2f",EIR7.1$EIR.low)
EIR7.1$EC.high<-sprintf("%0.0f",EIR7.1$EC.high)
EIR7.1$EIR.high<-sprintf("%0.2f",EIR7.1$EIR.high)

mvall7.2<-mvmeta(yori7.2~1, Sori7.2, method="reml")
summary(mvall7.2)

ER7.2<-ER.meta.2009(mvall7.2,3,59,2009)
ER7.2$ER.lower<-as.numeric(ER7.2$ER.lower)
ER7.2$ER.upper<-as.numeric(ER7.2$ER.upper)
ER7.2$ER<-as.numeric(ER7.2$ER)

EIR7.2<-EIR.2009(c(1:31),Hepatitis.B.agegroup11,"model7",c(3,11),(2+2*12)/12,10000)
EIR7.2<-data.frame(EIR7.2)
EIR7.2$EC.point<-sprintf("%0.0f",EIR7.2$EC.point)
EIR7.2$Population<-sprintf("%0.0f",EIR7.2$Population)
EIR7.2$EIR.point<-sprintf("%0.2f",EIR7.2$EIR.point)
EIR7.2$EC.low<-sprintf("%0.0f",EIR7.2$EC.low)
EIR7.2$EIR.low<-sprintf("%0.2f",EIR7.2$EIR.low)
EIR7.2$EC.high<-sprintf("%0.0f",EIR7.2$EC.high)
EIR7.2$EIR.high<-sprintf("%0.2f",EIR7.2$EIR.high)

Coefficient7.1<-matrix(0,nrow=32,ncol=6)
Coefficient7.1<-as.data.frame(Coefficient7.1)
colnames(Coefficient7.1)<-c("Age","PLADs","RR1 (95%CI)","RR2 (95%CI)","RR3 (95%CI)","Trend")
Coefficient7.1$Age<-c("10~14")
yori<-matrix(0,length(datalist),1,dimnames=list(regions, paste("beta",1,sep="")))
Sori<-vector("list", length(datalist)); names(Sori) <- regions
for (i in 1:31)
{model.number<-eval(parse(text=paste("model7.",i,sep="")))
Point1<-model.number$coefficients[2]
Se1<-summary(model.number)$se[2]
Point2<-model.number$coefficients[4]
Se2<-summary(model.number)$se[4]
length.coef<-length(model.number$coefficients)
Point3<-model.number$coefficients[length.coef-1]
Se3<-summary(model.number)$se[length.coef-1]
Coefficient7.1[i,2]<-Name[i]
Coefficient7.1[i,3]<-paste0(sprintf("%0.3f",Point1)," (",sprintf("%0.3f",Point1-1.96*Se1),
                            ", ",sprintf("%0.3f",Point1+1.96*Se1),")")
Coefficient7.1[i,4]<-paste0(sprintf("%0.3f",Point2)," (",sprintf("%0.3f",Point2-1.96*Se2),
                            ", ",sprintf("%0.3f",Point2+1.96*Se2),")") 
Coefficient7.1[i,5]<-paste0(sprintf("%0.3f",Point3)," (",sprintf("%0.3f",Point3-1.96*Se3),
                            ", ",sprintf("%0.3f",Point3+1.96*Se3),")")
Coefficient7.1[i,6]<-round(summary(model.number)$p.pv[length.coef],3)
yori[i,]<-as.data.frame(summary(model.number)$p.coeff)[c(4),1]
Sori[[i]]<-vcov(model.number)[c(4),c(4)]
}
mvall<-mvmeta(yori~1, Sori, method="reml")
mvall<-summary(mvall)
mvall7.1<-summary(mvall7.1)
Coefficient7.1[32,2]<-c("Chinese mainland")
Coefficient7.1[32,4]<-paste0(sprintf("%0.3f",mvall$coefficients[1])," (",
                             sprintf("%0.3f",mvall$coefficients[1]-1.96*mvall$coefficients[2]),
                             ", ",sprintf("%0.3f",mvall$coefficients[1]+1.96*mvall$coefficients[2]),")")
Coefficient7.1[32,3]<-paste0(sprintf("%0.3f",mvall7.1$coefficients[1,1])," (",
                             sprintf("%0.3f",mvall7.1$coefficients[1,1]-1.96*mvall7.1$coefficients[1,2]),
                             ", ",sprintf("%0.3f",mvall7.1$coefficients[1,1]+1.96*mvall7.1$coefficients[1,2]),")")
Coefficient7.1[32,5]<-paste0(sprintf("%0.3f",mvall7.1$coefficients[2,1])," (",
                             sprintf("%0.3f",mvall7.1$coefficients[2,1]-1.96*mvall7.1$coefficients[2,2]),
                             ", ",sprintf("%0.3f",mvall7.1$coefficients[2,1]+1.96*mvall7.1$coefficients[2,2]),")")
Coefficient7.1[32,6]<-round(mvall7.1$coefficients[2,4],3)

Coefficient7.2<-matrix(0,nrow=32,ncol=6)
Coefficient7.2<-as.data.frame(Coefficient7.2)
colnames(Coefficient7.2)<-c("Age","PLADs","RR1 (95%CI)","RR2 (95%CI)","RR3 (95%CI)","Trend")
Coefficient7.2$Age<-c("10~14")
yori<-matrix(0,length(datalist),1,dimnames=list(regions, paste("beta",1,sep="")))
Sori<-vector("list", length(datalist)); names(Sori) <- regions
for (i in 1:31)
{model.number<-eval(parse(text=paste("model7.",i,sep="")))
Point1<-model.number$coefficients[3]
Se1<-summary(model.number)$se[3]
Point2<-model.number$coefficients[4]
Se2<-summary(model.number)$se[4]
length.coef<-length(model.number$coefficients)
Point3<-model.number$coefficients[length.coef]
Se3<-summary(model.number)$se[length.coef]
Coefficient7.2[i,2]<-Name[i]
Coefficient7.2[i,3]<-paste0(sprintf("%0.3f",Point1)," (",sprintf("%0.3f",Point1-1.96*Se1),
                            ", ",sprintf("%0.3f",Point1+1.96*Se1),")")
Coefficient7.2[i,4]<-paste0(sprintf("%0.3f",Point2)," (",sprintf("%0.3f",Point2-1.96*Se2),
                            ", ",sprintf("%0.3f",Point2+1.96*Se2),")") 
Coefficient7.2[i,5]<-paste0(sprintf("%0.3f",Point3)," (",sprintf("%0.3f",Point3-1.96*Se3),
                            ", ",sprintf("%0.3f",Point3+1.96*Se3),")")
Coefficient7.2[i,6]<-round(summary(model.number)$p.pv[length.coef],3)
yori[i,]<-as.data.frame(summary(model.number)$p.coeff)[c(4),1]
Sori[[i]]<-vcov(model.number)[c(4),c(4)]
}
mvall<-mvmeta(yori~1, Sori, method="reml")
mvall<-summary(mvall)
mvall7.2<-summary(mvall7.2)
Coefficient7.2[32,2]<-c("Chinese mainland")
Coefficient7.2[32,4]<-paste0(sprintf("%0.3f",mvall$coefficients[1])," (",
                             sprintf("%0.3f",mvall$coefficients[1]-1.96*mvall$coefficients[2]),
                             ", ",sprintf("%0.3f",mvall$coefficients[1]+1.96*mvall$coefficients[2]),")")
Coefficient7.2[32,3]<-paste0(sprintf("%0.3f",mvall7.2$coefficients[1,1])," (",
                             sprintf("%0.3f",mvall7.2$coefficients[1,1]-1.96*mvall7.2$coefficients[1,2]),
                             ", ",sprintf("%0.3f",mvall7.2$coefficients[1,1]+1.96*mvall7.2$coefficients[1,2]),")")
Coefficient7.2[32,5]<-paste0(sprintf("%0.3f",mvall7.2$coefficients[2,1])," (",
                             sprintf("%0.3f",mvall7.2$coefficients[2,1]-1.96*mvall7.2$coefficients[2,2]),
                             ", ",sprintf("%0.3f",mvall7.2$coefficients[2,1]+1.96*mvall7.2$coefficients[2,2]),")")
Coefficient7.2[32,6]<-round(mvall7.2$coefficients[2,4],3)

mvall.7<-mvall
mvall7.1<-mvmeta(yori7.1~1, Sori7.1, method="reml")
mvall7.2<-mvmeta(yori7.2~1, Sori7.2, method="reml")
#######################Age group 12, both 2002 and 2009 intervention#######################
setwd("D:/~")
Hepatitis.B.agegroup12<-read.xlsx("agegroup12.xlsx")
Hepatitis.B.agegroup12$Month.factor<-factor(Hepatitis.B.agegroup12$Month)
Hepatitis.B.agegroup12$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup12$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup12$Intervention.2002<-factor(rep(c(rep(0,12*12),rep(1,12*3)),31))
Hepatitis.B.agegroup12$Intervention.2009<-factor(rep(c(rep(0,12*4),rep(0,10),rep(1,2+7*12),rep(0,12*3)),31))
Hepatitis.B.agegroup12$Interaction.2002<-rep(c(rep(0,12*12),c(0:(180-12*12-1))),31)
Hepatitis.B.agegroup12$Interaction.2009<-rep(c(rep(0,4*12),rep(0,10),c(0:(2+7*12-1)),rep(0,3*12)),31)
regions<-as.character(unique(Hepatitis.B.agegroup12$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup12[Hepatitis.B.agegroup12$Province==x, ])
seq(datalist)
yori8.1<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori8.1<-vector("list", length(datalist)); names(Sori8.1) <- regions
yori8.2<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori8.2<-vector("list", length(datalist)); names(Sori8.2) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+Intervention.2009+
              Time+ns(Temperature,df=3)+
              Holiday+Interaction.2002+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori8.1[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun-1),1]
Sori8.1[[i]] <- vcov(mfirst)[c(2,lengthfun-1),c(2,lengthfun-1)]
yori8.2[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(3,lengthfun),1]
Sori8.2[[i]] <- vcov(mfirst)[c(3,lengthfun),c(3,lengthfun)]
assign(paste("model8.",i,sep=""),mfirst)
}

mvall8.1<-mvmeta(yori8.1~1, Sori8.1, method="reml")
summary(mvall8.1)

ER8.1<-ER.meta.2002(mvall8.1,12)
ER8.1$ER.lower<-as.numeric(ER8.1$ER.lower)
ER8.1$ER.upper<-as.numeric(ER8.1$ER.upper)
ER8.1$ER<-as.numeric(ER8.1$ER)

EIR8.1<-EIR.2002(c(1:31),Hepatitis.B.agegroup12,"model8",c(2,10),3,10000)
EIR8.1<-data.frame(EIR8.1)
EIR8.1$EC.point<-sprintf("%0.0f",EIR8.1$EC.point)
EIR8.1$Population<-sprintf("%0.0f",EIR8.1$Population)
EIR8.1$EIR.point<-sprintf("%0.2f",EIR8.1$EIR.point)
EIR8.1$EC.low<-sprintf("%0.0f",EIR8.1$EC.low)
EIR8.1$EIR.low<-sprintf("%0.2f",EIR8.1$EIR.low)
EIR8.1$EC.high<-sprintf("%0.0f",EIR8.1$EC.high)
EIR8.1$EIR.high<-sprintf("%0.2f",EIR8.1$EIR.high)

mvall8.2<-mvmeta(yori8.2~1, Sori8.2, method="reml")
summary(mvall8.2)

ER8.2<-ER.meta.2009(mvall8.2,8,59,2009)
ER8.2$ER.lower<-as.numeric(ER8.2$ER.lower)
ER8.2$ER.upper<-as.numeric(ER8.2$ER.upper)
ER8.2$ER<-as.numeric(ER8.2$ER)

EIR8.2<-EIR.2009(c(1:31),Hepatitis.B.agegroup12,"model8",c(3,11),(2+7*12)/12,10000)
EIR8.2<-data.frame(EIR8.2)
EIR8.2$EC.point<-sprintf("%0.0f",EIR8.2$EC.point)
EIR8.2$Population<-sprintf("%0.0f",EIR8.2$Population)
EIR8.2$EIR.point<-sprintf("%0.2f",EIR8.2$EIR.point)
EIR8.2$EC.low<-sprintf("%0.0f",EIR8.2$EC.low)
EIR8.2$EIR.low<-sprintf("%0.2f",EIR8.2$EIR.low)
EIR8.2$EC.high<-sprintf("%0.0f",EIR8.2$EC.high)
EIR8.2$EIR.high<-sprintf("%0.2f",EIR8.2$EIR.high)

#Calculating the weights with proportion of population size and the standardised errors for each age group
Age.pop2<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup6[-which(Hepatitis.B.agegroup6$Order==1 |
                                                             Hepatitis.B.agegroup6$Order==2 |
                                                             Hepatitis.B.agegroup6$Order==6 |
                                                             Hepatitis.B.agegroup6$Order==9 |
                                                             Hepatitis.B.agegroup6$Order==26),],FUN=sum)[,2])
Age.pop3<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup7[-which(Hepatitis.B.agegroup7$Order==1 |
                                                             Hepatitis.B.agegroup7$Order==2 |
                                                             Hepatitis.B.agegroup7$Order==9 |
                                                             Hepatitis.B.agegroup7$Order==14 |
                                                             Hepatitis.B.agegroup7$Order==26),],FUN=sum)[,2])
Age.pop4<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup8[-which(Hepatitis.B.agegroup8$Order==1 |
                                                             Hepatitis.B.agegroup8$Order==2 |
                                                             Hepatitis.B.agegroup8$Order==9 |
                                                             Hepatitis.B.agegroup8$Order==14 |
                                                             Hepatitis.B.agegroup8$Order==26),],FUN=sum)[,2])
Age.pop5<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup9[-which(Hepatitis.B.agegroup9$Order==1 |
                                                             Hepatitis.B.agegroup9$Order==2 |
                                                             Hepatitis.B.agegroup9$Order==9 |
                                                             Hepatitis.B.agegroup9$Order==14 |
                                                             Hepatitis.B.agegroup9$Order==26),],FUN=sum)[,2])
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[-which(Hepatitis.B.agegroup10$Order==1 |
                                                              Hepatitis.B.agegroup10$Order==2 |
                                                              Hepatitis.B.agegroup10$Order==9 |
                                                              Hepatitis.B.agegroup10$Order==14 |
                                                              Hepatitis.B.agegroup10$Order==26),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,data=Hepatitis.B.agegroup11,FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,data=Hepatitis.B.agegroup12,FUN=sum)[,2])
Age.pop<-c(Age.pop2,Age.pop3,Age.pop4,Age.pop5,Age.pop6,Age.pop7,Age.pop8)
Prop.pop<-c(Age.pop2/sum(Age.pop),Age.pop3/sum(Age.pop),Age.pop4/sum(Age.pop),
            Age.pop5/sum(Age.pop),Age.pop6/sum(Age.pop),Age.pop7/sum(Age.pop),
            Age.pop8/sum(Age.pop))
Merge.coef<-Prop.pop[1]*mvall2$coefficients+Prop.pop[2]*mvall3$coefficients+
  Prop.pop[3]*mvall4$coefficients+Prop.pop[4]*mvall5.1$coefficients+
  Prop.pop[5]*mvall6.1$coefficients+Prop.pop[6]*mvall7.1$coefficients+
  Prop.pop[7]*mvall8.1$coefficients
Merge.vcov<-Prop.pop[1]^2*mvall2$vcov+Prop.pop[2]^2*mvall3$vcov+
  Prop.pop[3]^2*mvall4$vcov+Prop.pop[4]^2*mvall5.1$vcov+
  Prop.pop[5]^2*mvall6.1$vcov+Prop.pop[6]^2*mvall7.1$vcov+
  Prop.pop[7]^2*mvall8.1$vcov

#Age-standardized ERs for 2002 intervention
coef.meta<-as.matrix(t(Merge.coef))
time.meta<-as.matrix(cbind(rep(1,180+3*12),c(1:(180+3*12))-1))
cov.meta<-Merge.vcov
pre.ER.meta<-c()
ER.meta<-c()
se.ER.meta<-c()
ER.low.meta<-c()
ER.high.meta<-c()
for(l in 1:(180+3*12))
{pre.ER.meta[l]<-time.meta[l,]%*%coef.meta
ER.meta[l]<-(exp(pre.ER.meta[l])-1)*100
se.ER.meta[l]<-sqrt(time.meta[l,]%*%cov.meta%*%t(t(time.meta[l,]))) 
ER.low.meta[l]=(exp(pre.ER.meta[l]-1.96*se.ER.meta[l])-1)*100
ER.high.meta[l]=(exp(pre.ER.meta[l]+1.96*se.ER.meta[l])-1)*100
}
ER.whole<-ER.meta
ER.whole.low<-ER.low.meta
ER.whole.upper<-ER.high.meta
month.meta<-c(1:(180+3*12))
ER<-data.frame(Time=month.meta,ER=ER.meta,ER.lower=ER.low.meta,ER.upper=ER.high.meta)
ER$ER<-sprintf("%0.2f",ER$ER)
ER$ER.lower<-sprintf("%0.2f",ER$ER.lower)
ER$ER.upper<-sprintf("%0.2f",ER$ER.upper)

Coefficient8.1<-matrix(0,nrow=32,ncol=6)
Coefficient8.1<-as.data.frame(Coefficient8.1)
colnames(Coefficient8.1)<-c("Age","PLADs","RR1 (95%CI)","RR2 (95%CI)","RR3 (95%CI)","Trend")
Coefficient8.1$Age<-c("15~19")
yori<-matrix(0,length(datalist),1,dimnames=list(regions, paste("beta",1,sep="")))
Sori<-vector("list", length(datalist)); names(Sori) <- regions
for (i in 1:31)
{model.number<-eval(parse(text=paste("model8.",i,sep="")))
Point1<-model.number$coefficients[2]
Se1<-summary(model.number)$se[2]
Point2<-model.number$coefficients[4]
Se2<-summary(model.number)$se[4]
length.coef<-length(model.number$coefficients)
Point3<-model.number$coefficients[length.coef-1]
Se3<-summary(model.number)$se[length.coef-1]
Coefficient8.1[i,2]<-Name[i]
Coefficient8.1[i,3]<-paste0(sprintf("%0.3f",Point1)," (",sprintf("%0.3f",Point1-1.96*Se1),
                            ", ",sprintf("%0.3f",Point1+1.96*Se1),")")
Coefficient8.1[i,4]<-paste0(sprintf("%0.3f",Point2)," (",sprintf("%0.3f",Point2-1.96*Se2),
                            ", ",sprintf("%0.3f",Point2+1.96*Se2),")") 
Coefficient8.1[i,5]<-paste0(sprintf("%0.3f",Point3)," (",sprintf("%0.3f",Point3-1.96*Se3),
                            ", ",sprintf("%0.3f",Point3+1.96*Se3),")")
Coefficient8.1[i,6]<-round(summary(model.number)$p.pv[length.coef],3)
yori[i,]<-as.data.frame(summary(model.number)$p.coeff)[c(4),1]
Sori[[i]]<-vcov(model.number)[c(4),c(4)]
}
mvall<-mvmeta(yori~1, Sori, method="reml")
mvall<-summary(mvall)
mvall8.1<-summary(mvall8.1)
Coefficient8.1[32,2]<-c("Chinese mainland")
Coefficient8.1[32,4]<-paste0(sprintf("%0.3f",mvall$coefficients[1])," (",
                             sprintf("%0.3f",mvall$coefficients[1]-1.96*mvall$coefficients[2]),
                             ", ",sprintf("%0.3f",mvall$coefficients[1]+1.96*mvall$coefficients[2]),")")
Coefficient8.1[32,3]<-paste0(sprintf("%0.3f",mvall8.1$coefficients[1,1])," (",
                             sprintf("%0.3f",mvall8.1$coefficients[1,1]-1.96*mvall8.1$coefficients[1,2]),
                             ", ",sprintf("%0.3f",mvall8.1$coefficients[1,1]+1.96*mvall8.1$coefficients[1,2]),")")
Coefficient8.1[32,5]<-paste0(sprintf("%0.3f",mvall8.1$coefficients[2,1])," (",
                             sprintf("%0.3f",mvall8.1$coefficients[2,1]-1.96*mvall8.1$coefficients[2,2]),
                             ", ",sprintf("%0.3f",mvall8.1$coefficients[2,1]+1.96*mvall8.1$coefficients[2,2]),")")
Coefficient8.1[32,6]<-round(mvall8.1$coefficients[2,4],3)

Coefficient8.2<-matrix(0,nrow=32,ncol=6)
Coefficient8.2<-as.data.frame(Coefficient8.2)
colnames(Coefficient8.2)<-c("Age","PLADs","RR1 (95%CI)","RR2 (95%CI)","RR3 (95%CI)","Trend")
Coefficient8.2$Age<-c("15~19")
yori<-matrix(0,length(datalist),1,dimnames=list(regions, paste("beta",1,sep="")))
Sori<-vector("list", length(datalist)); names(Sori) <- regions
for (i in 1:31)
{model.number<-eval(parse(text=paste("model8.",i,sep="")))
Point1<-model.number$coefficients[3]
Se1<-summary(model.number)$se[3]
Point2<-model.number$coefficients[4]
Se2<-summary(model.number)$se[4]
length.coef<-length(model.number$coefficients)
Point3<-model.number$coefficients[length.coef]
Se3<-summary(model.number)$se[length.coef]
Coefficient8.2[i,2]<-Name[i]
Coefficient8.2[i,3]<-paste0(sprintf("%0.3f",Point1)," (",sprintf("%0.3f",Point1-1.96*Se1),
                            ", ",sprintf("%0.3f",Point1+1.96*Se1),")")
Coefficient8.2[i,4]<-paste0(sprintf("%0.3f",Point2)," (",sprintf("%0.3f",Point2-1.96*Se2),
                            ", ",sprintf("%0.3f",Point2+1.96*Se2),")") 
Coefficient8.2[i,5]<-paste0(sprintf("%0.3f",Point3)," (",sprintf("%0.3f",Point3-1.96*Se3),
                            ", ",sprintf("%0.3f",Point3+1.96*Se3),")")
Coefficient8.2[i,6]<-round(summary(model.number)$p.pv[length.coef],3)
yori[i,]<-as.data.frame(summary(model.number)$p.coeff)[c(4),1]
Sori[[i]]<-vcov(model.number)[c(4),c(4)]
}
mvall<-mvmeta(yori~1, Sori, method="reml")
mvall<-summary(mvall)
mvall8.2<-summary(mvall8.2)
Coefficient8.2[32,2]<-c("Chinese mainland")
Coefficient8.2[32,4]<-paste0(sprintf("%0.3f",mvall$coefficients[1])," (",
                             sprintf("%0.3f",mvall$coefficients[1]-1.96*mvall$coefficients[2]),
                             ", ",sprintf("%0.3f",mvall$coefficients[1]+1.96*mvall$coefficients[2]),")")
Coefficient8.2[32,3]<-paste0(sprintf("%0.3f",mvall8.2$coefficients[1,1])," (",
                             sprintf("%0.3f",mvall8.2$coefficients[1,1]-1.96*mvall8.2$coefficients[1,2]),
                             ", ",sprintf("%0.3f",mvall8.2$coefficients[1,1]+1.96*mvall8.2$coefficients[1,2]),")")
Coefficient8.2[32,5]<-paste0(sprintf("%0.3f",mvall8.2$coefficients[2,1])," (",
                             sprintf("%0.3f",mvall8.2$coefficients[2,1]-1.96*mvall8.2$coefficients[2,2]),
                             ", ",sprintf("%0.3f",mvall8.2$coefficients[2,1]+1.96*mvall8.2$coefficients[2,2]),")")
Coefficient8.2[32,6]<-round(mvall8.2$coefficients[2,4],3)

#Calculating the age-standardised coefficients
mvall.8<-mvall
Pre.coef<-Prop.pop[1]*mvall.2$coefficients[1]+Prop.pop[2]*mvall.3$coefficients[1]+
  Prop.pop[3]*mvall.4$coefficients[1]+Prop.pop[4]*mvall.5$coefficients[1]+
  Prop.pop[5]*mvall.6$coefficients[1]+Prop.pop[6]*mvall.7$coefficients[1]+
  Prop.pop[7]*mvall.8$coefficients[1]
Pre.vcov<-Prop.pop[1]^2*mvall.2$vcov+Prop.pop[2]^2*mvall.3$vcov+
  Prop.pop[3]^2*mvall.4$vcov+Prop.pop[4]^2*mvall.5$vcov+
  Prop.pop[5]^2*mvall.6$vcov+Prop.pop[6]^2*mvall.7$vcov+
  Prop.pop[7]^2*mvall.8$vcov
Pre.se<-sqrt(Pre.vcov)
Level.coef<-Prop.pop[1]*mvall2$coefficients[1]+Prop.pop[2]*mvall3$coefficients[1]+
  Prop.pop[3]*mvall4$coefficients[1]+Prop.pop[4]*mvall5.1$coefficients[1]+
  Prop.pop[5]*mvall6.1$coefficients[1]+Prop.pop[6]*mvall7.1$coefficients[1]+
  Prop.pop[7]*mvall8.1$coefficients[1]
Level.vcov<-Prop.pop[1]^2*mvall2$vcov[1,1]+Prop.pop[2]^2*mvall3$vcov[1,1]+
  Prop.pop[3]^2*mvall4$vcov[1,1]+Prop.pop[4]^2*mvall5.1$vcov[1,1]+
  Prop.pop[5]^2*mvall6.1$vcov[1,1]+Prop.pop[6]^2*mvall7.1$vcov[1,1]+
  Prop.pop[7]^2*mvall8.1$vcov[1,1]
Level.se<-sqrt(Level.vcov)
Post.coef<-Prop.pop[1]*mvall2$coefficients[2]+Prop.pop[2]*mvall3$coefficients[2]+
  Prop.pop[3]*mvall4$coefficients[2]+Prop.pop[4]*mvall5.1$coefficients[2]+
  Prop.pop[5]*mvall6.1$coefficients[2]+Prop.pop[6]*mvall7.1$coefficients[2]+
  Prop.pop[7]*mvall8.1$coefficients[2]
Post.vcov<-Prop.pop[1]^2*mvall2$vcov[2,2]+Prop.pop[2]^2*mvall3$vcov[2,2]+
  Prop.pop[3]^2*mvall4$vcov[2,2]+Prop.pop[4]^2*mvall5.1$vcov[2,2]+
  Prop.pop[5]^2*mvall6.1$vcov[2,2]+Prop.pop[6]^2*mvall7.1$vcov[2,2]+
  Prop.pop[7]^2*mvall8.1$vcov[2,2]
Post.se<-sqrt(Post.vcov)
Pre.2002<-paste0(sprintf("%0.3f",Pre.coef)," (",
                 sprintf("%0.3f",Pre.coef-1.96*Pre.se),
                 ", ",sprintf("%0.3f",Pre.coef+1.96*Pre.se),")")
Post.2002<-paste0(sprintf("%0.3f",Post.coef)," (",
                  sprintf("%0.3f",Post.coef-1.96*Post.se),
                  ", ",sprintf("%0.3f",Post.coef+1.96*Post.se),")")
Level.2002<-paste0(sprintf("%0.3f",Level.coef)," (",
                   sprintf("%0.3f",Level.coef-1.96*Level.se),
                   ", ",sprintf("%0.3f",Level.coef+1.96*Level.se),")")

mvall8.2<-mvmeta(yori8.2~1, Sori8.2, method="reml")
mvall8.1<-mvmeta(yori8.1~1, Sori8.1, method="reml")
#######################Age group 13, only 2009 intervention##############################
setwd("D:/~")
Hepatitis.B.agegroup13<-read.xlsx("Agegroup13.xlsx")
Hepatitis.B.agegroup13$Month.factor<-factor(Hepatitis.B.agegroup13$Month)
Hepatitis.B.agegroup13$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup13$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
#2009 intervention term
Hepatitis.B.agegroup13$Intervention.2009<-factor(rep(c(rep(0,9*12),rep(1,12*6)),31))
#2009 interaction term
Hepatitis.B.agegroup13$Interaction.2009<-rep(c(rep(0,9*12),c(0:(180-12*9-1))),31)
regions<-as.character(unique(Hepatitis.B.agegroup13$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup13[Hepatitis.B.agegroup13$Province==x, ])
seq(datalist)
yori9<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori9<-vector("list", length(datalist)); names(Sori9) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup13[which(Hepatitis.B.agegroup13$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2009+Time+
              ns(Temperature,df=3)+
              Holiday+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori9[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun),1]
Sori9[[i]] <- vcov(mfirst)[c(2,lengthfun),c(2,lengthfun)]
assign(paste("model9.",i,sep=""),mfirst)
}

mvall9<-mvmeta(yori9~1, Sori9, method="reml")
summary(mvall9)

ER9<-ER.meta.2009(mvall9,6,109,2014)
ER9$ER.lower<-as.numeric(ER9$ER.lower)
ER9$ER.upper<-as.numeric(ER9$ER.upper)
ER9$ER<-as.numeric(ER9$ER)

EIR9<-EIR.2009(c(1:31),Hepatitis.B.agegroup13,"model9",c(2,9),6,10000)
EIR9<-data.frame(EIR9)
EIR9$EC.point<-sprintf("%0.0f",EIR9$EC.point)
EIR9$Population<-sprintf("%0.0f",EIR9$Population)
EIR9$EIR.point<-sprintf("%0.2f",EIR9$EIR.point)
EIR9$EC.low<-sprintf("%0.0f",EIR9$EC.low)
EIR9$EIR.low<-sprintf("%0.2f",EIR9$EIR.low)
EIR9$EC.high<-sprintf("%0.0f",EIR9$EC.high)
EIR9$EIR.high<-sprintf("%0.2f",EIR9$EIR.high)

Coefficient9<-matrix(0,nrow=32,ncol=6)
Coefficient9<-as.data.frame(Coefficient9)
colnames(Coefficient9)<-c("Age","PLADs","RR1 (95%CI)","RR2 (95%CI)","RR3 (95%CI)","Trend")
Coefficient9$Age<-c("20~24")
yori<-matrix(0,length(datalist),1,dimnames=list(regions, paste("beta",1,sep="")))
Sori<-vector("list", length(datalist)); names(Sori) <- regions
for (i in 1:31)
{model.number<-eval(parse(text=paste("model9.",i,sep="")))
Point1<-model.number$coefficients[2]
Se1<-summary(model.number)$se[2]
Point2<-model.number$coefficients[3]
Se2<-summary(model.number)$se[3]
length.coef<-length(model.number$coefficients)
Point3<-model.number$coefficients[length.coef]
Se3<-summary(model.number)$se[length.coef]
Coefficient9[i,2]<-Name[i]
Coefficient9[i,3]<-paste0(sprintf("%0.3f",Point1)," (",sprintf("%0.3f",Point1-1.96*Se1),
                          ", ",sprintf("%0.3f",Point1+1.96*Se1),")")
Coefficient9[i,4]<-paste0(sprintf("%0.3f",Point2)," (",sprintf("%0.3f",Point2-1.96*Se2),
                          ", ",sprintf("%0.3f",Point2+1.96*Se2),")") 
Coefficient9[i,5]<-paste0(sprintf("%0.3f",Point3)," (",sprintf("%0.3f",Point3-1.96*Se3),
                          ", ",sprintf("%0.3f",Point3+1.96*Se3),")")
Coefficient9[i,6]<-round(summary(model.number)$p.pv[length.coef],3)
yori[i,]<-as.data.frame(summary(model.number)$p.coeff)[c(3),1]
Sori[[i]]<-vcov(model.number)[c(3),c(3)]
}
mvall<-mvmeta(yori~1, Sori, method="reml")
mvall<-summary(mvall)
mvall9<-summary(mvall9)
Coefficient9[32,2]<-c("Chinese mainland")
Coefficient9[32,4]<-paste0(sprintf("%0.3f",mvall$coefficients[1])," (",
                           sprintf("%0.3f",mvall$coefficients[1]-1.96*mvall$coefficients[2]),
                           ", ",sprintf("%0.3f",mvall$coefficients[1]+1.96*mvall$coefficients[2]),")")
Coefficient9[32,3]<-paste0(sprintf("%0.3f",mvall9$coefficients[1,1])," (",
                           sprintf("%0.3f",mvall9$coefficients[1,1]-1.96*mvall9$coefficients[1,2]),
                           ", ",sprintf("%0.3f",mvall9$coefficients[1,1]+1.96*mvall9$coefficients[1,2]),")")
Coefficient9[32,5]<-paste0(sprintf("%0.3f",mvall9$coefficients[2,1])," (",
                           sprintf("%0.3f",mvall9$coefficients[2,1]-1.96*mvall9$coefficients[2,2]),
                           ", ",sprintf("%0.3f",mvall9$coefficients[2,1]+1.96*mvall9$coefficients[2,2]),")")
Coefficient9[32,6]<-round(mvall9$coefficients[2,4],3)

mvall.9<-mvall
mvall9<-mvmeta(yori9~1, Sori9, method="reml")
#######################Age group 14, only 2009 intervention##############################
setwd("D:/~")
Hepatitis.B.agegroup14<-read.xlsx("Agegroup14.xlsx")
Hepatitis.B.agegroup14$Month.factor<-factor(Hepatitis.B.agegroup14$Month)
Hepatitis.B.agegroup14$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup14$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup14$Intervention.2009<-factor(rep(c(rep(0,14*12),rep(1,12*1)),31))
Hepatitis.B.agegroup14$Interaction.2009<-rep(c(rep(0,14*12),c(0:(180-12*14-1))),31)
regions<-as.character(unique(Hepatitis.B.agegroup14$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup14[Hepatitis.B.agegroup14$Province==x, ])
seq(datalist)
yori10<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori10<-vector("list", length(datalist)); names(Sori10) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup14[which(Hepatitis.B.agegroup14$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2009+Time+
              ns(Temperature,df=3)+
              Holiday+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori10[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun),1]
Sori10[[i]] <- vcov(mfirst)[c(2,lengthfun),c(2,lengthfun)]
assign(paste("model10.",i,sep=""),mfirst)
}

mvall10<-mvmeta(yori10~1, Sori10, method="reml")
summary(mvall10)

ER10<-ER.meta.2009(mvall10,1,169,2019)
ER10$ER.lower<-as.numeric(ER10$ER.lower)
ER10$ER.upper<-as.numeric(ER10$ER.upper)
ER10$ER<-as.numeric(ER10$ER)

EIR10<-EIR.2009(c(1:31),Hepatitis.B.agegroup14,"model10",c(2,9),1,10000)
EIR10<-data.frame(EIR10)
EIR10$EC.point<-sprintf("%0.0f",EIR10$EC.point)
EIR10$Population<-sprintf("%0.0f",EIR10$Population)
EIR10$EIR.point<-sprintf("%0.2f",EIR10$EIR.point)
EIR10$EC.low<-sprintf("%0.0f",EIR10$EC.low)
EIR10$EIR.low<-sprintf("%0.2f",EIR10$EIR.low)
EIR10$EC.high<-sprintf("%0.0f",EIR10$EC.high)
EIR10$EIR.high<-sprintf("%0.2f",EIR10$EIR.high)

Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[-which(Hepatitis.B.agegroup10$Order==1 |
                                                              Hepatitis.B.agegroup10$Order==2 |
                                                              Hepatitis.B.agegroup10$Order==9 |
                                                              Hepatitis.B.agegroup10$Order==26),],FUN=sum)[,2])
Age.pop9<-mean(aggregate(Population.book~Time,data=Hepatitis.B.agegroup13,FUN=sum)[,2])
Age.pop10<-mean(aggregate(Population.book~Time,data=Hepatitis.B.agegroup14,FUN=sum)[,2])
Age.pop.2009<-c(Age.pop6,Age.pop7,Age.pop8,Age.pop9,Age.pop10)
Prop.pop.2009<-c(Age.pop6/sum(Age.pop.2009),Age.pop7/sum(Age.pop.2009),
                 Age.pop8/sum(Age.pop.2009),Age.pop9/sum(Age.pop.2009),Age.pop10/sum(Age.pop.2009))
Merge.coef.2009<-Prop.pop.2009[1]*mvall6.2$coefficients+
  Prop.pop.2009[2]*mvall7.2$coefficients+Prop.pop.2009[3]*mvall8.2$coefficients+
  Prop.pop.2009[4]*mvall9$coefficients+Prop.pop.2009[5]*mvall10$coefficients
Merge.vcov.2009<-Prop.pop.2009[1]^2*mvall6.2$vcov+
  Prop.pop.2009[2]^2*mvall7.2$vcov+Prop.pop.2009[3]^2*mvall8.2$vcov+
  Prop.pop.2009[4]^2*mvall9$vcov+Prop.pop.2009[5]^2*mvall10$vcov

#Age-standardized ERs for 2009 intervention
coef.meta<-as.matrix(t(Merge.coef.2009))
time.meta<-as.matrix(cbind(rep(1,180-7*12),c((7*12):179)-7*12))
cov.meta<-Merge.vcov.2009
pre.ER.meta<-c()
ER.meta<-c()
se.ER.meta<-c()
ER.low.meta<-c()
ER.high.meta<-c()
for(l in 1:(180-7*12))
{pre.ER.meta[l]<-time.meta[l,]%*%coef.meta
ER.meta[l]<-(exp(pre.ER.meta[l])-1)*100
se.ER.meta[l]<-sqrt(time.meta[l,]%*%cov.meta%*%t(t(time.meta[l,]))) 
ER.low.meta[l]=(exp(pre.ER.meta[l]-1.96*se.ER.meta[l])-1)*100
ER.high.meta[l]=(exp(pre.ER.meta[l]+1.96*se.ER.meta[l])-1)*100
}
ER.whole<-ER.meta
ER.whole.low<-ER.low.meta
ER.whole.upper<-ER.high.meta
month.meta<-c((7*12+1):180)
ER<-data.frame(Time=month.meta,ER=ER.meta,ER.lower=ER.low.meta,ER.upper=ER.high.meta)
ER$ER<-sprintf("%0.2f",ER$ER)
ER$ER.lower<-sprintf("%0.2f",ER$ER.lower)
ER$ER.upper<-sprintf("%0.2f",ER$ER.upper)

Coefficient10<-matrix(0,nrow=32,ncol=6)
Coefficient10<-as.data.frame(Coefficient10)
colnames(Coefficient10)<-c("Age","PLADs","RR1 (95%CI)","RR2 (95%CI)","RR3 (95%CI)","Trend")
Coefficient10$Age<-c("25~29")
yori<-matrix(0,length(datalist),1,dimnames=list(regions, paste("beta",1,sep="")))
Sori<-vector("list", length(datalist)); names(Sori) <- regions
for (i in 1:31)
{model.number<-eval(parse(text=paste("model10.",i,sep="")))
Point1<-model.number$coefficients[2]
Se1<-summary(model.number)$se[2]
Point2<-model.number$coefficients[3]
Se2<-summary(model.number)$se[3]
length.coef<-length(model.number$coefficients)
Point3<-model.number$coefficients[length.coef]
Se3<-summary(model.number)$se[length.coef]
Coefficient10[i,2]<-Name[i]
Coefficient10[i,3]<-paste0(sprintf("%0.3f",Point1)," (",sprintf("%0.3f",Point1-1.96*Se1),
                           ", ",sprintf("%0.3f",Point1+1.96*Se1),")")
Coefficient10[i,4]<-paste0(sprintf("%0.3f",Point2)," (",sprintf("%0.3f",Point2-1.96*Se2),
                           ", ",sprintf("%0.3f",Point2+1.96*Se2),")") 
Coefficient10[i,5]<-paste0(sprintf("%0.3f",Point3)," (",sprintf("%0.3f",Point3-1.96*Se3),
                           ", ",sprintf("%0.3f",Point3+1.96*Se3),")")
Coefficient10[i,6]<-round(summary(model.number)$p.pv[length.coef],3)
yori[i,]<-as.data.frame(summary(model.number)$p.coeff)[c(3),1]
Sori[[i]]<-vcov(model.number)[c(3),c(3)]
}
mvall<-mvmeta(yori~1, Sori, method="reml")
mvall<-summary(mvall)
mvall10<-summary(mvall10)
Coefficient10[32,2]<-c("Chinese mainland")
Coefficient10[32,4]<-paste0(sprintf("%0.3f",mvall$coefficients[1])," (",
                            sprintf("%0.3f",mvall$coefficients[1]-1.96*mvall$coefficients[2]),
                            ", ",sprintf("%0.3f",mvall$coefficients[1]+1.96*mvall$coefficients[2]),")")
Coefficient10[32,3]<-paste0(sprintf("%0.3f",mvall10$coefficients[1,1])," (",
                            sprintf("%0.3f",mvall10$coefficients[1,1]-1.96*mvall10$coefficients[1,2]),
                            ", ",sprintf("%0.3f",mvall10$coefficients[1,1]+1.96*mvall10$coefficients[1,2]),")")
Coefficient10[32,5]<-paste0(sprintf("%0.3f",mvall10$coefficients[2,1])," (",
                            sprintf("%0.3f",mvall10$coefficients[2,1]-1.96*mvall10$coefficients[2,2]),
                            ", ",sprintf("%0.3f",mvall10$coefficients[2,1]+1.96*mvall10$coefficients[2,2]),")")
Coefficient10[32,6]<-round(mvall10$coefficients[2,4],3)

mvall.10<-mvall
Pre.coef<-Prop.pop.2009[1]*mvall.6$coefficients[1]+
  Prop.pop.2009[2]*mvall.7$coefficients[1]+Prop.pop.2009[3]*mvall.8$coefficients[1]+
  Prop.pop.2009[4]*mvall.9$coefficients[1]+Prop.pop.2009[5]*mvall.10$coefficients[1]
Pre.vcov<-Prop.pop.2009[1]^2*mvall.6$vcov+
  Prop.pop.2009[2]^2*mvall.7$vcov+Prop.pop.2009[3]^2*mvall.8$vcov+
  Prop.pop.2009[4]^2*mvall.9$vcov+Prop.pop.2009[5]^2*mvall.10$vcov
Pre.se<-sqrt(Pre.vcov)
Level.coef<-Prop.pop.2009[1]*mvall6.2$coefficients[1]+
  Prop.pop.2009[2]*mvall7.2$coefficients[1]+Prop.pop.2009[5]*mvall8.2$coefficients[1]+
  Prop.pop.2009[4]*mvall9$coefficients[1]+Prop.pop.2009[5]*mvall10$coefficients[1]
Level.vcov<-Prop.pop.2009[1]^2*mvall6.2$vcov[1,1]+
  Prop.pop.2009[2]^2*mvall7.2$vcov[1,1]+Prop.pop.2009[5]^2*mvall8.2$vcov[1,1]+
  Prop.pop.2009[4]^2*mvall9$vcov[1,1]+Prop.pop.2009[5]^2*mvall10$vcov[1,1]
Level.se<-sqrt(Level.vcov)
Post.coef<-Prop.pop.2009[1]*mvall6.2$coefficients[2]+
  Prop.pop.2009[2]*mvall7.2$coefficients[2]+Prop.pop.2009[3]*mvall8.2$coefficients[2]+
  Prop.pop.2009[4]*mvall9$coefficients[2]+Prop.pop.2009[5]*mvall10$coefficients[2]
Post.vcov<-Prop.pop.2009[1]^2*mvall6.2$vcov[2,2]+
  Prop.pop.2009[2]^2*mvall7.2$vcov[2,2]+Prop.pop.2009[3]^2*mvall8.2$vcov[2,2]+
  Prop.pop.2009[4]^2*mvall9$vcov[2,2]+Prop.pop.2009[5]^2*mvall10$vcov[2,2]
Post.se<-sqrt(Post.vcov)
Pre.2009<-paste0(sprintf("%0.3f",Pre.coef)," (",
                 sprintf("%0.3f",Pre.coef-1.96*Pre.se),
                 ", ",sprintf("%0.3f",Pre.coef+1.96*Pre.se),")")
Post.2009<-paste0(sprintf("%0.3f",Post.coef)," (",
                  sprintf("%0.3f",Post.coef-1.96*Post.se),
                  ", ",sprintf("%0.3f",Post.coef+1.96*Post.se),")")
Level.2009<-paste0(sprintf("%0.3f",Level.coef)," (",
                   sprintf("%0.3f",Level.coef-1.96*Level.se),
                   ", ",sprintf("%0.3f",Level.coef+1.96*Level.se),")")

mvall10<-mvmeta(yori10~1, Sori10, method="reml")
#######################Age-standardised EIR########################################
#2002 intervention
EIR.merge.2002<-function(Sequence,database,model,Location,intervention,n.cir){
  set.seed(202403)
  EMR.reserve<-matrix(0,nrow=31,ncol=8)
  EMR.reserve.diff<-c()
  MA.reserve<-matrix(0,nrow=n.cir,ncol=31)
  for(i in Sequence)
  {data.number<-database[which(database$Order==i),]
  model.number<-eval(parse(text=paste(model,".",i,sep="")))
  mod.number<-as.matrix(data.frame(model.matrix(model.number),offset=log(data.number$Population.book)))
  mod.number.0<-mod.number
  mod.number.0[,Location]<-0
  line.row<-which(data.number$Intervention.2002==1)
  length.row<-length(line.row)
  freq.number<-data.number$Case[line.row]
  pop.number<-data.number$Population.book[line.row]
  coef.number<-c(summary(model.number)$p.coeff,1)
  diff.number<-exp(mod.number%*% coef.number)-exp(mod.number.0%*%coef.number)
  diff.number<-diff.number[which(diff.number!=0)]
  EIR.number<-sum(diff.number)/intervention/(sum(pop.number)/length(pop.number))*100000
  EMR.reserve[i,1]<-sum(diff.number)
  EMR.reserve[i,2]<-intervention
  EMR.reserve[i,3]<-sum(pop.number)/length(pop.number)
  EMR.reserve[i,4]<-EIR.number
  coef.ini<-coef.number[Location]
  length.coef<-length(coef.ini)
  cov.number<-vcov(model.number)[Location,Location]
  eigen.number<-eigen(cov.number)
  r.norm<-rnorm(length.coef*n.cir)
  r.norm<-matrix(r.norm,n.cir)
  coef.sim<-coef.ini+eigen.number$vectors%*%diag(sqrt(eigen.number$values),length.coef)%*%t(r.norm)
  for(k in 1:n.cir){
    coef.number[Location[1]]<-coef.sim[1,k]
    coef.number[Location[2]]<-coef.sim[2,k]
    diff.number<-exp(mod.number%*%coef.number)-exp(mod.number.0%*%coef.number)
    diff.number<-diff.number[which(diff.number!=0)]
    EMR.reserve.diff[k]<-sum(diff.number)
  }
  MA.reserve[,i]<-EMR.reserve.diff
  EMR.reserve[i,5]<-sort(EMR.reserve.diff)[n.cir*0.025]
  EMR.reserve[i,7]<-sort(EMR.reserve.diff)[n.cir*0.975]
  EMR.reserve[i,6]<-EMR.reserve[i,5]/EMR.reserve[i,2]/EMR.reserve[i,3]*100000
  EMR.reserve[i,8]<-EMR.reserve[i,7]/EMR.reserve[i,2]/EMR.reserve[i,3]*100000
  }
  EMR.reserve[,2]<-intervention
  EMR.mainland.point<-sum(EMR.reserve[Sequence,1])/mean(EMR.reserve[Sequence,2])/sum(EMR.reserve[Sequence,3])*100000
  EMR.mainland<-rowSums(MA.reserve[,Sequence])/mean(EMR.reserve[Sequence,2])/sum(EMR.reserve[Sequence,3])*100000
  EC.mainland<-rowSums(MA.reserve[,Sequence])
  EIR.mainland<-data.frame(EC=EC.mainland,EIR=EMR.mainland,
                           Point.EC=sum(EMR.reserve[Sequence,1]),
                           Point.EIR=EMR.mainland.point)
  return(EIR.mainland)
}

#EIR for each age group
EIR.merge.2<-EIR.merge.2002(c(3:5,7,8,10:25,27:31),Hepatitis.B.agegroup6,"model2",c(2,9),13,10000)
EIR.merge.3<-EIR.merge.2002(c(3:8,10:13,15:25,27:31),Hepatitis.B.agegroup7,"model3",c(2,9),12,10000)
EIR.merge.4<-EIR.merge.2002(c(3:8,10:13,15:25,27:31),Hepatitis.B.agegroup8,"model4",c(2,9),11,10000)
EIR.merge.5.1<-EIR.merge.2002(c(3:8,10:13,15:25,27:31),Hepatitis.B.agegroup9,"model5",c(2,10),10,10000)
EIR.merge.6.1<-EIR.merge.2002(c(3:8,10:13,15:25,27:31),Hepatitis.B.agegroup10,"model6",c(2,10),9,10000)
EIR.merge.7.1<-EIR.merge.2002(c(1:31),Hepatitis.B.agegroup11,"model7",c(2,10),8,10000)
EIR.merge.8.1<-EIR.merge.2002(c(1:31),Hepatitis.B.agegroup12,"model8",c(2,10),3,10000)

#Calculating the weights with the proportion of population size for each age group
EIR.merge.2002<-Prop.pop[1]*EIR.merge.2$EIR+Prop.pop[2]*EIR.merge.3$EIR+
  Prop.pop[3]*EIR.merge.4$EIR+Prop.pop[4]*EIR.merge.5.1$EIR+
  Prop.pop[5]*EIR.merge.6.1$EIR+Prop.pop[6]*EIR.merge.7.1$EIR+
  Prop.pop[7]*EIR.merge.8.1$EIR
EIR.point.2002<-Prop.pop[1]*EIR.merge.2$Point.EIR[1]+Prop.pop[2]*EIR.merge.3$Point.EIR[1]+
  Prop.pop[3]*EIR.merge.4$Point.EIR[1]+Prop.pop[4]*EIR.merge.5.1$Point.EIR[1]+
  Prop.pop[5]*EIR.merge.6.1$Point.EIR[1]+Prop.pop[6]*EIR.merge.7.1$Point.EIR[1]+
  Prop.pop[7]*EIR.merge.8.1$Point.EIR[1]
#-24.66 (-25.01, -13.74)
EIR.2002.low<-sort(EIR.merge.2002)[n.cir*0.025]
EIR.2002.high<-sort(EIR.merge.2002)[n.cir*0.975]
EC.merge.2002<-EIR.merge.2$EC+EIR.merge.3$EC+
  EIR.merge.4$EC+EIR.merge.5.1$EC+
  EIR.merge.6.1$EC+EIR.merge.7.1$EC+
  EIR.merge.8.1$EC
EC.point.2002<-EIR.merge.2$Point.EC[1]+EIR.merge.3$Point.EC[1]+
  EIR.merge.4$Point.EC[1]+EIR.merge.5.1$Point.EC[1]+
  EIR.merge.6.1$Point.EC[1]+EIR.merge.7.1$Point.EC[1]+
  EIR.merge.8.1$Point.EC[1]
#-337277 (-340581, -10099)
EC.2002.low<-sort(EC.merge.2002)[n.cir*0.025]
EC.2002.high<-sort(EC.merge.2002)[n.cir*0.975]

#2009 intervention
EIR.merge.2009<-function(Sequence,database,model,Location,intervention,n.cir){
  set.seed(202403)
  EMR.reserve<-matrix(0,nrow=31,ncol=8)
  EMR.reserve.diff<-c()
  MA.reserve<-matrix(0,nrow=n.cir,ncol=31)
  for(i in Sequence)
  {data.number<-database[which(database$Order==i),]
  model.number<-eval(parse(text=paste(model,".",i,sep="")))
  mod.number<-as.matrix(data.frame(model.matrix(model.number),offset=log(data.number$Population.book)))
  mod.number.0<-mod.number
  mod.number.0[,Location]<-0
  line.row<-which(data.number$Intervention.2009==1)
  length.row<-length(line.row)
  freq.number<-data.number$Case[line.row]
  pop.number<-data.number$Population.book[line.row]
  coef.number<-c(summary(model.number)$p.coeff,1)
  diff.number<-exp(mod.number%*% coef.number)-exp(mod.number.0%*%coef.number)
  diff.number<-diff.number[which(diff.number!=0)]
  EIR.number<-sum(diff.number)/intervention/(sum(pop.number)/length(pop.number))*100000
  EMR.reserve[i,1]<-sum(diff.number)
  EMR.reserve[i,2]<-intervention
  EMR.reserve[i,3]<-sum(pop.number)/length(pop.number)
  EMR.reserve[i,4]<-EIR.number
  coef.ini<-coef.number[Location]
  length.coef<-length(coef.ini)
  cov.number<-vcov(model.number)[Location,Location]
  eigen.number<-eigen(cov.number)
  r.norm<-rnorm(length.coef*n.cir)
  r.norm<-matrix(r.norm,n.cir)
  coef.sim<-coef.ini+eigen.number$vectors%*%diag(sqrt(eigen.number$values),length.coef)%*%t(r.norm)
  for(k in 1:n.cir){
    coef.number[Location[1]]<-coef.sim[1,k]
    coef.number[Location[2]]<-coef.sim[2,k]
    diff.number<-exp(mod.number%*%coef.number)-exp(mod.number.0%*%coef.number)
    diff.number<-diff.number[which(diff.number!=0)]
    EMR.reserve.diff[k]<-sum(diff.number)
  }
  MA.reserve[,i]<-EMR.reserve.diff
  EMR.reserve[i,5]<-sort(EMR.reserve.diff)[n.cir*0.025]
  EMR.reserve[i,7]<-sort(EMR.reserve.diff)[n.cir*0.975]
  EMR.reserve[i,6]<-EMR.reserve[i,5]/EMR.reserve[i,2]/EMR.reserve[i,3]*100000
  EMR.reserve[i,8]<-EMR.reserve[i,7]/EMR.reserve[i,2]/EMR.reserve[i,3]*100000
  }
  EMR.mainland.point<-sum(EMR.reserve[Sequence,1])/mean(EMR.reserve[Sequence,2])/sum(EMR.reserve[Sequence,3])*100000
  EMR.mainland.point<-sum(EMR.reserve[Sequence,1])/mean(EMR.reserve[Sequence,2])/sum(EMR.reserve[Sequence,3])*100000
  EMR.mainland<-rowSums(MA.reserve[,Sequence])/mean(EMR.reserve[Sequence,2])/sum(EMR.reserve[Sequence,3])*100000
  EC.mainland<-rowSums(MA.reserve[,Sequence])
  EIR.mainland<-data.frame(EC=EC.mainland,EIR=EMR.mainland,
                           Point.EC=sum(EMR.reserve[Sequence,1]),
                           Point.EIR=EMR.mainland.point)
  return(EIR.mainland)
}

EIR.merge.6.2<-EIR.merge.2009(c(3:8,10:13,15:25,27:31),Hepatitis.B.agegroup10,"model6",c(3,11),(2+12)/12,10000)
EIR.merge.7.2<-EIR.merge.2009(c(1:31),Hepatitis.B.agegroup11,"model7",c(3,11),(2+12*2)/12,10000)
EIR.merge.8.2<-EIR.merge.2009(c(1:31),Hepatitis.B.agegroup12,"model8",c(3,11),(2+12*7)/12,10000)
EIR.merge.9<-EIR.merge.2009(c(1:31),Hepatitis.B.agegroup13,"model9",c(2,9),6,10000)
EIR.merge.10<-EIR.merge.2009(c(1:31),Hepatitis.B.agegroup14,"model10",c(2,9),1,10000)
EIR.merge.2009<-Prop.pop.2009[1]*EIR.merge.6.2$EIR+
  Prop.pop.2009[2]*EIR.merge.7.2$EIR+
  Prop.pop.2009[3]*EIR.merge.8.2$EIR+
  Prop.pop.2009[4]*EIR.merge.9$EIR+
  Prop.pop.2009[5]*EIR.merge.10$EIR
EIR.point.2009<-Prop.pop.2009[1]*EIR.merge.6.2$Point.EIR[1]+
  Prop.pop.2009[2]*EIR.merge.7.2$Point.EIR[1]+
  Prop.pop.2009[3]*EIR.merge.8.2$Point.EIR[1]+
  Prop.pop.2009[4]*EIR.merge.9$Point.EIR[1]+
  Prop.pop.2009[5]*EIR.merge.10$Point.EIR[1]
EIR.2009.low<-sort(EIR.merge.2009)[n.cir*0.025]
EIR.2009.high<-sort(EIR.merge.2009)[n.cir*0.975]
EC.merge.2009<-EIR.merge.6.2$EC+
  EIR.merge.7.2$EC+EIR.merge.8.2$EC+
  EIR.merge.9$EC+EIR.merge.10$EC
EC.point.2009<-EIR.merge.6.2$Point.EC[1]+
  EIR.merge.7.2$Point.EC[1]+EIR.merge.8.2$Point.EC[1]+
  EIR.merge.9$Point.EC[1]+EIR.merge.10$Point.EC[1]
EC.2009.low<-sort(EC.merge.2009)[n.cir*0.025]
EC.2009.high<-sort(EC.merge.2009)[n.cir*0.975]

#Saving the results
EIR<-data.frame(EC.point.2002=EC.point.2002,EC.low.2002=EC.2002.low,
                EC.high.2002=EC.2002.high,EIR.point.2002=EIR.point.2002,
                EIR.low.2002=EIR.2002.low,EIR.high.2002=EIR.2002.high,
                EC.point.2009=EC.point.2009,EC.low.2009=EC.2009.low,
                EC.high.2009=EC.2009.high,EIR.point.2009=EIR.point.2009,
                EIR.low.2009=EIR.2009.low,EIR.high.2009=EIR.2009.high)

EIR.merge.2002<-function(Sequence,database,model,Location,intervention,n.cir){
  set.seed(202403)
  EMR.reserve<-matrix(0,nrow=31,ncol=8)
  EMR.reserve.diff<-c()
  MA.reserve<-matrix(0,nrow=n.cir,ncol=31)
  for(i in Sequence)
  {data.number<-database[which(database$Order==i),]
  model.number<-eval(parse(text=paste(model,".",i,sep="")))
  mod.number<-as.matrix(data.frame(model.matrix(model.number),offset=log(data.number$Population.book)))
  mod.number.0<-mod.number
  mod.number.0[,Location]<-0
  line.row<-which(data.number$Intervention.2002==1)
  length.row<-length(line.row)
  freq.number<-data.number$Case[line.row]
  pop.number<-data.number$Population.book[line.row]
  #1对应人口
  coef.number<-c(summary(model.number)$p.coeff,1)
  diff.number<-exp(mod.number%*% coef.number)-exp(mod.number.0%*%coef.number)
  diff.number<-diff.number[which(diff.number!=0)]
  EIR.number<-sum(diff.number)/intervention/(sum(pop.number)/length(pop.number))*100000
  EMR.reserve[i,1]<-sum(diff.number)
  EMR.reserve[i,2]<-intervention
  EMR.reserve[i,3]<-sum(pop.number)/length(pop.number)
  EMR.reserve[i,4]<-EIR.number
  #区间计算
  coef.ini<-coef.number[Location]
  length.coef<-length(coef.ini)
  cov.number<-vcov(model.number)[Location,Location]
  eigen.number<-eigen(cov.number)
  r.norm<-rnorm(length.coef*n.cir)
  r.norm<-matrix(r.norm,n.cir)
  coef.sim<-coef.ini+eigen.number$vectors%*%diag(sqrt(eigen.number$values),length.coef)%*%t(r.norm)
  for(k in 1:n.cir){
    coef.number[Location[1]]<-coef.sim[1,k]
    coef.number[Location[2]]<-coef.sim[2,k]
    diff.number<-exp(mod.number%*%coef.number)-exp(mod.number.0%*%coef.number)
    diff.number<-diff.number[which(diff.number!=0)]
    EMR.reserve.diff[k]<-sum(diff.number)
  }
  MA.reserve[,i]<-EMR.reserve.diff
  EMR.reserve[i,5]<-sort(EMR.reserve.diff)[n.cir*0.025]
  EMR.reserve[i,7]<-sort(EMR.reserve.diff)[n.cir*0.975]
  EMR.reserve[i,6]<-EMR.reserve[i,5]/EMR.reserve[i,2]/EMR.reserve[i,3]*100000
  EMR.reserve[i,8]<-EMR.reserve[i,7]/EMR.reserve[i,2]/EMR.reserve[i,3]*100000
  }
  EMR.reserve[,2]<-intervention
  EMR.mainland.point<-sum(EMR.reserve[Sequence,1])/mean(EMR.reserve[Sequence,2])/sum(EMR.reserve[Sequence,3])*100000
  EMR.mainland<-rowSums(MA.reserve[,Sequence])/mean(EMR.reserve[Sequence,2])/sum(EMR.reserve[Sequence,3])*100000
  EC.mainland<-rowSums(MA.reserve[,Sequence])
  EIR.mainland<-data.frame(EC=EC.mainland,EIR=EMR.mainland,
                           Point.EC=sum(EMR.reserve[Sequence,1]),
                           Point.EIR=EMR.mainland.point)
  return(EIR.mainland)
}
EIR.merge.2009<-function(Sequence,database,model,Location,intervention,n.cir){
  set.seed(202403)
  EMR.reserve<-matrix(0,nrow=31,ncol=8)
  EMR.reserve.diff<-c()
  MA.reserve<-matrix(0,nrow=n.cir,ncol=31)
  for(i in Sequence)
  {data.number<-database[which(database$Order==i),]
  model.number<-eval(parse(text=paste(model,".",i,sep="")))
  mod.number<-as.matrix(data.frame(model.matrix(model.number),offset=log(data.number$Population.book)))
  mod.number.0<-mod.number
  mod.number.0[,Location]<-0
  line.row<-which(data.number$Intervention.2009==1)
  length.row<-length(line.row)
  freq.number<-data.number$Case[line.row]
  pop.number<-data.number$Population.book[line.row]
  #1对应人口
  coef.number<-c(summary(model.number)$p.coeff,1)
  diff.number<-exp(mod.number%*% coef.number)-exp(mod.number.0%*%coef.number)
  diff.number<-diff.number[which(diff.number!=0)]
  EIR.number<-sum(diff.number)/intervention/(sum(pop.number)/length(pop.number))*100000
  EMR.reserve[i,1]<-sum(diff.number)
  EMR.reserve[i,2]<-intervention
  EMR.reserve[i,3]<-sum(pop.number)/length(pop.number)
  EMR.reserve[i,4]<-EIR.number
  #区间计算
  coef.ini<-coef.number[Location]
  length.coef<-length(coef.ini)
  cov.number<-vcov(model.number)[Location,Location]
  eigen.number<-eigen(cov.number)
  r.norm<-rnorm(length.coef*n.cir)
  r.norm<-matrix(r.norm,n.cir)
  coef.sim<-coef.ini+eigen.number$vectors%*%diag(sqrt(eigen.number$values),length.coef)%*%t(r.norm)
  for(k in 1:n.cir){
    coef.number[Location[1]]<-coef.sim[1,k]
    coef.number[Location[2]]<-coef.sim[2,k]
    diff.number<-exp(mod.number%*%coef.number)-exp(mod.number.0%*%coef.number)
    diff.number<-diff.number[which(diff.number!=0)]
    EMR.reserve.diff[k]<-sum(diff.number)
  }
  MA.reserve[,i]<-EMR.reserve.diff
  EMR.reserve[i,5]<-sort(EMR.reserve.diff)[n.cir*0.025]
  EMR.reserve[i,7]<-sort(EMR.reserve.diff)[n.cir*0.975]
  EMR.reserve[i,6]<-EMR.reserve[i,5]/EMR.reserve[i,2]/EMR.reserve[i,3]*100000
  EMR.reserve[i,8]<-EMR.reserve[i,7]/EMR.reserve[i,2]/EMR.reserve[i,3]*100000
  }
  EMR.mainland.point<-sum(EMR.reserve[Sequence,1])/mean(EMR.reserve[Sequence,2])/sum(EMR.reserve[Sequence,3])*100000
  EMR.mainland.point<-sum(EMR.reserve[Sequence,1])/mean(EMR.reserve[Sequence,2])/sum(EMR.reserve[Sequence,3])*100000
  EMR.mainland<-rowSums(MA.reserve[,Sequence])/mean(EMR.reserve[Sequence,2])/sum(EMR.reserve[Sequence,3])*100000
  EC.mainland<-rowSums(MA.reserve[,Sequence])
  EIR.mainland<-data.frame(EC=EC.mainland,EIR=EMR.mainland,
                           Point.EC=sum(EMR.reserve[Sequence,1]),
                           Point.EIR=EMR.mainland.point)
  return(EIR.mainland)
}

#######################Region subgroup EIR, 2002###################################
#Region category
North<-c(1,2,3,4,5)
Northeast<-c(6,7,8)
East<-c(9,10,11,12,13,14,15)
Southcentral<-c(16,17,18,19,20,21)
Southwest<-c(22,23,24,25,26)
Northwest<-c(27,28,29,30,31)
#North
#Calculating the weights
Age.pop2<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup6[which(Hepatitis.B.agegroup6$Order==3|
                                                            Hepatitis.B.agegroup6$Order==4|
                                                            Hepatitis.B.agegroup6$Order==5),],FUN=sum)[,2])
Age.pop3<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup7[which(Hepatitis.B.agegroup7$Order==3|
                                                            Hepatitis.B.agegroup7$Order==4|
                                                            Hepatitis.B.agegroup7$Order==5),],FUN=sum)[,2])
Age.pop4<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup8[which(Hepatitis.B.agegroup8$Order==3|
                                                            Hepatitis.B.agegroup8$Order==4|
                                                            Hepatitis.B.agegroup8$Order==5),],FUN=sum)[,2])
Age.pop5<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup9[which(Hepatitis.B.agegroup9$Order==3|
                                                            Hepatitis.B.agegroup9$Order==4|
                                                            Hepatitis.B.agegroup9$Order==5),],FUN=sum)[,2])
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order==3|
                                                             Hepatitis.B.agegroup10$Order==4|
                                                             Hepatitis.B.agegroup10$Order==5),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order==1|
                                                             Hepatitis.B.agegroup11$Order==2|
                                                             Hepatitis.B.agegroup11$Order==3|
                                                             Hepatitis.B.agegroup11$Order==4|
                                                             Hepatitis.B.agegroup11$Order==5),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order==1|
                                                             Hepatitis.B.agegroup12$Order==2|
                                                             Hepatitis.B.agegroup12$Order==3|
                                                             Hepatitis.B.agegroup12$Order==4|
                                                             Hepatitis.B.agegroup12$Order==5),],FUN=sum)[,2])
Age.pop<-c(Age.pop2,Age.pop3,Age.pop4,Age.pop5,Age.pop6,Age.pop7,Age.pop8)
Prop.pop.north.2002<-c(Age.pop2/sum(Age.pop),Age.pop3/sum(Age.pop),Age.pop4/sum(Age.pop),
                       Age.pop5/sum(Age.pop),Age.pop6/sum(Age.pop),Age.pop7/sum(Age.pop),
                       Age.pop8/sum(Age.pop))
EIR.merge.2<-EIR.merge.2002(c(3:5),Hepatitis.B.agegroup6,"model2",c(2,9),13,10000)
EIR.merge.3<-EIR.merge.2002(c(3:5),Hepatitis.B.agegroup7,"model3",c(2,9),12,10000)
EIR.merge.4<-EIR.merge.2002(c(3:5),Hepatitis.B.agegroup8,"model4",c(2,9),11,10000)
EIR.merge.5.1<-EIR.merge.2002(c(3:5),Hepatitis.B.agegroup9,"model5",c(2,10),10,10000)
EIR.merge.6.1<-EIR.merge.2002(c(3:5),Hepatitis.B.agegroup10,"model6",c(2,10),9,10000)
EIR.merge.7.1<-EIR.merge.2002(c(1:5),Hepatitis.B.agegroup11,"model7",c(2,10),8,10000)
EIR.merge.8.1<-EIR.merge.2002(c(1:5),Hepatitis.B.agegroup12,"model8",c(2,10),3,10000)
EIR.merge.north.2002<-Prop.pop.north.2002[1]*EIR.merge.2$EIR+Prop.pop.north.2002[2]*EIR.merge.3$EIR+
  Prop.pop.north.2002[3]*EIR.merge.4$EIR+Prop.pop.north.2002[4]*EIR.merge.5.1$EIR+
  Prop.pop.north.2002[5]*EIR.merge.6.1$EIR+Prop.pop.north.2002[6]*EIR.merge.7.1$EIR+
  Prop.pop.north.2002[7]*EIR.merge.8.1$EIR
EIR.point.north.2002<-Prop.pop.north.2002[1]*EIR.merge.2$Point.EIR[1]+Prop.pop.north.2002[2]*EIR.merge.3$Point.EIR[1]+
  Prop.pop.north.2002[3]*EIR.merge.4$Point.EIR[1]+Prop.pop.north.2002[4]*EIR.merge.5.1$Point.EIR[1]+
  Prop.pop.north.2002[5]*EIR.merge.6.1$Point.EIR[1]+Prop.pop.north.2002[6]*EIR.merge.7.1$Point.EIR[1]+
  Prop.pop.north.2002[7]*EIR.merge.8.1$Point.EIR[1]
EIR.north.2002.low<-sort(EIR.merge.north.2002)[n.cir*0.025]
EIR.north.2002.high<-sort(EIR.merge.north.2002)[n.cir*0.975]
EC.merge.north.2002<-EIR.merge.2$EC+EIR.merge.3$EC+
  EIR.merge.4$EC+EIR.merge.5.1$EC+
  EIR.merge.6.1$EC+EIR.merge.7.1$EC+
  EIR.merge.8.1$EC
EC.point.north.2002<-EIR.merge.2$Point.EC[1]+EIR.merge.3$Point.EC[1]+
  EIR.merge.4$Point.EC[1]+EIR.merge.5.1$Point.EC[1]+
  EIR.merge.6.1$Point.EC[1]+EIR.merge.7.1$Point.EC[1]+
  EIR.merge.8.1$Point.EC[1]
EC.north.2002.low<-sort(EC.merge.north.2002)[n.cir*0.025]
EC.north.2002.high<-sort(EC.merge.north.2002)[n.cir*0.975]

#Northeast
Age.pop2<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup6[which(Hepatitis.B.agegroup6$Order==7|
                                                            Hepatitis.B.agegroup6$Order==8),],FUN=sum)[,2])
Age.pop3<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup7[which(Hepatitis.B.agegroup7$Order==6|
                                                            Hepatitis.B.agegroup7$Order==7|
                                                            Hepatitis.B.agegroup7$Order==8),],FUN=sum)[,2])
Age.pop4<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup8[which(Hepatitis.B.agegroup8$Order==6|
                                                            Hepatitis.B.agegroup8$Order==7|
                                                            Hepatitis.B.agegroup8$Order==8),],FUN=sum)[,2])
Age.pop5<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup9[which(Hepatitis.B.agegroup9$Order==6|
                                                            Hepatitis.B.agegroup9$Order==7|
                                                            Hepatitis.B.agegroup9$Order==8),],FUN=sum)[,2])
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order==6|
                                                             Hepatitis.B.agegroup10$Order==7|
                                                             Hepatitis.B.agegroup10$Order==8),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order==6|
                                                             Hepatitis.B.agegroup11$Order==7|
                                                             Hepatitis.B.agegroup11$Order==8),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order==6|
                                                             Hepatitis.B.agegroup12$Order==7|
                                                             Hepatitis.B.agegroup12$Order==8),],FUN=sum)[,2])
Age.pop<-c(Age.pop2,Age.pop3,Age.pop4,Age.pop5,Age.pop6,Age.pop7,Age.pop8)
Prop.pop.northeast.2002<-c(Age.pop2/sum(Age.pop),Age.pop3/sum(Age.pop),Age.pop4/sum(Age.pop),
                           Age.pop5/sum(Age.pop),Age.pop6/sum(Age.pop),Age.pop7/sum(Age.pop),
                           Age.pop8/sum(Age.pop))
EIR.merge.2<-EIR.merge.2002(c(7,8),Hepatitis.B.agegroup6,"model2",c(2,9),13,10000)
EIR.merge.3<-EIR.merge.2002(c(6:8),Hepatitis.B.agegroup7,"model3",c(2,9),12,10000)
EIR.merge.4<-EIR.merge.2002(c(6:8),Hepatitis.B.agegroup8,"model4",c(2,9),11,10000)
EIR.merge.5.1<-EIR.merge.2002(c(6:8),Hepatitis.B.agegroup9,"model5",c(2,10),10,10000)
EIR.merge.6.1<-EIR.merge.2002(c(6:8),Hepatitis.B.agegroup10,"model6",c(2,10),9,10000)
EIR.merge.7.1<-EIR.merge.2002(c(6:8),Hepatitis.B.agegroup11,"model7",c(2,10),8,10000)
EIR.merge.8.1<-EIR.merge.2002(c(6:8),Hepatitis.B.agegroup12,"model8",c(2,10),3,10000)
EIR.merge.northeast.2002<-Prop.pop.northeast.2002[1]*EIR.merge.2$EIR+Prop.pop.northeast.2002[2]*EIR.merge.3$EIR+
  Prop.pop.northeast.2002[3]*EIR.merge.4$EIR+Prop.pop.northeast.2002[4]*EIR.merge.5.1$EIR+
  Prop.pop.northeast.2002[5]*EIR.merge.6.1$EIR+Prop.pop.northeast.2002[6]*EIR.merge.7.1$EIR+
  Prop.pop.northeast.2002[7]*EIR.merge.8.1$EIR
EIR.point.northeast.2002<-Prop.pop.northeast.2002[1]*EIR.merge.2$Point.EIR[1]+Prop.pop.northeast.2002[2]*EIR.merge.3$Point.EIR[1]+
  Prop.pop.northeast.2002[3]*EIR.merge.4$Point.EIR[1]+Prop.pop.northeast.2002[4]*EIR.merge.5.1$Point.EIR[1]+
  Prop.pop.northeast.2002[5]*EIR.merge.6.1$Point.EIR[1]+Prop.pop.northeast.2002[6]*EIR.merge.7.1$Point.EIR[1]+
  Prop.pop.northeast.2002[7]*EIR.merge.8.1$Point.EIR[1]
EIR.northeast.2002.low<-sort(EIR.merge.northeast.2002)[n.cir*0.025]
EIR.northeast.2002.high<-sort(EIR.merge.northeast.2002)[n.cir*0.975]
EC.merge.northeast.2002<-EIR.merge.2$EC+EIR.merge.3$EC+
  EIR.merge.4$EC+EIR.merge.5.1$EC+
  EIR.merge.6.1$EC+EIR.merge.7.1$EC+
  EIR.merge.8.1$EC
EC.point.northeast.2002<-EIR.merge.2$Point.EC[1]+EIR.merge.3$Point.EC[1]+
  EIR.merge.4$Point.EC[1]+EIR.merge.5.1$Point.EC[1]+
  EIR.merge.6.1$Point.EC[1]+EIR.merge.7.1$Point.EC[1]+
  EIR.merge.8.1$Point.EC[1]
EC.northeast.2002.low<-sort(EC.merge.northeast.2002)[n.cir*0.025]
EC.northeast.2002.high<-sort(EC.merge.northeast.2002)[n.cir*0.975]

#East
Age.pop2<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup6[which(Hepatitis.B.agegroup6$Order>9&
                                                            Hepatitis.B.agegroup6$Order<16),],FUN=sum)[,2])
Age.pop3<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup7[which(Hepatitis.B.agegroup7$Order>9&
                                                            Hepatitis.B.agegroup7$Order<14|
                                                            Hepatitis.B.agegroup7$Order==15),],FUN=sum)[,2])
Age.pop4<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup8[which(Hepatitis.B.agegroup8$Order>9&
                                                            Hepatitis.B.agegroup8$Order<14|
                                                            Hepatitis.B.agegroup8$Order==15),],FUN=sum)[,2])
Age.pop5<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup9[which(Hepatitis.B.agegroup9$Order>9&
                                                            Hepatitis.B.agegroup9$Order<14|
                                                            Hepatitis.B.agegroup9$Order==15),],FUN=sum)[,2])
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order>9&
                                                             Hepatitis.B.agegroup10$Order<14|
                                                             Hepatitis.B.agegroup10$Order==15),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order>8&
                                                             Hepatitis.B.agegroup11$Order<16),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order>8&
                                                             Hepatitis.B.agegroup12$Order<16),],FUN=sum)[,2])
Age.pop<-c(Age.pop2,Age.pop3,Age.pop4,Age.pop5,Age.pop6,Age.pop7,Age.pop8)
Prop.pop.east.2002<-c(Age.pop2/sum(Age.pop),Age.pop3/sum(Age.pop),Age.pop4/sum(Age.pop),
                      Age.pop5/sum(Age.pop),Age.pop6/sum(Age.pop),Age.pop7/sum(Age.pop),
                      Age.pop8/sum(Age.pop))
EIR.merge.2<-EIR.merge.2002(c(10:15),Hepatitis.B.agegroup6,"model2",c(2,9),13,10000)
EIR.merge.3<-EIR.merge.2002(c(10:13,15),Hepatitis.B.agegroup7,"model3",c(2,9),12,10000)
EIR.merge.4<-EIR.merge.2002(c(10:13,15),Hepatitis.B.agegroup8,"model4",c(2,9),11,10000)
EIR.merge.5.1<-EIR.merge.2002(c(10:13,15),Hepatitis.B.agegroup9,"model5",c(2,10),10,10000)
EIR.merge.6.1<-EIR.merge.2002(c(10:13,15),Hepatitis.B.agegroup10,"model6",c(2,10),9,10000)
EIR.merge.7.1<-EIR.merge.2002(c(9:15),Hepatitis.B.agegroup11,"model7",c(2,10),8,10000)
EIR.merge.8.1<-EIR.merge.2002(c(9:15),Hepatitis.B.agegroup12,"model8",c(2,10),3,10000)
EIR.merge.east.2002<-Prop.pop.east.2002[1]*EIR.merge.2$EIR+Prop.pop.east.2002[2]*EIR.merge.3$EIR+
  Prop.pop.east.2002[3]*EIR.merge.4$EIR+Prop.pop.east.2002[4]*EIR.merge.5.1$EIR+
  Prop.pop.east.2002[5]*EIR.merge.6.1$EIR+Prop.pop.east.2002[6]*EIR.merge.7.1$EIR+
  Prop.pop.east.2002[7]*EIR.merge.8.1$EIR
EIR.point.east.2002<-Prop.pop.east.2002[1]*EIR.merge.2$Point.EIR[1]+Prop.pop.east.2002[2]*EIR.merge.3$Point.EIR[1]+
  Prop.pop.east.2002[3]*EIR.merge.4$Point.EIR[1]+Prop.pop.east.2002[4]*EIR.merge.5.1$Point.EIR[1]+
  Prop.pop.east.2002[5]*EIR.merge.6.1$Point.EIR[1]+Prop.pop.east.2002[6]*EIR.merge.7.1$Point.EIR[1]+
  Prop.pop.east.2002[7]*EIR.merge.8.1$Point.EIR[1]
EIR.east.2002.low<-sort(EIR.merge.east.2002)[n.cir*0.025]
EIR.east.2002.high<-sort(EIR.merge.east.2002)[n.cir*0.975]
EC.merge.east.2002<-EIR.merge.2$EC+EIR.merge.3$EC+
  EIR.merge.4$EC+EIR.merge.5.1$EC+
  EIR.merge.6.1$EC+EIR.merge.7.1$EC+
  EIR.merge.8.1$EC
EC.point.east.2002<-EIR.merge.2$Point.EC[1]+EIR.merge.3$Point.EC[1]+
  EIR.merge.4$Point.EC[1]+EIR.merge.5.1$Point.EC[1]+
  EIR.merge.6.1$Point.EC[1]+EIR.merge.7.1$Point.EC[1]+
  EIR.merge.8.1$Point.EC[1]
EC.east.2002.low<-sort(EC.merge.east.2002)[n.cir*0.025]
EC.east.2002.high<-sort(EC.merge.east.2002)[n.cir*0.975]

#South-centre
Age.pop2<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup6[which(Hepatitis.B.agegroup6$Order>15&
                                                            Hepatitis.B.agegroup6$Order<22),],FUN=sum)[,2])
Age.pop3<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup7[which(Hepatitis.B.agegroup7$Order>15&
                                                            Hepatitis.B.agegroup7$Order<22),],FUN=sum)[,2])
Age.pop4<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup8[which(Hepatitis.B.agegroup8$Order>15&
                                                            Hepatitis.B.agegroup8$Order<22),],FUN=sum)[,2])
Age.pop5<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup9[which(Hepatitis.B.agegroup9$Order>15&
                                                            Hepatitis.B.agegroup9$Order<22),],FUN=sum)[,2])
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order>15&
                                                             Hepatitis.B.agegroup10$Order<22),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order>15&
                                                             Hepatitis.B.agegroup11$Order<22),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order>15&
                                                             Hepatitis.B.agegroup12$Order<22),],FUN=sum)[,2])
Age.pop<-c(Age.pop2,Age.pop3,Age.pop4,Age.pop5,Age.pop6,Age.pop7,Age.pop8)
Prop.pop.southcentral.2002<-c(Age.pop2/sum(Age.pop),Age.pop3/sum(Age.pop),Age.pop4/sum(Age.pop),
                              Age.pop5/sum(Age.pop),Age.pop6/sum(Age.pop),Age.pop7/sum(Age.pop),
                              Age.pop8/sum(Age.pop))
EIR.merge.2<-EIR.merge.2002(c(16:21),Hepatitis.B.agegroup6,"model2",c(2,9),13,10000)
EIR.merge.3<-EIR.merge.2002(c(16:21),Hepatitis.B.agegroup7,"model3",c(2,9),12,10000)
EIR.merge.4<-EIR.merge.2002(c(16:21),Hepatitis.B.agegroup8,"model4",c(2,9),11,10000)
EIR.merge.5.1<-EIR.merge.2002(c(16:21),Hepatitis.B.agegroup9,"model5",c(2,10),10,10000)
EIR.merge.6.1<-EIR.merge.2002(c(16:21),Hepatitis.B.agegroup10,"model6",c(2,10),9,10000)
EIR.merge.7.1<-EIR.merge.2002(c(16:21),Hepatitis.B.agegroup11,"model7",c(2,10),8,10000)
EIR.merge.8.1<-EIR.merge.2002(c(16:21),Hepatitis.B.agegroup12,"model8",c(2,10),3,10000)
EIR.merge.southcentral.2002<-Prop.pop.southcentral.2002[1]*EIR.merge.2$EIR+Prop.pop.southcentral.2002[2]*EIR.merge.3$EIR+
  Prop.pop.southcentral.2002[3]*EIR.merge.4$EIR+Prop.pop.southcentral.2002[4]*EIR.merge.5.1$EIR+
  Prop.pop.southcentral.2002[5]*EIR.merge.6.1$EIR+Prop.pop.southcentral.2002[6]*EIR.merge.7.1$EIR+
  Prop.pop.southcentral.2002[7]*EIR.merge.8.1$EIR
EIR.point.southcentral.2002<-Prop.pop.southcentral.2002[1]*EIR.merge.2$Point.EIR[1]+Prop.pop.southcentral.2002[2]*EIR.merge.3$Point.EIR[1]+
  Prop.pop.southcentral.2002[3]*EIR.merge.4$Point.EIR[1]+Prop.pop.southcentral.2002[4]*EIR.merge.5.1$Point.EIR[1]+
  Prop.pop.southcentral.2002[5]*EIR.merge.6.1$Point.EIR[1]+Prop.pop.southcentral.2002[6]*EIR.merge.7.1$Point.EIR[1]+
  Prop.pop.southcentral.2002[7]*EIR.merge.8.1$Point.EIR[1]
EIR.southcentral.2002.low<-sort(EIR.merge.southcentral.2002)[n.cir*0.025]
EIR.southcentral.2002.high<-sort(EIR.merge.southcentral.2002)[n.cir*0.975]
EC.merge.southcentral.2002<-EIR.merge.2$EC+EIR.merge.3$EC+
  EIR.merge.4$EC+EIR.merge.5.1$EC+
  EIR.merge.6.1$EC+EIR.merge.7.1$EC+
  EIR.merge.8.1$EC
EC.point.southcentral.2002<-EIR.merge.2$Point.EC[1]+EIR.merge.3$Point.EC[1]+
  EIR.merge.4$Point.EC[1]+EIR.merge.5.1$Point.EC[1]+
  EIR.merge.6.1$Point.EC[1]+EIR.merge.7.1$Point.EC[1]+
  EIR.merge.8.1$Point.EC[1]
EC.southcentral.2002.low<-sort(EC.merge.southcentral.2002)[n.cir*0.025]
EC.southcentral.2002.high<-sort(EC.merge.southcentral.2002)[n.cir*0.975]

#Southwest
Age.pop2<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup6[which(Hepatitis.B.agegroup6$Order>21&
                                                            Hepatitis.B.agegroup6$Order<26),],FUN=sum)[,2])
Age.pop3<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup7[which(Hepatitis.B.agegroup7$Order>21&
                                                            Hepatitis.B.agegroup7$Order<26),],FUN=sum)[,2])
Age.pop4<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup8[which(Hepatitis.B.agegroup8$Order>21&
                                                            Hepatitis.B.agegroup8$Order<26),],FUN=sum)[,2])
Age.pop5<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup9[which(Hepatitis.B.agegroup9$Order>21&
                                                            Hepatitis.B.agegroup9$Order<26),],FUN=sum)[,2])
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order>21&
                                                             Hepatitis.B.agegroup10$Order<26),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order>21&
                                                             Hepatitis.B.agegroup11$Order<27),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order>21&
                                                             Hepatitis.B.agegroup12$Order<27),],FUN=sum)[,2])
Age.pop<-c(Age.pop2,Age.pop3,Age.pop4,Age.pop5,Age.pop6,Age.pop7,Age.pop8)
Prop.pop.southwest.2002<-c(Age.pop2/sum(Age.pop),Age.pop3/sum(Age.pop),Age.pop4/sum(Age.pop),
                           Age.pop5/sum(Age.pop),Age.pop6/sum(Age.pop),Age.pop7/sum(Age.pop),
                           Age.pop8/sum(Age.pop))
EIR.merge.2<-EIR.merge.2002(c(22:25),Hepatitis.B.agegroup6,"model2",c(2,9),13,10000)
EIR.merge.3<-EIR.merge.2002(c(22:25),Hepatitis.B.agegroup7,"model3",c(2,9),12,10000)
EIR.merge.4<-EIR.merge.2002(c(22:25),Hepatitis.B.agegroup8,"model4",c(2,9),11,10000)
EIR.merge.5.1<-EIR.merge.2002(c(22:25),Hepatitis.B.agegroup9,"model5",c(2,10),10,10000)
EIR.merge.6.1<-EIR.merge.2002(c(22:25),Hepatitis.B.agegroup10,"model6",c(2,10),9,10000)
EIR.merge.7.1<-EIR.merge.2002(c(22:26),Hepatitis.B.agegroup11,"model7",c(2,10),8,10000)
EIR.merge.8.1<-EIR.merge.2002(c(22:26),Hepatitis.B.agegroup12,"model8",c(2,10),3,10000)
EIR.merge.southwest.2002<-Prop.pop.southwest.2002[1]*EIR.merge.2$EIR+Prop.pop.southwest.2002[2]*EIR.merge.3$EIR+
  Prop.pop.southwest.2002[3]*EIR.merge.4$EIR+Prop.pop.southwest.2002[4]*EIR.merge.5.1$EIR+
  Prop.pop.southwest.2002[5]*EIR.merge.6.1$EIR+Prop.pop.southwest.2002[6]*EIR.merge.7.1$EIR+
  Prop.pop.southwest.2002[7]*EIR.merge.8.1$EIR
EIR.point.southwest.2002<-Prop.pop.southwest.2002[1]*EIR.merge.2$Point.EIR[1]+Prop.pop.southwest.2002[2]*EIR.merge.3$Point.EIR[1]+
  Prop.pop.southwest.2002[3]*EIR.merge.4$Point.EIR[1]+Prop.pop.southwest.2002[4]*EIR.merge.5.1$Point.EIR[1]+
  Prop.pop.southwest.2002[5]*EIR.merge.6.1$Point.EIR[1]+Prop.pop.southwest.2002[6]*EIR.merge.7.1$Point.EIR[1]+
  Prop.pop.southwest.2002[7]*EIR.merge.8.1$Point.EIR[1]
EIR.southwest.2002.low<-sort(EIR.merge.southwest.2002)[n.cir*0.025]
EIR.southwest.2002.high<-sort(EIR.merge.southwest.2002)[n.cir*0.975]
EC.merge.southwest.2002<-EIR.merge.2$EC+EIR.merge.3$EC+
  EIR.merge.4$EC+EIR.merge.5.1$EC+
  EIR.merge.6.1$EC+EIR.merge.7.1$EC+
  EIR.merge.8.1$EC
EC.point.southwest.2002<-EIR.merge.2$Point.EC[1]+EIR.merge.3$Point.EC[1]+
  EIR.merge.4$Point.EC[1]+EIR.merge.5.1$Point.EC[1]+
  EIR.merge.6.1$Point.EC[1]+EIR.merge.7.1$Point.EC[1]+
  EIR.merge.8.1$Point.EC[1]
EC.southwest.2002.low<-sort(EC.merge.southwest.2002)[n.cir*0.025]
EC.southwest.2002.high<-sort(EC.merge.southwest.2002)[n.cir*0.975]

#Northwest
Age.pop2<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup6[which(Hepatitis.B.agegroup6$Order>26&
                                                            Hepatitis.B.agegroup6$Order<32),],FUN=sum)[,2])
Age.pop3<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup7[which(Hepatitis.B.agegroup7$Order>26&
                                                            Hepatitis.B.agegroup7$Order<32),],FUN=sum)[,2])
Age.pop4<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup8[which(Hepatitis.B.agegroup8$Order>26&
                                                            Hepatitis.B.agegroup8$Order<32),],FUN=sum)[,2])
Age.pop5<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup9[which(Hepatitis.B.agegroup9$Order>26&
                                                            Hepatitis.B.agegroup9$Order<32),],FUN=sum)[,2])
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order>26&
                                                             Hepatitis.B.agegroup10$Order<32),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order>26&
                                                             Hepatitis.B.agegroup11$Order<32),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order>26&
                                                             Hepatitis.B.agegroup12$Order<32),],FUN=sum)[,2])
Age.pop<-c(Age.pop2,Age.pop3,Age.pop4,Age.pop5,Age.pop6,Age.pop7,Age.pop8)
Prop.pop.northwest.2002<-c(Age.pop2/sum(Age.pop),Age.pop3/sum(Age.pop),Age.pop4/sum(Age.pop),
                           Age.pop5/sum(Age.pop),Age.pop6/sum(Age.pop),Age.pop7/sum(Age.pop),
                           Age.pop8/sum(Age.pop))
EIR.merge.2<-EIR.merge.2002(c(27:31),Hepatitis.B.agegroup6,"model2",c(2,9),13,10000)
EIR.merge.3<-EIR.merge.2002(c(27:31),Hepatitis.B.agegroup7,"model3",c(2,9),12,10000)
EIR.merge.4<-EIR.merge.2002(c(27:31),Hepatitis.B.agegroup8,"model4",c(2,9),11,10000)
EIR.merge.5.1<-EIR.merge.2002(c(27:31),Hepatitis.B.agegroup9,"model5",c(2,10),10,10000)
EIR.merge.6.1<-EIR.merge.2002(c(27:31),Hepatitis.B.agegroup10,"model6",c(2,10),9,10000)
EIR.merge.7.1<-EIR.merge.2002(c(27:31),Hepatitis.B.agegroup11,"model7",c(2,10),8,10000)
EIR.merge.8.1<-EIR.merge.2002(c(27:31),Hepatitis.B.agegroup12,"model8",c(2,10),3,10000)
EIR.merge.northwest.2002<-Prop.pop.northwest.2002[1]*EIR.merge.2$EIR+Prop.pop.northwest.2002[2]*EIR.merge.3$EIR+
  Prop.pop.northwest.2002[3]*EIR.merge.4$EIR+Prop.pop.northwest.2002[4]*EIR.merge.5.1$EIR+
  Prop.pop.northwest.2002[5]*EIR.merge.6.1$EIR+Prop.pop.northwest.2002[6]*EIR.merge.7.1$EIR+
  Prop.pop.northwest.2002[7]*EIR.merge.8.1$EIR
EIR.point.northwest.2002<-Prop.pop.northwest.2002[1]*EIR.merge.2$Point.EIR[1]+Prop.pop.northwest.2002[2]*EIR.merge.3$Point.EIR[1]+
  Prop.pop.northwest.2002[3]*EIR.merge.4$Point.EIR[1]+Prop.pop.northwest.2002[4]*EIR.merge.5.1$Point.EIR[1]+
  Prop.pop.northwest.2002[5]*EIR.merge.6.1$Point.EIR[1]+Prop.pop.northwest.2002[6]*EIR.merge.7.1$Point.EIR[1]+
  Prop.pop.northwest.2002[7]*EIR.merge.8.1$Point.EIR[1]
EIR.northwest.2002.low<-sort(EIR.merge.northwest.2002)[n.cir*0.025]
EIR.northwest.2002.high<-sort(EIR.merge.northwest.2002)[n.cir*0.975]
EC.merge.northwest.2002<-EIR.merge.2$EC+EIR.merge.3$EC+
  EIR.merge.4$EC+EIR.merge.5.1$EC+
  EIR.merge.6.1$EC+EIR.merge.7.1$EC+
  EIR.merge.8.1$EC
EC.point.northwest.2002<-EIR.merge.2$Point.EC[1]+EIR.merge.3$Point.EC[1]+
  EIR.merge.4$Point.EC[1]+EIR.merge.5.1$Point.EC[1]+
  EIR.merge.6.1$Point.EC[1]+EIR.merge.7.1$Point.EC[1]+
  EIR.merge.8.1$Point.EC[1]
EC.northwest.2002.low<-sort(EC.merge.northwest.2002)[n.cir*0.025]
EC.northwest.2002.high<-sort(EC.merge.northwest.2002)[n.cir*0.975]

#######################HBsAg subgroup EIR, 2002###################################
#HBsAg category
Low<-c(1,2,3,4,5,6,7,9,10,15,16,25,27,28,31)
Middle<-c(8,11,12,17,18,22,23,24,29,30)
High<-c(13,14,19,20,21,26)
#Low HBsAg prevalence
Age.pop2<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup6[which(Hepatitis.B.agegroup6$Order>2&
                                                            Hepatitis.B.agegroup6$Order<6|
                                                            Hepatitis.B.agegroup6$Order==7|
                                                            Hepatitis.B.agegroup6$Order==10|
                                                            Hepatitis.B.agegroup6$Order==15|
                                                            Hepatitis.B.agegroup6$Order==16|
                                                            Hepatitis.B.agegroup6$Order==25|
                                                            Hepatitis.B.agegroup6$Order==27|
                                                            Hepatitis.B.agegroup6$Order==28|
                                                            Hepatitis.B.agegroup6$Order==31),],FUN=sum)[,2])
Age.pop3<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup7[which(Hepatitis.B.agegroup7$Order>2&
                                                            Hepatitis.B.agegroup7$Order<8|
                                                            Hepatitis.B.agegroup7$Order==10|
                                                            Hepatitis.B.agegroup7$Order==15|
                                                            Hepatitis.B.agegroup7$Order==16|
                                                            Hepatitis.B.agegroup7$Order==25|
                                                            Hepatitis.B.agegroup7$Order==27|
                                                            Hepatitis.B.agegroup7$Order==28|
                                                            Hepatitis.B.agegroup7$Order==31),],FUN=sum)[,2])
Age.pop4<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup8[which(Hepatitis.B.agegroup8$Order>2&
                                                            Hepatitis.B.agegroup8$Order<8|
                                                            Hepatitis.B.agegroup8$Order==10|
                                                            Hepatitis.B.agegroup8$Order==15|
                                                            Hepatitis.B.agegroup8$Order==16|
                                                            Hepatitis.B.agegroup8$Order==25|
                                                            Hepatitis.B.agegroup8$Order==27|
                                                            Hepatitis.B.agegroup8$Order==28|
                                                            Hepatitis.B.agegroup8$Order==31),],FUN=sum)[,2])
Age.pop5<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup9[which(Hepatitis.B.agegroup9$Order>2&
                                                            Hepatitis.B.agegroup9$Order<8|
                                                            Hepatitis.B.agegroup9$Order==10|
                                                            Hepatitis.B.agegroup9$Order==15|
                                                            Hepatitis.B.agegroup9$Order==16|
                                                            Hepatitis.B.agegroup9$Order==25|
                                                            Hepatitis.B.agegroup9$Order==27|
                                                            Hepatitis.B.agegroup9$Order==28|
                                                            Hepatitis.B.agegroup9$Order==31),],FUN=sum)[,2])
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order>2&
                                                             Hepatitis.B.agegroup10$Order<8|
                                                             Hepatitis.B.agegroup10$Order==10|
                                                             Hepatitis.B.agegroup10$Order==15|
                                                             Hepatitis.B.agegroup10$Order==16|
                                                             Hepatitis.B.agegroup10$Order==25|
                                                             Hepatitis.B.agegroup10$Order==27|
                                                             Hepatitis.B.agegroup10$Order==28|
                                                             Hepatitis.B.agegroup10$Order==31),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order>0&
                                                             Hepatitis.B.agegroup11$Order<8|
                                                             Hepatitis.B.agegroup11$Order==9|
                                                             Hepatitis.B.agegroup11$Order==10|
                                                             Hepatitis.B.agegroup11$Order==15|
                                                             Hepatitis.B.agegroup11$Order==16|
                                                             Hepatitis.B.agegroup11$Order==25|
                                                             Hepatitis.B.agegroup11$Order==27|
                                                             Hepatitis.B.agegroup11$Order==28|
                                                             Hepatitis.B.agegroup11$Order==31),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order>0&
                                                             Hepatitis.B.agegroup12$Order<8|
                                                             Hepatitis.B.agegroup12$Order==9|
                                                             Hepatitis.B.agegroup12$Order==10|
                                                             Hepatitis.B.agegroup12$Order==15|
                                                             Hepatitis.B.agegroup12$Order==16|
                                                             Hepatitis.B.agegroup12$Order==25|
                                                             Hepatitis.B.agegroup12$Order==27|
                                                             Hepatitis.B.agegroup12$Order==28|
                                                             Hepatitis.B.agegroup12$Order==31),],FUN=sum)[,2])
Age.pop<-c(Age.pop2,Age.pop3,Age.pop4,Age.pop5,Age.pop6,Age.pop7,Age.pop8)
Prop.pop.lowpr.2002<-c(Age.pop2/sum(Age.pop),Age.pop3/sum(Age.pop),Age.pop4/sum(Age.pop),
                       Age.pop5/sum(Age.pop),Age.pop6/sum(Age.pop),Age.pop7/sum(Age.pop),
                       Age.pop8/sum(Age.pop))
EIR.merge.2<-EIR.merge.2002(c(3,4,5,7,10,15,16,25,27,28,31),Hepatitis.B.agegroup6,"model2",c(2,9),13,10000)
EIR.merge.3<-EIR.merge.2002(c(3,4,5,6,7,10,15,16,25,27,28,31),Hepatitis.B.agegroup7,"model3",c(2,9),12,10000)
EIR.merge.4<-EIR.merge.2002(c(3,4,5,6,7,10,15,16,25,27,28,31),Hepatitis.B.agegroup8,"model4",c(2,9),11,10000)
EIR.merge.5.1<-EIR.merge.2002(c(3,4,5,6,7,10,15,16,25,27,28,31),Hepatitis.B.agegroup9,"model5",c(2,10),10,10000)
EIR.merge.6.1<-EIR.merge.2002(c(3,4,5,6,7,10,15,16,25,27,28,31),Hepatitis.B.agegroup10,"model6",c(2,10),9,10000)
EIR.merge.7.1<-EIR.merge.2002(c(1,2,3,4,5,6,7,9,10,15,16,25,27,28,31),Hepatitis.B.agegroup11,"model7",c(2,10),8,10000)
EIR.merge.8.1<-EIR.merge.2002(c(1,2,3,4,5,6,7,9,10,15,16,25,27,28,31),Hepatitis.B.agegroup12,"model8",c(2,10),3,10000)
EIR.merge.lowpr.2002<-Prop.pop.lowpr.2002[1]*EIR.merge.2$EIR+Prop.pop.lowpr.2002[2]*EIR.merge.3$EIR+
  Prop.pop.lowpr.2002[3]*EIR.merge.4$EIR+Prop.pop.lowpr.2002[4]*EIR.merge.5.1$EIR+
  Prop.pop.lowpr.2002[5]*EIR.merge.6.1$EIR+Prop.pop.lowpr.2002[6]*EIR.merge.7.1$EIR+
  Prop.pop.lowpr.2002[7]*EIR.merge.8.1$EIR
EIR.point.lowpr.2002<-Prop.pop.lowpr.2002[1]*EIR.merge.2$Point.EIR[1]+Prop.pop.lowpr.2002[2]*EIR.merge.3$Point.EIR[1]+
  Prop.pop.lowpr.2002[3]*EIR.merge.4$Point.EIR[1]+Prop.pop.lowpr.2002[4]*EIR.merge.5.1$Point.EIR[1]+
  Prop.pop.lowpr.2002[5]*EIR.merge.6.1$Point.EIR[1]+Prop.pop.lowpr.2002[6]*EIR.merge.7.1$Point.EIR[1]+
  Prop.pop.lowpr.2002[7]*EIR.merge.8.1$Point.EIR[1]
EIR.lowpr.2002.low<-sort(EIR.merge.lowpr.2002)[n.cir*0.025]
EIR.lowpr.2002.high<-sort(EIR.merge.lowpr.2002)[n.cir*0.975]
EC.merge.lowpr.2002<-EIR.merge.2$EC+EIR.merge.3$EC+
  EIR.merge.4$EC+EIR.merge.5.1$EC+
  EIR.merge.6.1$EC+EIR.merge.7.1$EC+
  EIR.merge.8.1$EC
EC.point.lowpr.2002<-EIR.merge.2$Point.EC[1]+EIR.merge.3$Point.EC[1]+
  EIR.merge.4$Point.EC[1]+EIR.merge.5.1$Point.EC[1]+
  EIR.merge.6.1$Point.EC[1]+EIR.merge.7.1$Point.EC[1]+
  EIR.merge.8.1$Point.EC[1]
EC.lowpr.2002.low<-sort(EC.merge.lowpr.2002)[n.cir*0.025]
EC.lowpr.2002.high<-sort(EC.merge.lowpr.2002)[n.cir*0.975]

#Middle HBsAg prevalence
Age.pop2<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup6[which(Hepatitis.B.agegroup6$Order==8|
                                                            Hepatitis.B.agegroup6$Order==11|
                                                            Hepatitis.B.agegroup6$Order==12|
                                                            Hepatitis.B.agegroup6$Order==17|
                                                            Hepatitis.B.agegroup6$Order==18|
                                                            Hepatitis.B.agegroup6$Order==22|
                                                            Hepatitis.B.agegroup6$Order==23|
                                                            Hepatitis.B.agegroup6$Order==24|
                                                            Hepatitis.B.agegroup6$Order==29|
                                                            Hepatitis.B.agegroup6$Order==30),],FUN=sum)[,2])
Age.pop3<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup7[which(Hepatitis.B.agegroup7$Order==8|
                                                            Hepatitis.B.agegroup7$Order==11|
                                                            Hepatitis.B.agegroup7$Order==12|
                                                            Hepatitis.B.agegroup7$Order==17|
                                                            Hepatitis.B.agegroup7$Order==18|
                                                            Hepatitis.B.agegroup7$Order==22|
                                                            Hepatitis.B.agegroup7$Order==23|
                                                            Hepatitis.B.agegroup7$Order==24|
                                                            Hepatitis.B.agegroup7$Order==29|
                                                            Hepatitis.B.agegroup7$Order==30),],FUN=sum)[,2])
Age.pop4<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup8[which(Hepatitis.B.agegroup8$Order==8|
                                                            Hepatitis.B.agegroup8$Order==11|
                                                            Hepatitis.B.agegroup8$Order==12|
                                                            Hepatitis.B.agegroup8$Order==17|
                                                            Hepatitis.B.agegroup8$Order==18|
                                                            Hepatitis.B.agegroup8$Order==22|
                                                            Hepatitis.B.agegroup8$Order==23|
                                                            Hepatitis.B.agegroup8$Order==24|
                                                            Hepatitis.B.agegroup8$Order==29|
                                                            Hepatitis.B.agegroup8$Order==30),],FUN=sum)[,2])
Age.pop5<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup9[which(Hepatitis.B.agegroup9$Order==8|
                                                            Hepatitis.B.agegroup9$Order==11|
                                                            Hepatitis.B.agegroup9$Order==12|
                                                            Hepatitis.B.agegroup9$Order==17|
                                                            Hepatitis.B.agegroup9$Order==18|
                                                            Hepatitis.B.agegroup9$Order==22|
                                                            Hepatitis.B.agegroup9$Order==23|
                                                            Hepatitis.B.agegroup9$Order==24|
                                                            Hepatitis.B.agegroup9$Order==29|
                                                            Hepatitis.B.agegroup9$Order==30),],FUN=sum)[,2])
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order==8|
                                                             Hepatitis.B.agegroup10$Order==11|
                                                             Hepatitis.B.agegroup10$Order==12|
                                                             Hepatitis.B.agegroup10$Order==17|
                                                             Hepatitis.B.agegroup10$Order==18|
                                                             Hepatitis.B.agegroup10$Order==22|
                                                             Hepatitis.B.agegroup10$Order==23|
                                                             Hepatitis.B.agegroup10$Order==24|
                                                             Hepatitis.B.agegroup10$Order==29|
                                                             Hepatitis.B.agegroup10$Order==30),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order==8|
                                                             Hepatitis.B.agegroup11$Order==11|
                                                             Hepatitis.B.agegroup11$Order==12|
                                                             Hepatitis.B.agegroup11$Order==17|
                                                             Hepatitis.B.agegroup11$Order==18|
                                                             Hepatitis.B.agegroup11$Order==22|
                                                             Hepatitis.B.agegroup11$Order==23|
                                                             Hepatitis.B.agegroup11$Order==24|
                                                             Hepatitis.B.agegroup11$Order==29|
                                                             Hepatitis.B.agegroup11$Order==30),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup11$Order==8|
                                                             Hepatitis.B.agegroup11$Order==11|
                                                             Hepatitis.B.agegroup11$Order==12|
                                                             Hepatitis.B.agegroup11$Order==17|
                                                             Hepatitis.B.agegroup11$Order==18|
                                                             Hepatitis.B.agegroup11$Order==22|
                                                             Hepatitis.B.agegroup11$Order==23|
                                                             Hepatitis.B.agegroup11$Order==24|
                                                             Hepatitis.B.agegroup11$Order==29|
                                                             Hepatitis.B.agegroup11$Order==30),],FUN=sum)[,2])
Age.pop<-c(Age.pop2,Age.pop3,Age.pop4,Age.pop5,Age.pop6,Age.pop7,Age.pop8)
Prop.pop.middlepr.2002<-c(Age.pop2/sum(Age.pop),Age.pop3/sum(Age.pop),Age.pop4/sum(Age.pop),
                          Age.pop5/sum(Age.pop),Age.pop6/sum(Age.pop),Age.pop7/sum(Age.pop),
                          Age.pop8/sum(Age.pop))
EIR.merge.2<-EIR.merge.2002(c(8,11,12,17,18,22,23,24,29,30),Hepatitis.B.agegroup6,"model2",c(2,9),13,10000)
EIR.merge.3<-EIR.merge.2002(c(8,11,12,17,18,22,23,24,29,30),Hepatitis.B.agegroup7,"model3",c(2,9),12,10000)
EIR.merge.4<-EIR.merge.2002(c(8,11,12,17,18,22,23,24,29,30),Hepatitis.B.agegroup8,"model4",c(2,9),11,10000)
EIR.merge.5.1<-EIR.merge.2002(c(8,11,12,17,18,22,23,24,29,30),Hepatitis.B.agegroup9,"model5",c(2,10),10,10000)
EIR.merge.6.1<-EIR.merge.2002(c(8,11,12,17,18,22,23,24,29,30),Hepatitis.B.agegroup10,"model6",c(2,10),9,10000)
EIR.merge.7.1<-EIR.merge.2002(c(8,11,12,17,18,22,23,24,29,30),Hepatitis.B.agegroup11,"model7",c(2,10),8,10000)
EIR.merge.8.1<-EIR.merge.2002(c(8,11,12,17,18,22,23,24,29,30),Hepatitis.B.agegroup12,"model8",c(2,10),3,10000)
EIR.merge.middlepr.2002<-Prop.pop.middlepr.2002[1]*EIR.merge.2$EIR+Prop.pop.middlepr.2002[2]*EIR.merge.3$EIR+
  Prop.pop.middlepr.2002[3]*EIR.merge.4$EIR+Prop.pop.middlepr.2002[4]*EIR.merge.5.1$EIR+
  Prop.pop.middlepr.2002[5]*EIR.merge.6.1$EIR+Prop.pop.middlepr.2002[6]*EIR.merge.7.1$EIR+
  Prop.pop.middlepr.2002[7]*EIR.merge.8.1$EIR
EIR.point.middlepr.2002<-Prop.pop.middlepr.2002[1]*EIR.merge.2$Point.EIR[1]+Prop.pop.middlepr.2002[2]*EIR.merge.3$Point.EIR[1]+
  Prop.pop.middlepr.2002[3]*EIR.merge.4$Point.EIR[1]+Prop.pop.middlepr.2002[4]*EIR.merge.5.1$Point.EIR[1]+
  Prop.pop.middlepr.2002[5]*EIR.merge.6.1$Point.EIR[1]+Prop.pop.middlepr.2002[6]*EIR.merge.7.1$Point.EIR[1]+
  Prop.pop.middlepr.2002[7]*EIR.merge.8.1$Point.EIR[1]
EIR.middlepr.2002.low<-sort(EIR.merge.middlepr.2002)[n.cir*0.025]
EIR.middlepr.2002.high<-sort(EIR.merge.middlepr.2002)[n.cir*0.975]
EC.merge.middlepr.2002<-EIR.merge.2$EC+EIR.merge.3$EC+
  EIR.merge.4$EC+EIR.merge.5.1$EC+
  EIR.merge.6.1$EC+EIR.merge.7.1$EC+
  EIR.merge.8.1$EC
EC.point.middlepr.2002<-EIR.merge.2$Point.EC[1]+EIR.merge.3$Point.EC[1]+
  EIR.merge.4$Point.EC[1]+EIR.merge.5.1$Point.EC[1]+
  EIR.merge.6.1$Point.EC[1]+EIR.merge.7.1$Point.EC[1]+
  EIR.merge.8.1$Point.EC[1]
EC.middlepr.2002.low<-sort(EC.merge.middlepr.2002)[n.cir*0.025]
EC.middlepr.2002.high<-sort(EC.merge.middlepr.2002)[n.cir*0.975]

#High HBsAg prevalence
Age.pop2<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup6[which(Hepatitis.B.agegroup6$Order==13|
                                                            Hepatitis.B.agegroup6$Order==14|
                                                            Hepatitis.B.agegroup6$Order==19|
                                                            Hepatitis.B.agegroup6$Order==20|
                                                            Hepatitis.B.agegroup6$Order==21),],FUN=sum)[,2])
Age.pop3<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup7[which(Hepatitis.B.agegroup7$Order==13|
                                                            Hepatitis.B.agegroup7$Order==19|
                                                            Hepatitis.B.agegroup7$Order==20|
                                                            Hepatitis.B.agegroup7$Order==21),],FUN=sum)[,2])
Age.pop4<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup8[which(Hepatitis.B.agegroup8$Order==13|
                                                            Hepatitis.B.agegroup8$Order==19|
                                                            Hepatitis.B.agegroup8$Order==20|
                                                            Hepatitis.B.agegroup8$Order==21),],FUN=sum)[,2])
Age.pop5<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup9[which(Hepatitis.B.agegroup9$Order==13|
                                                            Hepatitis.B.agegroup9$Order==19|
                                                            Hepatitis.B.agegroup9$Order==20|
                                                            Hepatitis.B.agegroup9$Order==21),],FUN=sum)[,2])
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order==13|
                                                             Hepatitis.B.agegroup10$Order==19|
                                                             Hepatitis.B.agegroup10$Order==20|
                                                             Hepatitis.B.agegroup10$Order==21),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order==13|
                                                             Hepatitis.B.agegroup11$Order==14|
                                                             Hepatitis.B.agegroup11$Order==19|
                                                             Hepatitis.B.agegroup11$Order==20|
                                                             Hepatitis.B.agegroup11$Order==21|
                                                             Hepatitis.B.agegroup11$Order==26),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order==13|
                                                             Hepatitis.B.agegroup12$Order==14|
                                                             Hepatitis.B.agegroup12$Order==19|
                                                             Hepatitis.B.agegroup12$Order==20|
                                                             Hepatitis.B.agegroup12$Order==21|
                                                             Hepatitis.B.agegroup12$Order==26),],FUN=sum)[,2])
Age.pop<-c(Age.pop2,Age.pop3,Age.pop4,Age.pop5,Age.pop6,Age.pop7,Age.pop8)
Prop.pop.highpr.2002<-c(Age.pop2/sum(Age.pop),Age.pop3/sum(Age.pop),Age.pop4/sum(Age.pop),
                        Age.pop5/sum(Age.pop),Age.pop6/sum(Age.pop),Age.pop7/sum(Age.pop),
                        Age.pop8/sum(Age.pop))
EIR.merge.2<-EIR.merge.2002(c(13,14,19,20,21),Hepatitis.B.agegroup6,"model2",c(2,9),13,10000)
EIR.merge.3<-EIR.merge.2002(c(13,19,20,21),Hepatitis.B.agegroup7,"model3",c(2,9),12,10000)
EIR.merge.4<-EIR.merge.2002(c(13,19,20,21),Hepatitis.B.agegroup8,"model4",c(2,9),11,10000)
EIR.merge.5.1<-EIR.merge.2002(c(13,19,20,21),Hepatitis.B.agegroup9,"model5",c(2,10),10,10000)
EIR.merge.6.1<-EIR.merge.2002(c(13,19,20,21),Hepatitis.B.agegroup10,"model6",c(2,10),9,10000)
EIR.merge.7.1<-EIR.merge.2002(c(13,14,19,20,21,26),Hepatitis.B.agegroup11,"model7",c(2,10),8,10000)
EIR.merge.8.1<-EIR.merge.2002(c(13,14,19,20,21,26),Hepatitis.B.agegroup12,"model8",c(2,10),3,10000)
EIR.merge.highpr.2002<-Prop.pop.highpr.2002[1]*EIR.merge.2$EIR+Prop.pop.highpr.2002[2]*EIR.merge.3$EIR+
  Prop.pop.highpr.2002[3]*EIR.merge.4$EIR+Prop.pop.highpr.2002[4]*EIR.merge.5.1$EIR+
  Prop.pop.highpr.2002[5]*EIR.merge.6.1$EIR+Prop.pop.highpr.2002[6]*EIR.merge.7.1$EIR+
  Prop.pop.highpr.2002[7]*EIR.merge.8.1$EIR
EIR.point.highpr.2002<-Prop.pop.highpr.2002[1]*EIR.merge.2$Point.EIR[1]+Prop.pop.highpr.2002[2]*EIR.merge.3$Point.EIR[1]+
  Prop.pop.highpr.2002[3]*EIR.merge.4$Point.EIR[1]+Prop.pop.highpr.2002[4]*EIR.merge.5.1$Point.EIR[1]+
  Prop.pop.highpr.2002[5]*EIR.merge.6.1$Point.EIR[1]+Prop.pop.highpr.2002[6]*EIR.merge.7.1$Point.EIR[1]+
  Prop.pop.highpr.2002[7]*EIR.merge.8.1$Point.EIR[1]
EIR.highpr.2002.low<-sort(EIR.merge.highpr.2002)[n.cir*0.025]
EIR.highpr.2002.high<-sort(EIR.merge.highpr.2002)[n.cir*0.975]
EC.merge.highpr.2002<-EIR.merge.2$EC+EIR.merge.3$EC+
  EIR.merge.4$EC+EIR.merge.5.1$EC+
  EIR.merge.6.1$EC+EIR.merge.7.1$EC+
  EIR.merge.8.1$EC
EC.point.highpr.2002<-EIR.merge.2$Point.EC[1]+EIR.merge.3$Point.EC[1]+
  EIR.merge.4$Point.EC[1]+EIR.merge.5.1$Point.EC[1]+
  EIR.merge.6.1$Point.EC[1]+EIR.merge.7.1$Point.EC[1]+
  EIR.merge.8.1$Point.EC[1]
EC.highpr.2002.low<-sort(EC.merge.highpr.2002)[n.cir*0.025]
EC.highpr.2002.high<-sort(EC.merge.highpr.2002)[n.cir*0.975]

#######################Urbanisation rate subgroup EIR, 2002###################################
#Urbanisation rate category
LowUR<-c(3,12,14,16,18,20,23:31)
HighUR<-c(1,2,4,5,6,7,8,9,10,11,13,15,17,19,21,22)
#Low urbanisation rate
Age.pop2<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup6[which(Hepatitis.B.agegroup6$Order==3|
                                                            Hepatitis.B.agegroup6$Order==12|
                                                            Hepatitis.B.agegroup6$Order==14|
                                                            Hepatitis.B.agegroup6$Order==16|
                                                            Hepatitis.B.agegroup6$Order==18|
                                                            Hepatitis.B.agegroup6$Order==20|
                                                            Hepatitis.B.agegroup6$Order>22&
                                                            Hepatitis.B.agegroup6$Order<26|
                                                            Hepatitis.B.agegroup6$Order>26),],FUN=sum)[,2])
Age.pop3<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup7[which(Hepatitis.B.agegroup7$Order==3|
                                                            Hepatitis.B.agegroup7$Order==12|
                                                            Hepatitis.B.agegroup7$Order==16|
                                                            Hepatitis.B.agegroup7$Order==18|
                                                            Hepatitis.B.agegroup7$Order==20|
                                                            Hepatitis.B.agegroup7$Order>22&
                                                            Hepatitis.B.agegroup7$Order<26|
                                                            Hepatitis.B.agegroup7$Order>26),],FUN=sum)[,2])
Age.pop4<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup8[which(Hepatitis.B.agegroup8$Order==3|
                                                            Hepatitis.B.agegroup8$Order==12|
                                                            Hepatitis.B.agegroup8$Order==16|
                                                            Hepatitis.B.agegroup8$Order==18|
                                                            Hepatitis.B.agegroup8$Order==20|
                                                            Hepatitis.B.agegroup8$Order>22&
                                                            Hepatitis.B.agegroup8$Order<26|
                                                            Hepatitis.B.agegroup8$Order>26),],FUN=sum)[,2])
Age.pop5<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup9[which(Hepatitis.B.agegroup9$Order==3|
                                                            Hepatitis.B.agegroup9$Order==12|
                                                            Hepatitis.B.agegroup9$Order==16|
                                                            Hepatitis.B.agegroup9$Order==18|
                                                            Hepatitis.B.agegroup9$Order==20|
                                                            Hepatitis.B.agegroup9$Order>22&
                                                            Hepatitis.B.agegroup9$Order<26|
                                                            Hepatitis.B.agegroup9$Order>26),],FUN=sum)[,2])
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order==3|
                                                             Hepatitis.B.agegroup10$Order==12|
                                                             Hepatitis.B.agegroup10$Order==16|
                                                             Hepatitis.B.agegroup10$Order==18|
                                                             Hepatitis.B.agegroup10$Order==20|
                                                             Hepatitis.B.agegroup10$Order>22&
                                                             Hepatitis.B.agegroup10$Order<26|
                                                             Hepatitis.B.agegroup10$Order>26),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup6$Order==3|
                                                             Hepatitis.B.agegroup6$Order==12|
                                                             Hepatitis.B.agegroup6$Order==14|
                                                             Hepatitis.B.agegroup6$Order==16|
                                                             Hepatitis.B.agegroup6$Order==18|
                                                             Hepatitis.B.agegroup6$Order==20|
                                                             Hepatitis.B.agegroup6$Order>22&
                                                             Hepatitis.B.agegroup6$Order<32),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup6$Order==3|
                                                             Hepatitis.B.agegroup6$Order==12|
                                                             Hepatitis.B.agegroup6$Order==14|
                                                             Hepatitis.B.agegroup6$Order==16|
                                                             Hepatitis.B.agegroup6$Order==18|
                                                             Hepatitis.B.agegroup6$Order==20|
                                                             Hepatitis.B.agegroup6$Order>22&
                                                             Hepatitis.B.agegroup6$Order<32),],FUN=sum)[,2])
Age.pop<-c(Age.pop2,Age.pop3,Age.pop4,Age.pop5,Age.pop6,Age.pop7,Age.pop8)
Prop.pop.lowur.2002<-c(Age.pop2/sum(Age.pop),Age.pop3/sum(Age.pop),Age.pop4/sum(Age.pop),
                       Age.pop5/sum(Age.pop),Age.pop6/sum(Age.pop),Age.pop7/sum(Age.pop),
                       Age.pop8/sum(Age.pop))
EIR.merge.2<-EIR.merge.2002(c(3,12,14,16,18,20,23:25,27:31),Hepatitis.B.agegroup6,"model2",c(2,9),13,10000)
EIR.merge.3<-EIR.merge.2002(c(3,12,16,18,20,23:25,27:31),Hepatitis.B.agegroup7,"model3",c(2,9),12,10000)
EIR.merge.4<-EIR.merge.2002(c(3,12,16,18,20,23:25,27:31),Hepatitis.B.agegroup8,"model4",c(2,9),11,10000)
EIR.merge.5.1<-EIR.merge.2002(c(3,12,16,18,20,23:25,27:31),Hepatitis.B.agegroup9,"model5",c(2,10),10,10000)
EIR.merge.6.1<-EIR.merge.2002(c(3,12,16,18,20,23:25,27:31),Hepatitis.B.agegroup10,"model6",c(2,10),9,10000)
EIR.merge.7.1<-EIR.merge.2002(c(3,12,14,16,18,20,23:31),Hepatitis.B.agegroup11,"model7",c(2,10),8,10000)
EIR.merge.8.1<-EIR.merge.2002(c(3,12,14,16,18,20,23:31),Hepatitis.B.agegroup12,"model8",c(2,10),3,10000)
EIR.merge.lowur.2002<-Prop.pop.lowur.2002[1]*EIR.merge.2$EIR+Prop.pop.lowur.2002[2]*EIR.merge.3$EIR+
  Prop.pop.lowur.2002[3]*EIR.merge.4$EIR+Prop.pop.lowur.2002[4]*EIR.merge.5.1$EIR+
  Prop.pop.lowur.2002[5]*EIR.merge.6.1$EIR+Prop.pop.lowur.2002[6]*EIR.merge.7.1$EIR+
  Prop.pop.lowur.2002[7]*EIR.merge.8.1$EIR
EIR.point.lowur.2002<-Prop.pop.lowur.2002[1]*EIR.merge.2$Point.EIR[1]+Prop.pop.lowur.2002[2]*EIR.merge.3$Point.EIR[1]+
  Prop.pop.lowur.2002[3]*EIR.merge.4$Point.EIR[1]+Prop.pop.lowur.2002[4]*EIR.merge.5.1$Point.EIR[1]+
  Prop.pop.lowur.2002[5]*EIR.merge.6.1$Point.EIR[1]+Prop.pop.lowur.2002[6]*EIR.merge.7.1$Point.EIR[1]+
  Prop.pop.lowur.2002[7]*EIR.merge.8.1$Point.EIR[1]
EIR.lowur.2002.low<-sort(EIR.merge.lowur.2002)[n.cir*0.025]
EIR.lowur.2002.high<-sort(EIR.merge.lowur.2002)[n.cir*0.975]
EC.merge.lowur.2002<-EIR.merge.2$EC+EIR.merge.3$EC+
  EIR.merge.4$EC+EIR.merge.5.1$EC+
  EIR.merge.6.1$EC+EIR.merge.7.1$EC+
  EIR.merge.8.1$EC
EC.point.lowur.2002<-EIR.merge.2$Point.EC[1]+EIR.merge.3$Point.EC[1]+
  EIR.merge.4$Point.EC[1]+EIR.merge.5.1$Point.EC[1]+
  EIR.merge.6.1$Point.EC[1]+EIR.merge.7.1$Point.EC[1]+
  EIR.merge.8.1$Point.EC[1]
EC.lowur.2002.low<-sort(EC.merge.lowur.2002)[n.cir*0.025]
EC.lowur.2002.high<-sort(EC.merge.lowur.2002)[n.cir*0.975]

#High urbanisation rate
Age.pop2<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup6[which(Hepatitis.B.agegroup6$Order==4|
                                                            Hepatitis.B.agegroup6$Order==5|
                                                            Hepatitis.B.agegroup6$Order==7|
                                                            Hepatitis.B.agegroup6$Order==8|
                                                            Hepatitis.B.agegroup6$Order==10|
                                                            Hepatitis.B.agegroup6$Order==11|
                                                            Hepatitis.B.agegroup6$Order==13|
                                                            Hepatitis.B.agegroup6$Order==15|
                                                            Hepatitis.B.agegroup6$Order==17|
                                                            Hepatitis.B.agegroup6$Order==19|
                                                            Hepatitis.B.agegroup6$Order==21|
                                                            Hepatitis.B.agegroup6$Order==22),],FUN=sum)[,2])
Age.pop3<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup7[which(Hepatitis.B.agegroup7$Order>3&
                                                            Hepatitis.B.agegroup7$Order<9|
                                                            Hepatitis.B.agegroup7$Order==10|
                                                            Hepatitis.B.agegroup7$Order==11|
                                                            Hepatitis.B.agegroup7$Order==13|
                                                            Hepatitis.B.agegroup7$Order==15|
                                                            Hepatitis.B.agegroup7$Order==17|
                                                            Hepatitis.B.agegroup7$Order==19|
                                                            Hepatitis.B.agegroup7$Order==21|
                                                            Hepatitis.B.agegroup7$Order==22),],FUN=sum)[,2])
Age.pop4<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup8[which(Hepatitis.B.agegroup8$Order>3&
                                                            Hepatitis.B.agegroup8$Order<9|
                                                            Hepatitis.B.agegroup8$Order==10|
                                                            Hepatitis.B.agegroup8$Order==11|
                                                            Hepatitis.B.agegroup8$Order==13|
                                                            Hepatitis.B.agegroup8$Order==15|
                                                            Hepatitis.B.agegroup8$Order==17|
                                                            Hepatitis.B.agegroup8$Order==19|
                                                            Hepatitis.B.agegroup8$Order==21|
                                                            Hepatitis.B.agegroup8$Order==22),],FUN=sum)[,2])
Age.pop5<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup9[which(Hepatitis.B.agegroup9$Order>3&
                                                            Hepatitis.B.agegroup9$Order<9|
                                                            Hepatitis.B.agegroup9$Order==10|
                                                            Hepatitis.B.agegroup9$Order==11|
                                                            Hepatitis.B.agegroup9$Order==13|
                                                            Hepatitis.B.agegroup9$Order==15|
                                                            Hepatitis.B.agegroup9$Order==17|
                                                            Hepatitis.B.agegroup9$Order==19|
                                                            Hepatitis.B.agegroup9$Order==21|
                                                            Hepatitis.B.agegroup9$Order==22),],FUN=sum)[,2])
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order>3&
                                                             Hepatitis.B.agegroup10$Order<9|
                                                             Hepatitis.B.agegroup10$Order==10|
                                                             Hepatitis.B.agegroup10$Order==11|
                                                             Hepatitis.B.agegroup10$Order==13|
                                                             Hepatitis.B.agegroup10$Order==15|
                                                             Hepatitis.B.agegroup10$Order==17|
                                                             Hepatitis.B.agegroup10$Order==19|
                                                             Hepatitis.B.agegroup10$Order==21|
                                                             Hepatitis.B.agegroup10$Order==22),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order==1|
                                                             Hepatitis.B.agegroup11$Order==2|
                                                             Hepatitis.B.agegroup11$Order>3&
                                                             Hepatitis.B.agegroup11$Order<12|
                                                             Hepatitis.B.agegroup11$Order==13|
                                                             Hepatitis.B.agegroup11$Order==15|
                                                             Hepatitis.B.agegroup11$Order==17|
                                                             Hepatitis.B.agegroup11$Order==19|
                                                             Hepatitis.B.agegroup11$Order==21|
                                                             Hepatitis.B.agegroup11$Order==22),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order==1|
                                                             Hepatitis.B.agegroup12$Order==2|
                                                             Hepatitis.B.agegroup12$Order>3&
                                                             Hepatitis.B.agegroup12$Order<12|
                                                             Hepatitis.B.agegroup12$Order==13|
                                                             Hepatitis.B.agegroup12$Order==15|
                                                             Hepatitis.B.agegroup12$Order==17|
                                                             Hepatitis.B.agegroup12$Order==19|
                                                             Hepatitis.B.agegroup12$Order==21|
                                                             Hepatitis.B.agegroup12$Order==22),],FUN=sum)[,2])
Age.pop<-c(Age.pop2,Age.pop3,Age.pop4,Age.pop5,Age.pop6,Age.pop7,Age.pop8)
Prop.pop.highur.2002<-c(Age.pop2/sum(Age.pop),Age.pop3/sum(Age.pop),Age.pop4/sum(Age.pop),
                        Age.pop5/sum(Age.pop),Age.pop6/sum(Age.pop),Age.pop7/sum(Age.pop),
                        Age.pop8/sum(Age.pop))
EIR.merge.2<-EIR.merge.2002(c(4,5,6,7,8,10,11,13,15,17,19,21,22),Hepatitis.B.agegroup6,"model2",c(2,9),13,10000)
EIR.merge.3<-EIR.merge.2002(c(4,5,6,7,8,10,11,13,15,17,19,21,22),Hepatitis.B.agegroup7,"model3",c(2,9),12,10000)
EIR.merge.4<-EIR.merge.2002(c(4,5,6,7,8,10,11,13,15,17,19,21,22),Hepatitis.B.agegroup8,"model4",c(2,9),11,10000)
EIR.merge.5.1<-EIR.merge.2002(c(4,5,6,7,8,10,11,13,15,17,19,21,22),Hepatitis.B.agegroup9,"model5",c(2,10),10,10000)
EIR.merge.6.1<-EIR.merge.2002(c(4,5,6,7,8,10,11,13,15,17,19,21,22),Hepatitis.B.agegroup10,"model6",c(2,10),9,10000)
EIR.merge.7.1<-EIR.merge.2002(c(1,2,4,5,6,7,8,9,10,11,13,15,17,19,21,22),Hepatitis.B.agegroup11,"model7",c(2,10),8,10000)
EIR.merge.8.1<-EIR.merge.2002(c(1,2,4,5,6,7,8,9,10,11,13,15,17,19,21,22),Hepatitis.B.agegroup12,"model8",c(2,10),3,10000)
EIR.merge.highur.2002<-Prop.pop.highur.2002[1]*EIR.merge.2$EIR+Prop.pop.highur.2002[2]*EIR.merge.3$EIR+
  Prop.pop.highur.2002[3]*EIR.merge.4$EIR+Prop.pop.highur.2002[4]*EIR.merge.5.1$EIR+
  Prop.pop.highur.2002[5]*EIR.merge.6.1$EIR+Prop.pop.highur.2002[6]*EIR.merge.7.1$EIR+
  Prop.pop.highur.2002[7]*EIR.merge.8.1$EIR
EIR.point.highur.2002<-Prop.pop.highur.2002[1]*EIR.merge.2$Point.EIR[1]+Prop.pop.highur.2002[2]*EIR.merge.3$Point.EIR[1]+
  Prop.pop.highur.2002[3]*EIR.merge.4$Point.EIR[1]+Prop.pop.highur.2002[4]*EIR.merge.5.1$Point.EIR[1]+
  Prop.pop.highur.2002[5]*EIR.merge.6.1$Point.EIR[1]+Prop.pop.highur.2002[6]*EIR.merge.7.1$Point.EIR[1]+
  Prop.pop.highur.2002[7]*EIR.merge.8.1$Point.EIR[1]
EIR.highur.2002.low<-sort(EIR.merge.highur.2002)[n.cir*0.025]
EIR.highur.2002.high<-sort(EIR.merge.highur.2002)[n.cir*0.975]
EC.merge.highur.2002<-EIR.merge.2$EC+EIR.merge.3$EC+
  EIR.merge.4$EC+EIR.merge.5.1$EC+
  EIR.merge.6.1$EC+EIR.merge.7.1$EC+
  EIR.merge.8.1$EC
EC.point.highur.2002<-EIR.merge.2$Point.EC[1]+EIR.merge.3$Point.EC[1]+
  EIR.merge.4$Point.EC[1]+EIR.merge.5.1$Point.EC[1]+
  EIR.merge.6.1$Point.EC[1]+EIR.merge.7.1$Point.EC[1]+
  EIR.merge.8.1$Point.EC[1]
EC.highur.2002.low<-sort(EC.merge.highur.2002)[n.cir*0.025]
EC.highur.2002.high<-sort(EC.merge.highur.2002)[n.cir*0.975]

#######################Hospitalisation bed density subgroup EIR, 2002###################################
#Hospitalisation bed density category
LowHB<-c(3,11:14,16,18:22,24:26,28)
HighHB<-c(1,2,4,5,6,7,8,9,10,15,17,23,27,29,30,31)
#Low hospitalisation bed density
Age.pop2<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup6[which(Hepatitis.B.agegroup6$Order==3|
                                                            Hepatitis.B.agegroup6$Order>10&
                                                            Hepatitis.B.agegroup6$Order<15|
                                                            Hepatitis.B.agegroup6$Order==16|
                                                            Hepatitis.B.agegroup6$Order>17&
                                                            Hepatitis.B.agegroup6$Order<23|
                                                            Hepatitis.B.agegroup6$Order==24|
                                                            Hepatitis.B.agegroup6$Order==25|
                                                            Hepatitis.B.agegroup6$Order==28),],FUN=sum)[,2])
Age.pop3<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup7[which(Hepatitis.B.agegroup7$Order==3|
                                                            Hepatitis.B.agegroup7$Order>10&
                                                            Hepatitis.B.agegroup7$Order<14|
                                                            Hepatitis.B.agegroup7$Order==16|
                                                            Hepatitis.B.agegroup7$Order>17&
                                                            Hepatitis.B.agegroup7$Order<23|
                                                            Hepatitis.B.agegroup7$Order==24|
                                                            Hepatitis.B.agegroup7$Order==25|
                                                            Hepatitis.B.agegroup7$Order==28),],FUN=sum)[,2])
Age.pop4<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup8[which(Hepatitis.B.agegroup8$Order==3|
                                                            Hepatitis.B.agegroup8$Order>10&
                                                            Hepatitis.B.agegroup8$Order<14|
                                                            Hepatitis.B.agegroup8$Order==16|
                                                            Hepatitis.B.agegroup8$Order>17&
                                                            Hepatitis.B.agegroup8$Order<23|
                                                            Hepatitis.B.agegroup8$Order==24|
                                                            Hepatitis.B.agegroup8$Order==25|
                                                            Hepatitis.B.agegroup8$Order==28),],FUN=sum)[,2])
Age.pop5<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup9[which(Hepatitis.B.agegroup9$Order==3|
                                                            Hepatitis.B.agegroup9$Order>10&
                                                            Hepatitis.B.agegroup9$Order<14|
                                                            Hepatitis.B.agegroup9$Order==16|
                                                            Hepatitis.B.agegroup9$Order>17&
                                                            Hepatitis.B.agegroup9$Order<23|
                                                            Hepatitis.B.agegroup9$Order==24|
                                                            Hepatitis.B.agegroup9$Order==25|
                                                            Hepatitis.B.agegroup9$Order==28),],FUN=sum)[,2])
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order==3|
                                                             Hepatitis.B.agegroup10$Order>10&
                                                             Hepatitis.B.agegroup10$Order<14|
                                                             Hepatitis.B.agegroup10$Order==16|
                                                             Hepatitis.B.agegroup10$Order>17&
                                                             Hepatitis.B.agegroup10$Order<23|
                                                             Hepatitis.B.agegroup10$Order==24|
                                                             Hepatitis.B.agegroup10$Order==25|
                                                             Hepatitis.B.agegroup10$Order==28),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order==3|
                                                             Hepatitis.B.agegroup11$Order>10&
                                                             Hepatitis.B.agegroup11$Order<15|
                                                             Hepatitis.B.agegroup11$Order==16|
                                                             Hepatitis.B.agegroup11$Order>17&
                                                             Hepatitis.B.agegroup11$Order<23|
                                                             Hepatitis.B.agegroup11$Order==24|
                                                             Hepatitis.B.agegroup11$Order==25|
                                                             Hepatitis.B.agegroup11$Order==26|
                                                             Hepatitis.B.agegroup11$Order==28),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order==3|
                                                             Hepatitis.B.agegroup12$Order>10&
                                                             Hepatitis.B.agegroup12$Order<15|
                                                             Hepatitis.B.agegroup12$Order==16|
                                                             Hepatitis.B.agegroup12$Order>17&
                                                             Hepatitis.B.agegroup12$Order<23|
                                                             Hepatitis.B.agegroup12$Order==24|
                                                             Hepatitis.B.agegroup12$Order==25|
                                                             Hepatitis.B.agegroup12$Order==26|
                                                             Hepatitis.B.agegroup12$Order==28),],FUN=sum)[,2])
Age.pop<-c(Age.pop2,Age.pop3,Age.pop4,Age.pop5,Age.pop6,Age.pop7,Age.pop8)
Prop.pop.lowhb.2002<-c(Age.pop2/sum(Age.pop),Age.pop3/sum(Age.pop),Age.pop4/sum(Age.pop),
                       Age.pop5/sum(Age.pop),Age.pop6/sum(Age.pop),Age.pop7/sum(Age.pop),
                       Age.pop8/sum(Age.pop))
EIR.merge.2<-EIR.merge.2002(c(3,11:14,16,18:22,24:25,28),Hepatitis.B.agegroup6,"model2",c(2,9),13,10000)
EIR.merge.3<-EIR.merge.2002(c(3,11:13,16,18:22,24:25,28),Hepatitis.B.agegroup7,"model3",c(2,9),12,10000)
EIR.merge.4<-EIR.merge.2002(c(3,11:13,16,18:22,24:25,28),Hepatitis.B.agegroup8,"model4",c(2,9),11,10000)
EIR.merge.5.1<-EIR.merge.2002(c(3,11:13,16,18:22,24:25,28),Hepatitis.B.agegroup9,"model5",c(2,10),10,10000)
EIR.merge.6.1<-EIR.merge.2002(c(3,11:13,16,18:22,24:25,28),Hepatitis.B.agegroup10,"model6",c(2,10),9,10000)
EIR.merge.7.1<-EIR.merge.2002(c(3,11:14,16,18:22,24:26,28),Hepatitis.B.agegroup11,"model7",c(2,10),8,10000)
EIR.merge.8.1<-EIR.merge.2002(c(3,11:14,16,18:22,24:26,28),Hepatitis.B.agegroup12,"model8",c(2,10),3,10000)
EIR.merge.lowhb.2002<-Prop.pop.lowhb.2002[1]*EIR.merge.2$EIR+Prop.pop.lowhb.2002[2]*EIR.merge.3$EIR+
  Prop.pop.lowhb.2002[3]*EIR.merge.4$EIR+Prop.pop.lowhb.2002[4]*EIR.merge.5.1$EIR+
  Prop.pop.lowhb.2002[5]*EIR.merge.6.1$EIR+Prop.pop.lowhb.2002[6]*EIR.merge.7.1$EIR+
  Prop.pop.lowhb.2002[7]*EIR.merge.8.1$EIR
EIR.point.lowhb.2002<-Prop.pop.lowhb.2002[1]*EIR.merge.2$Point.EIR[1]+Prop.pop.lowhb.2002[2]*EIR.merge.3$Point.EIR[1]+
  Prop.pop.lowhb.2002[3]*EIR.merge.4$Point.EIR[1]+Prop.pop.lowhb.2002[4]*EIR.merge.5.1$Point.EIR[1]+
  Prop.pop.lowhb.2002[5]*EIR.merge.6.1$Point.EIR[1]+Prop.pop.lowhb.2002[6]*EIR.merge.7.1$Point.EIR[1]+
  Prop.pop.lowhb.2002[7]*EIR.merge.8.1$Point.EIR[1]
EIR.lowhb.2002.low<-sort(EIR.merge.lowhb.2002)[n.cir*0.025]
EIR.lowhb.2002.high<-sort(EIR.merge.lowhb.2002)[n.cir*0.975]
EC.merge.lowhb.2002<-EIR.merge.2$EC+EIR.merge.3$EC+
  EIR.merge.4$EC+EIR.merge.5.1$EC+
  EIR.merge.6.1$EC+EIR.merge.7.1$EC+
  EIR.merge.8.1$EC
EC.point.lowhb.2002<-EIR.merge.2$Point.EC[1]+EIR.merge.3$Point.EC[1]+
  EIR.merge.4$Point.EC[1]+EIR.merge.5.1$Point.EC[1]+
  EIR.merge.6.1$Point.EC[1]+EIR.merge.7.1$Point.EC[1]+
  EIR.merge.8.1$Point.EC[1]
EC.lowhb.2002.low<-sort(EC.merge.lowhb.2002)[n.cir*0.025]
EC.lowhb.2002.high<-sort(EC.merge.lowhb.2002)[n.cir*0.975]

#High HBsAg prevalence
Age.pop2<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup6[which(Hepatitis.B.agegroup6$Order==4|
                                                            Hepatitis.B.agegroup6$Order==5|
                                                            Hepatitis.B.agegroup6$Order==7|
                                                            Hepatitis.B.agegroup6$Order==8|
                                                            Hepatitis.B.agegroup6$Order==10|
                                                            Hepatitis.B.agegroup6$Order==15|
                                                            Hepatitis.B.agegroup6$Order==17|
                                                            Hepatitis.B.agegroup6$Order==23|
                                                            Hepatitis.B.agegroup6$Order==27|
                                                            Hepatitis.B.agegroup6$Order==29|
                                                            Hepatitis.B.agegroup6$Order==30|
                                                            Hepatitis.B.agegroup6$Order==31),],FUN=sum)[,2])
Age.pop3<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup7[which(Hepatitis.B.agegroup7$Order>3&
                                                            Hepatitis.B.agegroup7$Order<9|
                                                            Hepatitis.B.agegroup7$Order==10|
                                                            Hepatitis.B.agegroup7$Order==15|
                                                            Hepatitis.B.agegroup7$Order==17|
                                                            Hepatitis.B.agegroup7$Order==23|
                                                            Hepatitis.B.agegroup7$Order==27|
                                                            Hepatitis.B.agegroup7$Order==29|
                                                            Hepatitis.B.agegroup7$Order==30|
                                                            Hepatitis.B.agegroup7$Order==31),],FUN=sum)[,2])
Age.pop4<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup8[which(Hepatitis.B.agegroup8$Order>3&
                                                            Hepatitis.B.agegroup8$Order<9|
                                                            Hepatitis.B.agegroup8$Order==10|
                                                            Hepatitis.B.agegroup8$Order==15|
                                                            Hepatitis.B.agegroup8$Order==17|
                                                            Hepatitis.B.agegroup8$Order==23|
                                                            Hepatitis.B.agegroup8$Order==27|
                                                            Hepatitis.B.agegroup8$Order==29|
                                                            Hepatitis.B.agegroup8$Order==30|
                                                            Hepatitis.B.agegroup8$Order==31),],FUN=sum)[,2])
Age.pop5<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup9[which(Hepatitis.B.agegroup9$Order>3&
                                                            Hepatitis.B.agegroup9$Order<9|
                                                            Hepatitis.B.agegroup9$Order==10|
                                                            Hepatitis.B.agegroup9$Order==15|
                                                            Hepatitis.B.agegroup9$Order==17|
                                                            Hepatitis.B.agegroup9$Order==23|
                                                            Hepatitis.B.agegroup9$Order==27|
                                                            Hepatitis.B.agegroup9$Order==29|
                                                            Hepatitis.B.agegroup9$Order==30|
                                                            Hepatitis.B.agegroup9$Order==31),],FUN=sum)[,2])
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order>3&
                                                             Hepatitis.B.agegroup10$Order<9|
                                                             Hepatitis.B.agegroup10$Order==10|
                                                             Hepatitis.B.agegroup10$Order==15|
                                                             Hepatitis.B.agegroup10$Order==17|
                                                             Hepatitis.B.agegroup10$Order==23|
                                                             Hepatitis.B.agegroup10$Order==27|
                                                             Hepatitis.B.agegroup10$Order==29|
                                                             Hepatitis.B.agegroup10$Order==30|
                                                             Hepatitis.B.agegroup10$Order==31),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order==1|
                                                             Hepatitis.B.agegroup11$Order==2|
                                                             Hepatitis.B.agegroup11$Order>3&
                                                             Hepatitis.B.agegroup11$Order<11|
                                                             Hepatitis.B.agegroup11$Order==15|
                                                             Hepatitis.B.agegroup11$Order==17|
                                                             Hepatitis.B.agegroup11$Order==23|
                                                             Hepatitis.B.agegroup11$Order==27|
                                                             Hepatitis.B.agegroup11$Order==29|
                                                             Hepatitis.B.agegroup11$Order==30|
                                                             Hepatitis.B.agegroup11$Order==31),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order==1|
                                                             Hepatitis.B.agegroup12$Order==2|
                                                             Hepatitis.B.agegroup12$Order>3&
                                                             Hepatitis.B.agegroup12$Order<11|
                                                             Hepatitis.B.agegroup12$Order==15|
                                                             Hepatitis.B.agegroup12$Order==17|
                                                             Hepatitis.B.agegroup12$Order==23|
                                                             Hepatitis.B.agegroup12$Order==27|
                                                             Hepatitis.B.agegroup12$Order==29|
                                                             Hepatitis.B.agegroup12$Order==30|
                                                             Hepatitis.B.agegroup12$Order==31),],FUN=sum)[,2])
Age.pop<-c(Age.pop2,Age.pop3,Age.pop4,Age.pop5,Age.pop6,Age.pop7,Age.pop8)
Prop.pop.highhb.2002<-c(Age.pop2/sum(Age.pop),Age.pop3/sum(Age.pop),Age.pop4/sum(Age.pop),
                        Age.pop5/sum(Age.pop),Age.pop6/sum(Age.pop),Age.pop7/sum(Age.pop),
                        Age.pop8/sum(Age.pop))
EIR.merge.2<-EIR.merge.2002(c(4,5,7,8,10,15,17,23,27,29,30,31),Hepatitis.B.agegroup6,"model2",c(2,9),13,10000)
EIR.merge.3<-EIR.merge.2002(c(4,5,6,7,8,10,15,17,23,27,29,30,31),Hepatitis.B.agegroup7,"model3",c(2,9),12,10000)
EIR.merge.4<-EIR.merge.2002(c(4,5,6,7,8,10,15,17,23,27,29,30,31),Hepatitis.B.agegroup8,"model4",c(2,9),11,10000)
EIR.merge.5.1<-EIR.merge.2002(c(4,5,6,7,8,10,15,17,23,27,29,30,31),Hepatitis.B.agegroup9,"model5",c(2,10),10,10000)
EIR.merge.6.1<-EIR.merge.2002(c(4,5,6,7,8,10,15,17,23,27,29,30,31),Hepatitis.B.agegroup10,"model6",c(2,10),9,10000)
EIR.merge.7.1<-EIR.merge.2002(c(1,2,4,5,6,7,8,9,10,15,17,23,27,29,30,31),Hepatitis.B.agegroup11,"model7",c(2,10),8,10000)
EIR.merge.8.1<-EIR.merge.2002(c(1,2,4,5,6,7,8,9,10,15,17,23,27,29,30,31),Hepatitis.B.agegroup12,"model8",c(2,10),3,10000)
EIR.merge.highhb.2002<-Prop.pop.highhb.2002[1]*EIR.merge.2$EIR+Prop.pop.highhb.2002[2]*EIR.merge.3$EIR+
  Prop.pop.highhb.2002[3]*EIR.merge.4$EIR+Prop.pop.highhb.2002[4]*EIR.merge.5.1$EIR+
  Prop.pop.highhb.2002[5]*EIR.merge.6.1$EIR+Prop.pop.highhb.2002[6]*EIR.merge.7.1$EIR+
  Prop.pop.highhb.2002[7]*EIR.merge.8.1$EIR
EIR.point.highhb.2002<-Prop.pop.highhb.2002[1]*EIR.merge.2$Point.EIR[1]+Prop.pop.highhb.2002[2]*EIR.merge.3$Point.EIR[1]+
  Prop.pop.highhb.2002[3]*EIR.merge.4$Point.EIR[1]+Prop.pop.highhb.2002[4]*EIR.merge.5.1$Point.EIR[1]+
  Prop.pop.highhb.2002[5]*EIR.merge.6.1$Point.EIR[1]+Prop.pop.highhb.2002[6]*EIR.merge.7.1$Point.EIR[1]+
  Prop.pop.highhb.2002[7]*EIR.merge.8.1$Point.EIR[1]
EIR.highhb.2002.low<-sort(EIR.merge.highhb.2002)[n.cir*0.025]
EIR.highhb.2002.high<-sort(EIR.merge.highhb.2002)[n.cir*0.975]
EC.merge.highhb.2002<-EIR.merge.2$EC+EIR.merge.3$EC+
  EIR.merge.4$EC+EIR.merge.5.1$EC+
  EIR.merge.6.1$EC+EIR.merge.7.1$EC+
  EIR.merge.8.1$EC
EC.point.highhb.2002<-EIR.merge.2$Point.EC[1]+EIR.merge.3$Point.EC[1]+
  EIR.merge.4$Point.EC[1]+EIR.merge.5.1$Point.EC[1]+
  EIR.merge.6.1$Point.EC[1]+EIR.merge.7.1$Point.EC[1]+
  EIR.merge.8.1$Point.EC[1]
EC.highhb.2002.low<-sort(EC.merge.highhb.2002)[n.cir*0.025]
EC.highhb.2002.high<-sort(EC.merge.highhb.2002)[n.cir*0.975]

#######################Illiteracy rate subgroup EIR, 2002###################################
#Illiteracy rate category
LowIR<-c(1:10,14,18,19,20,21,31)
HighIR<-c(11,12,13,15,16,17,22:30)
#Low Illiteracy rate
Age.pop2<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup6[which(Hepatitis.B.agegroup6$Order>2&
                                                            Hepatitis.B.agegroup6$Order<6|
                                                            Hepatitis.B.agegroup6$Order==7|
                                                            Hepatitis.B.agegroup6$Order==8|
                                                            Hepatitis.B.agegroup6$Order==10|
                                                            Hepatitis.B.agegroup6$Order==14|
                                                            Hepatitis.B.agegroup6$Order==18|
                                                            Hepatitis.B.agegroup6$Order==19|
                                                            Hepatitis.B.agegroup6$Order==20|
                                                            Hepatitis.B.agegroup6$Order==21|
                                                            Hepatitis.B.agegroup6$Order==31),],FUN=sum)[,2])
Age.pop3<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup7[which(Hepatitis.B.agegroup7$Order>2&
                                                            Hepatitis.B.agegroup7$Order<9|
                                                            Hepatitis.B.agegroup7$Order==10|
                                                            Hepatitis.B.agegroup7$Order==18|
                                                            Hepatitis.B.agegroup7$Order==19|
                                                            Hepatitis.B.agegroup7$Order==20|
                                                            Hepatitis.B.agegroup7$Order==21|
                                                            Hepatitis.B.agegroup7$Order==31),],FUN=sum)[,2])
Age.pop4<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup8[which(Hepatitis.B.agegroup8$Order>2&
                                                            Hepatitis.B.agegroup8$Order<9|
                                                            Hepatitis.B.agegroup8$Order==10|
                                                            Hepatitis.B.agegroup8$Order==18|
                                                            Hepatitis.B.agegroup8$Order==19|
                                                            Hepatitis.B.agegroup8$Order==20|
                                                            Hepatitis.B.agegroup8$Order==21|
                                                            Hepatitis.B.agegroup8$Order==31),],FUN=sum)[,2])
Age.pop5<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup9[which(Hepatitis.B.agegroup9$Order>2&
                                                            Hepatitis.B.agegroup9$Order<9|
                                                            Hepatitis.B.agegroup9$Order==10|
                                                            Hepatitis.B.agegroup9$Order==18|
                                                            Hepatitis.B.agegroup9$Order==19|
                                                            Hepatitis.B.agegroup9$Order==20|
                                                            Hepatitis.B.agegroup9$Order==21|
                                                            Hepatitis.B.agegroup9$Order==31),],FUN=sum)[,2])
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order>2&
                                                             Hepatitis.B.agegroup10$Order<9|
                                                             Hepatitis.B.agegroup10$Order==10|
                                                             Hepatitis.B.agegroup10$Order==18|
                                                             Hepatitis.B.agegroup10$Order==19|
                                                             Hepatitis.B.agegroup10$Order==20|
                                                             Hepatitis.B.agegroup10$Order==21|
                                                             Hepatitis.B.agegroup10$Order==31),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order>0&
                                                             Hepatitis.B.agegroup11$Order<11|
                                                             Hepatitis.B.agegroup11$Order==14|
                                                             Hepatitis.B.agegroup11$Order==18|
                                                             Hepatitis.B.agegroup11$Order==19|
                                                             Hepatitis.B.agegroup11$Order==20|
                                                             Hepatitis.B.agegroup11$Order==21|
                                                             Hepatitis.B.agegroup11$Order==31),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order>0&
                                                             Hepatitis.B.agegroup12$Order<11|
                                                             Hepatitis.B.agegroup12$Order==14|
                                                             Hepatitis.B.agegroup12$Order==18|
                                                             Hepatitis.B.agegroup12$Order==19|
                                                             Hepatitis.B.agegroup12$Order==20|
                                                             Hepatitis.B.agegroup12$Order==21|
                                                             Hepatitis.B.agegroup12$Order==31),],FUN=sum)[,2])
Age.pop<-c(Age.pop2,Age.pop3,Age.pop4,Age.pop5,Age.pop6,Age.pop7,Age.pop8)
Prop.pop.lowir.2002<-c(Age.pop2/sum(Age.pop),Age.pop3/sum(Age.pop),Age.pop4/sum(Age.pop),
                       Age.pop5/sum(Age.pop),Age.pop6/sum(Age.pop),Age.pop7/sum(Age.pop),
                       Age.pop8/sum(Age.pop))
EIR.merge.2<-EIR.merge.2002(c(3,4,5,7,8,10,14,18,19,20,21,31),Hepatitis.B.agegroup6,"model2",c(2,9),13,10000)
EIR.merge.3<-EIR.merge.2002(c(3:8,10,18,19,20,21,31),Hepatitis.B.agegroup7,"model3",c(2,9),12,10000)
EIR.merge.4<-EIR.merge.2002(c(3:8,10,18,19,20,21,31),Hepatitis.B.agegroup8,"model4",c(2,9),11,10000)
EIR.merge.5.1<-EIR.merge.2002(c(3:8,10,18,19,20,21,31),Hepatitis.B.agegroup9,"model5",c(2,10),10,10000)
EIR.merge.6.1<-EIR.merge.2002(c(3:8,10,18,19,20,21,31),Hepatitis.B.agegroup10,"model6",c(2,10),9,10000)
EIR.merge.7.1<-EIR.merge.2002(c(1:10,14,18,19,20,21,31),Hepatitis.B.agegroup11,"model7",c(2,10),8,10000)
EIR.merge.8.1<-EIR.merge.2002(c(1:10,14,18,19,20,21,31),Hepatitis.B.agegroup12,"model8",c(2,10),3,10000)
EIR.merge.lowir.2002<-Prop.pop.lowir.2002[1]*EIR.merge.2$EIR+Prop.pop.lowir.2002[2]*EIR.merge.3$EIR+
  Prop.pop.lowir.2002[3]*EIR.merge.4$EIR+Prop.pop.lowir.2002[4]*EIR.merge.5.1$EIR+
  Prop.pop.lowir.2002[5]*EIR.merge.6.1$EIR+Prop.pop.lowir.2002[6]*EIR.merge.7.1$EIR+
  Prop.pop.lowir.2002[7]*EIR.merge.8.1$EIR
EIR.point.lowir.2002<-Prop.pop.lowir.2002[1]*EIR.merge.2$Point.EIR[1]+Prop.pop.lowir.2002[2]*EIR.merge.3$Point.EIR[1]+
  Prop.pop.lowir.2002[3]*EIR.merge.4$Point.EIR[1]+Prop.pop.lowir.2002[4]*EIR.merge.5.1$Point.EIR[1]+
  Prop.pop.lowir.2002[5]*EIR.merge.6.1$Point.EIR[1]+Prop.pop.lowir.2002[6]*EIR.merge.7.1$Point.EIR[1]+
  Prop.pop.lowir.2002[7]*EIR.merge.8.1$Point.EIR[1]
EIR.lowir.2002.low<-sort(EIR.merge.lowir.2002)[n.cir*0.025]
EIR.lowir.2002.high<-sort(EIR.merge.lowir.2002)[n.cir*0.975]
EC.merge.lowir.2002<-EIR.merge.2$EC+EIR.merge.3$EC+
  EIR.merge.4$EC+EIR.merge.5.1$EC+
  EIR.merge.6.1$EC+EIR.merge.7.1$EC+
  EIR.merge.8.1$EC
EC.point.lowir.2002<-EIR.merge.2$Point.EC[1]+EIR.merge.3$Point.EC[1]+
  EIR.merge.4$Point.EC[1]+EIR.merge.5.1$Point.EC[1]+
  EIR.merge.6.1$Point.EC[1]+EIR.merge.7.1$Point.EC[1]+
  EIR.merge.8.1$Point.EC[1]
EC.lowir.2002.low<-sort(EC.merge.lowir.2002)[n.cir*0.025]
EC.lowir.2002.high<-sort(EC.merge.lowir.2002)[n.cir*0.975]

#High illiteracy rate
Age.pop2<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup6[which(Hepatitis.B.agegroup6$Order==11|
                                                            Hepatitis.B.agegroup6$Order==12|
                                                            Hepatitis.B.agegroup6$Order==13|
                                                            Hepatitis.B.agegroup6$Order==15|
                                                            Hepatitis.B.agegroup6$Order==16|
                                                            Hepatitis.B.agegroup6$Order==17|
                                                            Hepatitis.B.agegroup6$Order>21&
                                                            Hepatitis.B.agegroup6$Order<26|
                                                            Hepatitis.B.agegroup6$Order>26&
                                                            Hepatitis.B.agegroup6$Order<31),],FUN=sum)[,2])
Age.pop3<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup7[which(Hepatitis.B.agegroup7$Order==11|
                                                            Hepatitis.B.agegroup7$Order==12|
                                                            Hepatitis.B.agegroup7$Order==13|
                                                            Hepatitis.B.agegroup7$Order==15|
                                                            Hepatitis.B.agegroup7$Order==16|
                                                            Hepatitis.B.agegroup7$Order==17|
                                                            Hepatitis.B.agegroup7$Order>21&
                                                            Hepatitis.B.agegroup7$Order<26|
                                                            Hepatitis.B.agegroup7$Order>26&
                                                            Hepatitis.B.agegroup7$Order<31),],FUN=sum)[,2])
Age.pop4<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup8[which(Hepatitis.B.agegroup8$Order==11|
                                                            Hepatitis.B.agegroup8$Order==12|
                                                            Hepatitis.B.agegroup8$Order==13|
                                                            Hepatitis.B.agegroup8$Order==15|
                                                            Hepatitis.B.agegroup8$Order==16|
                                                            Hepatitis.B.agegroup8$Order==17|
                                                            Hepatitis.B.agegroup8$Order>21&
                                                            Hepatitis.B.agegroup8$Order<26|
                                                            Hepatitis.B.agegroup8$Order>26&
                                                            Hepatitis.B.agegroup8$Order<31),],FUN=sum)[,2])
Age.pop5<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup9[which(Hepatitis.B.agegroup9$Order==11|
                                                            Hepatitis.B.agegroup9$Order==12|
                                                            Hepatitis.B.agegroup9$Order==13|
                                                            Hepatitis.B.agegroup9$Order==15|
                                                            Hepatitis.B.agegroup9$Order==16|
                                                            Hepatitis.B.agegroup9$Order==17|
                                                            Hepatitis.B.agegroup9$Order>21&
                                                            Hepatitis.B.agegroup9$Order<26|
                                                            Hepatitis.B.agegroup9$Order>26&
                                                            Hepatitis.B.agegroup9$Order<31),],FUN=sum)[,2])
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order==11|
                                                             Hepatitis.B.agegroup10$Order==12|
                                                             Hepatitis.B.agegroup10$Order==13|
                                                             Hepatitis.B.agegroup10$Order==15|
                                                             Hepatitis.B.agegroup10$Order==16|
                                                             Hepatitis.B.agegroup10$Order==17|
                                                             Hepatitis.B.agegroup10$Order>21&
                                                             Hepatitis.B.agegroup10$Order<26|
                                                             Hepatitis.B.agegroup10$Order>26&
                                                             Hepatitis.B.agegroup10$Order<31),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order==11|
                                                             Hepatitis.B.agegroup11$Order==12|
                                                             Hepatitis.B.agegroup11$Order==13|
                                                             Hepatitis.B.agegroup11$Order==15|
                                                             Hepatitis.B.agegroup11$Order==16|
                                                             Hepatitis.B.agegroup11$Order==17|
                                                             Hepatitis.B.agegroup11$Order>21&
                                                             Hepatitis.B.agegroup11$Order<31),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order==11|
                                                             Hepatitis.B.agegroup12$Order==12|
                                                             Hepatitis.B.agegroup12$Order==13|
                                                             Hepatitis.B.agegroup12$Order==15|
                                                             Hepatitis.B.agegroup12$Order==16|
                                                             Hepatitis.B.agegroup12$Order==17|
                                                             Hepatitis.B.agegroup12$Order>21&
                                                             Hepatitis.B.agegroup12$Order<31),],FUN=sum)[,2])
Age.pop<-c(Age.pop2,Age.pop3,Age.pop4,Age.pop5,Age.pop6,Age.pop7,Age.pop8)
Prop.pop.highir.2002<-c(Age.pop2/sum(Age.pop),Age.pop3/sum(Age.pop),Age.pop4/sum(Age.pop),
                        Age.pop5/sum(Age.pop),Age.pop6/sum(Age.pop),Age.pop7/sum(Age.pop),
                        Age.pop8/sum(Age.pop))
EIR.merge.2<-EIR.merge.2002(c(11,12,13,15,16,17,22:25,27:30),Hepatitis.B.agegroup6,"model2",c(2,9),13,10000)
EIR.merge.3<-EIR.merge.2002(c(11,12,13,15,16,17,22:25,27:30),Hepatitis.B.agegroup7,"model3",c(2,9),12,10000)
EIR.merge.4<-EIR.merge.2002(c(11,12,13,15,16,17,22:25,27:30),Hepatitis.B.agegroup8,"model4",c(2,9),11,10000)
EIR.merge.5.1<-EIR.merge.2002(c(11,12,13,15,16,17,22:25,27:30),Hepatitis.B.agegroup9,"model5",c(2,10),10,10000)
EIR.merge.6.1<-EIR.merge.2002(c(11,12,13,15,16,17,22:25,27:30),Hepatitis.B.agegroup10,"model6",c(2,10),9,10000)
EIR.merge.7.1<-EIR.merge.2002(c(11,12,13,15,16,17,22:30),Hepatitis.B.agegroup11,"model7",c(2,10),8,10000)
EIR.merge.8.1<-EIR.merge.2002(c(11,12,13,15,16,17,22:30),Hepatitis.B.agegroup12,"model8",c(2,10),3,10000)
EIR.merge.highir.2002<-Prop.pop.highir.2002[1]*EIR.merge.2$EIR+Prop.pop.highir.2002[2]*EIR.merge.3$EIR+
  Prop.pop.highir.2002[3]*EIR.merge.4$EIR+Prop.pop.highir.2002[4]*EIR.merge.5.1$EIR+
  Prop.pop.highir.2002[5]*EIR.merge.6.1$EIR+Prop.pop.highir.2002[6]*EIR.merge.7.1$EIR+
  Prop.pop.highir.2002[7]*EIR.merge.8.1$EIR
EIR.point.highir.2002<-Prop.pop.highir.2002[1]*EIR.merge.2$Point.EIR[1]+Prop.pop.highir.2002[2]*EIR.merge.3$Point.EIR[1]+
  Prop.pop.highir.2002[3]*EIR.merge.4$Point.EIR[1]+Prop.pop.highir.2002[4]*EIR.merge.5.1$Point.EIR[1]+
  Prop.pop.highir.2002[5]*EIR.merge.6.1$Point.EIR[1]+Prop.pop.highir.2002[6]*EIR.merge.7.1$Point.EIR[1]+
  Prop.pop.highir.2002[7]*EIR.merge.8.1$Point.EIR[1]
EIR.highir.2002.low<-sort(EIR.merge.highir.2002)[n.cir*0.025]
EIR.highir.2002.high<-sort(EIR.merge.highir.2002)[n.cir*0.975]
EC.merge.highir.2002<-EIR.merge.2$EC+EIR.merge.3$EC+
  EIR.merge.4$EC+EIR.merge.5.1$EC+
  EIR.merge.6.1$EC+EIR.merge.7.1$EC+
  EIR.merge.8.1$EC
EC.point.highir.2002<-EIR.merge.2$Point.EC[1]+EIR.merge.3$Point.EC[1]+
  EIR.merge.4$Point.EC[1]+EIR.merge.5.1$Point.EC[1]+
  EIR.merge.6.1$Point.EC[1]+EIR.merge.7.1$Point.EC[1]+
  EIR.merge.8.1$Point.EC[1]
EC.highir.2002.low<-sort(EC.merge.highir.2002)[n.cir*0.025]
EC.highir.2002.high<-sort(EC.merge.highir.2002)[n.cir*0.975]

#######################Region subgroup EIR, 2009###################################
#Region category
North<-c(1,2,3,4,5)
Northeast<-c(6,7,8)
East<-c(9,10,11,12,13,14,15)
Southcentral<-c(16,17,18,19,20,21)
Southwest<-c(22,23,24,25,26)
Northwest<-c(27,28,29,30,31)
#North
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order==3|
                                                             Hepatitis.B.agegroup10$Order==4|
                                                             Hepatitis.B.agegroup10$Order==5),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order>0&
                                                             Hepatitis.B.agegroup11$Order<6),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order>0&
                                                             Hepatitis.B.agegroup12$Order<6),],FUN=sum)[,2])
Age.pop9<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup13[which(Hepatitis.B.agegroup13$Order>0&
                                                             Hepatitis.B.agegroup13$Order<6),],FUN=sum)[,2])
Age.pop10<-mean(aggregate(Population.book~Time,
                          data=Hepatitis.B.agegroup14[which(Hepatitis.B.agegroup14$Order>0&
                                                              Hepatitis.B.agegroup14$Order<6),],FUN=sum)[,2])
Age.pop.2009<-c(Age.pop6,Age.pop7,Age.pop8,Age.pop9,Age.pop10)
Prop.pop.north.2009<-c(Age.pop6/sum(Age.pop.2009),Age.pop7/sum(Age.pop.2009),
                       Age.pop8/sum(Age.pop.2009),Age.pop9/sum(Age.pop.2009),Age.pop10/sum(Age.pop.2009))
EIR.merge.6.2<-EIR.merge.2009(c(3,4,5),Hepatitis.B.agegroup10,"model6",c(3,11),(2+12)/12,10000)
EIR.merge.7.2<-EIR.merge.2009(c(1,2,3,4,5),Hepatitis.B.agegroup11,"model7",c(3,11),(2+2*12)/12,10000)
EIR.merge.8.2<-EIR.merge.2009(c(1,2,3,4,5),Hepatitis.B.agegroup12,"model8",c(3,11),(2+7*12)/12,10000)
EIR.merge.9<-EIR.merge.2009(c(1,2,3,4,5),Hepatitis.B.agegroup13,"model9",c(2,9),6,10000)
EIR.merge.10<-EIR.merge.2009(c(1,2,3,4,5),Hepatitis.B.agegroup14,"model10",c(2,9),1,10000)
EIR.merge.north.2009<-Prop.pop.north.2009[1]*EIR.merge.6.2$EIR+
  Prop.pop.north.2009[2]*EIR.merge.7.2$EIR+Prop.pop.north.2009[3]*EIR.merge.8.2$EIR+
  Prop.pop.north.2009[4]*EIR.merge.9$EIR+Prop.pop.north.2009[5]*EIR.merge.10$EIR
EIR.point.north.2009<-Prop.pop.north.2009[1]*EIR.merge.6.2$Point.EIR[1]+
  Prop.pop.north.2009[2]*EIR.merge.7.2$Point.EIR[1]+Prop.pop.north.2009[3]*EIR.merge.8.2$Point.EIR[1]+
  Prop.pop.north.2009[4]*EIR.merge.9$Point.EIR[1]+Prop.pop.north.2009[5]*EIR.merge.10$Point.EIR[1]
EIR.north.2009.low<-sort(EIR.merge.north.2009)[n.cir*0.025]
EIR.north.2009.high<-sort(EIR.merge.north.2009)[n.cir*0.975]
EC.merge.north.2009<-EIR.merge.6.2$EC+
  EIR.merge.7.2$EC+EIR.merge.8.2$EC+
  EIR.merge.9$EC+EIR.merge.10$EC
EC.point.north.2009<-EIR.merge.6.2$Point.EC[1]+
  EIR.merge.7.2$Point.EC[1]+EIR.merge.8.2$Point.EC[1]+
  EIR.merge.9$Point.EC[1]+EIR.merge.10$Point.EC[1]
EC.north.2009.low<-sort(EC.merge.north.2009)[n.cir*0.025]
EC.north.2009.high<-sort(EC.merge.north.2009)[n.cir*0.975]

#Northeast
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order>5&
                                                             Hepatitis.B.agegroup10$Order<9),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order>5&
                                                             Hepatitis.B.agegroup11$Order<9),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order>5&
                                                             Hepatitis.B.agegroup12$Order<9),],FUN=sum)[,2])
Age.pop9<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup13[which(Hepatitis.B.agegroup13$Order>5&
                                                             Hepatitis.B.agegroup13$Order<9),],FUN=sum)[,2])
Age.pop10<-mean(aggregate(Population.book~Time,
                          data=Hepatitis.B.agegroup14[which(Hepatitis.B.agegroup14$Order>5&
                                                              Hepatitis.B.agegroup14$Order<9),],FUN=sum)[,2])
Age.pop.2009<-c(Age.pop6,Age.pop7,Age.pop8,Age.pop9,Age.pop10)
Prop.pop.northeast.2009<-c(Age.pop6/sum(Age.pop.2009),Age.pop7/sum(Age.pop.2009),
                           Age.pop8/sum(Age.pop.2009),Age.pop9/sum(Age.pop.2009),Age.pop10/sum(Age.pop.2009))
EIR.merge.6.2<-EIR.merge.2009(c(6:8),Hepatitis.B.agegroup10,"model6",c(3,11),(2+12)/12,10000)
EIR.merge.7.2<-EIR.merge.2009(c(6:8),Hepatitis.B.agegroup11,"model7",c(3,11),(2+2*12)/12,10000)
EIR.merge.8.2<-EIR.merge.2009(c(6:8),Hepatitis.B.agegroup12,"model8",c(3,11),(2+7*12)/12,10000)
EIR.merge.9<-EIR.merge.2009(c(6:8),Hepatitis.B.agegroup13,"model9",c(2,9),6,10000)
EIR.merge.10<-EIR.merge.2009(c(6:8),Hepatitis.B.agegroup14,"model10",c(2,9),1,10000)
EIR.merge.northeast.2009<-Prop.pop.northeast.2009[1]*EIR.merge.6.2$EIR+
  Prop.pop.northeast.2009[2]*EIR.merge.7.2$EIR+Prop.pop.northeast.2009[3]*EIR.merge.8.2$EIR+
  Prop.pop.northeast.2009[4]*EIR.merge.9$EIR+Prop.pop.northeast.2009[5]*EIR.merge.10$EIR
EIR.point.northeast.2009<-Prop.pop.northeast.2009[1]*EIR.merge.6.2$Point.EIR[1]+
  Prop.pop.northeast.2009[2]*EIR.merge.7.2$Point.EIR[1]+Prop.pop.northeast.2009[3]*EIR.merge.8.2$Point.EIR[1]+
  Prop.pop.northeast.2009[4]*EIR.merge.9$Point.EIR[1]+Prop.pop.northeast.2009[5]*EIR.merge.10$Point.EIR[1]
EIR.northeast.2009.low<-sort(EIR.merge.northeast.2009)[n.cir*0.025]
EIR.northeast.2009.high<-sort(EIR.merge.northeast.2009)[n.cir*0.975]
EC.merge.northeast.2009<-EIR.merge.6.2$EC+
  EIR.merge.7.2$EC+EIR.merge.8.2$EC+
  EIR.merge.9$EC+EIR.merge.10$EC
EC.point.northeast.2009<-EIR.merge.6.2$Point.EC[1]+
  EIR.merge.7.2$Point.EC[1]+EIR.merge.8.2$Point.EC[1]+
  EIR.merge.9$Point.EC[1]+EIR.merge.10$Point.EC[1]
EC.northeast.2009.low<-sort(EC.merge.northeast.2009)[n.cir*0.025]
EC.northeast.2009.high<-sort(EC.merge.northeast.2009)[n.cir*0.975]

#East
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order>9&
                                                             Hepatitis.B.agegroup10$Order<16),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order>8&
                                                             Hepatitis.B.agegroup11$Order<16),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order>8&
                                                             Hepatitis.B.agegroup12$Order<16),],FUN=sum)[,2])
Age.pop9<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup13[which(Hepatitis.B.agegroup13$Order>8&
                                                             Hepatitis.B.agegroup13$Order<16),],FUN=sum)[,2])
Age.pop10<-mean(aggregate(Population.book~Time,
                          data=Hepatitis.B.agegroup14[which(Hepatitis.B.agegroup14$Order>8&
                                                              Hepatitis.B.agegroup14$Order<16),],FUN=sum)[,2])
Age.pop.2009<-c(Age.pop6,Age.pop7,Age.pop8,Age.pop9,Age.pop10)
Prop.pop.east.2009<-c(Age.pop6/sum(Age.pop.2009),Age.pop7/sum(Age.pop.2009),
                      Age.pop8/sum(Age.pop.2009),Age.pop9/sum(Age.pop.2009),Age.pop10/sum(Age.pop.2009))
EIR.merge.6.2<-EIR.merge.2009(c(10:15),Hepatitis.B.agegroup10,"model6",c(3,11),(2+12)/12,10000)
EIR.merge.7.2<-EIR.merge.2009(c(9:15),Hepatitis.B.agegroup11,"model7",c(3,11),(2+2*12)/12,10000)
EIR.merge.8.2<-EIR.merge.2009(c(9:15),Hepatitis.B.agegroup12,"model8",c(3,11),(2+7*12)/12,10000)
EIR.merge.9<-EIR.merge.2009(c(9:15),Hepatitis.B.agegroup13,"model9",c(2,9),6,10000)
EIR.merge.10<-EIR.merge.2009(c(9:15),Hepatitis.B.agegroup14,"model10",c(2,9),1,10000)
EIR.merge.east.2009<-Prop.pop.east.2009[1]*EIR.merge.6.2$EIR+
  Prop.pop.east.2009[2]*EIR.merge.7.2$EIR+Prop.pop.east.2009[3]*EIR.merge.8.2$EIR+
  Prop.pop.east.2009[4]*EIR.merge.9$EIR+Prop.pop.east.2009[5]*EIR.merge.10$EIR
EIR.point.east.2009<-Prop.pop.east.2009[1]*EIR.merge.6.2$Point.EIR[1]+
  Prop.pop.east.2009[2]*EIR.merge.7.2$Point.EIR[1]+Prop.pop.east.2009[3]*EIR.merge.8.2$Point.EIR[1]+
  Prop.pop.east.2009[4]*EIR.merge.9$Point.EIR[1]+Prop.pop.east.2009[5]*EIR.merge.10$Point.EIR[1]
EIR.east.2009.low<-sort(EIR.merge.east.2009)[n.cir*0.025]
EIR.east.2009.high<-sort(EIR.merge.east.2009)[n.cir*0.975]
EC.merge.east.2009<-EIR.merge.6.2$EC+
  EIR.merge.7.2$EC+EIR.merge.8.2$EC+
  EIR.merge.9$EC+EIR.merge.10$EC
EC.point.east.2009<-EIR.merge.6.2$Point.EC[1]+
  EIR.merge.7.2$Point.EC[1]+EIR.merge.8.2$Point.EC[1]+
  EIR.merge.9$Point.EC[1]+EIR.merge.10$Point.EC[1]
EC.east.2009.low<-sort(EC.merge.east.2009)[n.cir*0.025]
EC.east.2009.high<-sort(EC.merge.east.2009)[n.cir*0.975]

#South-centre
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order>15&
                                                             Hepatitis.B.agegroup10$Order<22),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order>15&
                                                             Hepatitis.B.agegroup11$Order<22),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order>15&
                                                             Hepatitis.B.agegroup12$Order<22),],FUN=sum)[,2])
Age.pop9<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup13[which(Hepatitis.B.agegroup13$Order>15&
                                                             Hepatitis.B.agegroup13$Order<22),],FUN=sum)[,2])
Age.pop10<-mean(aggregate(Population.book~Time,
                          data=Hepatitis.B.agegroup14[which(Hepatitis.B.agegroup14$Order>15&
                                                              Hepatitis.B.agegroup14$Order<22),],FUN=sum)[,2])
Age.pop.2009<-c(Age.pop6,Age.pop7,Age.pop8,Age.pop9,Age.pop10)
Prop.pop.southcentral.2009<-c(Age.pop6/sum(Age.pop.2009),Age.pop7/sum(Age.pop.2009),
                              Age.pop8/sum(Age.pop.2009),Age.pop9/sum(Age.pop.2009),Age.pop10/sum(Age.pop.2009))
EIR.merge.6.2<-EIR.merge.2009(c(16:21),Hepatitis.B.agegroup10,"model6",c(3,11),(2+12)/12,10000)
EIR.merge.7.2<-EIR.merge.2009(c(16:21),Hepatitis.B.agegroup11,"model7",c(3,11),(2+2*12)/12,10000)
EIR.merge.8.2<-EIR.merge.2009(c(16:21),Hepatitis.B.agegroup12,"model8",c(3,11),(2+7*12)/12,10000)
EIR.merge.9<-EIR.merge.2009(c(16:21),Hepatitis.B.agegroup13,"model9",c(2,9),6,10000)
EIR.merge.10<-EIR.merge.2009(c(16:21),Hepatitis.B.agegroup14,"model10",c(2,9),1,10000)
EIR.merge.southcentral.2009<-Prop.pop.southcentral.2009[1]*EIR.merge.6.2$EIR+
  Prop.pop.southcentral.2009[2]*EIR.merge.7.2$EIR+Prop.pop.southcentral.2009[3]*EIR.merge.8.2$EIR+
  Prop.pop.southcentral.2009[4]*EIR.merge.9$EIR+Prop.pop.southcentral.2009[5]*EIR.merge.10$EIR
EIR.point.southcentral.2009<-Prop.pop.southcentral.2009[1]*EIR.merge.6.2$Point.EIR[1]+
  Prop.pop.southcentral.2009[2]*EIR.merge.7.2$Point.EIR[1]+Prop.pop.southcentral.2009[3]*EIR.merge.8.2$Point.EIR[1]+
  Prop.pop.southcentral.2009[4]*EIR.merge.9$Point.EIR[1]+Prop.pop.southcentral.2009[5]*EIR.merge.10$Point.EIR[1]
EIR.southcentral.2009.low<-sort(EIR.merge.southcentral.2009)[n.cir*0.025]
EIR.southcentral.2009.high<-sort(EIR.merge.southcentral.2009)[n.cir*0.975]
EC.merge.southcentral.2009<-EIR.merge.6.2$EC+
  EIR.merge.7.2$EC+EIR.merge.8.2$EC+
  EIR.merge.9$EC+EIR.merge.10$EC
EC.point.southcentral.2009<-EIR.merge.6.2$Point.EC[1]+
  EIR.merge.7.2$Point.EC[1]+EIR.merge.8.2$Point.EC[1]+
  EIR.merge.9$Point.EC[1]+EIR.merge.10$Point.EC[1]
EC.southcentral.2009.low<-sort(EC.merge.southcentral.2009)[n.cir*0.025]
EC.southcentral.2009.high<-sort(EC.merge.southcentral.2009)[n.cir*0.975]

#Southwest
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order>21&
                                                             Hepatitis.B.agegroup10$Order<26),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order>21&
                                                             Hepatitis.B.agegroup11$Order<27),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order>21&
                                                             Hepatitis.B.agegroup12$Order<27),],FUN=sum)[,2])
Age.pop9<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup13[which(Hepatitis.B.agegroup13$Order>21&
                                                             Hepatitis.B.agegroup13$Order<27),],FUN=sum)[,2])
Age.pop10<-mean(aggregate(Population.book~Time,
                          data=Hepatitis.B.agegroup14[which(Hepatitis.B.agegroup14$Order>21&
                                                              Hepatitis.B.agegroup14$Order<27),],FUN=sum)[,2])
Age.pop.2009<-c(Age.pop6,Age.pop7,Age.pop8,Age.pop9,Age.pop10)
Prop.pop.southwest.2009<-c(Age.pop6/sum(Age.pop.2009),Age.pop7/sum(Age.pop.2009),
                           Age.pop8/sum(Age.pop.2009),Age.pop9/sum(Age.pop.2009),Age.pop10/sum(Age.pop.2009))
EIR.merge.6.2<-EIR.merge.2009(c(22:25),Hepatitis.B.agegroup10,"model6",c(3,11),(2+12)/12,10000)
EIR.merge.7.2<-EIR.merge.2009(c(22:26),Hepatitis.B.agegroup11,"model7",c(3,11),(2+2*12)/12,10000)
EIR.merge.8.2<-EIR.merge.2009(c(22:26),Hepatitis.B.agegroup12,"model8",c(3,11),(2+7*12)/12,10000)
EIR.merge.9<-EIR.merge.2009(c(22:26),Hepatitis.B.agegroup13,"model9",c(2,9),6,10000)
EIR.merge.10<-EIR.merge.2009(c(22:26),Hepatitis.B.agegroup14,"model10",c(2,9),1,10000)
EIR.merge.southwest.2009<-Prop.pop.southwest.2009[1]*EIR.merge.6.2$EIR+
  Prop.pop.southwest.2009[2]*EIR.merge.7.2$EIR+Prop.pop.southwest.2009[3]*EIR.merge.8.2$EIR+
  Prop.pop.southwest.2009[4]*EIR.merge.9$EIR+Prop.pop.southwest.2009[5]*EIR.merge.10$EIR
EIR.point.southwest.2009<-Prop.pop.southwest.2009[1]*EIR.merge.6.2$Point.EIR[1]+
  Prop.pop.southwest.2009[2]*EIR.merge.7.2$Point.EIR[1]+Prop.pop.southwest.2009[3]*EIR.merge.8.2$Point.EIR[1]+
  Prop.pop.southwest.2009[4]*EIR.merge.9$Point.EIR[1]+Prop.pop.southwest.2009[5]*EIR.merge.10$Point.EIR[1]
EIR.southwest.2009.low<-sort(EIR.merge.southwest.2009)[n.cir*0.025]
EIR.southwest.2009.high<-sort(EIR.merge.southwest.2009)[n.cir*0.975]
EC.merge.southwest.2009<-EIR.merge.6.2$EC+
  EIR.merge.7.2$EC+EIR.merge.8.2$EC+
  EIR.merge.9$EC+EIR.merge.10$EC
EC.point.southwest.2009<-EIR.merge.6.2$Point.EC[1]+
  EIR.merge.7.2$Point.EC[1]+EIR.merge.8.2$Point.EC[1]+
  EIR.merge.9$Point.EC[1]+EIR.merge.10$Point.EC[1]
EC.southwest.2009.low<-sort(EC.merge.southwest.2009)[n.cir*0.025]
EC.southwest.2009.high<-sort(EC.merge.southwest.2009)[n.cir*0.975]

#Northwest
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order>26),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order>26),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order>26),],FUN=sum)[,2])
Age.pop9<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup13[which(Hepatitis.B.agegroup13$Order>26),],FUN=sum)[,2])
Age.pop10<-mean(aggregate(Population.book~Time,
                          data=Hepatitis.B.agegroup14[which(Hepatitis.B.agegroup14$Order>26),],FUN=sum)[,2])
Age.pop.2009<-c(Age.pop6,Age.pop7,Age.pop8,Age.pop9,Age.pop10)
Prop.pop.northwest.2009<-c(Age.pop6/sum(Age.pop.2009),Age.pop7/sum(Age.pop.2009),
                           Age.pop8/sum(Age.pop.2009),Age.pop9/sum(Age.pop.2009),Age.pop10/sum(Age.pop.2009))
EIR.merge.6.2<-EIR.merge.2009(c(27:31),Hepatitis.B.agegroup10,"model6",c(3,11),(2+12)/12,10000)
EIR.merge.7.2<-EIR.merge.2009(c(27:31),Hepatitis.B.agegroup11,"model7",c(3,11),(2+2*12)/12,10000)
EIR.merge.8.2<-EIR.merge.2009(c(27:31),Hepatitis.B.agegroup12,"model8",c(3,11),(2+7*12)/12,10000)
EIR.merge.9<-EIR.merge.2009(c(27:31),Hepatitis.B.agegroup13,"model9",c(2,9),6,10000)
EIR.merge.10<-EIR.merge.2009(c(27:31),Hepatitis.B.agegroup14,"model10",c(2,9),1,10000)
EIR.merge.northwest.2009<-Prop.pop.northwest.2009[1]*EIR.merge.6.2$EIR+
  Prop.pop.northwest.2009[2]*EIR.merge.7.2$EIR+Prop.pop.northwest.2009[3]*EIR.merge.8.2$EIR+
  Prop.pop.northwest.2009[4]*EIR.merge.9$EIR+Prop.pop.northwest.2009[5]*EIR.merge.10$EIR
EIR.point.northwest.2009<-Prop.pop.northwest.2009[1]*EIR.merge.6.2$Point.EIR[1]+
  Prop.pop.northwest.2009[2]*EIR.merge.7.2$Point.EIR[1]+Prop.pop.northwest.2009[3]*EIR.merge.8.2$Point.EIR[1]+
  Prop.pop.northwest.2009[4]*EIR.merge.9$Point.EIR[1]+Prop.pop.northwest.2009[5]*EIR.merge.10$Point.EIR[1]
EIR.northwest.2009.low<-sort(EIR.merge.northwest.2009)[n.cir*0.025]
EIR.northwest.2009.high<-sort(EIR.merge.northwest.2009)[n.cir*0.975]
EC.merge.northwest.2009<-EIR.merge.6.2$EC+
  EIR.merge.7.2$EC+EIR.merge.8.2$EC+
  EIR.merge.9$EC+EIR.merge.10$EC
EC.point.northwest.2009<-EIR.merge.6.2$Point.EC[1]+
  EIR.merge.7.2$Point.EC[1]+EIR.merge.8.2$Point.EC[1]+
  EIR.merge.9$Point.EC[1]+EIR.merge.10$Point.EC[1]
EC.northwest.2009.low<-sort(EC.merge.northwest.2009)[n.cir*0.025]
EC.northwest.2009.high<-sort(EC.merge.northwest.2009)[n.cir*0.975]

#######################HBsAg subgroup EIR, 2009###################################
#Prevalence of HBsAg category
Low<-c(1,2,3,4,5,6,7,9,10,15,16,25,27,28,31)
Middle<-c(8,11,12,17,18,22,23,24,29,30)
High<-c(13,14,19,20,21,26)
#Low HBsAg prevalence
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order>2&
                                                             Hepatitis.B.agegroup10$Order<8|
                                                             Hepatitis.B.agegroup10$Order==10|
                                                             Hepatitis.B.agegroup10$Order==15|
                                                             Hepatitis.B.agegroup10$Order==16|
                                                             Hepatitis.B.agegroup10$Order==25|
                                                             Hepatitis.B.agegroup10$Order==27|
                                                             Hepatitis.B.agegroup10$Order==28|
                                                             Hepatitis.B.agegroup10$Order==31),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order>0&
                                                             Hepatitis.B.agegroup11$Order<8|
                                                             Hepatitis.B.agegroup11$Order==9|
                                                             Hepatitis.B.agegroup11$Order==10|
                                                             Hepatitis.B.agegroup11$Order==15|
                                                             Hepatitis.B.agegroup11$Order==16|
                                                             Hepatitis.B.agegroup11$Order==25|
                                                             Hepatitis.B.agegroup11$Order==27|
                                                             Hepatitis.B.agegroup11$Order==28|
                                                             Hepatitis.B.agegroup11$Order==31),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order>0&
                                                             Hepatitis.B.agegroup12$Order<8|
                                                             Hepatitis.B.agegroup12$Order==9|
                                                             Hepatitis.B.agegroup12$Order==10|
                                                             Hepatitis.B.agegroup12$Order==15|
                                                             Hepatitis.B.agegroup12$Order==16|
                                                             Hepatitis.B.agegroup12$Order==25|
                                                             Hepatitis.B.agegroup12$Order==27|
                                                             Hepatitis.B.agegroup12$Order==28|
                                                             Hepatitis.B.agegroup12$Order==31),],FUN=sum)[,2])
Age.pop9<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup13[which(Hepatitis.B.agegroup13$Order>0&
                                                             Hepatitis.B.agegroup13$Order<8|
                                                             Hepatitis.B.agegroup13$Order==9|
                                                             Hepatitis.B.agegroup13$Order==10|
                                                             Hepatitis.B.agegroup13$Order==15|
                                                             Hepatitis.B.agegroup13$Order==16|
                                                             Hepatitis.B.agegroup13$Order==25|
                                                             Hepatitis.B.agegroup13$Order==27|
                                                             Hepatitis.B.agegroup13$Order==28|
                                                             Hepatitis.B.agegroup13$Order==31),],FUN=sum)[,2])
Age.pop10<-mean(aggregate(Population.book~Time,
                          data=Hepatitis.B.agegroup14[which(Hepatitis.B.agegroup14$Order>0&
                                                              Hepatitis.B.agegroup14$Order<8|
                                                              Hepatitis.B.agegroup14$Order==9|
                                                              Hepatitis.B.agegroup14$Order==10|
                                                              Hepatitis.B.agegroup14$Order==15|
                                                              Hepatitis.B.agegroup14$Order==16|
                                                              Hepatitis.B.agegroup14$Order==25|
                                                              Hepatitis.B.agegroup14$Order==27|
                                                              Hepatitis.B.agegroup14$Order==28|
                                                              Hepatitis.B.agegroup14$Order==31),],FUN=sum)[,2])
Age.pop.2009<-c(Age.pop6,Age.pop7,Age.pop8,Age.pop9,Age.pop10)
Prop.pop.lowpr.2009<-c(Age.pop6/sum(Age.pop.2009),Age.pop7/sum(Age.pop.2009),
                       Age.pop8/sum(Age.pop.2009),Age.pop9/sum(Age.pop.2009),Age.pop10/sum(Age.pop.2009))
EIR.merge.6.2<-EIR.merge.2009(c(3,4,5,6,7,10,15,16,25,27,28,31),Hepatitis.B.agegroup10,"model6",c(3,11),(2+12)/12,10000)
EIR.merge.7.2<-EIR.merge.2009(c(1,2,3,4,5,6,7,9,10,15,16,25,27,28,31),Hepatitis.B.agegroup11,"model7",c(3,11),(2+2*12)/12,10000)
EIR.merge.8.2<-EIR.merge.2009(c(1,2,3,4,5,6,7,9,10,15,16,25,27,28,31),Hepatitis.B.agegroup12,"model8",c(3,11),(2+7*12)/12,10000)
EIR.merge.9<-EIR.merge.2009(c(1,2,3,4,5,6,7,9,10,15,16,25,27,28,31),Hepatitis.B.agegroup13,"model9",c(2,9),6,10000)
EIR.merge.10<-EIR.merge.2009(c(1,2,3,4,5,6,7,9,10,15,16,25,27,28,31),Hepatitis.B.agegroup14,"model10",c(2,9),1,10000)
EIR.merge.lowpr.2009<-Prop.pop.lowpr.2009[1]*EIR.merge.6.2$EIR+
  Prop.pop.lowpr.2009[2]*EIR.merge.7.2$EIR+Prop.pop.lowpr.2009[3]*EIR.merge.8.2$EIR+
  Prop.pop.lowpr.2009[4]*EIR.merge.9$EIR+Prop.pop.lowpr.2009[5]*EIR.merge.10$EIR
EIR.point.lowpr.2009<-Prop.pop.lowpr.2009[1]*EIR.merge.6.2$Point.EIR[1]+
  Prop.pop.lowpr.2009[2]*EIR.merge.7.2$Point.EIR[1]+Prop.pop.lowpr.2009[3]*EIR.merge.8.2$Point.EIR[1]+
  Prop.pop.lowpr.2009[4]*EIR.merge.9$Point.EIR[1]+Prop.pop.lowpr.2009[5]*EIR.merge.10$Point.EIR[1]
EIR.lowpr.2009.low<-sort(EIR.merge.lowpr.2009)[n.cir*0.025]
EIR.lowpr.2009.high<-sort(EIR.merge.lowpr.2009)[n.cir*0.975]
EC.merge.lowpr.2009<-EIR.merge.6.2$EC+
  EIR.merge.7.2$EC+EIR.merge.8.2$EC+
  EIR.merge.9$EC+EIR.merge.10$EC
EC.point.lowpr.2009<-EIR.merge.6.2$Point.EC[1]+
  EIR.merge.7.2$Point.EC[1]+EIR.merge.8.2$Point.EC[1]+
  EIR.merge.9$Point.EC[1]+EIR.merge.10$Point.EC[1]
EC.lowpr.2009.low<-sort(EC.merge.lowpr.2009)[n.cir*0.025]
EC.lowpr.2009.high<-sort(EC.merge.lowpr.2009)[n.cir*0.975]

#Middle HBsAg prevalence
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order==8|
                                                             Hepatitis.B.agegroup10$Order==11|
                                                             Hepatitis.B.agegroup10$Order==12|
                                                             Hepatitis.B.agegroup10$Order==17|
                                                             Hepatitis.B.agegroup10$Order==18|
                                                             Hepatitis.B.agegroup10$Order==22|
                                                             Hepatitis.B.agegroup10$Order==23|
                                                             Hepatitis.B.agegroup10$Order==24|
                                                             Hepatitis.B.agegroup10$Order==29|
                                                             Hepatitis.B.agegroup10$Order==30),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order==8|
                                                             Hepatitis.B.agegroup11$Order==11|
                                                             Hepatitis.B.agegroup11$Order==12|
                                                             Hepatitis.B.agegroup11$Order==17|
                                                             Hepatitis.B.agegroup11$Order==18|
                                                             Hepatitis.B.agegroup11$Order==22|
                                                             Hepatitis.B.agegroup11$Order==23|
                                                             Hepatitis.B.agegroup11$Order==24|
                                                             Hepatitis.B.agegroup11$Order==29|
                                                             Hepatitis.B.agegroup11$Order==30),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order==8|
                                                             Hepatitis.B.agegroup12$Order==11|
                                                             Hepatitis.B.agegroup12$Order==12|
                                                             Hepatitis.B.agegroup12$Order==17|
                                                             Hepatitis.B.agegroup12$Order==18|
                                                             Hepatitis.B.agegroup12$Order==22|
                                                             Hepatitis.B.agegroup12$Order==23|
                                                             Hepatitis.B.agegroup12$Order==24|
                                                             Hepatitis.B.agegroup12$Order==29|
                                                             Hepatitis.B.agegroup12$Order==30),],FUN=sum)[,2])
Age.pop9<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup13[which(Hepatitis.B.agegroup13$Order==8|
                                                             Hepatitis.B.agegroup13$Order==11|
                                                             Hepatitis.B.agegroup13$Order==12|
                                                             Hepatitis.B.agegroup13$Order==17|
                                                             Hepatitis.B.agegroup13$Order==18|
                                                             Hepatitis.B.agegroup13$Order==22|
                                                             Hepatitis.B.agegroup13$Order==23|
                                                             Hepatitis.B.agegroup13$Order==24|
                                                             Hepatitis.B.agegroup13$Order==29|
                                                             Hepatitis.B.agegroup13$Order==30),],FUN=sum)[,2])
Age.pop10<-mean(aggregate(Population.book~Time,
                          data=Hepatitis.B.agegroup14[which(Hepatitis.B.agegroup14$Order==8|
                                                              Hepatitis.B.agegroup14$Order==11|
                                                              Hepatitis.B.agegroup14$Order==12|
                                                              Hepatitis.B.agegroup14$Order==17|
                                                              Hepatitis.B.agegroup14$Order==18|
                                                              Hepatitis.B.agegroup14$Order==22|
                                                              Hepatitis.B.agegroup14$Order==23|
                                                              Hepatitis.B.agegroup14$Order==24|
                                                              Hepatitis.B.agegroup14$Order==29|
                                                              Hepatitis.B.agegroup14$Order==30),],FUN=sum)[,2])
Age.pop.2009<-c(Age.pop6,Age.pop7,Age.pop8,Age.pop9,Age.pop10)
Prop.pop.middlepr.2009<-c(Age.pop6/sum(Age.pop.2009),Age.pop7/sum(Age.pop.2009),
                          Age.pop8/sum(Age.pop.2009),Age.pop9/sum(Age.pop.2009),Age.pop10/sum(Age.pop.2009))
EIR.merge.6.2<-EIR.merge.2009(c(8,11,12,17,18,22,23,24,29,30),Hepatitis.B.agegroup10,"model6",c(3,11),(2+12)/12,10000)
EIR.merge.7.2<-EIR.merge.2009(c(8,11,12,17,18,22,23,24,29,30),Hepatitis.B.agegroup11,"model7",c(3,11),(2+2*12)/12,10000)
EIR.merge.8.2<-EIR.merge.2009(c(8,11,12,17,18,22,23,24,29,30),Hepatitis.B.agegroup12,"model8",c(3,11),(2+7*12)/12,10000)
EIR.merge.9<-EIR.merge.2009(c(8,11,12,17,18,22,23,24,29,30),Hepatitis.B.agegroup13,"model9",c(2,9),6,10000)
EIR.merge.10<-EIR.merge.2009(c(8,11,12,17,18,22,23,24,29,30),Hepatitis.B.agegroup14,"model10",c(2,9),1,10000)
EIR.merge.middlepr.2009<-Prop.pop.middlepr.2009[1]*EIR.merge.6.2$EIR+
  Prop.pop.middlepr.2009[2]*EIR.merge.7.2$EIR+Prop.pop.middlepr.2009[3]*EIR.merge.8.2$EIR+
  Prop.pop.middlepr.2009[4]*EIR.merge.9$EIR+Prop.pop.middlepr.2009[5]*EIR.merge.10$EIR
EIR.point.middlepr.2009<-Prop.pop.middlepr.2009[1]*EIR.merge.6.2$Point.EIR[1]+
  Prop.pop.middlepr.2009[2]*EIR.merge.7.2$Point.EIR[1]+Prop.pop.middlepr.2009[3]*EIR.merge.8.2$Point.EIR[1]+
  Prop.pop.middlepr.2009[4]*EIR.merge.9$Point.EIR[1]+Prop.pop.middlepr.2009[5]*EIR.merge.10$Point.EIR[1]
EIR.middlepr.2009.low<-sort(EIR.merge.middlepr.2009)[n.cir*0.025]
EIR.middlepr.2009.high<-sort(EIR.merge.middlepr.2009)[n.cir*0.975]
EC.merge.middlepr.2009<-EIR.merge.6.2$EC+
  EIR.merge.7.2$EC+EIR.merge.8.2$EC+
  EIR.merge.9$EC+EIR.merge.10$EC
EC.point.middlepr.2009<-EIR.merge.6.2$Point.EC[1]+
  EIR.merge.7.2$Point.EC[1]+EIR.merge.8.2$Point.EC[1]+
  EIR.merge.9$Point.EC[1]+EIR.merge.10$Point.EC[1]
EC.middlepr.2009.low<-sort(EC.merge.middlepr.2009)[n.cir*0.025]
EC.middlepr.2009.high<-sort(EC.merge.middlepr.2009)[n.cir*0.975]

#High HBsAg prevalence
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order==13|
                                                             Hepatitis.B.agegroup10$Order==14|
                                                             Hepatitis.B.agegroup10$Order==19|
                                                             Hepatitis.B.agegroup10$Order==20|
                                                             Hepatitis.B.agegroup10$Order==21),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order==13|
                                                             Hepatitis.B.agegroup11$Order==14|
                                                             Hepatitis.B.agegroup11$Order==19|
                                                             Hepatitis.B.agegroup11$Order==20|
                                                             Hepatitis.B.agegroup11$Order==21|
                                                             Hepatitis.B.agegroup11$Order==26),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order==13|
                                                             Hepatitis.B.agegroup12$Order==14|
                                                             Hepatitis.B.agegroup12$Order==19|
                                                             Hepatitis.B.agegroup12$Order==20|
                                                             Hepatitis.B.agegroup12$Order==21|
                                                             Hepatitis.B.agegroup12$Order==26),],FUN=sum)[,2])
Age.pop9<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup13[which(Hepatitis.B.agegroup13$Order==13|
                                                             Hepatitis.B.agegroup13$Order==14|
                                                             Hepatitis.B.agegroup13$Order==19|
                                                             Hepatitis.B.agegroup13$Order==20|
                                                             Hepatitis.B.agegroup13$Order==21|
                                                             Hepatitis.B.agegroup13$Order==26),],FUN=sum)[,2])
Age.pop10<-mean(aggregate(Population.book~Time,
                          data=Hepatitis.B.agegroup14[which(Hepatitis.B.agegroup14$Order==13|
                                                              Hepatitis.B.agegroup14$Order==14|
                                                              Hepatitis.B.agegroup14$Order==19|
                                                              Hepatitis.B.agegroup14$Order==20|
                                                              Hepatitis.B.agegroup14$Order==21|
                                                              Hepatitis.B.agegroup14$Order==26),],FUN=sum)[,2])
Age.pop.2009<-c(Age.pop6,Age.pop7,Age.pop8,Age.pop9,Age.pop10)
Prop.pop.highpr.2009<-c(Age.pop6/sum(Age.pop.2009),Age.pop7/sum(Age.pop.2009),
                        Age.pop8/sum(Age.pop.2009),Age.pop9/sum(Age.pop.2009),Age.pop10/sum(Age.pop.2009))
EIR.merge.6.2<-EIR.merge.2009(c(13,14,19,20,21),Hepatitis.B.agegroup10,"model6",c(3,11),(2+12)/12,10000)
EIR.merge.7.2<-EIR.merge.2009(c(13,14,19,20,21,26),Hepatitis.B.agegroup11,"model7",c(3,11),(2+2*12)/14,10000)
EIR.merge.8.2<-EIR.merge.2009(c(13,14,19,20,21,26),Hepatitis.B.agegroup12,"model8",c(3,11),(2+7*12)/14,10000)
EIR.merge.9<-EIR.merge.2009(c(13,14,19,20,21,26),Hepatitis.B.agegroup13,"model9",c(2,9),6,10000)
EIR.merge.10<-EIR.merge.2009(c(13,14,19,20,21,26),Hepatitis.B.agegroup14,"model10",c(2,9),1,10000)
EIR.merge.highpr.2009<-Prop.pop.highpr.2009[1]*EIR.merge.6.2$EIR+
  Prop.pop.highpr.2009[2]*EIR.merge.7.2$EIR+Prop.pop.highpr.2009[3]*EIR.merge.8.2$EIR+
  Prop.pop.highpr.2009[4]*EIR.merge.9$EIR+Prop.pop.highpr.2009[5]*EIR.merge.10$EIR
EIR.point.highpr.2009<-Prop.pop.highpr.2009[1]*EIR.merge.6.2$Point.EIR[1]+
  Prop.pop.highpr.2009[2]*EIR.merge.7.2$Point.EIR[1]+Prop.pop.highpr.2009[3]*EIR.merge.8.2$Point.EIR[1]+
  Prop.pop.highpr.2009[4]*EIR.merge.9$Point.EIR[1]+Prop.pop.highpr.2009[5]*EIR.merge.10$Point.EIR[1]
EIR.highpr.2009.low<-sort(EIR.merge.highpr.2009)[n.cir*0.025]
EIR.highpr.2009.high<-sort(EIR.merge.highpr.2009)[n.cir*0.975]
EC.merge.highpr.2009<-EIR.merge.6.2$EC+
  EIR.merge.7.2$EC+EIR.merge.8.2$EC+
  EIR.merge.9$EC+EIR.merge.10$EC
EC.point.highpr.2009<-EIR.merge.6.2$Point.EC[1]+
  EIR.merge.7.2$Point.EC[1]+EIR.merge.8.2$Point.EC[1]+
  EIR.merge.9$Point.EC[1]+EIR.merge.10$Point.EC[1]
EC.highpr.2009.low<-sort(EC.merge.highpr.2009)[n.cir*0.025]
EC.highpr.2009.high<-sort(EC.merge.highpr.2009)[n.cir*0.975]

#######################Urbanisation rate subgroup EIR, 2009###################################
#Urbanisation rate category
LowUR<-c(3,12,14,16,18,20,23:31)
HighUR<-c(1,2,4,5,6,7,8,9,10,11,13,15,17,19,21,22)
#Low urbanisation
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order==3|
                                                             Hepatitis.B.agegroup10$Order==12|
                                                             Hepatitis.B.agegroup10$Order==14|
                                                             Hepatitis.B.agegroup10$Order==16|
                                                             Hepatitis.B.agegroup10$Order==18|
                                                             Hepatitis.B.agegroup10$Order==20|
                                                             Hepatitis.B.agegroup10$Order>22&
                                                             Hepatitis.B.agegroup10$Order<26|
                                                             Hepatitis.B.agegroup10$Order>26),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order==3|
                                                             Hepatitis.B.agegroup11$Order==12|
                                                             Hepatitis.B.agegroup11$Order==14|
                                                             Hepatitis.B.agegroup11$Order==16|
                                                             Hepatitis.B.agegroup11$Order==18|
                                                             Hepatitis.B.agegroup11$Order==20|
                                                             Hepatitis.B.agegroup11$Order>22),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order==3|
                                                             Hepatitis.B.agegroup12$Order==12|
                                                             Hepatitis.B.agegroup12$Order==14|
                                                             Hepatitis.B.agegroup12$Order==16|
                                                             Hepatitis.B.agegroup12$Order==18|
                                                             Hepatitis.B.agegroup12$Order==20|
                                                             Hepatitis.B.agegroup12$Order>22),],FUN=sum)[,2])
Age.pop9<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup13[which(Hepatitis.B.agegroup13$Order==3|
                                                             Hepatitis.B.agegroup13$Order==12|
                                                             Hepatitis.B.agegroup13$Order==14|
                                                             Hepatitis.B.agegroup13$Order==16|
                                                             Hepatitis.B.agegroup13$Order==18|
                                                             Hepatitis.B.agegroup13$Order==20|
                                                             Hepatitis.B.agegroup13$Order>22),],FUN=sum)[,2])
Age.pop10<-mean(aggregate(Population.book~Time,
                          data=Hepatitis.B.agegroup14[which(Hepatitis.B.agegroup14$Order==3|
                                                              Hepatitis.B.agegroup14$Order==12|
                                                              Hepatitis.B.agegroup14$Order==14|
                                                              Hepatitis.B.agegroup14$Order==16|
                                                              Hepatitis.B.agegroup14$Order==18|
                                                              Hepatitis.B.agegroup14$Order==20|
                                                              Hepatitis.B.agegroup14$Order>22),],FUN=sum)[,2])
Age.pop.2009<-c(Age.pop6,Age.pop7,Age.pop8,Age.pop9,Age.pop10)
Prop.pop.lowur.2009<-c(Age.pop6/sum(Age.pop.2009),Age.pop7/sum(Age.pop.2009),
                       Age.pop8/sum(Age.pop.2009),Age.pop9/sum(Age.pop.2009),Age.pop10/sum(Age.pop.2009))
EIR.merge.6.2<-EIR.merge.2009(c(3,12,14,16,18,20,23:25,27:31),Hepatitis.B.agegroup10,"model6",c(3,11),(2+12)/12,10000)
EIR.merge.7.2<-EIR.merge.2009(c(3,12,14,16,18,20,23:31),Hepatitis.B.agegroup11,"model7",c(3,11),(2+2*12)/12,10000)
EIR.merge.8.2<-EIR.merge.2009(c(3,12,14,16,18,20,23:31),Hepatitis.B.agegroup12,"model8",c(3,11),(2+7*12)/12,10000)
EIR.merge.9<-EIR.merge.2009(c(3,12,14,16,18,20,23:31),Hepatitis.B.agegroup13,"model9",c(2,9),6,10000)
EIR.merge.10<-EIR.merge.2009(c(3,12,14,16,18,20,23:31),Hepatitis.B.agegroup14,"model10",c(2,9),1,10000)
EIR.merge.lowur.2009<-Prop.pop.lowur.2009[1]*EIR.merge.6.2$EIR+
  Prop.pop.lowur.2009[2]*EIR.merge.7.2$EIR+Prop.pop.lowur.2009[3]*EIR.merge.8.2$EIR+
  Prop.pop.lowur.2009[4]*EIR.merge.9$EIR+Prop.pop.lowur.2009[5]*EIR.merge.10$EIR
EIR.point.lowur.2009<-Prop.pop.lowur.2009[1]*EIR.merge.6.2$Point.EIR[1]+
  Prop.pop.lowur.2009[2]*EIR.merge.7.2$Point.EIR[1]+Prop.pop.lowur.2009[3]*EIR.merge.8.2$Point.EIR[1]+
  Prop.pop.lowur.2009[4]*EIR.merge.9$Point.EIR[1]+Prop.pop.lowur.2009[5]*EIR.merge.10$Point.EIR[1]
EIR.lowur.2009.low<-sort(EIR.merge.lowur.2009)[n.cir*0.025]
EIR.lowur.2009.high<-sort(EIR.merge.lowur.2009)[n.cir*0.975]
EC.merge.lowur.2009<-EIR.merge.6.2$EC+
  EIR.merge.7.2$EC+EIR.merge.8.2$EC+
  EIR.merge.9$EC+EIR.merge.10$EC
EC.point.lowur.2009<-EIR.merge.6.2$Point.EC[1]+
  EIR.merge.7.2$Point.EC[1]+EIR.merge.8.2$Point.EC[1]+
  EIR.merge.9$Point.EC[1]+EIR.merge.10$Point.EC[1]
EC.lowur.2009.low<-sort(EC.merge.lowur.2009)[n.cir*0.025]
EC.lowur.2009.high<-sort(EC.merge.lowur.2009)[n.cir*0.975]

#High urbanisation rate
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order>3&
                                                             Hepatitis.B.agegroup10$Order<9|
                                                             Hepatitis.B.agegroup10$Order==10|
                                                             Hepatitis.B.agegroup10$Order==11|
                                                             Hepatitis.B.agegroup10$Order==13|
                                                             Hepatitis.B.agegroup10$Order==15|
                                                             Hepatitis.B.agegroup10$Order==17|
                                                             Hepatitis.B.agegroup10$Order==19|
                                                             Hepatitis.B.agegroup10$Order==21|
                                                             Hepatitis.B.agegroup10$Order==22),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order<3|
                                                             Hepatitis.B.agegroup11$Order>3&
                                                             Hepatitis.B.agegroup11$Order<12|
                                                             Hepatitis.B.agegroup11$Order==13|
                                                             Hepatitis.B.agegroup11$Order==15|
                                                             Hepatitis.B.agegroup11$Order==17|
                                                             Hepatitis.B.agegroup11$Order==19|
                                                             Hepatitis.B.agegroup11$Order==21|
                                                             Hepatitis.B.agegroup11$Order==22),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order<3|
                                                             Hepatitis.B.agegroup12$Order>3&
                                                             Hepatitis.B.agegroup12$Order<12|
                                                             Hepatitis.B.agegroup12$Order==13|
                                                             Hepatitis.B.agegroup12$Order==15|
                                                             Hepatitis.B.agegroup12$Order==17|
                                                             Hepatitis.B.agegroup12$Order==19|
                                                             Hepatitis.B.agegroup12$Order==21|
                                                             Hepatitis.B.agegroup12$Order==22),],FUN=sum)[,2])
Age.pop9<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup13[which(Hepatitis.B.agegroup13$Order<3|
                                                             Hepatitis.B.agegroup13$Order>3&
                                                             Hepatitis.B.agegroup13$Order<12|
                                                             Hepatitis.B.agegroup13$Order==13|
                                                             Hepatitis.B.agegroup13$Order==15|
                                                             Hepatitis.B.agegroup13$Order==17|
                                                             Hepatitis.B.agegroup13$Order==19|
                                                             Hepatitis.B.agegroup13$Order==21|
                                                             Hepatitis.B.agegroup13$Order==22),],FUN=sum)[,2])
Age.pop10<-mean(aggregate(Population.book~Time,
                          data=Hepatitis.B.agegroup14[which(Hepatitis.B.agegroup14$Order<3|
                                                              Hepatitis.B.agegroup14$Order>3&
                                                              Hepatitis.B.agegroup14$Order<12|
                                                              Hepatitis.B.agegroup14$Order==13|
                                                              Hepatitis.B.agegroup14$Order==15|
                                                              Hepatitis.B.agegroup14$Order==17|
                                                              Hepatitis.B.agegroup14$Order==19|
                                                              Hepatitis.B.agegroup14$Order==21|
                                                              Hepatitis.B.agegroup14$Order==22),],FUN=sum)[,2])
Age.pop.2009<-c(Age.pop6,Age.pop7,Age.pop8,Age.pop9,Age.pop10)
Prop.pop.highur.2009<-c(Age.pop6/sum(Age.pop.2009),Age.pop7/sum(Age.pop.2009),
                        Age.pop8/sum(Age.pop.2009),Age.pop9/sum(Age.pop.2009),Age.pop10/sum(Age.pop.2009))
EIR.merge.6.2<-EIR.merge.2009(c(4,5,6,7,8,10,11,13,15,17,19,21,22),Hepatitis.B.agegroup10,"model6",c(3,11),(2+12)/12,10000)
EIR.merge.7.2<-EIR.merge.2009(c(1,2,4,5,6,7,8,9,10,11,13,15,17,19,21,22),Hepatitis.B.agegroup11,"model7",c(3,11),(2+2*12)/12,10000)
EIR.merge.8.2<-EIR.merge.2009(c(1,2,4,5,6,7,8,9,10,11,13,15,17,19,21,22),Hepatitis.B.agegroup12,"model8",c(3,11),(2+7*12)/12,10000)
EIR.merge.9<-EIR.merge.2009(c(1,2,4,5,6,7,8,9,10,11,13,15,17,19,21,22),Hepatitis.B.agegroup13,"model9",c(2,9),6,10000)
EIR.merge.10<-EIR.merge.2009(c(1,2,4,5,6,7,8,9,10,11,13,15,17,19,21,22),Hepatitis.B.agegroup14,"model10",c(2,9),1,10000)
EIR.merge.highur.2009<-Prop.pop.highur.2009[1]*EIR.merge.6.2$EIR+
  Prop.pop.highur.2009[2]*EIR.merge.7.2$EIR+Prop.pop.highur.2009[3]*EIR.merge.8.2$EIR+
  Prop.pop.highur.2009[4]*EIR.merge.9$EIR+Prop.pop.highur.2009[5]*EIR.merge.10$EIR
EIR.point.highur.2009<-Prop.pop.highur.2009[1]*EIR.merge.6.2$Point.EIR[1]+
  Prop.pop.highur.2009[2]*EIR.merge.7.2$Point.EIR[1]+Prop.pop.highur.2009[3]*EIR.merge.8.2$Point.EIR[1]+
  Prop.pop.highur.2009[4]*EIR.merge.9$Point.EIR[1]+Prop.pop.highur.2009[5]*EIR.merge.10$Point.EIR[1]
EIR.highur.2009.low<-sort(EIR.merge.highur.2009)[n.cir*0.025]
EIR.highur.2009.high<-sort(EIR.merge.highur.2009)[n.cir*0.975]
EC.merge.highur.2009<-EIR.merge.6.2$EC+
  EIR.merge.7.2$EC+EIR.merge.8.2$EC+
  EIR.merge.9$EC+EIR.merge.10$EC
EC.point.highur.2009<-EIR.merge.6.2$Point.EC[1]+
  EIR.merge.7.2$Point.EC[1]+EIR.merge.8.2$Point.EC[1]+
  EIR.merge.9$Point.EC[1]+EIR.merge.10$Point.EC[1]
EC.highur.2009.low<-sort(EC.merge.highur.2009)[n.cir*0.025]
EC.highur.2009.high<-sort(EC.merge.highur.2009)[n.cir*0.975]

#######################Hospitalisation rate subgroup EIR, 2009###################################
#Hospitalisation bed density
LowHB<-c(3,11:14,16,18:22,24:26,28)
HighHB<-c(1,2,4,5,6,7,8,9,10,15,17,23,27,29,30,31)
#Low hospitalisation bed density
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order==3|
                                                             Hepatitis.B.agegroup10$Order>10&
                                                             Hepatitis.B.agegroup10$Order<15|
                                                             Hepatitis.B.agegroup10$Order==16|
                                                             Hepatitis.B.agegroup10$Order>17&
                                                             Hepatitis.B.agegroup10$Order<23|
                                                             Hepatitis.B.agegroup10$Order==24|
                                                             Hepatitis.B.agegroup10$Order==25|
                                                             Hepatitis.B.agegroup10$Order==28),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order==3|
                                                             Hepatitis.B.agegroup11$Order>10&
                                                             Hepatitis.B.agegroup11$Order<15|
                                                             Hepatitis.B.agegroup11$Order==16|
                                                             Hepatitis.B.agegroup11$Order>17&
                                                             Hepatitis.B.agegroup11$Order<23|
                                                             Hepatitis.B.agegroup11$Order==24|
                                                             Hepatitis.B.agegroup11$Order==25|
                                                             Hepatitis.B.agegroup11$Order==26|
                                                             Hepatitis.B.agegroup11$Order==28),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order==3|
                                                             Hepatitis.B.agegroup12$Order>10&
                                                             Hepatitis.B.agegroup12$Order<15|
                                                             Hepatitis.B.agegroup12$Order==16|
                                                             Hepatitis.B.agegroup12$Order>17&
                                                             Hepatitis.B.agegroup12$Order<23|
                                                             Hepatitis.B.agegroup12$Order==24|
                                                             Hepatitis.B.agegroup12$Order==25|
                                                             Hepatitis.B.agegroup12$Order==26|
                                                             Hepatitis.B.agegroup12$Order==28),],FUN=sum)[,2])
Age.pop9<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup13[which(Hepatitis.B.agegroup13$Order==3|
                                                             Hepatitis.B.agegroup13$Order>10&
                                                             Hepatitis.B.agegroup13$Order<15|
                                                             Hepatitis.B.agegroup13$Order==16|
                                                             Hepatitis.B.agegroup13$Order>17&
                                                             Hepatitis.B.agegroup13$Order<23|
                                                             Hepatitis.B.agegroup13$Order==24|
                                                             Hepatitis.B.agegroup13$Order==25|
                                                             Hepatitis.B.agegroup13$Order==26|
                                                             Hepatitis.B.agegroup13$Order==28),],FUN=sum)[,2])
Age.pop10<-mean(aggregate(Population.book~Time,
                          data=Hepatitis.B.agegroup14[which(Hepatitis.B.agegroup14$Order==3|
                                                              Hepatitis.B.agegroup14$Order>10&
                                                              Hepatitis.B.agegroup14$Order<15|
                                                              Hepatitis.B.agegroup14$Order==16|
                                                              Hepatitis.B.agegroup14$Order>17&
                                                              Hepatitis.B.agegroup14$Order<23|
                                                              Hepatitis.B.agegroup14$Order==24|
                                                              Hepatitis.B.agegroup14$Order==25|
                                                              Hepatitis.B.agegroup14$Order==26|
                                                              Hepatitis.B.agegroup14$Order==28),],FUN=sum)[,2])
Age.pop.2009<-c(Age.pop6,Age.pop7,Age.pop8,Age.pop9,Age.pop10)
Prop.pop.lowhb.2009<-c(Age.pop6/sum(Age.pop.2009),Age.pop7/sum(Age.pop.2009),
                       Age.pop8/sum(Age.pop.2009),Age.pop9/sum(Age.pop.2009),Age.pop10/sum(Age.pop.2009))
EIR.merge.6.2<-EIR.merge.2009(c(3,11:14,16,18:22,24:25,28),Hepatitis.B.agegroup10,"model6",c(3,11),(2+12)/12,10000)
EIR.merge.7.2<-EIR.merge.2009(c(3,11:14,16,18:22,24:26,28),Hepatitis.B.agegroup11,"model7",c(3,11),(2+2*12)/12,10000)
EIR.merge.8.2<-EIR.merge.2009(c(3,11:14,16,18:22,24:26,28),Hepatitis.B.agegroup12,"model8",c(3,11),(2+7*12)/12,10000)
EIR.merge.9<-EIR.merge.2009(c(3,11:14,16,18:22,24:26,28),Hepatitis.B.agegroup13,"model9",c(2,9),6,10000)
EIR.merge.10<-EIR.merge.2009(c(3,11:14,16,18:22,24:26,28),Hepatitis.B.agegroup14,"model10",c(2,9),1,10000)
EIR.merge.lowhb.2009<-Prop.pop.lowhb.2009[1]*EIR.merge.6.2$EIR+
  Prop.pop.lowhb.2009[2]*EIR.merge.7.2$EIR+Prop.pop.lowhb.2009[3]*EIR.merge.8.2$EIR+
  Prop.pop.lowhb.2009[4]*EIR.merge.9$EIR+Prop.pop.lowhb.2009[5]*EIR.merge.10$EIR
EIR.point.lowhb.2009<-Prop.pop.lowhb.2009[1]*EIR.merge.6.2$Point.EIR[1]+
  Prop.pop.lowhb.2009[2]*EIR.merge.7.2$Point.EIR[1]+Prop.pop.lowhb.2009[3]*EIR.merge.8.2$Point.EIR[1]+
  Prop.pop.lowhb.2009[4]*EIR.merge.9$Point.EIR[1]+Prop.pop.lowhb.2009[5]*EIR.merge.10$Point.EIR[1]
EIR.lowhb.2009.low<-sort(EIR.merge.lowhb.2009)[n.cir*0.025]
EIR.lowhb.2009.high<-sort(EIR.merge.lowhb.2009)[n.cir*0.975]
EC.merge.lowhb.2009<-EIR.merge.6.2$EC+
  EIR.merge.7.2$EC+EIR.merge.8.2$EC+
  EIR.merge.9$EC+EIR.merge.10$EC
EC.point.lowhb.2009<-EIR.merge.6.2$Point.EC[1]+
  EIR.merge.7.2$Point.EC[1]+EIR.merge.8.2$Point.EC[1]+
  EIR.merge.9$Point.EC[1]+EIR.merge.10$Point.EC[1]
EC.lowhb.2009.low<-sort(EC.merge.lowhb.2009)[n.cir*0.025]
EC.lowhb.2009.high<-sort(EC.merge.lowhb.2009)[n.cir*0.975]

#High hospitalisation bed density
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order>3&
                                                             Hepatitis.B.agegroup10$Order<9|
                                                             Hepatitis.B.agegroup10$Order==10|
                                                             Hepatitis.B.agegroup10$Order==15|
                                                             Hepatitis.B.agegroup10$Order==17|
                                                             Hepatitis.B.agegroup10$Order==23|
                                                             Hepatitis.B.agegroup10$Order==27|
                                                             Hepatitis.B.agegroup10$Order==29|
                                                             Hepatitis.B.agegroup10$Order==30|
                                                             Hepatitis.B.agegroup10$Order==31),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup10$Order==1|
                                                             Hepatitis.B.agegroup10$Order==2|
                                                             Hepatitis.B.agegroup10$Order>3&
                                                             Hepatitis.B.agegroup10$Order<10|
                                                             Hepatitis.B.agegroup10$Order==15|
                                                             Hepatitis.B.agegroup10$Order==17|
                                                             Hepatitis.B.agegroup10$Order==23|
                                                             Hepatitis.B.agegroup10$Order==27|
                                                             Hepatitis.B.agegroup10$Order==29|
                                                             Hepatitis.B.agegroup10$Order==30|
                                                             Hepatitis.B.agegroup10$Order==31),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup10$Order==1|
                                                             Hepatitis.B.agegroup10$Order==2|
                                                             Hepatitis.B.agegroup10$Order>3&
                                                             Hepatitis.B.agegroup10$Order<10|
                                                             Hepatitis.B.agegroup10$Order==15|
                                                             Hepatitis.B.agegroup10$Order==17|
                                                             Hepatitis.B.agegroup10$Order==23|
                                                             Hepatitis.B.agegroup10$Order==27|
                                                             Hepatitis.B.agegroup10$Order==29|
                                                             Hepatitis.B.agegroup10$Order==30|
                                                             Hepatitis.B.agegroup10$Order==31),],FUN=sum)[,2])
Age.pop9<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup13[which(Hepatitis.B.agegroup10$Order==1|
                                                             Hepatitis.B.agegroup10$Order==2|
                                                             Hepatitis.B.agegroup10$Order>3&
                                                             Hepatitis.B.agegroup10$Order<10|
                                                             Hepatitis.B.agegroup10$Order==15|
                                                             Hepatitis.B.agegroup10$Order==17|
                                                             Hepatitis.B.agegroup10$Order==23|
                                                             Hepatitis.B.agegroup10$Order==27|
                                                             Hepatitis.B.agegroup10$Order==29|
                                                             Hepatitis.B.agegroup10$Order==30|
                                                             Hepatitis.B.agegroup10$Order==31),],FUN=sum)[,2])
Age.pop10<-mean(aggregate(Population.book~Time,
                          data=Hepatitis.B.agegroup14[which(Hepatitis.B.agegroup10$Order==1|
                                                              Hepatitis.B.agegroup10$Order==2|
                                                              Hepatitis.B.agegroup10$Order>3&
                                                              Hepatitis.B.agegroup10$Order<10|
                                                              Hepatitis.B.agegroup10$Order==15|
                                                              Hepatitis.B.agegroup10$Order==17|
                                                              Hepatitis.B.agegroup10$Order==23|
                                                              Hepatitis.B.agegroup10$Order==27|
                                                              Hepatitis.B.agegroup10$Order==29|
                                                              Hepatitis.B.agegroup10$Order==30|
                                                              Hepatitis.B.agegroup10$Order==31),],FUN=sum)[,2])
Age.pop.2009<-c(Age.pop6,Age.pop7,Age.pop8,Age.pop9,Age.pop10)
Prop.pop.highhb.2009<-c(Age.pop6/sum(Age.pop.2009),Age.pop7/sum(Age.pop.2009),
                        Age.pop8/sum(Age.pop.2009),Age.pop9/sum(Age.pop.2009),Age.pop10/sum(Age.pop.2009))
EIR.merge.6.2<-EIR.merge.2009(c(4,5,6,7,8,10,15,17,23,27,29,30,31),Hepatitis.B.agegroup10,"model6",c(3,11),(2+12)/12,10000)
EIR.merge.7.2<-EIR.merge.2009(c(1,2,4,5,6,7,8,9,10,15,17,23,27,29,30,31),Hepatitis.B.agegroup11,"model7",c(3,11),(2+2*12)/12,10000)
EIR.merge.8.2<-EIR.merge.2009(c(1,2,4,5,6,7,8,9,10,15,17,23,27,29,30,31),Hepatitis.B.agegroup12,"model8",c(3,11),(2+7*12)/12,10000)
EIR.merge.9<-EIR.merge.2009(c(1,2,4,5,6,7,8,9,10,15,17,23,27,29,30,31),Hepatitis.B.agegroup13,"model9",c(2,9),6,10000)
EIR.merge.10<-EIR.merge.2009(c(1,2,4,5,6,7,8,9,10,15,17,23,27,29,30,31),Hepatitis.B.agegroup14,"model10",c(2,9),1,10000)
EIR.merge.highhb.2009<-Prop.pop.highhb.2009[1]*EIR.merge.6.2$EIR+
  Prop.pop.highhb.2009[2]*EIR.merge.7.2$EIR+Prop.pop.highhb.2009[3]*EIR.merge.8.2$EIR+
  Prop.pop.highhb.2009[4]*EIR.merge.9$EIR+Prop.pop.highhb.2009[5]*EIR.merge.10$EIR
EIR.point.highhb.2009<-Prop.pop.highhb.2009[1]*EIR.merge.6.2$Point.EIR[1]+
  Prop.pop.highhb.2009[2]*EIR.merge.7.2$Point.EIR[1]+Prop.pop.highhb.2009[3]*EIR.merge.8.2$Point.EIR[1]+
  Prop.pop.highhb.2009[4]*EIR.merge.9$Point.EIR[1]+Prop.pop.highhb.2009[5]*EIR.merge.10$Point.EIR[1]
EIR.highhb.2009.low<-sort(EIR.merge.highhb.2009)[n.cir*0.025]
EIR.highhb.2009.high<-sort(EIR.merge.highhb.2009)[n.cir*0.975]
EC.merge.highhb.2009<-EIR.merge.6.2$EC+
  EIR.merge.7.2$EC+EIR.merge.8.2$EC+
  EIR.merge.9$EC+EIR.merge.10$EC
EC.point.highhb.2009<-EIR.merge.6.2$Point.EC[1]+
  EIR.merge.7.2$Point.EC[1]+EIR.merge.8.2$Point.EC[1]+
  EIR.merge.9$Point.EC[1]+EIR.merge.10$Point.EC[1]
EC.highhb.2009.low<-sort(EC.merge.highhb.2009)[n.cir*0.025]
EC.highhb.2009.high<-sort(EC.merge.highhb.2009)[n.cir*0.975]

#######################Illiteracy rate subgroup EIR, 2009###################################
#Illiteracy rate category
LowIR<-c(1:10,14,18,19,20,21,31)
HighIR<-c(11,12,13,15,16,17,22:30)
#Low illiteracy rate
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order>2&
                                                             Hepatitis.B.agegroup10$Order<9|
                                                             Hepatitis.B.agegroup10$Order==10|
                                                             Hepatitis.B.agegroup10$Order==14|
                                                             Hepatitis.B.agegroup10$Order==18|
                                                             Hepatitis.B.agegroup10$Order==19|
                                                             Hepatitis.B.agegroup10$Order==20|
                                                             Hepatitis.B.agegroup10$Order==21|
                                                             Hepatitis.B.agegroup10$Order==31),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order<11|
                                                             Hepatitis.B.agegroup11$Order==14|
                                                             Hepatitis.B.agegroup11$Order==18|
                                                             Hepatitis.B.agegroup11$Order==19|
                                                             Hepatitis.B.agegroup11$Order==20|
                                                             Hepatitis.B.agegroup11$Order==21|
                                                             Hepatitis.B.agegroup11$Order==31),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order<11|
                                                             Hepatitis.B.agegroup12$Order==14|
                                                             Hepatitis.B.agegroup12$Order==18|
                                                             Hepatitis.B.agegroup12$Order==19|
                                                             Hepatitis.B.agegroup12$Order==20|
                                                             Hepatitis.B.agegroup12$Order==21|
                                                             Hepatitis.B.agegroup12$Order==31),],FUN=sum)[,2])
Age.pop9<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup13[which(Hepatitis.B.agegroup13$Order<11|
                                                             Hepatitis.B.agegroup13$Order==14|
                                                             Hepatitis.B.agegroup13$Order==18|
                                                             Hepatitis.B.agegroup13$Order==19|
                                                             Hepatitis.B.agegroup13$Order==20|
                                                             Hepatitis.B.agegroup13$Order==21|
                                                             Hepatitis.B.agegroup13$Order==31),],FUN=sum)[,2])
Age.pop10<-mean(aggregate(Population.book~Time,
                          data=Hepatitis.B.agegroup14[which(Hepatitis.B.agegroup14$Order<11|
                                                              Hepatitis.B.agegroup14$Order==14|
                                                              Hepatitis.B.agegroup14$Order==18|
                                                              Hepatitis.B.agegroup14$Order==19|
                                                              Hepatitis.B.agegroup14$Order==20|
                                                              Hepatitis.B.agegroup14$Order==21|
                                                              Hepatitis.B.agegroup14$Order==31),],FUN=sum)[,2])
Age.pop.2009<-c(Age.pop6,Age.pop7,Age.pop8,Age.pop9,Age.pop10)
Prop.pop.lowir.2009<-c(Age.pop6/sum(Age.pop.2009),Age.pop7/sum(Age.pop.2009),
                       Age.pop8/sum(Age.pop.2009),Age.pop9/sum(Age.pop.2009),Age.pop10/sum(Age.pop.2009))
EIR.merge.6.2<-EIR.merge.2009(c(3:8,10,14,18,19,20,21,31),Hepatitis.B.agegroup10,"model6",c(3,11),(2+12)/12,10000)
EIR.merge.7.2<-EIR.merge.2009(c(1:10,14,18,19,20,21,31),Hepatitis.B.agegroup11,"model7",c(3,11),(2+2*12)/12,10000)
EIR.merge.8.2<-EIR.merge.2009(c(1:10,14,18,19,20,21,31),Hepatitis.B.agegroup12,"model8",c(3,11),(2+7*12)/12,10000)
EIR.merge.9<-EIR.merge.2009(c(1:10,14,18,19,20,21,31),Hepatitis.B.agegroup13,"model9",c(2,9),6,10000)
EIR.merge.10<-EIR.merge.2009(c(1:10,14,18,19,20,21,31),Hepatitis.B.agegroup14,"model10",c(2,9),1,10000)
EIR.merge.lowir.2009<-Prop.pop.lowir.2009[1]*EIR.merge.6.2$EIR+
  Prop.pop.lowir.2009[2]*EIR.merge.7.2$EIR+Prop.pop.lowir.2009[3]*EIR.merge.8.2$EIR+
  Prop.pop.lowir.2009[4]*EIR.merge.9$EIR+Prop.pop.lowir.2009[5]*EIR.merge.10$EIR
EIR.point.lowir.2009<-Prop.pop.lowir.2009[1]*EIR.merge.6.2$Point.EIR[1]+
  Prop.pop.lowir.2009[2]*EIR.merge.7.2$Point.EIR[1]+Prop.pop.lowir.2009[3]*EIR.merge.8.2$Point.EIR[1]+
  Prop.pop.lowir.2009[4]*EIR.merge.9$Point.EIR[1]+Prop.pop.lowir.2009[5]*EIR.merge.10$Point.EIR[1]
EIR.lowir.2009.low<-sort(EIR.merge.lowir.2009)[n.cir*0.025]
EIR.lowir.2009.high<-sort(EIR.merge.lowir.2009)[n.cir*0.975]
EC.merge.lowir.2009<-EIR.merge.6.2$EC+
  EIR.merge.7.2$EC+EIR.merge.8.2$EC+
  EIR.merge.9$EC+EIR.merge.10$EC
EC.point.lowir.2009<-EIR.merge.6.2$Point.EC[1]+
  EIR.merge.7.2$Point.EC[1]+EIR.merge.8.2$Point.EC[1]+
  EIR.merge.9$Point.EC[1]+EIR.merge.10$Point.EC[1]
EC.lowir.2009.low<-sort(EC.merge.lowir.2009)[n.cir*0.025]
EC.lowir.2009.high<-sort(EC.merge.lowir.2009)[n.cir*0.975]

#High illiteracy rate
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order==11|
                                                             Hepatitis.B.agegroup10$Order==12|
                                                             Hepatitis.B.agegroup10$Order==13|
                                                             Hepatitis.B.agegroup10$Order==15|
                                                             Hepatitis.B.agegroup10$Order==16|
                                                             Hepatitis.B.agegroup10$Order==17|
                                                             Hepatitis.B.agegroup10$Order>21&
                                                             Hepatitis.B.agegroup10$Order<26|
                                                             Hepatitis.B.agegroup10$Order>26|
                                                             Hepatitis.B.agegroup10$Order<31),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order==11|
                                                             Hepatitis.B.agegroup11$Order==12|
                                                             Hepatitis.B.agegroup11$Order==13|
                                                             Hepatitis.B.agegroup11$Order==15|
                                                             Hepatitis.B.agegroup11$Order==16|
                                                             Hepatitis.B.agegroup11$Order==17|
                                                             Hepatitis.B.agegroup11$Order>21&
                                                             Hepatitis.B.agegroup11$Order<31),],FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order==11|
                                                             Hepatitis.B.agegroup12$Order==12|
                                                             Hepatitis.B.agegroup12$Order==13|
                                                             Hepatitis.B.agegroup12$Order==15|
                                                             Hepatitis.B.agegroup12$Order==16|
                                                             Hepatitis.B.agegroup12$Order==17|
                                                             Hepatitis.B.agegroup12$Order>21&
                                                             Hepatitis.B.agegroup12$Order<31),],FUN=sum)[,2])
Age.pop9<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup13[which(Hepatitis.B.agegroup13$Order==11|
                                                             Hepatitis.B.agegroup13$Order==12|
                                                             Hepatitis.B.agegroup13$Order==13|
                                                             Hepatitis.B.agegroup13$Order==15|
                                                             Hepatitis.B.agegroup13$Order==16|
                                                             Hepatitis.B.agegroup13$Order==17|
                                                             Hepatitis.B.agegroup13$Order>21&
                                                             Hepatitis.B.agegroup13$Order<31),],FUN=sum)[,2])
Age.pop10<-mean(aggregate(Population.book~Time,
                          data=Hepatitis.B.agegroup14[which(Hepatitis.B.agegroup14$Order==11|
                                                              Hepatitis.B.agegroup14$Order==12|
                                                              Hepatitis.B.agegroup14$Order==13|
                                                              Hepatitis.B.agegroup14$Order==15|
                                                              Hepatitis.B.agegroup14$Order==16|
                                                              Hepatitis.B.agegroup14$Order==17|
                                                              Hepatitis.B.agegroup14$Order>21&
                                                              Hepatitis.B.agegroup14$Order<31),],FUN=sum)[,2])
Age.pop.2009<-c(Age.pop6,Age.pop7,Age.pop8,Age.pop9,Age.pop10)
Prop.pop.highir.2009<-c(Age.pop6/sum(Age.pop.2009),Age.pop7/sum(Age.pop.2009),
                        Age.pop8/sum(Age.pop.2009),Age.pop9/sum(Age.pop.2009),Age.pop10/sum(Age.pop.2009))
EIR.merge.6.2<-EIR.merge.2009(c(11,12,13,15,16,17,22:25,27:30),Hepatitis.B.agegroup10,"model6",c(3,11),(2+12)/12,10000)
EIR.merge.7.2<-EIR.merge.2009(c(11,12,13,15,16,17,22:30),Hepatitis.B.agegroup11,"model7",c(3,11),(2+2*12)/12,10000)
EIR.merge.8.2<-EIR.merge.2009(c(11,12,13,15,16,17,22:30),Hepatitis.B.agegroup12,"model8",c(3,11),(2+7*12)/12,10000)
EIR.merge.9<-EIR.merge.2009(c(11,12,13,15,16,17,22:30),Hepatitis.B.agegroup13,"model9",c(2,9),6,10000)
EIR.merge.10<-EIR.merge.2009(c(11,12,13,15,16,17,22:30),Hepatitis.B.agegroup14,"model10",c(2,9),1,10000)
EIR.merge.highir.2009<-Prop.pop.highir.2009[1]*EIR.merge.6.2$EIR+
  Prop.pop.highir.2009[2]*EIR.merge.7.2$EIR+Prop.pop.highir.2009[3]*EIR.merge.8.2$EIR+
  Prop.pop.highir.2009[4]*EIR.merge.9$EIR+Prop.pop.highir.2009[5]*EIR.merge.10$EIR
EIR.point.highir.2009<-Prop.pop.highir.2009[1]*EIR.merge.6.2$Point.EIR[1]+
  Prop.pop.highir.2009[2]*EIR.merge.7.2$Point.EIR[1]+Prop.pop.highir.2009[3]*EIR.merge.8.2$Point.EIR[1]+
  Prop.pop.highir.2009[4]*EIR.merge.9$Point.EIR[1]+Prop.pop.highir.2009[5]*EIR.merge.10$Point.EIR[1]
EIR.highir.2009.low<-sort(EIR.merge.highir.2009)[n.cir*0.025]
EIR.highir.2009.high<-sort(EIR.merge.highir.2009)[n.cir*0.975]
EC.merge.highir.2009<-EIR.merge.6.2$EC+
  EIR.merge.7.2$EC+EIR.merge.8.2$EC+
  EIR.merge.9$EC+EIR.merge.10$EC
EC.point.highir.2009<-EIR.merge.6.2$Point.EC[1]+
  EIR.merge.7.2$Point.EC[1]+EIR.merge.8.2$Point.EC[1]+
  EIR.merge.9$Point.EC[1]+EIR.merge.10$Point.EC[1]
EC.highir.2009.low<-sort(EC.merge.highir.2009)[n.cir*0.025]
EC.highir.2009.high<-sort(EC.merge.highir.2009)[n.cir*0.975]

#######################Pooling the results of EIR and EC########################################
EIR<-data.frame(EC.point.2002=EC.point.2002,EC.low.2002=EC.2002.low,
                EC.high.2002=EC.2002.high,EIR.point.2002=EIR.point.2002,
                EIR.low.2002=EIR.2002.low,EIR.high.2002=EIR.2002.high,
                EC.point.2009=EC.point.2009,EC.low.2009=EC.2009.low,
                EC.high.2009=EC.2009.high,EIR.point.2009=EIR.point.2009,
                EIR.low.2009=EIR.2009.low,EIR.high.2009=EIR.2009.high)
line1<-c(EC.point.north.2002,EC.north.2002.low,EC.north.2002.high,
         EIR.point.north.2002,EIR.north.2002.low,EIR.north.2002.high,
         EC.point.north.2009,EC.north.2009.low,EC.north.2009.high,
         EIR.point.north.2009,EIR.north.2009.low,EIR.north.2009.high)
line2<-c(EC.point.northeast.2002,EC.northeast.2002.low,EC.northeast.2002.high,
         EIR.point.northeast.2002,EIR.northeast.2002.low,EIR.northeast.2002.high,
         EC.point.northeast.2009,EC.northeast.2009.low,EC.northeast.2009.high,
         EIR.point.northeast.2009,EIR.northeast.2009.low,EIR.northeast.2009.high)
line3<-c(EC.point.east.2002,EC.east.2002.low,EC.east.2002.high,
         EIR.point.east.2002,EIR.east.2002.low,EIR.east.2002.high,
         EC.point.east.2009,EC.east.2009.low,EC.east.2009.high,
         EIR.point.east.2009,EIR.east.2009.low,EIR.east.2009.high)
line4<-c(EC.point.southcentral.2002,EC.southcentral.2002.low,EC.southcentral.2002.high,
         EIR.point.southcentral.2002,EIR.southcentral.2002.low,EIR.southcentral.2002.high,
         EC.point.southcentral.2009,EC.southcentral.2009.low,EC.southcentral.2009.high,
         EIR.point.southcentral.2009,EIR.southcentral.2009.low,EIR.southcentral.2009.high)
line5<-c(EC.point.southwest.2002,EC.southwest.2002.low,EC.southwest.2002.high,
         EIR.point.southwest.2002,EIR.southwest.2002.low,EIR.southwest.2002.high,
         EC.point.southwest.2009,EC.southwest.2009.low,EC.southwest.2009.high,
         EIR.point.southwest.2009,EIR.southwest.2009.low,EIR.southwest.2009.high)
line6<-c(EC.point.northwest.2002,EC.northwest.2002.low,EC.northwest.2002.high,
         EIR.point.northwest.2002,EIR.northwest.2002.low,EIR.northwest.2002.high,
         EC.point.northwest.2009,EC.northwest.2009.low,EC.northwest.2009.high,
         EIR.point.northwest.2009,EIR.northwest.2009.low,EIR.northwest.2009.high)
line7<-c(EC.point.lowpr.2002,EC.lowpr.2002.low,EC.lowpr.2002.high,
         EIR.point.lowpr.2002,EIR.lowpr.2002.low,EIR.lowpr.2002.high,
         EC.point.lowpr.2009,EC.lowpr.2009.low,EC.lowpr.2009.high,
         EIR.point.lowpr.2009,EIR.lowpr.2009.low,EIR.lowpr.2009.high)
line8<-c(EC.point.middlepr.2002,EC.middlepr.2002.low,EC.middlepr.2002.high,
         EIR.point.middlepr.2002,EIR.middlepr.2002.low,EIR.middlepr.2002.high,
         EC.point.middlepr.2009,EC.middlepr.2009.low,EC.middlepr.2009.high,
         EIR.point.middlepr.2009,EIR.middlepr.2009.low,EIR.middlepr.2009.high)
line9<-c(EC.point.highpr.2002,EC.highpr.2002.low,EC.highpr.2002.high,
         EIR.point.highpr.2002,EIR.highpr.2002.low,EIR.highpr.2002.high,
         EC.point.highpr.2009,EC.highpr.2009.low,EC.highpr.2009.high,
         EIR.point.highpr.2009,EIR.highpr.2009.low,EIR.highpr.2009.high)
line10<-c(EC.point.lowur.2002,EC.lowur.2002.low,EC.lowur.2002.high,
          EIR.point.lowur.2002,EIR.lowur.2002.low,EIR.lowur.2002.high,
          EC.point.lowur.2009,EC.lowur.2009.low,EC.lowur.2009.high,
          EIR.point.lowur.2009,EIR.lowur.2009.low,EIR.lowur.2009.high)
line11<-c(EC.point.highur.2002,EC.highur.2002.low,EC.highur.2002.high,
          EIR.point.highur.2002,EIR.highur.2002.low,EIR.highur.2002.high,
          EC.point.highur.2009,EC.highur.2009.low,EC.highur.2009.high,
          EIR.point.highur.2009,EIR.highur.2009.low,EIR.highur.2009.high)
line12<-c(EC.point.lowhb.2002,EC.lowhb.2002.low,EC.lowhb.2002.high,
          EIR.point.lowhb.2002,EIR.lowhb.2002.low,EIR.lowhb.2002.high,
          EC.point.lowhb.2009,EC.lowhb.2009.low,EC.lowhb.2009.high,
          EIR.point.lowhb.2009,EIR.lowhb.2009.low,EIR.lowhb.2009.high)
line13<-c(EC.point.highhb.2002,EC.highhb.2002.low,EC.highhb.2002.high,
          EIR.point.highhb.2002,EIR.highhb.2002.low,EIR.highhb.2002.high,
          EC.point.highhb.2009,EC.highhb.2009.low,EC.highhb.2009.high,
          EIR.point.highhb.2009,EIR.highhb.2009.low,EIR.highhb.2009.high)
line14<-c(EC.point.lowir.2002,EC.lowir.2002.low,EC.lowir.2002.high,
          EIR.point.lowir.2002,EIR.lowir.2002.low,EIR.lowir.2002.high,
          EC.point.lowir.2009,EC.lowir.2009.low,EC.lowir.2009.high,
          EIR.point.lowir.2009,EIR.lowir.2009.low,EIR.lowir.2009.high)
line15<-c(EC.point.highir.2002,EC.highir.2002.low,EC.highir.2002.high,
          EIR.point.highir.2002,EIR.highir.2002.low,EIR.highir.2002.high,
          EC.point.highir.2009,EC.highir.2009.low,EC.highir.2009.high,
          EIR.point.highir.2009,EIR.highir.2009.low,EIR.highir.2009.high)
EIR<-rbind(EIR,line1,line2,line3,line4,line5,line6,line7,line8,line9,line10,
           line11,line12,line13,line14,line15)
row.names(EIR)<-c("Overall","Northern","Northeast","East","Southcentral",
                  "Southwest","Northwest","Low prevalence","Middle prevalence",
                  "High prevalence","Low urbanization rate","High urbanization rate",
                  "Low hospitalization bed","High hospitalization bed",
                  "Low illiteracy rate","High illiteracy rate")
#######################Pooling province level beta, 2002###################
#Beijing, Tianjin, Shanghai, Tibet
Para.single<-as.data.frame(matrix(ncol=6,nrow=31,0))
yori.specific<-matrix(0,length(datalist),2,dimnames=list(regions, paste("beta",c(1,2),sep="")))
Sori.specific<-vector("list",length(datalist))
names(Sori.specific)<-regions
for(i in c(1,2,9,26)){
  y7<-yori7.1[i,]
  y8<-yori8.1[i,]
  s7<-Sori7.1[i]
  s8<-Sori8.1[i]
  Pop7<-mean(Hepatitis.B.agegroup11$Population.book[
    which(Hepatitis.B.agegroup11$Order==i)])
  Pop8<-mean(Hepatitis.B.agegroup12$Population.book[
    which(Hepatitis.B.agegroup12$Order==i)])
  Pop<-c(Pop7,Pop8)
  Prop1.1<-c(Pop[1]/sum(Pop),Pop[2]/sum(Pop))
  model.number1<-eval(parse(text=paste("model7.",i,sep="")))
  model.number2<-eval(parse(text=paste("model8.",i,sep="")))
  yori.specific[i,]<-y7*Prop1.1[1]+y8*Prop1.1[2]
  Sori.specific[[i]]<-Prop1.1[1]^2*s7[[1]]+s8[[1]]*Prop1.1[2]^2
  #Level change
  Point1<-matrix(nrow=2,ncol=1,c(model.number1$coefficients[2],
                                 model.number2$coefficients[2]))
  Para.single[i,1]<-Prop1.1%*%Point1
  Se1<-matrix(nrow=2,ncol=1,c(summary(model.number1)$se[2],
                              summary(model.number2)$se[2]))
  Para.single[i,2]<-sqrt(Prop1.1^2%*%Se1^2)
  #Pre-intervention slope
  Point2<-matrix(nrow=2,ncol=1,c(model.number1$coefficients[4],
                                 model.number2$coefficients[4]))
  Para.single[i,3]<-Prop1.1%*%Point2
  Se2<-matrix(nrow=2,ncol=1,c(summary(model.number1)$se[4],
                              summary(model.number2)$se[4]))
  Para.single[i,4]<-sqrt(Prop1.1^2%*%Se2^2)
  #Post-intervention slope
  length.coef1<-length(model.number1$coefficients)
  Point3<-matrix(nrow=2,ncol=1,c(model.number1$coefficients[length.coef1-1],
                                 model.number2$coefficients[length.coef1-1]))
  Para.single[i,5]<-Prop1.1%*%Point3
  Se3<-matrix(nrow=2,ncol=1,c(summary(model.number1)$se[length.coef1-1],
                              summary(model.number2)$se[length.coef1-1]))
  Para.single[i,6]<-sqrt(Prop1.1^2%*%Se3^2)
}

#Liaoning
for(i in c(6)){
  y3<-yori3[i,]
  y4<-yori4[i,]
  y5<-yori5.1[i,]
  y6<-yori6.1[i,]
  y7<-yori7.1[i,]
  y8<-yori8.1[i,]
  s3<-Sori3[i]
  s4<-Sori4[i]
  s5<-Sori5.1[i]
  s6<-Sori6.1[i]
  s7<-Sori7.1[i]
  s8<-Sori8.1[i]
  Pop3<-mean(Hepatitis.B.agegroup7$Population.book[
    which(Hepatitis.B.agegroup7$Order==i)])
  Pop4<-mean(Hepatitis.B.agegroup8$Population.book[
    which(Hepatitis.B.agegroup8$Order==i)])
  Pop5<-mean(Hepatitis.B.agegroup9$Population.book[
    which(Hepatitis.B.agegroup9$Order==i)])
  Pop6<-mean(Hepatitis.B.agegroup10$Population.book[
    which(Hepatitis.B.agegroup10$Order==i)])
  Pop7<-mean(Hepatitis.B.agegroup11$Population.book[
    which(Hepatitis.B.agegroup11$Order==i)])
  Pop8<-mean(Hepatitis.B.agegroup12$Population.book[
    which(Hepatitis.B.agegroup12$Order==i)])
  Pop<-c(Pop3,Pop4,Pop5,Pop6,Pop7,Pop8)
  Prop1.2<-c(Pop[1]/sum(Pop),Pop[2]/sum(Pop),
             Pop[3]/sum(Pop),Pop[4]/sum(Pop),
             Pop[5]/sum(Pop),Pop[6]/sum(Pop))
  model.number1<-eval(parse(text=paste("model3.",i,sep="")))
  model.number2<-eval(parse(text=paste("model4.",i,sep="")))
  model.number3<-eval(parse(text=paste("model5.",i,sep="")))
  model.number4<-eval(parse(text=paste("model6.",i,sep="")))
  model.number5<-eval(parse(text=paste("model7.",i,sep="")))
  model.number6<-eval(parse(text=paste("model8.",i,sep="")))
  yori.specific[i,]<-y3*Prop1.2[1]+y4*Prop1.2[2]+
    y5*Prop1.2[3]+y6*Prop1.2[4]+
    y7*Prop1.2[5]+y8*Prop1.2[6]
  Sori.specific[[i]]<-Prop1.2[1]^2*s3[[1]]+s4[[1]]*Prop1.2[2]^2+
    Prop1.2[3]^2*s5[[1]]+s6[[1]]*Prop1.2[4]^2+
    Prop1.2[5]^2*s7[[1]]+s8[[1]]*Prop1.2[6]^2
  Point1<-matrix(nrow=6,ncol=1,c(model.number1$coefficients[2],
                                 model.number2$coefficients[2],
                                 model.number3$coefficients[2],
                                 model.number4$coefficients[2],
                                 model.number5$coefficients[2],
                                 model.number6$coefficients[2]))
  Para.single[i,1]<-Prop1.2%*%Point1
  Se1<-matrix(nrow=6,ncol=1,c(summary(model.number1)$se[2],
                              summary(model.number2)$se[2],
                              summary(model.number3)$se[2],
                              summary(model.number4)$se[2],
                              summary(model.number5)$se[2],
                              summary(model.number6)$se[2]))
  Para.single[i,2]<-sqrt(Prop1.2^2%*%Se1^2)
  Point2<-matrix(nrow=6,ncol=1,c(model.number1$coefficients[3],
                                 model.number2$coefficients[3],
                                 model.number3$coefficients[4],
                                 model.number4$coefficients[4],
                                 model.number5$coefficients[4],
                                 model.number6$coefficients[4]))
  Para.single[i,3]<-Prop1.2%*%Point2
  Se2<-matrix(nrow=6,ncol=1,c(summary(model.number1)$se[3],
                              summary(model.number2)$se[3],
                              summary(model.number3)$se[4],
                              summary(model.number4)$se[4],
                              summary(model.number5)$se[4],
                              summary(model.number6)$se[4]))
  Para.single[i,4]<-sqrt(Prop1.2^2%*%Se2^2)
  length.coef1<-length(model.number1$coefficients)
  length.coef2<-length(model.number3$coefficients)
  Point3<-matrix(nrow=6,ncol=1,c(model.number1$coefficients[length.coef1],
                                 model.number2$coefficients[length.coef1],
                                 model.number3$coefficients[length.coef2-1],
                                 model.number4$coefficients[length.coef2-1],
                                 model.number5$coefficients[length.coef2-1],
                                 model.number6$coefficients[length.coef2-1]))
  Para.single[i,5]<-Prop1.2%*%Point3
  Se3<-matrix(nrow=6,ncol=1,c(summary(model.number1)$se[length.coef1],
                              summary(model.number2)$se[length.coef1],
                              summary(model.number3)$se[length.coef2-1],
                              summary(model.number4)$se[length.coef2-1],
                              summary(model.number5)$se[length.coef2-1],
                              summary(model.number6)$se[length.coef2-1]))
  Para.single[i,6]<-sqrt(Prop1.2^2%*%Se3^2)
}

#Jiangxi
for(i in c(14)){
  y2<-yori2[i,]
  y7<-yori7.1[i,]
  y8<-yori8.1[i,]
  s2<-Sori2[i]
  s7<-Sori7.1[i]
  s8<-Sori8.1[i]
  Pop2<-mean(Hepatitis.B.agegroup6$Population.book[
    which(Hepatitis.B.agegroup6$Order==i)])
  Pop7<-mean(Hepatitis.B.agegroup11$Population.book[
    which(Hepatitis.B.agegroup11$Order==i)])
  Pop8<-mean(Hepatitis.B.agegroup12$Population.book[
    which(Hepatitis.B.agegroup12$Order==i)])
  Pop<-c(Pop2,Pop7,Pop8)
  Prop1.3<-c(Pop[1]/sum(Pop),Pop[2]/sum(Pop),Pop[3]/sum(Pop))
  model.number1<-eval(parse(text=paste("model2.",i,sep="")))
  model.number2<-eval(parse(text=paste("model7.",i,sep="")))
  model.number3<-eval(parse(text=paste("model8.",i,sep="")))
  yori.specific[i,]<-y2*Prop1.3[1]+y7*Prop1.3[2]+y8*Prop1.3[3]
  Sori.specific[[i]]<-Prop1.3[1]^2*s2[[1]]+s7[[1]]*Prop1.3[2]^2+
    Prop1.3[3]^2*s8[[1]]
  Point1<-matrix(nrow=3,ncol=1,c(model.number1$coefficients[2],
                                 model.number2$coefficients[2],
                                 model.number3$coefficients[2]))
  Para.single[i,1]<-Prop1.3%*%Point1
  Se1<-matrix(nrow=3,ncol=1,c(summary(model.number1)$se[2],
                              summary(model.number2)$se[2],
                              summary(model.number3)$se[2]))
  Para.single[i,2]<-sqrt(Prop1.3^2%*%Se1^2)
  Point2<-matrix(nrow=3,ncol=1,c(model.number1$coefficients[3],
                                 model.number2$coefficients[4],
                                 model.number3$coefficients[4]))
  Para.single[i,3]<-Prop1.3%*%Point2
  Se2<-matrix(nrow=3,ncol=1,c(summary(model.number1)$se[3],
                              summary(model.number2)$se[4],
                              summary(model.number3)$se[4]))
  Para.single[i,4]<-sqrt(Prop1.3^2%*%Se2^2)
  length.coef1<-length(model.number1$coefficients)
  length.coef2<-length(model.number2$coefficients)
  Point3<-matrix(nrow=3,ncol=1,c(model.number1$coefficients[length.coef1],
                                 model.number2$coefficients[length.coef2-1],
                                 model.number3$coefficients[length.coef2-1]))
  Para.single[i,5]<-Prop1.3%*%Point3
  Se3<-matrix(nrow=3,ncol=1,c(summary(model.number1)$se[length.coef1],
                              summary(model.number2)$se[length.coef2-1],
                              summary(model.number3)$se[length.coef2-1]))
  Para.single[i,6]<-sqrt(Prop1.3^2%*%Se3^2)
}

#Other provinces
for(i in c(3:5,7,8,10:13,15:25,27:31)){
  y2<-yori2[i,]
  y3<-yori3[i,]
  y4<-yori4[i,]
  y5<-yori5.1[i,]
  y6<-yori6.1[i,]
  y7<-yori7.1[i,]
  y8<-yori8.1[i,]
  s2<-Sori2[i]
  s3<-Sori3[i]
  s4<-Sori4[i]
  s5<-Sori5.1[i]
  s6<-Sori6.1[i]
  s7<-Sori7.1[i]
  s8<-Sori8.1[i]
  Pop2<-mean(Hepatitis.B.agegroup6$Population.book[
    which(Hepatitis.B.agegroup6$Order==i)])
  Pop3<-mean(Hepatitis.B.agegroup7$Population.book[
    which(Hepatitis.B.agegroup7$Order==i)])
  Pop4<-mean(Hepatitis.B.agegroup8$Population.book[
    which(Hepatitis.B.agegroup8$Order==i)])
  Pop5<-mean(Hepatitis.B.agegroup9$Population.book[
    which(Hepatitis.B.agegroup9$Order==i)])
  Pop6<-mean(Hepatitis.B.agegroup10$Population.book[
    which(Hepatitis.B.agegroup10$Order==i)])
  Pop7<-mean(Hepatitis.B.agegroup11$Population.book[
    which(Hepatitis.B.agegroup11$Order==i)])
  Pop8<-mean(Hepatitis.B.agegroup12$Population.book[
    which(Hepatitis.B.agegroup12$Order==i)])
  Pop<-c(Pop2,Pop3,Pop4,Pop5,Pop6,Pop7,Pop8)
  Prop1.4<-c(Pop[1]/sum(Pop),Pop[2]/sum(Pop),
             Pop[3]/sum(Pop),Pop[4]/sum(Pop),
             Pop[5]/sum(Pop),Pop[6]/sum(Pop),Pop[7]/sum(Pop))
  model.number1<-eval(parse(text=paste("model2.",i,sep="")))
  model.number2<-eval(parse(text=paste("model3.",i,sep="")))
  model.number3<-eval(parse(text=paste("model4.",i,sep="")))
  model.number4<-eval(parse(text=paste("model5.",i,sep="")))
  model.number5<-eval(parse(text=paste("model6.",i,sep="")))
  model.number6<-eval(parse(text=paste("model7.",i,sep="")))
  model.number7<-eval(parse(text=paste("model8.",i,sep="")))
  yori.specific[i,]<-y2*Prop1.4[1]+y3*Prop1.4[2]+y4*Prop1.4[3]+
    y5*Prop1.4[4]+y6*Prop1.4[5]+
    y7*Prop1.4[6]+y8*Prop1.4[7]
  Sori.specific[[i]]<-Prop1.4[1]^2*s2[[1]]+Prop1.4[2]^2*s3[[1]]+
    s4[[1]]*Prop1.4[3]^2+
    Prop1.4[4]^2*s5[[1]]+s6[[1]]*Prop1.4[5]^2+
    Prop1.4[6]^2*s7[[1]]+s8[[1]]*Prop1.4[7]^2
  Point1<-matrix(nrow=7,ncol=1,c(model.number1$coefficients[2],
                                 model.number2$coefficients[2],
                                 model.number3$coefficients[2],
                                 model.number4$coefficients[2],
                                 model.number5$coefficients[2],
                                 model.number6$coefficients[2],
                                 model.number7$coefficients[2]))
  Para.single[i,1]<-Prop1.4%*%Point1
  Se1<-matrix(nrow=7,ncol=1,c(summary(model.number1)$se[2],
                              summary(model.number2)$se[2],
                              summary(model.number3)$se[2],
                              summary(model.number4)$se[2],
                              summary(model.number5)$se[2],
                              summary(model.number6)$se[2],
                              summary(model.number6)$se[2]))
  Para.single[i,2]<-sqrt(Prop1.4^2%*%Se1^2)
  Point2<-matrix(nrow=7,ncol=1,c(model.number1$coefficients[3],
                                 model.number2$coefficients[3],
                                 model.number3$coefficients[3],
                                 model.number4$coefficients[4],
                                 model.number5$coefficients[4],
                                 model.number6$coefficients[4],
                                 model.number7$coefficients[4]))
  Para.single[i,3]<-Prop1.4%*%Point2
  Se2<-matrix(nrow=7,ncol=1,c(summary(model.number1)$se[3],
                              summary(model.number2)$se[3],
                              summary(model.number3)$se[3],
                              summary(model.number4)$se[4],
                              summary(model.number5)$se[4],
                              summary(model.number6)$se[4],
                              summary(model.number7)$se[4]))
  Para.single[i,4]<-sqrt(Prop1.4^2%*%Se2^2)
  length.coef1<-length(model.number1$coefficients)
  length.coef2<-length(model.number4$coefficients)
  Point3<-matrix(nrow=7,ncol=1,c(model.number1$coefficients[length.coef1],
                                 model.number2$coefficients[length.coef1],
                                 model.number3$coefficients[length.coef1],
                                 model.number4$coefficients[length.coef2-1],
                                 model.number5$coefficients[length.coef2-1],
                                 model.number6$coefficients[length.coef2-1],
                                 model.number7$coefficients[length.coef2-1]))
  Para.single[i,5]<-Prop1.4%*%Point3
  Se3<-matrix(nrow=7,ncol=1,c(summary(model.number1)$se[length.coef1],
                              summary(model.number2)$se[length.coef1],
                              summary(model.number3)$se[length.coef1],
                              summary(model.number4)$se[length.coef2-1],
                              summary(model.number5)$se[length.coef2-1],
                              summary(model.number6)$se[length.coef2-1],
                              summary(model.number7)$se[length.coef2-1]))
  Para.single[i,6]<-sqrt(Prop1.4^2%*%Se3^2)
}
Para<-as.data.frame(matrix(ncol=4,nrow=31,0))
Para$V1<-paste0(sprintf("%0.3f",Para.single$V1)," (",
                sprintf("%0.3f",Para.single$V1-1.96*Para.single$V2),
                ", ",sprintf("%0.3f",Para.single$V1+1.96*Para.single$V2),")")
Para$V2<-paste0(sprintf("%0.3f",Para.single$V3)," (",
                sprintf("%0.3f",Para.single$V3-1.96*Para.single$V4),
                ", ",sprintf("%0.3f",Para.single$V3+1.96*Para.single$V4),")")
Para$V3<-paste0(sprintf("%0.2f",Para.single$V5)," (",
                sprintf("%0.2f",Para.single$V5-1.96*Para.single$V6),
                ", ",sprintf("%0.2f",Para.single$V5+1.96*Para.single$V6),")")
Para$V4<-paste0(sprintf("%0.3f",Para.single$V5)," (",
                sprintf("%0.3f",Para.single$V5-1.96*Para.single$V6),
                ", ",sprintf("%0.3f",Para.single$V5+1.96*Para.single$V6),")")

colnames(Para)<-c("Beta1","Beta2","Beta3-0.01","Beta3-0.001")
row.names(Para)<-Name

#######################Province level ER, 2002###########################
ER.province<-as.data.frame(matrix(0,ncol=3,nrow=31))
for(i in c(1:31)){
  coef.meta<-as.matrix(t(yori.specific[i,]))
  time.meta<-as.matrix(cbind(rep(1,(180+3*12)),c(1:(180+3*12))))
  cov.meta<-Sori.specific[[i]]
  pre.ER.meta<-c()
  ER.meta<-c()
  se.ER.meta<-c()
  ER.low.meta<-c()
  ER.high.meta<-c()
  for(l in 1:(180+3*12))
  {pre.ER.meta[l]<-t(time.meta[l,])%*%t(coef.meta)
  ER.meta[l]<-(exp(pre.ER.meta[l])-1)*100
  se.ER.meta[l]<-sqrt(time.meta[l,]%*%cov.meta%*%t(t(time.meta[l,]))) 
  ER.low.meta[l]=(exp(pre.ER.meta[l]-1.96*se.ER.meta[l])-1)*100
  ER.high.meta[l]=(exp(pre.ER.meta[l]+1.96*se.ER.meta[l])-1)*100
  }
  ER.province[i,1]<-mean(ER.meta)
  ER.province[i,2]<-mean(ER.low.meta)
  ER.province[i,3]<-mean(ER.high.meta)
}
colnames(ER.province)<-c("ER.point","ER.low","ER.high")
row.names(ER.province)<-Name
ER.province$ER<-paste0(sprintf("%0.2f",ER.province$ER.point)," (",
                       sprintf("%0.2f",ER.province$ER.low),", ",
                       sprintf("%0.2f",ER.province$ER.high),")")

#######################Province level EIR, 2002###########################
EIR.specific.2002<-function(i,database,model,Location,intervention,n.cir){
  set.seed(202403)
  EMR.reserve<-matrix(0,nrow=31,ncol=8)
  EMR.reserve.diff<-c()
  MA.reserve<-matrix(0,nrow=n.cir,ncol=31)
  data.number<-database[which(database$Order==i),]
  model.number<-eval(parse(text=paste(model,".",i,sep="")))
  mod.number<-as.matrix(data.frame(model.matrix(model.number),offset=log(data.number$Population.book)))
  mod.number.0<-mod.number
  mod.number.0[,Location]<-0
  line.row<-which(data.number$Intervention.2002==1)
  length.row<-length(line.row)
  freq.number<-data.number$Case[line.row]
  pop.number<-data.number$Population.book[line.row]
  coef.number<-c(summary(model.number)$p.coeff,1)
  diff.number<-exp(mod.number%*% coef.number)-exp(mod.number.0%*%coef.number)
  diff.number<-diff.number[which(diff.number!=0)]
  EIR.number<-sum(diff.number)/intervention/(sum(pop.number)/length(pop.number))*100000
  EMR.reserve[i,1]<-sum(diff.number)
  EMR.reserve[i,2]<-intervention
  EMR.reserve[i,3]<-sum(pop.number)/length(pop.number)
  EMR.reserve[i,4]<-EIR.number
  coef.ini<-coef.number[Location]
  length.coef<-length(coef.ini)
  cov.number<-vcov(model.number)[Location,Location]
  eigen.number<-eigen(cov.number)
  r.norm<-rnorm(length.coef*n.cir)
  r.norm<-matrix(r.norm,n.cir)
  coef.sim<-coef.ini+eigen.number$vectors%*%diag(sqrt(eigen.number$values),length.coef)%*%t(r.norm)
  for(k in 1:n.cir){
    coef.number[Location[1]]<-coef.sim[1,k]
    coef.number[Location[2]]<-coef.sim[2,k]
    diff.number<-exp(mod.number%*%coef.number)-exp(mod.number.0%*%coef.number)
    diff.number<-diff.number[which(diff.number!=0)]
    EMR.reserve.diff[k]<-sum(diff.number)
  }
  MA.reserve[,i]<-EMR.reserve.diff
  EMR.reserve[i,5]<-sort(EMR.reserve.diff)[n.cir*0.025]
  EMR.reserve[i,7]<-sort(EMR.reserve.diff)[n.cir*0.975]
  EMR.reserve[i,6]<-EMR.reserve[i,5]/EMR.reserve[i,2]/EMR.reserve[i,3]*100000
  EMR.reserve[i,8]<-EMR.reserve[i,7]/EMR.reserve[i,2]/EMR.reserve[i,3]*100000
  EMR.reserve[,2]<-intervention
  EMR.mainland.point<-sum(EMR.reserve[i,1])/mean(EMR.reserve[i,2])/sum(EMR.reserve[i,3])*100000
  EMR.mainland<-MA.reserve[,i]/mean(EMR.reserve[i,2])/sum(EMR.reserve[i,3])*100000
  EC.mainland<-MA.reserve[,i]
  EIR.mainland<-data.frame(EC=EC.mainland,EIR=EMR.mainland,
                           Point.EC=sum(EMR.reserve[i,1]),
                           Point.EIR=EMR.mainland.point)
  return(EIR.mainland)
}

EC.specific<-as.data.frame(matrix(0,nrow=31,ncol=6))
EIR.specific<-as.data.frame(matrix(0,nrow=31,ncol=6))
for(i in c(1,2,9,26)){
  EIR.specific.7.1<-EIR.specific.2002(c(i),Hepatitis.B.agegroup11,"model7",c(2,10),8,10000)
  EIR.specific.8.1<-EIR.specific.2002(c(i),Hepatitis.B.agegroup12,"model8",c(2,10),3,10000)
  EIR.specific.2002.1<-Prop1.1[1]*EIR.specific.7.1$EIR+Prop1.1[2]*EIR.specific.8.1$EIR
  EIR.specific[i,1]<-Prop1.1[1]*EIR.specific.7.1$Point.EIR[1]+Prop1.1[2]*EIR.specific.8.1$Point.EIR[1]
  EIR.specific[i,2]<-sort(EIR.specific.2002.1)[n.cir*0.025]
  EIR.specific[i,3]<-sort(EIR.specific.2002.1)[n.cir*0.975]
  EC.specific.2002.1<-EIR.specific.7.1$EC+EIR.specific.8.1$EC
  EC.specific[i,1]<-EIR.specific.7.1$Point.EC[1]+EIR.specific.8.1$Point.EC[1]
  EC.specific[i,2]<-sort(EC.specific.2002.1)[n.cir*0.025]
  EC.specific[i,3]<-sort(EC.specific.2002.1)[n.cir*0.975]
}

for(i in c(6)){
  EIR.specific.3<-EIR.specific.2002(c(i),Hepatitis.B.agegroup7,"model3",c(2,9),12,10000)
  EIR.specific.4<-EIR.specific.2002(c(i),Hepatitis.B.agegroup8,"model4",c(2,9),11,10000)
  EIR.specific.5.1<-EIR.specific.2002(c(i),Hepatitis.B.agegroup9,"model5",c(2,10),10,10000)
  EIR.specific.6.1<-EIR.specific.2002(c(i),Hepatitis.B.agegroup10,"model6",c(2,10),9,10000)
  EIR.specific.7.1<-EIR.specific.2002(c(i),Hepatitis.B.agegroup11,"model7",c(2,10),8,10000)
  EIR.specific.8.1<-EIR.specific.2002(c(i),Hepatitis.B.agegroup12,"model8",c(2,10),3,10000)
  EIR.specific.2002.2<-Prop1.2[1]*EIR.specific.3$EIR+Prop1.2[2]*EIR.specific.4$EIR+
    Prop1.2[3]*EIR.specific.5.1$EIR+Prop1.2[4]*EIR.specific.6.1$EIR+
    Prop1.2[5]*EIR.specific.7.1$EIR+Prop1.2[6]*EIR.specific.8.1$EIR
  EIR.specific[i,1]<-Prop1.2[1]*EIR.specific.3$Point.EIR[1]+Prop1.2[2]*EIR.specific.4$Point.EIR[1]+
    Prop1.2[3]*EIR.specific.5.1$Point.EIR[1]+Prop1.2[4]*EIR.specific.6.1$Point.EIR[1]+
    Prop1.2[5]*EIR.specific.7.1$Point.EIR[1]+Prop1.2[6]*EIR.specific.8.1$Point.EIR[1]
  EIR.specific[i,2]<-sort(EIR.specific.2002.2)[n.cir*0.025]
  EIR.specific[i,3]<-sort(EIR.specific.2002.2)[n.cir*0.975]
  EC.specific.2002.2<-EIR.specific.3$EC+EIR.specific.4$EC+
    EIR.specific.5.1$EC+EIR.specific.6.1$EC+
    EIR.specific.7.1$EC+EIR.specific.8.1$EC
  EC.specific[i,1]<-EIR.specific.3$Point.EC[1]+EIR.specific.4$Point.EC[1]+
    EIR.specific.5.1$Point.EC[1]+EIR.specific.6.1$Point.EC[1]+
    EIR.specific.7.1$Point.EC[1]+EIR.specific.8.1$Point.EC[1]
  EC.specific[i,2]<-sort(EC.specific.2002.2)[n.cir*0.025]
  EC.specific[i,3]<-sort(EC.specific.2002.2)[n.cir*0.975]
}

for(i in c(14)){
  EIR.specific.2<-EIR.specific.2002(c(i),Hepatitis.B.agegroup6,"model2",c(2,9),13,10000)
  EIR.specific.7.1<-EIR.specific.2002(c(i),Hepatitis.B.agegroup11,"model7",c(2,10),8,10000)
  EIR.specific.8.1<-EIR.specific.2002(c(i),Hepatitis.B.agegroup12,"model8",c(2,10),3,10000)
  EIR.specific.2002.3<-Prop1.3[1]*EIR.specific.2$EIR+
    Prop1.3[2]*EIR.specific.7.1$EIR+Prop1.3[3]*EIR.specific.8.1$EIR
  EIR.specific[i,1]<-Prop1.3[1]*EIR.specific.2$Point.EIR[1]+
    Prop1.3[2]*EIR.specific.7.1$Point.EIR[1]+Prop1.3[3]*EIR.specific.8.1$Point.EIR[1]
  EIR.specific[i,2]<-sort(EIR.specific.2002.3)[n.cir*0.025]
  EIR.specific[i,3]<-sort(EIR.specific.2002.3)[n.cir*0.975]
  EC.specific.2002.3<-EIR.specific.2$EC+
    EIR.specific.7.1$EC+EIR.specific.8.1$EC
  EC.specific[i,1]<-EIR.specific.2$Point.EC[1]+
    EIR.specific.7.1$Point.EC[1]+EIR.specific.8.1$Point.EC[1]
  EC.specific[i,2]<-sort(EC.specific.2002.3)[n.cir*0.025]
  EC.specific[i,3]<-sort(EC.specific.2002.3)[n.cir*0.975]
}

for(i in c(3:5,7,8,10:13,15:25,27:31)){
  EIR.specific.2<-EIR.specific.2002(c(i),Hepatitis.B.agegroup6,"model2",c(2,9),13,10000)
  EIR.specific.3<-EIR.specific.2002(c(i),Hepatitis.B.agegroup7,"model3",c(2,9),12,10000)
  EIR.specific.4<-EIR.specific.2002(c(i),Hepatitis.B.agegroup8,"model4",c(2,9),11,10000)
  EIR.specific.5.1<-EIR.specific.2002(c(i),Hepatitis.B.agegroup9,"model5",c(2,10),10,10000)
  EIR.specific.6.1<-EIR.specific.2002(c(i),Hepatitis.B.agegroup10,"model6",c(2,10),9,10000)
  EIR.specific.7.1<-EIR.specific.2002(c(i),Hepatitis.B.agegroup11,"model7",c(2,10),8,10000)
  EIR.specific.8.1<-EIR.specific.2002(c(i),Hepatitis.B.agegroup12,"model8",c(2,10),3,10000)
  EIR.specific.2002.4<-Prop1.4[1]*EIR.specific.2$EIR+
    Prop1.4[2]*EIR.specific.3$EIR+Prop1.4[3]*EIR.specific.4$EIR+
    Prop1.4[4]*EIR.specific.5.1$EIR+Prop1.4[5]*EIR.specific.6.1$EIR+
    Prop1.4[6]*EIR.specific.7.1$EIR+Prop1.4[7]*EIR.specific.8.1$EIR
  EIR.specific[i,1]<-Prop1.4[1]*EIR.specific.2$Point.EIR[1]+
    Prop1.4[2]*EIR.specific.3$Point.EIR[1]+Prop1.4[3]*EIR.specific.4$Point.EIR[1]+
    Prop1.4[4]*EIR.specific.5.1$Point.EIR[1]+Prop1.4[5]*EIR.specific.6.1$Point.EIR[1]+
    Prop1.4[6]*EIR.specific.7.1$Point.EIR[1]+Prop1.4[7]*EIR.specific.8.1$Point.EIR[1]
  EIR.specific[i,2]<-sort(EIR.specific.2002.4)[n.cir*0.025]
  EIR.specific[i,3]<-sort(EIR.specific.2002.4)[n.cir*0.975]
  EC.specific.2002.4<-EIR.specific.2$EC+
    EIR.specific.3$EC+EIR.specific.4$EC+
    EIR.specific.5.1$EC+EIR.specific.6.1$EC+
    EIR.specific.7.1$EC+EIR.specific.8.1$EC
  EC.specific[i,1]<-EIR.specific.2$Point.EC[1]+
    EIR.specific.3$Point.EC[1]+EIR.specific.4$Point.EC[1]+
    EIR.specific.5.1$Point.EC[1]+EIR.specific.6.1$Point.EC[1]+
    EIR.specific.7.1$Point.EC[1]+EIR.specific.8.1$Point.EC[1]
  EC.specific[i,2]<-sort(EC.specific.2002.4)[n.cir*0.025]
  EC.specific[i,3]<-sort(EC.specific.2002.4)[n.cir*0.975]
}

#######################Pooling province level beta, 2009###################
#Beijing, Tianjin, Shanghai, Tibet
Para.single<-as.data.frame(matrix(ncol=6,nrow=31,0))
yori.specific<-matrix(0,length(datalist),2,dimnames=list(regions, paste("beta",c(1,2),sep="")))
Sori.specific<-vector("list",length(datalist))
names(Sori.specific)<-regions
for(i in c(1,2,9,26)){
  y7<-yori7.2[i,]
  y8<-yori8.2[i,]
  y9<-yori9[i,]
  y10<-yori10[i,]
  s7<-Sori7.2[i]
  s8<-Sori8.2[i]
  s9<-Sori9[i]
  s10<-Sori10[i]
  Pop7<-mean(Hepatitis.B.agegroup11$Population.book[
    which(Hepatitis.B.agegroup11$Order==i)])
  Pop8<-mean(Hepatitis.B.agegroup12$Population.book[
    which(Hepatitis.B.agegroup12$Order==i)])
  Pop9<-mean(Hepatitis.B.agegroup13$Population.book[
    which(Hepatitis.B.agegroup13$Order==i)])
  Pop10<-mean(Hepatitis.B.agegroup14$Population.book[
    which(Hepatitis.B.agegroup14$Order==i)])
  Pop<-c(Pop7,Pop8,Pop9,Pop10)
  Prop2.1<-c(Pop[1]/sum(Pop),Pop[2]/sum(Pop),
             Pop[3]/sum(Pop),Pop[4]/sum(Pop))
  model.number1<-eval(parse(text=paste("model7.",i,sep="")))
  model.number2<-eval(parse(text=paste("model8.",i,sep="")))
  model.number3<-eval(parse(text=paste("model9.",i,sep="")))
  model.number4<-eval(parse(text=paste("model10.",i,sep="")))
  yori.specific[i,]<-y7*Prop2.1[1]+y8*Prop2.1[2]+
    y9*Prop2.1[3]+y10*Prop2.1[4]
  Sori.specific[[i]]<-Prop2.1[1]^2*s7[[1]]+s8[[1]]*Prop2.1[2]^2+
    Prop2.1[3]^2*s9[[1]]+s10[[1]]*Prop2.1[4]^2
  Point1<-matrix(nrow=4,ncol=1,c(model.number1$coefficients[3],
                                 model.number2$coefficients[3],
                                 model.number3$coefficients[2],
                                 model.number4$coefficients[2]))
  Para.single[i,1]<-Prop2.1%*%Point1
  Se1<-matrix(nrow=4,ncol=1,c(summary(model.number1)$se[3],
                              summary(model.number2)$se[3],
                              summary(model.number3)$se[2],
                              summary(model.number4)$se[2]))
  Para.single[i,2]<-sqrt(Prop2.1^2%*%Se1^2)
  Point2<-matrix(nrow=4,ncol=1,c(model.number1$coefficients[4],
                                 model.number2$coefficients[4],
                                 model.number3$coefficients[3],
                                 model.number4$coefficients[3]))
  Para.single[i,3]<-Prop2.1%*%Point2
  Se2<-matrix(nrow=4,ncol=1,c(summary(model.number1)$se[4],
                              summary(model.number2)$se[4],
                              summary(model.number3)$se[3],
                              summary(model.number4)$se[3]))
  Para.single[i,4]<-sqrt(Prop2.1^2%*%Se2^2)
  length.coef1<-length(model.number1$coefficients)
  length.coef2<-length(model.number3$coefficients)
  Point3<-matrix(nrow=4,ncol=1,c(model.number1$coefficients[length.coef1],
                                 model.number2$coefficients[length.coef1],
                                 model.number3$coefficients[length.coef2],
                                 model.number4$coefficients[length.coef2]))
  Para.single[i,5]<-Prop2.1%*%Point3
  Se3<-matrix(nrow=4,ncol=1,c(summary(model.number1)$se[length.coef1],
                              summary(model.number2)$se[length.coef1],
                              summary(model.number3)$se[length.coef2],
                              summary(model.number4)$se[length.coef2]))
  Para.single[i,6]<-sqrt(Prop2.1^2%*%Se3^2)
}

#Others
for(i in c(3:8,10:25,27:31)){
  y6<-yori6.2[i,]
  y7<-yori7.2[i,]
  y8<-yori8.2[i,]
  y9<-yori9[i,]
  y10<-yori10[i,]
  s6<-Sori6.2[i]
  s7<-Sori7.2[i]
  s8<-Sori8.2[i]
  s9<-Sori9[i]
  s10<-Sori10[i]
  Pop6<-mean(Hepatitis.B.agegroup10$Population.book[
    which(Hepatitis.B.agegroup10$Order==i)])
  Pop7<-mean(Hepatitis.B.agegroup11$Population.book[
    which(Hepatitis.B.agegroup11$Order==i)])
  Pop8<-mean(Hepatitis.B.agegroup12$Population.book[
    which(Hepatitis.B.agegroup12$Order==i)])
  Pop9<-mean(Hepatitis.B.agegroup13$Population.book[
    which(Hepatitis.B.agegroup13$Order==i)])
  Pop10<-mean(Hepatitis.B.agegroup14$Population.book[
    which(Hepatitis.B.agegroup14$Order==i)])
  Pop<-c(Pop6,Pop7,Pop8,Pop9,Pop10)
  Prop2.2<-c(Pop[1]/sum(Pop),Pop[2]/sum(Pop),Pop[3]/sum(Pop),
             Pop[4]/sum(Pop),Pop[5]/sum(Pop))
  model.number1<-eval(parse(text=paste("model6.",i,sep="")))
  model.number2<-eval(parse(text=paste("model7.",i,sep="")))
  model.number3<-eval(parse(text=paste("model8.",i,sep="")))
  model.number4<-eval(parse(text=paste("model9.",i,sep="")))
  model.number5<-eval(parse(text=paste("model10.",i,sep="")))
  yori.specific[i,]<-y6*Prop2.2[1]+y7*Prop2.2[2]+y8*Prop2.2[3]+
    y9*Prop2.2[4]+y10*Prop2.2[5]
  Sori.specific[[i]]<-Prop2.2[1]^2*s6[[1]]+s7[[1]]*Prop2.2[2]^2+
    Prop2.2[3]^2*s8[[1]]+Prop2.2[4]^2*s9[[1]]+s10[[1]]*Prop2.2[5]^2
  Point1<-matrix(nrow=5,ncol=1,c(model.number1$coefficients[3],
                                 model.number2$coefficients[3],
                                 model.number3$coefficients[3],
                                 model.number4$coefficients[2],
                                 model.number5$coefficients[2]))
  Para.single[i,1]<-Prop2.2%*%Point1
  Se1<-matrix(nrow=5,ncol=1,c(summary(model.number1)$se[3],
                              summary(model.number2)$se[3],
                              summary(model.number3)$se[3],
                              summary(model.number4)$se[2],
                              summary(model.number5)$se[2]))
  Para.single[i,2]<-sqrt(Prop2.2^2%*%Se1^2)
  Point2<-matrix(nrow=5,ncol=1,c(model.number1$coefficients[4],
                                 model.number2$coefficients[4],
                                 model.number3$coefficients[4],
                                 model.number4$coefficients[3],
                                 model.number5$coefficients[3]))
  Para.single[i,3]<-Prop2.2%*%Point2
  Se2<-matrix(nrow=5,ncol=1,c(summary(model.number1)$se[4],
                              summary(model.number2)$se[4],
                              summary(model.number3)$se[4],
                              summary(model.number4)$se[3],
                              summary(model.number5)$se[3]))
  Para.single[i,4]<-sqrt(Prop2.2^2%*%Se2^2)
  length.coef1<-length(model.number1$coefficients)
  length.coef2<-length(model.number4$coefficients)
  Point3<-matrix(nrow=5,ncol=1,c(model.number1$coefficients[length.coef1],
                                 model.number2$coefficients[length.coef1],
                                 model.number3$coefficients[length.coef1],
                                 model.number4$coefficients[length.coef2],
                                 model.number5$coefficients[length.coef2]))
  Para.single[i,5]<-Prop2.2%*%Point3
  Se3<-matrix(nrow=5,ncol=1,c(summary(model.number1)$se[length.coef1],
                              summary(model.number2)$se[length.coef1],
                              summary(model.number3)$se[length.coef1],
                              summary(model.number4)$se[length.coef2],
                              summary(model.number5)$se[length.coef2]))
  Para.single[i,6]<-sqrt(Prop2.2^2%*%Se3^2)
}

#Other provinces
for(i in c(3,4,6:8,10,13,16:20,22:24,27:31)){
  y6<-yori6.2[i,]
  y7<-yori7.2[i,]
  y8<-yori8.2[i,]
  y9<-yori9[i,]
  y10<-yori10[i,]
  s6<-Sori6.2[i]
  s7<-Sori7.2[i]
  s8<-Sori8.2[i]
  s9<-Sori9[i]
  s10<-Sori10[i]
  Pop6<-mean(Hepatitis.B.agegroup10$Population.book[
    which(Hepatitis.B.agegroup10$Order==i)])
  Pop7<-mean(Hepatitis.B.agegroup11$Population.book[
    which(Hepatitis.B.agegroup11$Order==i)])
  Pop8<-mean(Hepatitis.B.agegroup12$Population.book[
    which(Hepatitis.B.agegroup12$Order==i)])
  Pop9<-mean(Hepatitis.B.agegroup13$Population.book[
    which(Hepatitis.B.agegroup13$Order==i)])
  Pop10<-mean(Hepatitis.B.agegroup14$Population.book[
    which(Hepatitis.B.agegroup14$Order==i)])
  Pop<-c(Pop6,Pop7,Pop8,Pop9,Pop10)
  Prop2.2<-c(Pop[1]/sum(Pop),Pop[2]/sum(Pop),Pop[3]/sum(Pop),
             Pop[4]/sum(Pop),Pop[5]/sum(Pop))
  model.number1<-eval(parse(text=paste("model6.",i,sep="")))
  model.number2<-eval(parse(text=paste("model7.",i,sep="")))
  model.number3<-eval(parse(text=paste("model8.",i,sep="")))
  model.number4<-eval(parse(text=paste("model9.",i,sep="")))
  model.number5<-eval(parse(text=paste("model10.",i,sep="")))
  yori.specific[i,]<-y6*Prop2.2[1]+y7*Prop2.2[2]+y8*Prop2.2[3]+
    y9*Prop2.2[4]+y10*Prop2.2[5]
  Sori.specific[[i]]<-Prop2.2[1]^2*s6[[1]]+s7[[1]]*Prop2.2[2]^2+
    Prop2.2[3]^2*s8[[1]]+Prop2.2[4]^2*s9[[1]]+s10[[1]]*Prop2.2[5]^2
  Point1<-matrix(nrow=5,ncol=1,c(model.number1$coefficients[3],
                                 model.number2$coefficients[3],
                                 model.number3$coefficients[3],
                                 model.number4$coefficients[2],
                                 model.number5$coefficients[2]))
  Para.single[i,1]<-Prop2.2%*%Point1
  Se1<-matrix(nrow=5,ncol=1,c(summary(model.number1)$se[3],
                              summary(model.number2)$se[3],
                              summary(model.number3)$se[3],
                              summary(model.number4)$se[2],
                              summary(model.number5)$se[2]))
  Para.single[i,2]<-sqrt(Prop2.2^2%*%Se1^2)
  Point2<-matrix(nrow=5,ncol=1,c(model.number1$coefficients[4],
                                 model.number2$coefficients[4],
                                 model.number3$coefficients[4],
                                 model.number4$coefficients[3],
                                 model.number5$coefficients[3]))
  Para.single[i,3]<-Prop2.2%*%Point2
  Se2<-matrix(nrow=5,ncol=1,c(summary(model.number1)$se[4],
                              summary(model.number2)$se[4],
                              summary(model.number3)$se[4],
                              summary(model.number4)$se[3],
                              summary(model.number5)$se[3]))
  Para.single[i,4]<-sqrt(Prop2.2^2%*%Se2^2)
  length.coef1<-length(model.number1$coefficients)
  length.coef2<-length(model.number4$coefficients)
  Point3<-matrix(nrow=5,ncol=1,c(model.number1$coefficients[length.coef1],
                                 model.number2$coefficients[length.coef1],
                                 model.number3$coefficients[length.coef1],
                                 model.number4$coefficients[length.coef2],
                                 model.number5$coefficients[length.coef2]))
  Para.single[i,5]<-Prop2.2%*%Point3
  Se3<-matrix(nrow=5,ncol=1,c(summary(model.number1)$se[length.coef1],
                              summary(model.number2)$se[length.coef1],
                              summary(model.number3)$se[length.coef1],
                              summary(model.number4)$se[length.coef2],
                              summary(model.number5)$se[length.coef2]))
  Para.single[i,6]<-sqrt(Prop2.2^2%*%Se3^2)
}
Para<-as.data.frame(matrix(ncol=4,nrow=31,0))
Para$V1<-paste0(sprintf("%0.3f",Para.single$V1)," (",
                sprintf("%0.3f",Para.single$V1-1.96*Para.single$V2),
                ", ",sprintf("%0.3f",Para.single$V1+1.96*Para.single$V2),")")
Para$V2<-paste0(sprintf("%0.3f",Para.single$V3)," (",
                sprintf("%0.3f",Para.single$V3-1.96*Para.single$V4),
                ", ",sprintf("%0.3f",Para.single$V3+1.96*Para.single$V4),")")
Para$V3<-paste0(sprintf("%0.2f",Para.single$V5)," (",
                sprintf("%0.2f",Para.single$V5-1.96*Para.single$V6),
                ", ",sprintf("%0.2f",Para.single$V5+1.96*Para.single$V6),")")
Para$V4<-paste0(sprintf("%0.3f",Para.single$V5)," (",
                sprintf("%0.3f",Para.single$V5-1.96*Para.single$V6),
                ", ",sprintf("%0.3f",Para.single$V5+1.96*Para.single$V6),")")
colnames(Para)<-c("Beta1","Beta2","Beta3-0.01","Beta3-0.001")
row.names(Para)<-Name

#######################Province level ER, 2002###########################
ER.province<-as.data.frame(matrix(0,ncol=3,nrow=31))
for(i in c(1:31)){
  coef.meta<-as.matrix(t(yori.specific[i,]))
  time.meta<-as.matrix(cbind(rep(1,180-7*12),c((7*12):179)-7*12))
  cov.meta<-Sori.specific[[i]]
  pre.ER.meta<-c()
  ER.meta<-c()
  se.ER.meta<-c()
  ER.low.meta<-c()
  ER.high.meta<-c()
  for(l in 1:(180-7*12))
  {pre.ER.meta[l]<-t(time.meta[l,])%*%t(coef.meta)
  ER.meta[l]<-(exp(pre.ER.meta[l])-1)*100
  se.ER.meta[l]<-sqrt(time.meta[l,]%*%cov.meta%*%t(t(time.meta[l,]))) 
  ER.low.meta[l]=(exp(pre.ER.meta[l]-1.96*se.ER.meta[l])-1)*100
  ER.high.meta[l]=(exp(pre.ER.meta[l]+1.96*se.ER.meta[l])-1)*100
  }
  ER.province[i,1]<-mean(ER.meta)
  ER.province[i,2]<-mean(ER.low.meta)
  ER.province[i,3]<-mean(ER.high.meta)
}
colnames(ER.province)<-c("ER.point","ER.low","ER.high")
row.names(ER.province)<-Name
ER.province$ER<-paste0(sprintf("%0.2f",ER.province$ER.point)," (",
                       sprintf("%0.2f",ER.province$ER.low),", ",
                       sprintf("%0.2f",ER.province$ER.high),")")
#######################Province level EIR, 2009###########################
EIR.specific.2009<-function(i,database,model,Location,intervention,n.cir){
  set.seed(202403)
  EMR.reserve<-matrix(0,nrow=31,ncol=8)
  EMR.reserve.diff<-c()
  MA.reserve<-matrix(0,nrow=n.cir,ncol=31)
  data.number<-database[which(database$Order==i),]
  model.number<-eval(parse(text=paste(model,".",i,sep="")))
  mod.number<-as.matrix(data.frame(model.matrix(model.number),offset=log(data.number$Population.book)))
  mod.number.0<-mod.number
  mod.number.0[,Location]<-0
  line.row<-which(data.number$Intervention.2009==1)
  length.row<-length(line.row)
  freq.number<-data.number$Case[line.row]
  pop.number<-data.number$Population.book[line.row]
  coef.number<-c(summary(model.number)$p.coeff,1)
  diff.number<-exp(mod.number%*% coef.number)-exp(mod.number.0%*%coef.number)
  diff.number<-diff.number[which(diff.number!=0)]
  EIR.number<-sum(diff.number)/intervention/(sum(pop.number)/length(pop.number))*100000
  EMR.reserve[i,1]<-sum(diff.number)
  EMR.reserve[i,2]<-intervention
  EMR.reserve[i,3]<-sum(pop.number)/length(pop.number)
  EMR.reserve[i,4]<-EIR.number
  coef.ini<-coef.number[Location]
  length.coef<-length(coef.ini)
  cov.number<-vcov(model.number)[Location,Location]
  eigen.number<-eigen(cov.number)
  r.norm<-rnorm(length.coef*n.cir)
  r.norm<-matrix(r.norm,n.cir)
  coef.sim<-coef.ini+eigen.number$vectors%*%diag(sqrt(eigen.number$values),length.coef)%*%t(r.norm)
  for(k in 1:n.cir){
    coef.number[Location[1]]<-coef.sim[1,k]
    coef.number[Location[2]]<-coef.sim[2,k]
    diff.number<-exp(mod.number%*%coef.number)-exp(mod.number.0%*%coef.number)
    diff.number<-diff.number[which(diff.number!=0)]
    EMR.reserve.diff[k]<-sum(diff.number)
  }
  MA.reserve[,i]<-EMR.reserve.diff
  EMR.reserve[i,5]<-sort(EMR.reserve.diff)[n.cir*0.025]
  EMR.reserve[i,7]<-sort(EMR.reserve.diff)[n.cir*0.975]
  EMR.reserve[i,6]<-EMR.reserve[i,5]/EMR.reserve[i,2]/EMR.reserve[i,3]*100000
  EMR.reserve[i,8]<-EMR.reserve[i,7]/EMR.reserve[i,2]/EMR.reserve[i,3]*100000
  EMR.reserve[,2]<-intervention
  EMR.mainland.point<-sum(EMR.reserve[i,1])/mean(EMR.reserve[i,2])/sum(EMR.reserve[i,3])*100000
  EMR.mainland<-MA.reserve[,i]/mean(EMR.reserve[i,2])/sum(EMR.reserve[i,3])*100000
  EC.mainland<-MA.reserve[,i]
  EIR.mainland<-data.frame(EC=EC.mainland,EIR=EMR.mainland,
                           Point.EC=sum(EMR.reserve[i,1]),
                           Point.EIR=EMR.mainland.point)
  return(EIR.mainland)
}

for(i in c(1,2,9,26)){
  EIR.specific.7.2<-EIR.specific.2009(c(i),Hepatitis.B.agegroup11,"model7",c(3,11),(2+2*12)/12,10000)
  EIR.specific.8.2<-EIR.specific.2009(c(i),Hepatitis.B.agegroup12,"model8",c(3,11),(2+7*12)/12,10000)
  EIR.specific.9<-EIR.specific.2009(c(i),Hepatitis.B.agegroup13,"model9",c(2,9),6,10000)
  EIR.specific.10<-EIR.specific.2009(c(i),Hepatitis.B.agegroup14,"model10",c(2,9),1,10000)
  EIR.specific.2009.1<-Prop2.1[1]*EIR.specific.7.2$EIR+Prop2.1[2]*EIR.specific.8.2$EIR+
    Prop2.1[3]*EIR.specific.9$EIR+Prop2.1[4]*EIR.specific.10$EIR
  EIR.specific[i,4]<-Prop2.1[1]*EIR.specific.7.2$Point.EIR[1]+Prop2.1[2]*EIR.specific.8.2$Point.EIR[1]+
    Prop2.1[3]*EIR.specific.9$Point.EIR[1]+Prop2.1[4]*EIR.specific.10$Point.EIR[1]
  EIR.specific[i,5]<-sort(EIR.specific.2009.1)[n.cir*0.025]
  EIR.specific[i,6]<-sort(EIR.specific.2009.1)[n.cir*0.975]
  EC.specific.2009.1<-EIR.specific.7.2$EC+EIR.specific.8.2$EC+
    EIR.specific.9$EC+EIR.specific.10$EC
  EC.specific[i,4]<-EIR.specific.7.2$Point.EC[1]+EIR.specific.8.2$Point.EC[1]+
    EIR.specific.9$Point.EC[1]+EIR.specific.10$Point.EC[1]
  EC.specific[i,5]<-sort(EC.specific.2009.1)[n.cir*0.025]
  EC.specific[i,6]<-sort(EC.specific.2009.1)[n.cir*0.975]
}

for(i in c(3:8,10:25,27:31)){
  EIR.specific.6.2<-EIR.specific.2009(c(i),Hepatitis.B.agegroup10,"model6",c(3,11),(2+12)/12,10000)
  EIR.specific.7.2<-EIR.specific.2009(c(i),Hepatitis.B.agegroup11,"model7",c(3,11),(2+2*12)/12,10000)
  EIR.specific.8.2<-EIR.specific.2009(c(i),Hepatitis.B.agegroup12,"model8",c(3,11),(2+7*12)/12,10000)
  EIR.specific.9<-EIR.specific.2009(c(i),Hepatitis.B.agegroup13,"model9",c(2,9),6,10000)
  EIR.specific.10<-EIR.specific.2009(c(i),Hepatitis.B.agegroup14,"model10",c(2,9),1,10000)
  EIR.specific.2009.2<-Prop2.2[1]*EIR.specific.6.2$EIR+
    Prop2.2[2]*EIR.specific.7.2$EIR+Prop2.2[3]*EIR.specific.8.2$EIR+
    Prop2.2[4]*EIR.specific.9$EIR+Prop2.2[5]*EIR.specific.10$EIR
  EIR.specific[i,4]<-Prop2.2[1]*EIR.specific.6.2$Point.EIR[1]+
    Prop2.2[2]*EIR.specific.7.2$Point.EIR[1]+Prop2.2[3]*EIR.specific.8.2$Point.EIR[1]+
    Prop2.2[4]*EIR.specific.9$Point.EIR[1]+Prop2.2[5]*EIR.specific.10$Point.EIR[1]
  EIR.specific[i,5]<-sort(EIR.specific.2009.2)[n.cir*0.025]
  EIR.specific[i,6]<-sort(EIR.specific.2009.2)[n.cir*0.975]
  EC.specific.2009.2<-EIR.specific.6.2$EC+
    EIR.specific.7.2$EC+EIR.specific.8.2$EC+
    EIR.specific.9$EC+EIR.specific.10$EC
  EC.specific[i,4]<-EIR.specific.6.2$Point.EC[1]+
    EIR.specific.7.2$Point.EC[1]+EIR.specific.8.2$Point.EC[1]+
    EIR.specific.9$Point.EC[1]+EIR.specific.10$Point.EC[1]
  EC.specific[i,5]<-sort(EC.specific.2009.2)[n.cir*0.025]
  EC.specific[i,6]<-sort(EC.specific.2009.2)[n.cir*0.975]
}
colnames(EIR.specific)<-c("Point.2002","Low.2002","High.2002",
                          "Point.2009","Low.2009","High.2009")
colnames(EC.specific)<-c("Point.2002","Low.2002","High.2002",
                         "Point.2009","Low.2009","High.2009")
row.names(EIR.specific)<-Name
row.names(EC.specific)<-Name
for(i in 1:6)
{EIR.specific[,i]<-sprintf("%0.2f",EIR.specific[,i])
EC.specific[,i]<-sprintf("%0.0f",EC.specific[,i]) 
}

#######################Sensitivity analysis 1#################################
#######################Change the start point of interventions################
#For the intervention in 2009, December were selected
#December, 2009
#Age group 6
Hepatitis.B.agegroup6<-read.xlsx("agegroup6.xlsx")
Hepatitis.B.agegroup6$Month.factor<-factor(Hepatitis.B.agegroup6$Month)
Hepatitis.B.agegroup6$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup6$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup6$Intervention.2002<-factor(rep(c(rep(0,12*2),rep(1,13*12)),31))
Hepatitis.B.agegroup6$Interaction.2002<-rep(c(rep(0,2*12),c(0:(180-2*12-1))),31)
regions<-as.character(unique(Hepatitis.B.agegroup6$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup6[Hepatitis.B.agegroup6$Province==x, ])
seq(datalist)
yori2<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori2<-vector("list", length(datalist)); names(Sori2) <- regionsfor(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup6[which(Hepatitis.B.agegroup6$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+Time+
              ns(Temperature,df=3)+
              Holiday+Interaction.2002,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori2[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun),1]
Sori2[[i]] <- vcov(mfirst)[c(2,lengthfun),c(2,lengthfun)]
assign(paste("model2.",i,sep=""),mfirst)
}
mvall2<-mvmeta(yori2[-c(1,2,6,9,26),]~1, Sori2[-c(1,2,6,9,26)], method="reml")
summary(mvall2)

#Age group 7
Hepatitis.B.agegroup7<-read.xlsx("agegroup7.xlsx")
Hepatitis.B.agegroup7$Month.factor<-factor(Hepatitis.B.agegroup7$Month)
Hepatitis.B.agegroup7$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup7$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup7$Intervention.2002<-factor(rep(c(rep(0,12*3),rep(1,12*12)),31))
Hepatitis.B.agegroup7$Interaction.2002<-rep(c(rep(0,3*12),c(0:(180-3*12-1))),31)
regions<-as.character(unique(Hepatitis.B.agegroup7$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup7[Hepatitis.B.agegroup7$Province==x, ])
seq(datalist)
yori3<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori3<-vector("list", length(datalist)); names(Sori3) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup7[which(Hepatitis.B.agegroup7$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+Time+
              ns(Temperature,df=3)+
              Holiday+Interaction.2002,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori3[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun),1]
Sori3[[i]] <- vcov(mfirst)[c(2,lengthfun),c(2,lengthfun)]
assign(paste("model3.",i,sep=""),mfirst)
}
mvall3<-mvmeta(yori3[-c(1,2,9,14,26),]~1, Sori3[-c(1,2,9,14,26)], method="reml")
summary(mvall3)

#Age group 8
Hepatitis.B.agegroup8<-read.xlsx("agegroup8.xlsx")
Hepatitis.B.agegroup8$Month.factor<-factor(Hepatitis.B.agegroup8$Month)
Hepatitis.B.agegroup8$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup8$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup8$Intervention.2002<-factor(rep(c(rep(0,12*4),rep(1,12*11)),31))
Hepatitis.B.agegroup8$Interaction.2002<-rep(c(rep(0,4*12),c(0:(180-4*12-1))),31)
regions<-as.character(unique(Hepatitis.B.agegroup8$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup8[Hepatitis.B.agegroup8$Province==x, ])
seq(datalist)
yori4<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori4<-vector("list", length(datalist)); names(Sori4) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup8[which(Hepatitis.B.agegroup8$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+Time+
              ns(Temperature,df=3)+
              Holiday+Interaction.2002,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori4[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun),1]
Sori4[[i]] <- vcov(mfirst)[c(2,lengthfun),c(2,lengthfun)]
assign(paste("model4.",i,sep=""),mfirst)
}
mvall4<-mvmeta(yori4[-c(1,2,9,14,26),]~1, Sori4[-c(1,2,9,14,26)], method="reml")
summary(mvall4)

#Age group 9
Hepatitis.B.agegroup9<-read.xlsx("agegroup9.xlsx")
Hepatitis.B.agegroup9$Month.factor<-factor(Hepatitis.B.agegroup9$Month)
Hepatitis.B.agegroup9$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup9$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup9$Intervention.2002<-factor(rep(c(rep(0,12*5),rep(1,12*10)),31))
Hepatitis.B.agegroup9$Intervention.2009<-factor(rep(c(rep(0,12*4),rep(0,11),rep(1,1),rep(0,12*10)),31))
Hepatitis.B.agegroup9$Interaction.2002<-rep(c(rep(0,5*12),c(0:(180-5*12-1))),31)
regions<-as.character(unique(Hepatitis.B.agegroup9$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup9[Hepatitis.B.agegroup9$Province==x, ])
seq(datalist)
yori5.1<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori5.1<-vector("list", length(datalist)); names(Sori5.1) <- regions
yori5.2<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori5.2<-vector("list", length(datalist)); names(Sori5.2) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup9[which(Hepatitis.B.agegroup9$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+#Intervention.2009+
              Time+ns(Temperature,df=3)+
              Holiday+Interaction.2002+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori5.1[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun-1),1]
Sori5.1[[i]] <- vcov(mfirst)[c(2,lengthfun-1),c(2,lengthfun-1)]
assign(paste("model5.",i,sep=""),mfirst)
}
mvall5.1<-mvmeta(yori5.1[-c(1,2,9,14,26),]~1, Sori5.1[-c(1,2,9,14,26)], method="reml")
summary(mvall5.1)

#Age group 10
Hepatitis.B.agegroup10<-read.xlsx("agegroup10.xlsx")
Hepatitis.B.agegroup10$Month.factor<-factor(Hepatitis.B.agegroup10$Month)
Hepatitis.B.agegroup10$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup10$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup10$Intervention.2002<-factor(rep(c(rep(0,12*6),rep(1,12*9)),31))
Hepatitis.B.agegroup10$Intervention.2009<-factor(rep(c(rep(0,12*4),rep(0,11),rep(1,1+12),rep(0,12*9)),31))
Hepatitis.B.agegroup10$Interaction.2002<-rep(c(rep(0,6*12),c(0:(180-6*12-1))),31)
Hepatitis.B.agegroup10$Interaction.2009<-rep(c(rep(0,4*12),rep(0,11),c(0:(1+12-1)),rep(0,9*12)),31)
regions<-as.character(unique(Hepatitis.B.agegroup10$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup10[Hepatitis.B.agegroup10$Province==x, ])
seq(datalist)
yori6.1<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori6.1<-vector("list", length(datalist)); names(Sori6.1) <- regions
yori6.2<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori6.2<-vector("list", length(datalist)); names(Sori6.2) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+Intervention.2009+
              Time+ns(Temperature,df=3)+
              Holiday+Interaction.2002+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori6.1[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun-1),1]
Sori6.1[[i]] <- vcov(mfirst)[c(2,lengthfun-1),c(2,lengthfun-1)]
yori6.2[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(3,lengthfun),1]
Sori6.2[[i]] <- vcov(mfirst)[c(3,lengthfun),c(3,lengthfun)]
assign(paste("model6.",i,sep=""),mfirst)
}
mvall6.1<-mvmeta(yori6.1[-c(1,2,9,14,26),]~1, Sori6.1[-c(1,2,9,14,26)], method="reml")
summary(mvall6.1)
mvall6.2<-mvmeta(yori6.2[-c(1,2,9,14,26),]~1, Sori6.2[-c(1,2,9,14,26)], method="reml")
summary(mvall6.2)

#Age group 11
Hepatitis.B.agegroup11<-read.xlsx("agegroup11.xlsx")
Hepatitis.B.agegroup11$Month.factor<-factor(Hepatitis.B.agegroup11$Month)
Hepatitis.B.agegroup11$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup11$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup11$Intervention.2002<-factor(rep(c(rep(0,12*7),rep(1,12*8)),31))
Hepatitis.B.agegroup11$Intervention.2009<-factor(rep(c(rep(0,12*4),rep(0,11),rep(1,1+2*12),rep(0,12*8)),31))
Hepatitis.B.agegroup11$Interaction.2002<-rep(c(rep(0,7*12),c(0:(180-7*12-1))),31)
Hepatitis.B.agegroup11$Interaction.2009<-rep(c(rep(0,4*12),rep(0,11),c(0:(1+2*12-1)),rep(0,8*12)),31)
regions<-as.character(unique(Hepatitis.B.agegroup11$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup11[Hepatitis.B.agegroup11$Province==x, ])
seq(datalist)
yori7.1<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori7.1<-vector("list", length(datalist)); names(Sori7.1) <- regions
yori7.2<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori7.2<-vector("list", length(datalist)); names(Sori7.2) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+Intervention.2009+
              Time+ns(Temperature,df=3)+
              Holiday+Interaction.2002+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori7.1[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun-1),1]
Sori7.1[[i]] <- vcov(mfirst)[c(2,lengthfun-1),c(2,lengthfun-1)]
yori7.2[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(3,lengthfun),1]
Sori7.2[[i]] <- vcov(mfirst)[c(3,lengthfun),c(3,lengthfun)]
assign(paste("model7.",i,sep=""),mfirst)
}
mvall7.1<-mvmeta(yori7.1~1, Sori7.1, method="reml")
summary(mvall7.1)
mvall7.2<-mvmeta(yori7.2~1, Sori7.2, method="reml")
summary(mvall7.2)

#Age group 12
Hepatitis.B.agegroup12<-read.xlsx("agegroup12.xlsx")
Hepatitis.B.agegroup12$Month.factor<-factor(Hepatitis.B.agegroup12$Month)
Hepatitis.B.agegroup12$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup12$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup12$Intervention.2002<-factor(rep(c(rep(0,12*12),rep(1,12*3)),31))
Hepatitis.B.agegroup12$Intervention.2009<-factor(rep(c(rep(0,12*4),rep(0,11),rep(1,1+7*12),rep(0,12*3)),31))
Hepatitis.B.agegroup12$Interaction.2002<-rep(c(rep(0,12*12),c(0:(180-12*12-1))),31)
Hepatitis.B.agegroup12$Interaction.2009<-rep(c(rep(0,4*12),rep(0,11),c(0:(1+7*12-1)),rep(0,3*12)),31)
regions<-as.character(unique(Hepatitis.B.agegroup12$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup12[Hepatitis.B.agegroup12$Province==x, ])
seq(datalist)
yori8.1<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori8.1<-vector("list", length(datalist)); names(Sori8.1) <- regions
yori8.2<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori8.2<-vector("list", length(datalist)); names(Sori8.2) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+Intervention.2009+
              Time+ns(Temperature,df=3)+
              Holiday+Interaction.2002+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori8.1[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun-1),1]
Sori8.1[[i]] <- vcov(mfirst)[c(2,lengthfun-1),c(2,lengthfun-1)]
yori8.2[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(3,lengthfun),1]
Sori8.2[[i]] <- vcov(mfirst)[c(3,lengthfun),c(3,lengthfun)]
assign(paste("model8.",i,sep=""),mfirst)
}
mvall8.1<-mvmeta(yori8.1~1, Sori8.1, method="reml")
summary(mvall8.1)
mvall8.2<-mvmeta(yori8.2~1, Sori8.2, method="reml")
summary(mvall8.2)

#Calculate the weights
Age.pop2<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup6[-which(Hepatitis.B.agegroup6$Order==1 |
                                                             Hepatitis.B.agegroup6$Order==2 |
                                                             Hepatitis.B.agegroup6$Order==6 |
                                                             Hepatitis.B.agegroup6$Order==9 |
                                                             Hepatitis.B.agegroup6$Order==26),],FUN=sum)[,2])
Age.pop3<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup7[-which(Hepatitis.B.agegroup7$Order==1 |
                                                             Hepatitis.B.agegroup7$Order==2 |
                                                             Hepatitis.B.agegroup7$Order==9 |
                                                             Hepatitis.B.agegroup7$Order==14 |
                                                             Hepatitis.B.agegroup7$Order==26),],FUN=sum)[,2])
Age.pop4<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup8[-which(Hepatitis.B.agegroup8$Order==1 |
                                                             Hepatitis.B.agegroup8$Order==2 |
                                                             Hepatitis.B.agegroup8$Order==9 |
                                                             Hepatitis.B.agegroup8$Order==14 |
                                                             Hepatitis.B.agegroup8$Order==26),],FUN=sum)[,2])
Age.pop5<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup9[-which(Hepatitis.B.agegroup9$Order==1 |
                                                             Hepatitis.B.agegroup9$Order==2 |
                                                             Hepatitis.B.agegroup9$Order==9 |
                                                             Hepatitis.B.agegroup9$Order==14 |
                                                             Hepatitis.B.agegroup9$Order==26),],FUN=sum)[,2])
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[-which(Hepatitis.B.agegroup10$Order==1 |
                                                              Hepatitis.B.agegroup10$Order==2 |
                                                              Hepatitis.B.agegroup10$Order==9 |
                                                              Hepatitis.B.agegroup10$Order==14 |
                                                              Hepatitis.B.agegroup10$Order==26),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,data=Hepatitis.B.agegroup11,FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,data=Hepatitis.B.agegroup12,FUN=sum)[,2])
Age.pop<-c(Age.pop2,Age.pop3,Age.pop4,Age.pop5,Age.pop6,Age.pop7,Age.pop8)
Prop.pop<-c(Age.pop2/sum(Age.pop),Age.pop3/sum(Age.pop),Age.pop4/sum(Age.pop),
            Age.pop5/sum(Age.pop),Age.pop6/sum(Age.pop),Age.pop7/sum(Age.pop),
            Age.pop8/sum(Age.pop))
Merge.coef<-Prop.pop[1]*mvall2$coefficients+Prop.pop[2]*mvall3$coefficients+
  Prop.pop[3]*mvall4$coefficients+Prop.pop[4]*mvall5.1$coefficients+
  Prop.pop[5]*mvall6.1$coefficients+Prop.pop[6]*mvall7.1$coefficients+
  Prop.pop[7]*mvall8.1$coefficients
Merge.vcov<-Prop.pop[1]^2*mvall2$vcov+Prop.pop[2]^2*mvall3$vcov+
  Prop.pop[3]^2*mvall4$vcov+Prop.pop[4]^2*mvall5.1$vcov+
  Prop.pop[5]^2*mvall6.1$vcov+Prop.pop[6]^2*mvall7.1$vcov+
  Prop.pop[7]^2*mvall8.1$vcov

#Age-standardized ER
coef.meta<-as.matrix(t(Merge.coef))
time.meta<-as.matrix(cbind(rep(1,180+3*12),c(1:(180+3*12))-1))
cov.meta<-Merge.vcov
pre.ER.meta<-c()
ER.meta<-c()
se.ER.meta<-c()
ER.low.meta<-c()
ER.high.meta<-c()
for(l in 1:(180+3*12))
{pre.ER.meta[l]<-time.meta[l,]%*%coef.meta
ER.meta[l]<-(exp(pre.ER.meta[l])-1)*100
se.ER.meta[l]<-sqrt(time.meta[l,]%*%cov.meta%*%t(t(time.meta[l,]))) 
ER.low.meta[l]=(exp(pre.ER.meta[l]-1.96*se.ER.meta[l])-1)*100
ER.high.meta[l]=(exp(pre.ER.meta[l]+1.96*se.ER.meta[l])-1)*100
}
ER.whole<-ER.meta
ER.whole.low<-ER.low.meta
ER.whole.upper<-ER.high.meta
month.meta<-c(1:(180+3*12))
ER<-data.frame(Time=month.meta,ER=ER.meta,ER.lower=ER.low.meta,ER.upper=ER.high.meta)
ER$ER<-sprintf("%0.2f",ER$ER)
ER$ER.lower<-sprintf("%0.2f",ER$ER.lower)
ER$ER.upper<-sprintf("%0.2f",ER$ER.upper)

#Age group 13
Hepatitis.B.agegroup13<-read.xlsx("Agegroup13.xlsx")
Hepatitis.B.agegroup13$Month.factor<-factor(Hepatitis.B.agegroup13$Month)
Hepatitis.B.agegroup13$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup13$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup13$Intervention.2009<-factor(rep(c(rep(0,9*12),rep(1,12*6)),31))
Hepatitis.B.agegroup13$Interaction.2009<-rep(c(rep(0,9*12),c(0:(180-12*9-1))),31)
regions<-as.character(unique(Hepatitis.B.agegroup13$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup13[Hepatitis.B.agegroup13$Province==x, ])
seq(datalist)
yori9<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori9<-vector("list", length(datalist)); names(Sori9) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup13[which(Hepatitis.B.agegroup13$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2009+Time+
              ns(Temperature,df=3)+
              Holiday+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori9[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun),1]
Sori9[[i]] <- vcov(mfirst)[c(2,lengthfun),c(2,lengthfun)]
assign(paste("model9.",i,sep=""),mfirst)
}
mvall9<-mvmeta(yori9~1, Sori9, method="reml")
summary(mvall9)

#Age group 14
Hepatitis.B.agegroup14<-read.xlsx("Agegroup14.xlsx")
Hepatitis.B.agegroup14$Month.factor<-factor(Hepatitis.B.agegroup14$Month)
Hepatitis.B.agegroup14$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup14$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup14$Intervention.2009<-factor(rep(c(rep(0,14*12),rep(1,12*1)),31))
Hepatitis.B.agegroup14$Interaction.2009<-rep(c(rep(0,14*12),c(0:(180-12*14-1))),31)
regions<-as.character(unique(Hepatitis.B.agegroup14$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup14[Hepatitis.B.agegroup14$Province==x, ])
seq(datalist)
yori10<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori10<-vector("list", length(datalist)); names(Sori10) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup14[which(Hepatitis.B.agegroup14$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2009+Time+
              ns(Temperature,df=3)+
              Holiday+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori10[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun),1]
Sori10[[i]] <- vcov(mfirst)[c(2,lengthfun),c(2,lengthfun)]
assign(paste("model10.",i,sep=""),mfirst)
}
mvall10<-mvmeta(yori10~1, Sori10, method="reml")
summary(mvall10)

#Calculate the weights
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[-which(Hepatitis.B.agegroup10$Order==1 |
                                                              Hepatitis.B.agegroup10$Order==2 |
                                                              Hepatitis.B.agegroup10$Order==9 |
                                                              Hepatitis.B.agegroup10$Order==26),],FUN=sum)[,2])
Age.pop9<-mean(aggregate(Population.book~Time,data=Hepatitis.B.agegroup13,FUN=sum)[,2])
Age.pop10<-mean(aggregate(Population.book~Time,data=Hepatitis.B.agegroup14,FUN=sum)[,2])
Age.pop.2009<-c(Age.pop6,Age.pop7,Age.pop8,Age.pop9,Age.pop10)
Prop.pop.2009<-c(Age.pop6/sum(Age.pop.2009),Age.pop7/sum(Age.pop.2009),
                 Age.pop8/sum(Age.pop.2009),Age.pop9/sum(Age.pop.2009),Age.pop10/sum(Age.pop.2009))
Merge.coef.2009<-Prop.pop.2009[1]*mvall6.2$coefficients+
  Prop.pop.2009[2]*mvall7.2$coefficients+Prop.pop.2009[3]*mvall8.2$coefficients+
  Prop.pop.2009[4]*mvall9$coefficients+Prop.pop.2009[5]*mvall10$coefficients
Merge.vcov.2009<-Prop.pop.2009[1]^2*mvall6.2$vcov+
  Prop.pop.2009[2]^2*mvall7.2$vcov+Prop.pop.2009[3]^2*mvall8.2$vcov+
  Prop.pop.2009[4]^2*mvall9$vcov+Prop.pop.2009[5]^2*mvall10$vcov

#Age-standardized ER
coef.meta<-as.matrix(t(Merge.coef.2009))
time.meta<-as.matrix(cbind(rep(1,180-7*12),c((7*12):179)-7*12))
cov.meta<-Merge.vcov.2009
pre.ER.meta<-c()
ER.meta<-c()
se.ER.meta<-c()
ER.low.meta<-c()
ER.high.meta<-c()
for(l in 1:(180-7*12))
{pre.ER.meta[l]<-time.meta[l,]%*%coef.meta
ER.meta[l]<-(exp(pre.ER.meta[l])-1)*100
se.ER.meta[l]<-sqrt(time.meta[l,]%*%cov.meta%*%t(t(time.meta[l,]))) 
ER.low.meta[l]=(exp(pre.ER.meta[l]-1.96*se.ER.meta[l])-1)*100
ER.high.meta[l]=(exp(pre.ER.meta[l]+1.96*se.ER.meta[l])-1)*100
}
ER.whole<-ER.meta
ER.whole.low<-ER.low.meta
ER.whole.upper<-ER.high.meta
month.meta<-c((7*12+1):180)
ER<-data.frame(Time=month.meta,ER=ER.meta,ER.lower=ER.low.meta,ER.upper=ER.high.meta)
ER$ER<-sprintf("%0.2f",ER$ER)
ER$ER.lower<-sprintf("%0.2f",ER$ER.lower)
ER$ER.upper<-sprintf("%0.2f",ER$ER.upper)

#######################The random sampling of under-reporting distribution#############
set.seed(20250717)#set.seed(202507)
#A left-truncated normal distribution (μ=25%, σ=5%, and truncated point=0) provides a close approximation to
#a normal distribution (μ=25%, σ=5%).
Inflation<-c(rnorm(180,25,5))/100
setwd("D:/~")
#Age group 6
Hepatitis.B.agegroup6<-read.xlsx("agegroup6.xlsx")
Hepatitis.B.agegroup6$Month.factor<-factor(Hepatitis.B.agegroup6$Month)
Hepatitis.B.agegroup6$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup6$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup6$Intervention.2002<-factor(rep(c(rep(0,12*2),rep(1,13*12)),31))
Hepatitis.B.agegroup6$Interaction.2002<-rep(c(rep(0,2*12),c(0:(180-2*12-1))),31)
Hepatitis.B.agegroup6$Case.updated<-round(Hepatitis.B.agegroup6$Case*(1+Inflation))
regions<-as.character(unique(Hepatitis.B.agegroup6$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup6[Hepatitis.B.agegroup6$Province==x, ])
seq(datalist)
yori2<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori2<-vector("list", length(datalist)); names(Sori2) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup6[which(Hepatitis.B.agegroup6$Order==i),]
mfirst<-gam(Case.updated~ offset(log(Population.book))+Intervention.2002+Time+
              ns(Temperature,df=3)+
              Holiday+Interaction.2002,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori2[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun),1]
Sori2[[i]] <- vcov(mfirst)[c(2,lengthfun),c(2,lengthfun)]
assign(paste("model2.",i,sep=""),mfirst)
}
mvall2<-mvmeta(yori2[-c(1,2,6,9,26),]~1, Sori2[-c(1,2,6,9,26)], method="reml")
summary(mvall2)

#Age group 7
Hepatitis.B.agegroup7<-read.xlsx("agegroup7.xlsx")
Hepatitis.B.agegroup7$Month.factor<-factor(Hepatitis.B.agegroup7$Month)
Hepatitis.B.agegroup7$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup7$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup7$Intervention.2002<-factor(rep(c(rep(0,12*3),rep(1,12*12)),31))
Hepatitis.B.agegroup7$Interaction.2002<-rep(c(rep(0,3*12),c(0:(180-3*12-1))),31)
Hepatitis.B.agegroup7$Case.updated<-round(Hepatitis.B.agegroup7$Case*(1+Inflation))
regions<-as.character(unique(Hepatitis.B.agegroup7$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup7[Hepatitis.B.agegroup7$Province==x, ])
seq(datalist)
yori3<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori3<-vector("list", length(datalist)); names(Sori3) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup7[which(Hepatitis.B.agegroup7$Order==i),]
mfirst<-gam(Case.updated~ offset(log(Population.book))+Intervention.2002+Time+
              ns(Temperature,df=3)+
              Holiday+Interaction.2002,#+Holiday+GDP+Birth_Rate
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori3[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun),1]
Sori3[[i]] <- vcov(mfirst)[c(2,lengthfun),c(2,lengthfun)]
assign(paste("model3.",i,sep=""),mfirst)
}
mvall3<-mvmeta(yori3[-c(1,2,9,14,26),]~1, Sori3[-c(1,2,9,14,26)], method="reml")
summary(mvall3)

#Age group 8
Hepatitis.B.agegroup8<-read.xlsx("agegroup8.xlsx")
Hepatitis.B.agegroup8$Month.factor<-factor(Hepatitis.B.agegroup8$Month)
Hepatitis.B.agegroup8$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup8$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup8$Intervention.2002<-factor(rep(c(rep(0,12*4),rep(1,12*11)),31))
Hepatitis.B.agegroup8$Interaction.2002<-rep(c(rep(0,4*12),c(0:(180-4*12-1))),31)
Hepatitis.B.agegroup8$Case.updated<-round(Hepatitis.B.agegroup8$Case*(1+Inflation))
regions<-as.character(unique(Hepatitis.B.agegroup8$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup8[Hepatitis.B.agegroup8$Province==x, ])
seq(datalist)
yori4<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori4<-vector("list", length(datalist)); names(Sori4) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup8[which(Hepatitis.B.agegroup8$Order==i),]
mfirst<-gam(Case.updated~ offset(log(Population.book))+Intervention.2002+Time+
              ns(Temperature,df=3)+
              Holiday+Interaction.2002,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori4[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun),1]
Sori4[[i]] <- vcov(mfirst)[c(2,lengthfun),c(2,lengthfun)]
assign(paste("model4.",i,sep=""),mfirst)
}
mvall4<-mvmeta(yori4[-c(1,2,9,14,26),]~1, Sori4[-c(1,2,9,14,26)], method="reml")
summary(mvall4)

#Age group 9
setwd("D:/研究生/乙肝干预/全国研究")
Hepatitis.B.agegroup9<-read.xlsx("agegroup9.xlsx")
Hepatitis.B.agegroup9$Month.factor<-factor(Hepatitis.B.agegroup9$Month)
Hepatitis.B.agegroup9$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup9$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup9$Intervention.2002<-factor(rep(c(rep(0,12*5),rep(1,12*10)),31))
Hepatitis.B.agegroup9$Intervention.2009<-factor(rep(c(rep(0,12*4),rep(0,10),rep(1,2),rep(0,12*10)),31))
Hepatitis.B.agegroup9$Interaction.2002<-rep(c(rep(0,5*12),c(0:(180-5*12-1))),31)
Hepatitis.B.agegroup9$Interaction.2009<-rep(c(rep(0,4*12),rep(0,10),c(0:(2-1)),rep(0,10*12)),31)
Hepatitis.B.agegroup9$Case.updated<-round(Hepatitis.B.agegroup9$Case*(1+Inflation))
regions<-as.character(unique(Hepatitis.B.agegroup9$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup9[Hepatitis.B.agegroup9$Province==x, ])
seq(datalist)
yori5.1<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori5.1<-vector("list", length(datalist)); names(Sori5.1) <- regions
yori5.2<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori5.2<-vector("list", length(datalist)); names(Sori5.2) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup9[which(Hepatitis.B.agegroup9$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+#Intervention.2009+
              Time+ns(Temperature,df=3)+
              Holiday+Interaction.2002+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori5.1[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun-1),1]
Sori5.1[[i]] <- vcov(mfirst)[c(2,lengthfun-1),c(2,lengthfun-1)]
yori5.2[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(3,lengthfun),1]
Sori5.2[[i]] <- vcov(mfirst)[c(3,lengthfun),c(3,lengthfun)]
assign(paste("model5.",i,sep=""),mfirst)
}
mvall5.1<-mvmeta(yori5.1[-c(1,2,9,14,26),]~1, Sori5.1[-c(1,2,9,14,26)], method="reml")
summary(mvall5.1)

#Age group 10
setwd("D:/研究生/乙肝干预/全国研究")
Hepatitis.B.agegroup10<-read.xlsx("agegroup10.xlsx")
Hepatitis.B.agegroup10$Month.factor<-factor(Hepatitis.B.agegroup10$Month)
Hepatitis.B.agegroup10$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup10$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup10$Intervention.2002<-factor(rep(c(rep(0,12*6),rep(1,12*9)),31))
Hepatitis.B.agegroup10$Intervention.2009<-factor(rep(c(rep(0,12*4),rep(0,10),rep(1,2+12),rep(0,12*9)),31))
Hepatitis.B.agegroup10$Interaction.2002<-rep(c(rep(0,6*12),c(0:(180-6*12-1))),31)
Hepatitis.B.agegroup10$Interaction.2009<-rep(c(rep(0,4*12),rep(0,10),c(0:(2+12-1)),rep(0,9*12)),31)
Hepatitis.B.agegroup10$Case.updated<-round(Hepatitis.B.agegroup10$Case*(1+Inflation))
regions<-as.character(unique(Hepatitis.B.agegroup10$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup10[Hepatitis.B.agegroup10$Province==x, ])
seq(datalist)
yori6.1<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori6.1<-vector("list", length(datalist)); names(Sori6.1) <- regions
yori6.2<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori6.2<-vector("list", length(datalist)); names(Sori6.2) <- regions
Percent<-as.data.frame(matrix(0,nrow=31,ncol=6))
Percent[,1]<-regions
year.pre<-rep(6,31)*12
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order==i),]
mfirst<-gam(Case.updated~ offset(log(Population.book))+Intervention.2002+Intervention.2009+
              Time+ns(Temperature,df=3)+
              Holiday+Interaction.2002+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori6.1[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun-1),1]
Sori6.1[[i]] <- vcov(mfirst)[c(2,lengthfun-1),c(2,lengthfun-1)]
yori6.2[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(3,lengthfun),1]
Sori6.2[[i]] <- vcov(mfirst)[c(3,lengthfun),c(3,lengthfun)]
assign(paste("model6.",i,sep=""),mfirst)
}
mvall6.1<-mvmeta(yori6.1[-c(1,2,9,14,26),]~1, Sori6.1[-c(1,2,9,14,26)], method="reml")
summary(mvall6.1)
mvall6.2<-mvmeta(yori6.2[-c(1,2,9,14,26),]~1, Sori6.2[-c(1,2,9,14,26)], method="reml")
summary(mvall6.2)

#Age group 11
Hepatitis.B.agegroup11<-read.xlsx("agegroup11.xlsx")
Hepatitis.B.agegroup11$Month.factor<-factor(Hepatitis.B.agegroup11$Month)
Hepatitis.B.agegroup11$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup11$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup11$Intervention.2002<-factor(rep(c(rep(0,12*7),rep(1,12*8)),31))
Hepatitis.B.agegroup11$Intervention.2009<-factor(rep(c(rep(0,12*4),rep(0,10),rep(1,2+2*12),rep(0,12*8)),31))
Hepatitis.B.agegroup11$Interaction.2002<-rep(c(rep(0,7*12),c(0:(180-7*12-1))),31)
Hepatitis.B.agegroup11$Interaction.2009<-rep(c(rep(0,4*12),rep(0,10),c(0:(2+2*12-1)),rep(0,8*12)),31)
Hepatitis.B.agegroup11$Case.updated<-round(Hepatitis.B.agegroup11$Case*(1+Inflation))
regions<-as.character(unique(Hepatitis.B.agegroup11$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup11[Hepatitis.B.agegroup11$Province==x, ])
seq(datalist)
yori7.1<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori7.1<-vector("list", length(datalist)); names(Sori7.1) <- regions
yori7.2<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori7.2<-vector("list", length(datalist)); names(Sori7.2) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order==i),]
mfirst<-gam(Case.updated~ offset(log(Population.book))+Intervention.2002+Intervention.2009+
              Time+ns(Temperature,df=3)+
              Holiday+Interaction.2002+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori7.1[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun-1),1]
Sori7.1[[i]] <- vcov(mfirst)[c(2,lengthfun-1),c(2,lengthfun-1)]
yori7.2[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(3,lengthfun),1]
Sori7.2[[i]] <- vcov(mfirst)[c(3,lengthfun),c(3,lengthfun)]
assign(paste("model7.",i,sep=""),mfirst)
}
mvall7.1<-mvmeta(yori7.1~1, Sori7.1, method="reml")
summary(mvall7.1)
mvall7.2<-mvmeta(yori7.2~1, Sori7.2, method="reml")
summary(mvall7.2)

#Age group 12
Hepatitis.B.agegroup12<-read.xlsx("agegroup12.xlsx")
Hepatitis.B.agegroup12$Month.factor<-factor(Hepatitis.B.agegroup12$Month)
Hepatitis.B.agegroup12$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup12$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup12$Intervention.2002<-factor(rep(c(rep(0,12*12),rep(1,12*3)),31))
Hepatitis.B.agegroup12$Intervention.2009<-factor(rep(c(rep(0,12*4),rep(0,10),rep(1,2+7*12),rep(0,12*3)),31))
Hepatitis.B.agegroup12$Interaction.2002<-rep(c(rep(0,12*12),c(0:(180-12*12-1))),31)
Hepatitis.B.agegroup12$Interaction.2009<-rep(c(rep(0,4*12),rep(0,10),c(0:(2+7*12-1)),rep(0,3*12)),31)
Hepatitis.B.agegroup12$Case.updated<-round(Hepatitis.B.agegroup12$Case*(1+Inflation))
regions<-as.character(unique(Hepatitis.B.agegroup12$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup12[Hepatitis.B.agegroup12$Province==x, ])
seq(datalist)
yori8.1<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori8.1<-vector("list", length(datalist)); names(Sori8.1) <- regions
yori8.2<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori8.2<-vector("list", length(datalist)); names(Sori8.2) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order==i),]
mfirst<-gam(Case.updated~ offset(log(Population.book))+Intervention.2002+Intervention.2009+
              Time+ns(Temperature,df=3)+
              Holiday+Interaction.2002+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori8.1[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun-1),1]
Sori8.1[[i]] <- vcov(mfirst)[c(2,lengthfun-1),c(2,lengthfun-1)]
yori8.2[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(3,lengthfun),1]
Sori8.2[[i]] <- vcov(mfirst)[c(3,lengthfun),c(3,lengthfun)]
assign(paste("model8.",i,sep=""),mfirst)
}
mvall8.1<-mvmeta(yori8.1~1, Sori8.1, method="reml")
summary(mvall8.1)
mvall8.2<-mvmeta(yori8.2~1, Sori8.2, method="reml")
summary(mvall8.2)

#Calculate the weights
Age.pop2<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup6[-which(Hepatitis.B.agegroup6$Order==1 |
                                                             Hepatitis.B.agegroup6$Order==2 |
                                                             Hepatitis.B.agegroup6$Order==6 |
                                                             Hepatitis.B.agegroup6$Order==9 |
                                                             Hepatitis.B.agegroup6$Order==26),],FUN=sum)[,2])
Age.pop3<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup7[-which(Hepatitis.B.agegroup7$Order==1 |
                                                             Hepatitis.B.agegroup7$Order==2 |
                                                             Hepatitis.B.agegroup7$Order==9 |
                                                             Hepatitis.B.agegroup7$Order==14 |
                                                             Hepatitis.B.agegroup7$Order==26),],FUN=sum)[,2])
Age.pop4<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup8[-which(Hepatitis.B.agegroup8$Order==1 |
                                                             Hepatitis.B.agegroup8$Order==2 |
                                                             Hepatitis.B.agegroup8$Order==9 |
                                                             Hepatitis.B.agegroup8$Order==14 |
                                                             Hepatitis.B.agegroup8$Order==26),],FUN=sum)[,2])
Age.pop5<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup9[-which(Hepatitis.B.agegroup9$Order==1 |
                                                             Hepatitis.B.agegroup9$Order==2 |
                                                             Hepatitis.B.agegroup9$Order==9 |
                                                             Hepatitis.B.agegroup9$Order==14 |
                                                             Hepatitis.B.agegroup9$Order==26),],FUN=sum)[,2])
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[-which(Hepatitis.B.agegroup10$Order==1 |
                                                              Hepatitis.B.agegroup10$Order==2 |
                                                              Hepatitis.B.agegroup10$Order==9 |
                                                              Hepatitis.B.agegroup10$Order==14 |
                                                              Hepatitis.B.agegroup10$Order==26),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,data=Hepatitis.B.agegroup11,FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,data=Hepatitis.B.agegroup12,FUN=sum)[,2])
Age.pop<-c(Age.pop2,Age.pop3,Age.pop4,Age.pop5,Age.pop6,Age.pop7,Age.pop8)
Prop.pop<-c(Age.pop2/sum(Age.pop),Age.pop3/sum(Age.pop),Age.pop4/sum(Age.pop),
            Age.pop5/sum(Age.pop),Age.pop6/sum(Age.pop),Age.pop7/sum(Age.pop),
            Age.pop8/sum(Age.pop))
Merge.coef<-Prop.pop[1]*mvall2$coefficients+Prop.pop[2]*mvall3$coefficients+
  Prop.pop[3]*mvall4$coefficients+Prop.pop[4]*mvall5.1$coefficients+
  Prop.pop[5]*mvall6.1$coefficients+Prop.pop[6]*mvall7.1$coefficients+
  Prop.pop[7]*mvall8.1$coefficients
Merge.vcov<-Prop.pop[1]^2*mvall2$vcov+Prop.pop[2]^2*mvall3$vcov+
  Prop.pop[3]^2*mvall4$vcov+Prop.pop[4]^2*mvall5.1$vcov+
  Prop.pop[5]^2*mvall6.1$vcov+Prop.pop[6]^2*mvall7.1$vcov+
  Prop.pop[7]^2*mvall8.1$vcov

#Age-standardized ER
coef.meta<-as.matrix(t(Merge.coef))
time.meta<-as.matrix(cbind(rep(1,180+3*12),c(1:(180+3*12))-1))
cov.meta<-Merge.vcov
pre.ER.meta<-c()
ER.meta<-c()
se.ER.meta<-c()
ER.low.meta<-c()
ER.high.meta<-c()
for(l in 1:(180+3*12))
{pre.ER.meta[l]<-time.meta[l,]%*%coef.meta
ER.meta[l]<-(exp(pre.ER.meta[l])-1)*100
se.ER.meta[l]<-sqrt(time.meta[l,]%*%cov.meta%*%t(t(time.meta[l,]))) 
ER.low.meta[l]=(exp(pre.ER.meta[l]-1.96*se.ER.meta[l])-1)*100
ER.high.meta[l]=(exp(pre.ER.meta[l]+1.96*se.ER.meta[l])-1)*100
}
ER.whole<-ER.meta
ER.whole.low<-ER.low.meta
ER.whole.upper<-ER.high.meta
month.meta<-c(1:(180+3*12))
ER<-data.frame(Time=month.meta,ER=ER.meta,ER.lower=ER.low.meta,ER.upper=ER.high.meta)
ER$ER<-sprintf("%0.2f",ER$ER)
ER$ER.lower<-sprintf("%0.2f",ER$ER.lower)
ER$ER.upper<-sprintf("%0.2f",ER$ER.upper)

#Age group 13
Hepatitis.B.agegroup13<-read.xlsx("Agegroup13.xlsx")
Hepatitis.B.agegroup13$Month.factor<-factor(Hepatitis.B.agegroup13$Month)
Hepatitis.B.agegroup13$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup13$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup13$Intervention.2009<-factor(rep(c(rep(0,9*12),rep(1,12*6)),31))
Hepatitis.B.agegroup13$Interaction.2009<-rep(c(rep(0,9*12),c(0:(180-12*9-1))),31)
Hepatitis.B.agegroup13$Case.updated<-round(Hepatitis.B.agegroup13$Case*(1+Inflation))
regions<-as.character(unique(Hepatitis.B.agegroup13$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup13[Hepatitis.B.agegroup13$Province==x, ])
seq(datalist)
yori9<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori9<-vector("list", length(datalist)); names(Sori9) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup13[which(Hepatitis.B.agegroup13$Order==i),]
mfirst<-gam(Case.updated~ offset(log(Population.book))+Intervention.2009+Time+
              ns(Temperature,df=3)+
              Holiday+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori9[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun),1]
Sori9[[i]] <- vcov(mfirst)[c(2,lengthfun),c(2,lengthfun)]
assign(paste("model9.",i,sep=""),mfirst)
}
mvall9<-mvmeta(yori9~1, Sori9, method="reml")
summary(mvall9)

#Age group 14
Hepatitis.B.agegroup14<-read.xlsx("Agegroup14.xlsx")
Hepatitis.B.agegroup14$Month.factor<-factor(Hepatitis.B.agegroup14$Month)
Hepatitis.B.agegroup14$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup14$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup14$Intervention.2009<-factor(rep(c(rep(0,14*12),rep(1,12*1)),31))
Hepatitis.B.agegroup14$Interaction.2009<-rep(c(rep(0,14*12),c(0:(180-12*14-1))),31)
Hepatitis.B.agegroup14$Case.updated<-round(Hepatitis.B.agegroup14$Case*(1+Inflation))
regions<-as.character(unique(Hepatitis.B.agegroup14$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup14[Hepatitis.B.agegroup14$Province==x, ])
seq(datalist)
yori10<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori10<-vector("list", length(datalist)); names(Sori10) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup14[which(Hepatitis.B.agegroup14$Order==i),]
mfirst<-gam(Case.updated~ offset(log(Population.book))+Intervention.2009+Time+
              ns(Temperature,df=3)+
              Holiday+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori10[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun),1]
Sori10[[i]] <- vcov(mfirst)[c(2,lengthfun),c(2,lengthfun)]
assign(paste("model10.",i,sep=""),mfirst)
}
mvall10<-mvmeta(yori10~1, Sori10, method="reml")
summary(mvall10)

#Calculate the weights
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[-which(Hepatitis.B.agegroup10$Order==1 |
                                                              Hepatitis.B.agegroup10$Order==2 |
                                                              Hepatitis.B.agegroup10$Order==9 |
                                                              Hepatitis.B.agegroup10$Order==26),],FUN=sum)[,2])
Age.pop9<-mean(aggregate(Population.book~Time,data=Hepatitis.B.agegroup13,FUN=sum)[,2])
Age.pop10<-mean(aggregate(Population.book~Time,data=Hepatitis.B.agegroup14,FUN=sum)[,2])
Age.pop.2009<-c(Age.pop6,Age.pop7,Age.pop8,Age.pop9,Age.pop10)
Prop.pop.2009<-c(Age.pop6/sum(Age.pop.2009),Age.pop7/sum(Age.pop.2009),Age.pop8/sum(Age.pop.2009),
                 Age.pop9/sum(Age.pop.2009),Age.pop10/sum(Age.pop.2009))
Merge.coef.2009<-Prop.pop.2009[1]*mvall6.2$coefficients+
  Prop.pop.2009[2]*mvall7.2$coefficients+Prop.pop.2009[3]*mvall8.2$coefficients+
  Prop.pop.2009[4]*mvall9$coefficients+Prop.pop.2009[5]*mvall10$coefficients
Merge.vcov.2009<-Prop.pop.2009[1]^2*mvall6.2$vcov+
  Prop.pop.2009[2]^2*mvall7.2$vcov+Prop.pop.2009[3]^2*mvall8.2$vcov+
  Prop.pop.2009[4]^2*mvall9$vcov+Prop.pop.2009[5]^2*mvall10$vcov

#Age-standardized ER
coef.meta<-as.matrix(t(Merge.coef.2009))
time.meta<-as.matrix(cbind(rep(1,180-7*12),c((7*12):179)-7*12))
cov.meta<-Merge.vcov.2009
pre.ER.meta<-c()
ER.meta<-c()
se.ER.meta<-c()
ER.low.meta<-c()
ER.high.meta<-c()
for(l in 1:(180-7*12))
{pre.ER.meta[l]<-time.meta[l,]%*%coef.meta
ER.meta[l]<-(exp(pre.ER.meta[l])-1)*100
se.ER.meta[l]<-sqrt(time.meta[l,]%*%cov.meta%*%t(t(time.meta[l,]))) 
ER.low.meta[l]=(exp(pre.ER.meta[l]-1.96*se.ER.meta[l])-1)*100
ER.high.meta[l]=(exp(pre.ER.meta[l]+1.96*se.ER.meta[l])-1)*100
}
ER.whole<-ER.meta
ER.whole.low<-ER.low.meta
ER.whole.upper<-ER.high.meta
month.meta<-c((7*12+1):180)
ER<-data.frame(Time=month.meta,ER=ER.meta,ER.lower=ER.low.meta,ER.upper=ER.high.meta)
ER$ER<-sprintf("%0.2f",ER$ER)
ER$ER.lower<-sprintf("%0.2f",ER$ER.lower)
ER$ER.upper<-sprintf("%0.2f",ER$ER.upper)

#######################The consideration of individuals aged 4y#################
#Age group 5
Hepatitis.B.agegroup5<-read.xlsx("Agegroup5.xlsx")
Hepatitis.B.agegroup5$Month.factor<-factor(Hepatitis.B.agegroup5$Month)
Hepatitis.B.agegroup5$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup5$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup5$Intervention.2002<-factor(rep(c(rep(0,12*1),rep(1,14*12)),31))
Hepatitis.B.agegroup5$Interaction.2002<-rep(c(rep(0,1*12),c(0:(180-1*12-1))),31)
regions<-as.character(unique(Hepatitis.B.agegroup5$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup5[Hepatitis.B.agegroup5$Province==x, ])
seq(datalist)
yori1<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori1<-vector("list", length(datalist)); names(Sori2) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup5[which(Hepatitis.B.agegroup5$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+Time+
              ns(Temperature,df=3)+
              Holiday+Interaction.2002,
            family=quasipoisson(link='log'),data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori1[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun),1]
Sori1[[i]] <- vcov(mfirst)[c(2,lengthfun),c(2,lengthfun)]
assign(paste("model1.",i,sep=""),mfirst)
}
mvall1<-mvmeta(yori1[-c(1,2,6,9,26,31),]~1, Sori1[-c(1,2,6,9,26,31)], method="reml")
summary(mvall1)

#Age group 6
Hepatitis.B.agegroup6<-read.xlsx("agegroup6.xlsx")
Hepatitis.B.agegroup6$Month.factor<-factor(Hepatitis.B.agegroup6$Month)
Hepatitis.B.agegroup6$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup6$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup6$Intervention.2002<-factor(rep(c(rep(0,12*2),rep(1,13*12)),31))
Hepatitis.B.agegroup6$Interaction.2002<-rep(c(rep(0,2*12),c(0:(180-2*12-1))),31)
regions<-as.character(unique(Hepatitis.B.agegroup6$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup6[Hepatitis.B.agegroup6$Province==x, ])
seq(datalist)
yori2<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori2<-vector("list", length(datalist)); names(Sori2) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup6[which(Hepatitis.B.agegroup6$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+Time+
              ns(Temperature,df=3)+
              Holiday+Interaction.2002,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori2[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun),1]
Sori2[[i]] <- vcov(mfirst)[c(2,lengthfun),c(2,lengthfun)]
assign(paste("model2.",i,sep=""),mfirst)
}
mvall2<-mvmeta(yori2[-c(1,2,6,9,26),]~1, Sori2[-c(1,2,6,9,26)], method="reml")
summary(mvall2)

#Age group 7
Hepatitis.B.agegroup7<-read.xlsx("agegroup7.xlsx")
Hepatitis.B.agegroup7$Month.factor<-factor(Hepatitis.B.agegroup7$Month)
Hepatitis.B.agegroup7$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup7$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup7$Intervention.2002<-factor(rep(c(rep(0,12*3),rep(1,12*12)),31))
Hepatitis.B.agegroup7$Interaction.2002<-rep(c(rep(0,3*12),c(0:(180-3*12-1))),31)
regions<-as.character(unique(Hepatitis.B.agegroup7$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup7[Hepatitis.B.agegroup7$Province==x, ])
seq(datalist)
yori3<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori3<-vector("list", length(datalist)); names(Sori3) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup7[which(Hepatitis.B.agegroup7$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+Time+
              ns(Temperature,df=3)+
              Holiday+Interaction.2002,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori3[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun),1]
Sori3[[i]] <- vcov(mfirst)[c(2,lengthfun),c(2,lengthfun)]
assign(paste("model3.",i,sep=""),mfirst)
}
mvall3<-mvmeta(yori3[-c(1,2,9,14,26),]~1, Sori3[-c(1,2,9,14,26)], method="reml")
summary(mvall3)

#Age group 8
Hepatitis.B.agegroup8<-read.xlsx("agegroup8.xlsx")
Hepatitis.B.agegroup8$Month.factor<-factor(Hepatitis.B.agegroup8$Month)
Hepatitis.B.agegroup8$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup8$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup8$Intervention.2002<-factor(rep(c(rep(0,12*4),rep(1,12*11)),31))
Hepatitis.B.agegroup8$Interaction.2002<-rep(c(rep(0,4*12),c(0:(180-4*12-1))),31)
regions<-as.character(unique(Hepatitis.B.agegroup8$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup8[Hepatitis.B.agegroup8$Province==x, ])
seq(datalist)
yori4<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori4<-vector("list", length(datalist)); names(Sori4) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup8[which(Hepatitis.B.agegroup8$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+Time+
              ns(Temperature,df=3)+
              Holiday+Interaction.2002,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori4[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun),1]
Sori4[[i]] <- vcov(mfirst)[c(2,lengthfun),c(2,lengthfun)]
assign(paste("model4.",i,sep=""),mfirst)
}
mvall4<-mvmeta(yori4[-c(1,2,9,14,26),]~1, Sori4[-c(1,2,9,14,26)], method="reml")
summary(mvall4)

#Age group 9
Hepatitis.B.agegroup9<-read.xlsx("agegroup9.xlsx")
Hepatitis.B.agegroup9$Month.factor<-factor(Hepatitis.B.agegroup9$Month)
Hepatitis.B.agegroup9$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup9$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup9$Intervention.2002<-factor(rep(c(rep(0,12*5),rep(1,12*10)),31))
Hepatitis.B.agegroup9$Intervention.2009<-factor(rep(c(rep(0,12*4),rep(0,10),rep(1,2),rep(0,12*10)),31))
Hepatitis.B.agegroup9$Interaction.2002<-rep(c(rep(0,5*12),c(0:(180-5*12-1))),31)
Hepatitis.B.agegroup9$Interaction.2009<-rep(c(rep(0,4*12),rep(0,10),c(0:(2-1)),rep(0,10*12)),31)
regions<-as.character(unique(Hepatitis.B.agegroup9$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup9[Hepatitis.B.agegroup9$Province==x, ])
seq(datalist)
yori5.1<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori5.1<-vector("list", length(datalist)); names(Sori5.1) <- regions
yori5.2<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori5.2<-vector("list", length(datalist)); names(Sori5.2) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup9[which(Hepatitis.B.agegroup9$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+Intervention.2009+
              Time+ns(Temperature,df=3)+
              Holiday+Interaction.2002+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori5.1[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun-1),1]
Sori5.1[[i]] <- vcov(mfirst)[c(2,lengthfun-1),c(2,lengthfun-1)]
yori5.2[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(3,lengthfun),1]
Sori5.2[[i]] <- vcov(mfirst)[c(3,lengthfun),c(3,lengthfun)]
assign(paste("model5.",i,sep=""),mfirst)
}
mvall5.1<-mvmeta(yori5.1[-c(1,2,9,14,26),]~1, Sori5.1[-c(1,2,9,14,26)], method="reml")
summary(mvall5.1)
mvall5.2<-mvmeta(yori5.2[-c(1,2,5,9,11,12,14,15,21,25,26),]~1, Sori5.2[-c(1,2,5,9,11,12,14,15,21,25,26)], method="reml")
summary(mvall5.2)

#Age group 10
Hepatitis.B.agegroup10<-read.xlsx("agegroup10.xlsx")
Hepatitis.B.agegroup10$Month.factor<-factor(Hepatitis.B.agegroup10$Month)
Hepatitis.B.agegroup10$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup10$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup10$Intervention.2002<-factor(rep(c(rep(0,12*6),rep(1,12*9)),31))
Hepatitis.B.agegroup10$Intervention.2009<-factor(rep(c(rep(0,12*4),rep(0,10),rep(1,2+12),rep(0,12*9)),31))
Hepatitis.B.agegroup10$Interaction.2002<-rep(c(rep(0,6*12),c(0:(180-6*12-1))),31)
Hepatitis.B.agegroup10$Interaction.2009<-rep(c(rep(0,4*12),rep(0,10),c(0:(2+12-1)),rep(0,9*12)),31)
regions<-as.character(unique(Hepatitis.B.agegroup10$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup10[Hepatitis.B.agegroup10$Province==x, ])
seq(datalist)
yori6.1<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori6.1<-vector("list", length(datalist)); names(Sori6.1) <- regions
yori6.2<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori6.2<-vector("list", length(datalist)); names(Sori6.2) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup10[which(Hepatitis.B.agegroup10$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+Intervention.2009+
              Time+ns(Temperature,df=3)+
              Holiday+Interaction.2002+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori6.1[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun-1),1]
Sori6.1[[i]] <- vcov(mfirst)[c(2,lengthfun-1),c(2,lengthfun-1)]
yori6.2[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(3,lengthfun),1]
Sori6.2[[i]] <- vcov(mfirst)[c(3,lengthfun),c(3,lengthfun)]
assign(paste("model6.",i,sep=""),mfirst)
}
mvall6.1<-mvmeta(yori6.1[-c(1,2,9,14,26),]~1, Sori6.1[-c(1,2,9,14,26)], method="reml")
summary(mvall6.1)
mvall6.2<-mvmeta(yori6.2[-c(1,2,9,14,26),]~1, Sori6.2[-c(1,2,9,14,26)], method="reml")
summary(mvall6.2)

#Age group 11
Hepatitis.B.agegroup11<-read.xlsx("agegroup11.xlsx")
Hepatitis.B.agegroup11$Month.factor<-factor(Hepatitis.B.agegroup11$Month)
Hepatitis.B.agegroup11$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup11$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup11$Intervention.2002<-factor(rep(c(rep(0,12*7),rep(1,12*8)),31))
Hepatitis.B.agegroup11$Intervention.2009<-factor(rep(c(rep(0,12*4),rep(0,10),rep(1,2+2*12),rep(0,12*8)),31))
Hepatitis.B.agegroup11$Interaction.2002<-rep(c(rep(0,7*12),c(0:(180-7*12-1))),31)
Hepatitis.B.agegroup11$Interaction.2009<-rep(c(rep(0,4*12),rep(0,10),c(0:(2+2*12-1)),rep(0,8*12)),31)
regions<-as.character(unique(Hepatitis.B.agegroup11$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup11[Hepatitis.B.agegroup11$Province==x, ])
seq(datalist)
yori7.1<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori7.1<-vector("list", length(datalist)); names(Sori7.1) <- regions
yori7.2<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori7.2<-vector("list", length(datalist)); names(Sori7.2) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup11[which(Hepatitis.B.agegroup11$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+Intervention.2009+
              Time+ns(Temperature,df=3)+
              Holiday+Interaction.2002+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori7.1[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun-1),1]
Sori7.1[[i]] <- vcov(mfirst)[c(2,lengthfun-1),c(2,lengthfun-1)]
yori7.2[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(3,lengthfun),1]
Sori7.2[[i]] <- vcov(mfirst)[c(3,lengthfun),c(3,lengthfun)]
assign(paste("model7.",i,sep=""),mfirst)
}
mvall7.1<-mvmeta(yori7.1~1, Sori7.1, method="reml")
summary(mvall7.1)
mvall7.2<-mvmeta(yori7.2~1, Sori7.2, method="reml")
summary(mvall7.2)

#Age group 12
Hepatitis.B.agegroup12<-read.xlsx("agegroup12.xlsx")
Hepatitis.B.agegroup12$Month.factor<-factor(Hepatitis.B.agegroup12$Month)
Hepatitis.B.agegroup12$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup12$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup12$Intervention.2002<-factor(rep(c(rep(0,12*12),rep(1,12*3)),31))
Hepatitis.B.agegroup12$Intervention.2009<-factor(rep(c(rep(0,12*4),rep(0,10),rep(1,2+7*12),rep(0,12*3)),31))
Hepatitis.B.agegroup12$Interaction.2002<-rep(c(rep(0,12*12),c(0:(180-12*12-1))),31)
Hepatitis.B.agegroup12$Interaction.2009<-rep(c(rep(0,4*12),rep(0,10),c(0:(2+7*12-1)),rep(0,3*12)),31)
regions<-as.character(unique(Hepatitis.B.agegroup12$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup12[Hepatitis.B.agegroup12$Province==x, ])
seq(datalist)
yori8.1<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori8.1<-vector("list", length(datalist)); names(Sori8.1) <- regions
yori8.2<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori8.2<-vector("list", length(datalist)); names(Sori8.2) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup12[which(Hepatitis.B.agegroup12$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2002+Intervention.2009+
              Time+ns(Temperature,df=3)+
              Holiday+Interaction.2002+Interaction.2009,
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori8.1[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun-1),1]
Sori8.1[[i]] <- vcov(mfirst)[c(2,lengthfun-1),c(2,lengthfun-1)]
yori8.2[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(3,lengthfun),1]
Sori8.2[[i]] <- vcov(mfirst)[c(3,lengthfun),c(3,lengthfun)]
assign(paste("model8.",i,sep=""),mfirst)
}
mvall8.1<-mvmeta(yori8.1~1, Sori8.1, method="reml")
summary(mvall8.1)
mvall8.2<-mvmeta(yori8.2~1, Sori8.2, method="reml")
summary(mvall8.2)

#Calculate the weights
Age.pop1<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup5[-which(Hepatitis.B.agegroup5$Order==1 |
                                                             Hepatitis.B.agegroup5$Order==2 |
                                                             Hepatitis.B.agegroup5$Order==6 |
                                                             Hepatitis.B.agegroup5$Order==9 |
                                                             Hepatitis.B.agegroup5$Order==26|
                                                             Hepatitis.B.agegroup5$Order==31),],FUN=sum)[,2])
Age.pop2<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup6[-which(Hepatitis.B.agegroup6$Order==1 |
                                                             Hepatitis.B.agegroup6$Order==2 |
                                                             Hepatitis.B.agegroup6$Order==6 |
                                                             Hepatitis.B.agegroup6$Order==9 |
                                                             Hepatitis.B.agegroup6$Order==26),],FUN=sum)[,2])
Age.pop3<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup7[-which(Hepatitis.B.agegroup7$Order==1 |
                                                             Hepatitis.B.agegroup7$Order==2 |
                                                             Hepatitis.B.agegroup7$Order==9 |
                                                             Hepatitis.B.agegroup7$Order==14 |
                                                             Hepatitis.B.agegroup7$Order==26),],FUN=sum)[,2])
Age.pop4<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup8[-which(Hepatitis.B.agegroup8$Order==1 |
                                                             Hepatitis.B.agegroup8$Order==2 |
                                                             Hepatitis.B.agegroup8$Order==9 |
                                                             Hepatitis.B.agegroup8$Order==14 |
                                                             Hepatitis.B.agegroup8$Order==26),],FUN=sum)[,2])
Age.pop5<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup9[-which(Hepatitis.B.agegroup9$Order==1 |
                                                             Hepatitis.B.agegroup9$Order==2 |
                                                             Hepatitis.B.agegroup9$Order==9 |
                                                             Hepatitis.B.agegroup9$Order==14 |
                                                             Hepatitis.B.agegroup9$Order==26),],FUN=sum)[,2])
Age.pop6<-mean(aggregate(Population.book~Time,
                         data=Hepatitis.B.agegroup10[-which(Hepatitis.B.agegroup10$Order==1 |
                                                              Hepatitis.B.agegroup10$Order==2 |
                                                              Hepatitis.B.agegroup10$Order==9 |
                                                              Hepatitis.B.agegroup10$Order==14 |
                                                              Hepatitis.B.agegroup10$Order==26),],FUN=sum)[,2])
Age.pop7<-mean(aggregate(Population.book~Time,data=Hepatitis.B.agegroup11,FUN=sum)[,2])
Age.pop8<-mean(aggregate(Population.book~Time,data=Hepatitis.B.agegroup12,FUN=sum)[,2])
Age.pop<-c(Age.pop1,Age.pop2,Age.pop3,Age.pop4,Age.pop5,Age.pop6,Age.pop7,Age.pop8)
Prop.pop<-c(Age.pop1/sum(Age.pop),Age.pop2/sum(Age.pop),Age.pop3/sum(Age.pop),
            Age.pop4/sum(Age.pop),Age.pop5/sum(Age.pop),Age.pop6/sum(Age.pop),
            Age.pop7/sum(Age.pop),Age.pop8/sum(Age.pop))
Merge.coef<-Prop.pop[1]*mvall1$coefficients+Prop.pop[2]*mvall2$coefficients+
  Prop.pop[3]*mvall3$coefficients+Prop.pop[4]*mvall4$coefficients+
  Prop.pop[5]*mvall5.1$coefficients+Prop.pop[6]*mvall6.1$coefficients+
  Prop.pop[7]*mvall7.1$coefficients+Prop.pop[8]*mvall8.1$coefficients
Merge.vcov<-Prop.pop[1]^2*mvall1$vcov+Prop.pop[2]^2*mvall2$vcov+
  Prop.pop[3]^2*mvall3$vcov+Prop.pop[4]^2*mvall4$vcov+
  Prop.pop[5]^2*mvall5.1$vcov+Prop.pop[6]^2*mvall6.1$vcov+
  Prop.pop[7]^2*mvall7.1$vcov+Prop.pop[8]^2*mvall8.1$vcov

coef.meta<-as.matrix(t(Merge.coef))
time.meta<-as.matrix(cbind(rep(1,180+3*12),c(1:(180+3*12))-1))
cov.meta<-Merge.vcov
pre.ER.meta<-c()
ER.meta<-c()
se.ER.meta<-c()
ER.low.meta<-c()
ER.high.meta<-c()
for(l in 1:(180+3*12))
{pre.ER.meta[l]<-time.meta[l,]%*%coef.meta
ER.meta[l]<-(exp(pre.ER.meta[l])-1)*100
se.ER.meta[l]<-sqrt(time.meta[l,]%*%cov.meta%*%t(t(time.meta[l,]))) 
ER.low.meta[l]=(exp(pre.ER.meta[l]-1.96*se.ER.meta[l])-1)*100
ER.high.meta[l]=(exp(pre.ER.meta[l]+1.96*se.ER.meta[l])-1)*100
}
ER.whole<-ER.meta
ER.whole.low<-ER.low.meta
ER.whole.upper<-ER.high.meta
month.meta<-c(1:(180+3*12))
ER<-data.frame(Time=month.meta,ER=ER.meta,ER.lower=ER.low.meta,ER.upper=ER.high.meta)
ER$ER<-sprintf("%0.2f",ER$ER)
ER$ER.lower<-sprintf("%0.2f",ER$ER.lower)
ER$ER.upper<-sprintf("%0.2f",ER$ER.upper)

#Age group 13
Hepatitis.B.agegroup13<-read.xlsx("Agegroup13.xlsx")
Hepatitis.B.agegroup13$Month.factor<-factor(Hepatitis.B.agegroup13$Month)
Hepatitis.B.agegroup13$Time<-rep(c(0:179),31)
Hepatitis.B.agegroup13$Holiday<-factor(rep(rep(c(1,1,0,0,0,0,2,2,0,0,0,0),15),31))
Hepatitis.B.agegroup13$Intervention.2009<-factor(rep(c(rep(0,9*12),rep(1,12*6)),31))
Hepatitis.B.agegroup13$Interaction.2009<-rep(c(rep(0,9*12),c(0:(180-12*9-1))),31)
regions<-as.character(unique(Hepatitis.B.agegroup13$Province))
datalist<-lapply(regions, function(x) Hepatitis.B.agegroup13[Hepatitis.B.agegroup13$Province==x, ])
seq(datalist)
yori9<-matrix(0, length(datalist), 2, dimnames=list(regions, paste("beta", seq(2,3), sep="")))
Sori9<-vector("list", length(datalist)); names(Sori9) <- regions
for(i in c(1:31))
{sub.data<-Hepatitis.B.agegroup13[which(Hepatitis.B.agegroup13$Order==i),]
mfirst<-gam(Case~ offset(log(Population.book))+Intervention.2009+Time+
              ns(Temperature,df=3)+
              #ns(Rain,df=4)+ns(Wind,df=4)+
              #ns(Sunshine,df=4)+
              Holiday+Interaction.2009,#+Holiday+GDP+Birth_Rate
            family=quasipoisson(link='log'), data=sub.data)
lengthfun<-length(summary(mfirst)$p.coeff)
yori9[i,c(1:2)] <- as.data.frame(summary(mfirst)$p.coeff)[c(2,lengthfun),1]
Sori9[[i]] <- vcov(mfirst)[c(2,lengthfun),c(2,lengthfun)]
assign(paste("model9.",i,sep=""),mfirst)
}
mvall9<-mvmeta(yori9~1, Sori9, method="reml")
summary(mvall9)
