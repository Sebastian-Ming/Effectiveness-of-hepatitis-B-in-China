library(data.table)
library(tidyverse)
library(mvmeta)
library(patchwork)


# -------------------------
# Function Definitions
# -------------------------


# Calculate heterogeneity index I² and Q-test p-value from mvmeta object
get_I2 <- function(x){
  qtest_result <- mvmeta::qtest(x)
  I2 <- as.numeric((qtest_result$Q[1]-qtest_result$df[1])/qtest_result$Q[1])
  I2 <- paste0(sprintf("%.1f",I2*100), "%")
  p.value <- as.numeric(qtest_result$pvalue[1])
  p.value <- sprintf("%.3f",p.value)
  data.frame(I2=I2, p.value=p.value)
}

# Perform meta-regression and output I² and related results
# c: indices of rows to be excluded; num: index of the age group
do_mere <- function(information){
  num <- information[1]
  c <- information[2]
  if (c != " ") {
    groups_used <- eval(parse(text = paste0("groups[-", c, ",]")))
    eval(parse(text = paste0("yori_used <- yori", num, "[-", c, ",]")))
    eval(parse(text = paste0("Sori_used <- Sori", num, "[-", c, "]")))
  }
  else {
    groups_used <- groups
    yori_used <- eval(parse(text = paste0("yori",num)))
    Sori_used <- eval(parse(text = paste0("Sori",num)))
  }
  groups_used <- data.frame(lapply(groups_used ,as.factor))
  
  mere <-list()
  for(i in 2:6){mere[[i-1]] <- mvmeta(yori_used ~ groups_used[[i]], Sori_used, method="reml")}
  mere[[6]] <-  mvmeta(yori_used ~ 1, Sori_used, method="reml")  # meta regression
  mere %>%
    lapply(get_I2) %>%
    rbindlist() %>%
    mutate(term = c(names(groups)[-1],"bind"),
           age_group = num)
}

# Based on meta-regression results, predict the effect estimates 
# under different levels of covariates
do_meta <- function(information){
  num <- information[1]
  c <- information[2]
  if (c != " ") {
    groups_used <- eval(parse(text = paste0("groups[-", c, ",]")))
    eval(parse(text = paste0("yori_used <- yori", num, "[-", c, ",]")))
    eval(parse(text = paste0("Sori_used <- Sori", num, "[-", c, "]")))
  }
  else {
    groups_used <- groups
    yori_used <- eval(parse(text = paste0("yori",num)))
    Sori_used <- eval(parse(text = paste0("Sori",num)))
  }
  groups_used <- data.frame(lapply(groups_used ,as.factor))
  
  mere <- predict_list <- list()
  for(i in 2:6){
    levels_used <- groups_used[[i]]
    mere[[i-1]] <- mvmeta(yori_used ~ levels_used, Sori_used, method="reml")
    levels <- mere[[i-1]]$xlevels
    newdata <- data.frame(levels_used=factor(levels[[1]])) 
    predict_list[[i-1]] <- predict(mere[[i-1]],newdata,vcov = T)
  }
  predict_list
}

#setwd("C:/~")

# ------------------------------
# Data Loading and Preprocessing
# ------------------------------

groups <- read.csv("groups.csv")

# Select relevant columns for analysis
groups <- groups %>% dplyr::select(Name, HBsAgB, Geo, Ill, Bed, Ur)

# Age groups to be used
num <-c(1,2,3,4,5.1,5.2,6.1,6.2,7.1,7.2,8.1,8.2,9,10)

# For each age group, specify rows (provinces) to exclude according to previous criteria
# Rows are specified as character strings of R vectors, e.g., "c(1,2,6,9,26)"
delete_rows <- c("c(1,2,6,9,26)", "c(1,2,6,9,26)", "c(1,2,9,14,26)", "c(1,2,9,14,26)", "c(1,2,9,14,26)", "c(1,2,9,14,26)", "c(1,2,9,14,26)", "c(1,2,9,14,26)", 
                 " "," ", " "," "," "," ")

num1 <- c(1,2,3,4,5.1,6.1,7.1,8.1) # Age groups for 2002 NIP
num2 <- c(5.2,6.2,7.2,8.2,9,10) # Age groups for 2009 catch-up vaccine

delete_rows1 <- c("c(1,2,6,9,26)", "c(1,2,6,9,26)", "c(1,2,9,14,26)", "c(1,2,9,14,26)", "c(1,2,9,14,26)", "c(1,2,9,14,26)", " "," ")
delete_rows2 <- c("c(1,2,9,14,26)", "c(1,2,9,14,26)",  " "," ", " "," ")

# Intervention years corresponding to groups in each set
intervention_years1 <- c(14,13,12,11,10,9,8,3)
intervention_years2 <- c(1,2,3,8,6,1)

information <- cbind(num, delete_rows)
information1 <- cbind(num1, delete_rows1, intervention_years1)
information2 <- cbind(num2, delete_rows2, intervention_years2)


# Define comparison group indices and grouping numbers for analysis
compare1 <- c(1,1,2,1,1,1,1,1,2,2,2,2,3,3,3,4,4,5,1,1,1)
compare2 <- c(2,3,3,2,3,4,5,6,3,4,5,6,4,5,6,5,6,6,2,2,2)
group_num <- c(1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,4,5)
compare <- cbind(compare1, compare2, group_num)

# ----------------------------------------------------
# Calculate population proportions by age group & province
# ----------------------------------------------------

population1_list <- list()
for(i in 5:14){
  eval(parse(text=paste0("population1_list[[",i,"]] <- Hepatitis.B.agegroup", i, " %>% group_by(Province) %>% summarise(sum(Population.book)) %>% mutate(age_group=",i,")")))
}

population11 <- rbindlist(population1_list) %>% filter(age_group %in% c(6:12)) 
population12 <- rbindlist(population1_list) %>% filter(age_group %in% c(10:14)) 
population11 <- pivot_wider(population11,names_from = Province,values_from = `sum(Population.book)`)
population12 <- pivot_wider(population12,names_from = Province,values_from = `sum(Population.book)`)
population11_allages <- apply(population11,2,sum)
population12_allages <- apply(population11,2,sum)


population11_groups <- population11 %>% dplyr::select(all_of(Name)) %>% as.matrix()
population12_groups <- population12 %>% dplyr::select(all_of(Name)) %>% as.matrix()

for(i in 1:31){
  population11[,i+1]<-population11[,(i+1)]/population11_allages[i+1]
  population12[,i+1]<-population12[,(i+1)]/population12_allages[i+1]
}

population11 <- population11 %>% dplyr::select(all_of(Name)) %>% as.matrix()
population12 <- population12 %>% dplyr::select(all_of(Name)) %>% as.matrix()

population11_groups <- cbind(as.data.frame(t(population11_groups)),groups[,-1])
population12_groups <- cbind(as.data.frame(t(population12_groups)),groups[,-1])

# Calculate weights
population_weight1 <- apply(population11_groups,2,sum)[1:7]/sum(apply(population11_groups,2,sum)[1:7])
population_weight2 <- apply(population12_groups,2,sum)[1:5]/sum(apply(population12_groups,2,sum)[1:5])

population11_groups_HBsAgB <- population11_groups %>% group_by(HBsAgB) %>% summarise(sum(V1),sum(V2),sum(V3),sum(V4),sum(V5),sum(V6),sum(V7))
population11_groups_Geo <- population11_groups %>% group_by(Geo) %>% summarise(sum(V1),sum(V2),sum(V3),sum(V4),sum(V5),sum(V6),sum(V7))
population11_groups_Ill <- population11_groups %>% group_by(Ill) %>% summarise(sum(V1),sum(V2),sum(V3),sum(V4),sum(V5),sum(V6),sum(V7))
population11_groups_Bed <- population11_groups %>% group_by(Bed) %>% summarise(sum(V1),sum(V2),sum(V3),sum(V4),sum(V5),sum(V6),sum(V7))
population11_groups_Ur <- population11_groups %>% group_by(Ur) %>% summarise(sum(V1),sum(V2),sum(V3),sum(V4),sum(V5),sum(V6),sum(V7))
population11_groups <- rbind(population11_groups_HBsAgB[,-1],population11_groups_Geo[,-1],population11_groups_Ill[,-1],population11_groups_Bed[,-1],population11_groups_Ur[,-1])
population11_groups$sum <- apply(population11_groups,1,sum)
population11_groups_weight <- apply(population11_groups[,1:7],2,function(x) x/population11_groups$sum)

population12_groups_HBsAgB <- population12_groups %>% group_by(HBsAgB) %>% summarise(sum(V1),sum(V2),sum(V3),sum(V4),sum(V5))
population12_groups_Geo <- population12_groups %>% group_by(Geo) %>% summarise(sum(V1),sum(V2),sum(V3),sum(V4),sum(V5))
population12_groups_Ill <- population12_groups %>% group_by(Ill) %>% summarise(sum(V1),sum(V2),sum(V3),sum(V4),sum(V5))
population12_groups_Bed <- population12_groups %>% group_by(Bed) %>% summarise(sum(V1),sum(V2),sum(V3),sum(V4),sum(V5))
population12_groups_Ur <- population12_groups %>% group_by(Ur) %>% summarise(sum(V1),sum(V2),sum(V3),sum(V4),sum(V5))
population12_groups <- rbind(population12_groups_HBsAgB[,-1],population12_groups_Geo[,-1],population12_groups_Ill[,-1],population12_groups_Bed[,-1],population12_groups_Ur[,-1])
population12_groups$sum <- apply(population12_groups,1,sum)
population12_groups_weight <- apply(population12_groups[,1:5],2,function(x) x/population12_groups$sum)

# ----------------------------------
# Construct weighted age group datasets
# ----------------------------------
compare_weight1 <- cbind(as.data.frame(population11_groups_weight), levels=c(1,2,3,0,1,2,3,4,5,1,2,1,2,1,2),groups=c(1,1,1,2,2,2,2,2,2,3,3,4,4,5,5))
compare_weight2 <- cbind(as.data.frame(population12_groups_weight), levels=c(1,2,3,0,1,2,3,4,5,1,2,1,2,1,2),groups=c(1,1,1,2,2,2,2,2,2,3,3,4,4,5,5))

# ------------------------
# Meta-regression and subgroup analyses
# ------------------------

# Run meta-regression for each row of information
mere_result_list <- apply(information[-1,], 1, do_mere)
mere_result <- rbindlist(mere_result_list) 

# Use meta-regression results to predict coefficients at different covariate levels
apply(information[-1,], 1, do_meta) -> model_list 


read.csv("groups_used_plot.csv")->groups_used_plot

groups_used_plot <- groups_used_plot %>% dplyr::select(Name,HBsAgB,Geo,Ill,Bed,Ur)
groups_used_plot <- groups_used_plot %>% mutate(HBsAgB = factor(HBsAgB, levels = c("Low", "Middle", "High")),
                                                Geo = factor(Geo, levels=c("South-centre", "East", "South West", "North", "North West", "North East")),
                                                Ill = factor(Ill, levels=c("Low", "High")),
                                                Bed = factor(Bed, levels=c("Low", "High")),
                                                Ur = factor(Ur, levels=c("Low", "High")),
)


compare_weight1$levels <- c("Low","Middle","High","South-centre", "East", "South West", "North", "North West", "North East", "Low", "High", "Low", "High", "Low", "High" )
compare_weight2$levels <- c("Low","Middle","High","South-centre", "East", "South West", "North", "North West", "North East", "Low", "High", "Low", "High", "Low", "High" )

# --------------------------
# Calculate the age-standardize coefficients and variance-covariance matrices for subgroups
# --------------------------

# For 2002 intervention
coef1 <-vcov1 <- list()

for (j in 2:6){
  coef1[[j-1]]<- vcov1[[j-1]] <-list()
  levels <- levels(groups_used_plot[,j]) 
  for (level in levels){
    num_level <- which(levels == level)
    coef1[[j-1]][[level]] <- model_list[[1]][[j-1]][[num_level]]$fit * compare_weight1[which(compare_weight1$levels==level & compare_weight1$groups==j-1),1] 
    vcov1[[j-1]][[level]] <- model_list[[1]][[j-1]][[num_level]]$vcov * compare_weight1[which(compare_weight1$levels==level & compare_weight1$groups==j-1),1]  * compare_weight1[which(compare_weight1$levels==level & compare_weight1$groups==j-1),1] 
  }
}
for (i in 2:7){
  num <- c(2,3,4,6,8,10)[i-1]
  for (j in 2:6){
    levels <- levels(groups_used_plot[,j]) 
    for (level in levels){
      num_level <- which(levels == level)
      coef1[[j-1]][[level]] <- coef1[[j-1]][[level]] + model_list[[num]][[j-1]][[num_level]]$fit * compare_weight1[which(compare_weight1$levels==level & compare_weight1$groups==j-1),i]  
      vcov1[[j-1]][[level]] <- vcov1[[j-1]][[level]] + model_list[[num]][[j-1]][[num_level]]$vcov * compare_weight1[which(compare_weight1$levels==level & compare_weight1$groups==j-1),i] * compare_weight1[which(compare_weight1$levels==level & compare_weight1$groups==j-1),i] 
    }
  }
}

# For 2009 intervention
coef2 <-vcov2 <- list()

for (j in 2:6){
  coef2[[j-1]]<- vcov2[[j-1]] <-list()
  levels <- levels(groups_used_plot[,j]) 
  for (level in levels){
    num_level <- which(levels == level)
    coef2[[j-1]][[level]] <- model_list[[7]][[j-1]][[num_level]]$fit * compare_weight2[which(compare_weight2$levels==level & compare_weight2$groups==j-1),1] 
    vcov2[[j-1]][[level]] <- model_list[[7]][[j-1]][[num_level]]$vcov * compare_weight2[which(compare_weight2$levels==level & compare_weight2$groups==j-1),1] * compare_weight2[which(compare_weight2$levels==level & compare_weight2$groups==j-1),1]
  }
}
for (i in 2:5){
  num <- c(9,11,12,13)[i-1]
  for (j in 2:6){
    levels <- levels(groups_used_plot[,j]) 
    for (level in levels){
      num_level <- which(levels == level)
      coef2[[j-1]][[level]] <- coef2[[j-1]][[level]] + model_list[[num]][[j-1]][[num_level]]$fit * compare_weight2[which(compare_weight2$levels==level & compare_weight2$groups==j-1),i] 
      vcov2[[j-1]][[level]] <- vcov2[[j-1]][[level]] + model_list[[num]][[j-1]][[num_level]]$vcov * compare_weight2[which(compare_weight2$levels==level & compare_weight2$groups==j-1),i] * compare_weight2[which(compare_weight2$levels==level & compare_weight2$groups==j-1),i]
    }
  }
}

# --------------------------
# Age-standardize ER values for subgroups
# --------------------------

plot1<-list()

for (j in 2:6){
  sER_bind1<-data.frame(NULL)
  levels <- levels(groups_used_plot[,j]) 
  num_of_level <- 0
  for (level in levels){
    num_of_level <- num_of_level + 1
    time <- as.matrix(cbind(rep(1,204),c(0:203)))
    spre_ER1 <- time %*% matrix(coef1[[j-1]][[num_of_level]], 2)
    sER_point1 <- (exp(spre_ER1)-1)*100
    sER_se1 <- sqrt(diag(time %*% matrix(vcov1[[j-1]][[num_of_level]],2) %*% t(time)))
    sER_low1 <- (exp(spre_ER1-1.96*sER_se1)-1)*100
    sER_high1 <- (exp(spre_ER1+1.96*sER_se1)-1)*100
    month <- c(0:203)   
    sER1 <- data.frame(Time=month,ER=sER_point1,ER.lower=sER_low1,ER.upper=sER_high1)
    sER1$ER<-sprintf("%0.2f",sER1$ER) %>% as.numeric()
    sER1$ER.lower<-sprintf("%0.2f",sER1$ER.lower) %>% as.numeric()
    sER1$ER.upper<-sprintf("%0.2f",sER1$ER.upper) %>% as.numeric()
    sER1$Years <- as.integer(sER1$Time/12)+1
    sER1$level <- rep(level, 204) %>% as.factor()
    sER_bind1 <- rbind(sER_bind1, sER1)
  }
  
  if(j==2)  sER_bind1$level <- factor(sER_bind1$level, levels = c("Low", "Middle", "High"))
  else if (j==3) sER_bind1$level <- factor(sER_bind1$level, levels = c("South-centre", "East", "South West", "North", "North West", "North East"))
  else sER_bind1$level <- factor(sER_bind1$level, levels = c("Low", "High"))
  
  # sER_bind1 <- sER_bind1 %>% filter(Time %in% c(11,107,203)) %>% mutate(CI = paste0(ER, " (", ER.lower, ", ", ER.upper, ")"))

  if(j==3){
    plot1[[2]] <- ggplot(sER_bind1,aes(x=Time,color=level,group=level)) +
      geom_ribbon(aes(ymin=ER.lower,ymax=ER.upper,group=level,fill=level),alpha=0.2,linetype = 2,color = NA )+
      geom_line(aes(y=ER,color=level,group=level),linewidth=0.5)+
      scale_x_continuous(breaks=seq(from=0,to=215,by=12), labels = 0:17)+
      scale_y_continuous(limits=c(-100,800),breaks=c(-100,0,100,200,300,400,500,600,700,800))+
      xlab("Years after the implementation of the NIP")+ylab("Excess risk (%)")+
      theme_bw()+
      theme(
        #panel.background = element_rect(color='black'),
        plot.margin = margin(10,5,5,5,unit = 'mm'))+
      scale_color_manual(values = c("#931e18","#247d3f","#20235b","#04a3bd", "#f0be3d",  "#da7901"))+
      scale_fill_manual(values = c("#931e18","#247d3f","#20235b","#04a3bd", "#f0be3d",  "#da7901"))+
      geom_hline(yintercept=0,lty=2)+
      labs(color = NULL, fill = NULL) + 
      theme(
        text = element_text(family = "Times New Roman"),  
        legend.title = element_text(size = 14,family = "Times New Roman"),
        legend.text = element_text(size = 14,family = "Times New Roman"),
        axis.title = element_text(size = 18,family = "Times New Roman"),
        axis.text = element_text(size = 18, family = "Times New Roman", color = "black"), 
        plot.title = element_text(size = 14,family = "Times New Roman"),
        plot.subtitle = element_text(size = 14,family = "Times New Roman"),
        plot.tag = element_text(size = 24, family = "Times New Roman"),
        plot.caption = element_text(size = 14,family = "Times New Roman")
      )
  }
  else{
    plot1[[j-1]] <- ggplot(sER_bind1,aes(x=Time,color=level,group=level)) +
      geom_ribbon(aes(ymin=ER.lower,ymax=ER.upper,group=level,fill=level),alpha=0.2,linetype = 2,color = NA )+
      geom_line(aes(y=ER,color=level,group=level),linewidth=0.5)+
      scale_x_continuous(breaks=seq(from=0,to=215,by=12), labels = 0:17)+
      scale_y_continuous(limits=c(-100,100),breaks=c(-100,-50,0,50,100))+
      xlab("Years after the implementation of the NIP")+ylab("Excess risk (%)")+
      theme_bw()+
      theme(
      plot.margin = margin(10,5,5,5,unit = 'mm'))+
      scale_color_manual(values = c("#931e18","#247d3f","#20235b","#04a3bd", "#f0be3d",  "#da7901"))+
      scale_fill_manual(values = c("#931e18","#247d3f","#20235b","#04a3bd", "#f0be3d",  "#da7901"))+
      geom_hline(yintercept=0,lty=2)+
      labs(color = NULL, fill = NULL) + 
      theme(
        text = element_text(family = "Times New Roman"),  
        legend.title = element_text(size = 14,family = "Times New Roman"),
        legend.text = element_text(size = 14,family = "Times New Roman"),
        axis.title = element_text(size = 18,family = "Times New Roman"),
        axis.text = element_text(size = 18, family = "Times New Roman", color = "black"),  
        plot.title = element_text(size = 14,family = "Times New Roman"),
        plot.subtitle = element_text(size = 14,family = "Times New Roman"),
        plot.tag = element_text(size = 24, family = "Times New Roman"),
        plot.caption = element_text(size = 14,family = "Times New Roman")
      )
  }
  if (j == 2) {
    print(plot)
  } else {
    print(plot, new = TRUE)
  }
}

plot2<-list()

for (j in 2:6){
  sER_bind2<-data.frame(NULL)
  levels <- levels(groups_used_plot[,j]) 
  num_of_level <- 0
  for (level in levels){
    num_of_level <- num_of_level +1
    time <- as.matrix(cbind(rep(1,84),c(0:83)))
    spre_ER2 <- time %*% matrix(coef2[[j-1]][[num_of_level]], 2)
    sER_point2 <- (exp(spre_ER2)-1)*100
    sER_se2 <- sqrt(diag(time %*% matrix(vcov2[[j-1]][[num_of_level]],2) %*% t(time)))
    sER_low2 <- (exp(spre_ER2-1.96*sER_se2)-1)*100
    sER_high2 <- (exp(spre_ER2+1.96*sER_se2)-1)*100
    month <- c(0:83)
    sER2 <- data.frame(Time=month,ER=sER_point2,ER.lower=sER_low2,ER.upper=sER_high2)
    sER2$ER<-sprintf("%0.2f",sER2$ER) %>% as.numeric()
    sER2$ER.lower<-sprintf("%0.2f",sER2$ER.lower) %>% as.numeric()
    sER2$ER.upper<-sprintf("%0.2f",sER2$ER.upper) %>% as.numeric()
    sER2$Years <- as.integer(sER2$Time/12)+1
    sER2$level <- rep(level, 84) %>% as.factor()
    sER_bind2 <- rbind(sER_bind2, sER2)
  }
  if(j==2)  sER_bind2$level <- factor(sER_bind2$level, levels = c("Low", "Middle", "High"))
  else if (j==3) sER_bind2$level <- factor(sER_bind2$level, levels = c("South-centre", "East", "South West", "North", "North West", "North East"))
  else sER_bind2$level <- factor(sER_bind2$level, levels = c("Low", "High"))
  
  #sER_bind2 <- sER_bind2 %>% filter(Time %in% c(11,47,83)) %>% mutate(CI = paste0(ER, " (", ER.lower, ", ", ER.upper, ")"))
  #write.csv(sER_bind2,paste0("2009亚组ER值",j-1,".csv"),row.names = F)
  
  plot2[[j-1]] <- ggplot(sER_bind2,aes(x=Time,color=level,group=level)) +
    geom_ribbon(aes(ymin=ER.lower,ymax=ER.upper,group=level,fill=level),alpha=0.2,linetype = 2,color = NA )+
    geom_line(aes(y=ER,color=level,group=level),linewidth=0.5)+
    scale_x_continuous(breaks=seq(from=0,to=95,by=12), labels = 0:7)+
    scale_y_continuous(limits=c(-100,100),breaks=c(-100,-50,0,50,100))+
    xlab("Years after the catch-up vaccination")+ylab("Excess risk (%)")+
    theme_bw()+
    theme(
      plot.margin = margin(10,5,5,5,unit = 'mm'))+
    scale_color_manual(values = c("#931e18","#247d3f","#20235b","#04a3bd", "#f0be3d",  "#da7901"))+
    scale_fill_manual(values = c("#931e18","#247d3f","#20235b","#04a3bd", "#f0be3d",  "#da7901"))+
    geom_hline(yintercept=0,lty=2)+
    labs(color = NULL, fill = NULL)+ 
    theme(
      text = element_text(family = "Times New Roman"),  
      legend.title = element_text(size = 14,family = "Times New Roman"),
      legend.text = element_text(size = 14,family = "Times New Roman"),
      axis.title = element_text(size = 18,family = "Times New Roman"),
      axis.text = element_text(size = 18,color="black",family = "Times New Roman"),
      plot.title = element_text(size = 14,family = "Times New Roman"),
      plot.subtitle = element_text(size = 14,family = "Times New Roman"),
      plot.tag = element_text(size = 24, family = "Times New Roman"),
      plot.caption = element_text(size = 14,family = "Times New Roman")
      
    )
  #ggsave(plot,filename=paste0("2009年干预-",names(groups)[j],".png"), width = 8, height = 5)http://127.0.0.1:44273/graphics/plot_zoom_png?width=1200&height=900
}

# Geographic region plots
Geo_plot <- (plot1[[2]] | plot2[[2]]) + 
  plot_annotation(
    tag_levels = list(c("A", "B")),
    theme = theme(
      text = element_text(family = "Times New Roman"),
      plot.tag = element_text(size = 24,  family = "Times New Roman"),
      axis.title = element_text(family = "Times New Roman"),
      axis.ticks = element_line(color = "black"), 
      axis.text = element_text(size = 18, family = "Times New Roman", color = "black"),  # 设置坐标轴文本的颜色为黑色
      legend.title = element_text(family = "Times New Roman"),
      legend.text = element_text(family = "Times New Roman"),
      strip.text = element_text(family = "Times New Roman")
    )
  )
#ggsave("Geo_plot.png", Geo_plot, height=8, width=24, limitsize = FALSE, units = "in", dpi = 300)
print(Geo_plot)

# HBsAg plots
HBsAg_plot <- (plot1[[1]] | plot2[[1]]) + 
  plot_annotation(
    tag_levels = list(c("A", "B")),
    theme = theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(size = 16,  family = "Times New Roman"),
      plot.subtitle = element_text(size = 16, family = "Times New Roman"),
      plot.tag = element_text(size = 24,  family = "Times New Roman"),
      axis.title = element_text(size = 16,family = "Times New Roman"),
      axis.text = element_text(size = 18, family = "Times New Roman", color = "black"), 
      legend.title = element_text(size = 16,family = "Times New Roman"),
      legend.text = element_text(size = 16,family = "Times New Roman"),
      strip.text = element_text(size = 16,family = "Times New Roman")
    )
  )
#ggsave("HBsAg_plot.png", HBsAg_plot, height=8, width=24, limitsize = FALSE, units = "in", dpi = 300)
print(HBsAg_plot)

# Other plots
plot11 <- list()
for(i in 3:5) plot11[[i-2]] <- plot1[[i]]

plot21 <- list()
for(i in 3:5) plot21[[i-2]] <- plot2[[i]]

first_column <- wrap_plots(plot11, ncol = 1)
second_column <- wrap_plots(plot21, ncol = 1)

combined_plot <- (first_column | second_column) + 
  plot_annotation(
    tag_levels = list(c("A", "C", "E", "B", "D", "F")),
    theme = theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(size = 16, face = "bold", family = "Times New Roman"),
      plot.subtitle = element_text(size = 14, family = "Times New Roman"),
      plot.tag = element_text(size = 24,  family = "Times New Roman"),
      axis.title = element_text(family = "Times New Roman"),
      axis.text = element_text(size = 18, family = "Times New Roman", color = "black"),  
      legend.title = element_text(family = "Times New Roman"),
      legend.text = element_text(family = "Times New Roman"),
      strip.text = element_text(family = "Times New Roman")
    )
  )
ggsave("combined_plot.png", combined_plot, height=24, width=24, limitsize = FALSE, units = "in", dpi = 300)
print(combined_plot)



#### Use Hotelling’s T² test (based on F-statistic) to assess differences in ER across subgroups####
num1 <- c(1,2,3,4,5.1,6.1,7.1,8.1)
num2 <- c(5.2,6.2,7.2,8.2,9,10)

#Compute the coef first
#2002 NIP
coef11 <-vcov11 <- list()
coef11 <- yori2*population11[1,] 
for (i in 3:8){
  coef11 <- coef11 + eval(parse(text=paste0("yori",num1[i],"*population11[i-1,]")))
}

#2009 catch-up vaccine
coef21 <-vcov21 <- list()
coef21 <- yori6.2*population12[1,] 
for (i in 3:6){
  coef21 <- coef21 + eval(parse(text=paste0("yori",num2[i-1],"*population12[i-1,]")))
}

coef_groups1 <- cbind(coef11, groups[,-1])
coef_groups2 <- cbind(coef21, groups[,-1])
sta_for_Ftest1 <- sta_for_Ftest2 <- list()

sta_for_Ftest1$HBsAgB <- coef_groups1 %>% group_by(HBsAgB) %>% summarise(n(),mean(beta2),mean(beta3),var(beta2),var(beta3),cov(beta2,beta3)) 
sta_for_Ftest1$Geo <- coef_groups1 %>% group_by(Geo) %>% summarise(n(),mean(beta2),mean(beta3),var(beta2),var(beta3),cov(beta2,beta3)) 
sta_for_Ftest1$Ill <- coef_groups1 %>% group_by(Ill) %>% summarise(n(),mean(beta2),mean(beta3),var(beta2),var(beta3),cov(beta2,beta3)) 
sta_for_Ftest1$Bed <- coef_groups1 %>% group_by(Bed) %>% summarise(n(),mean(beta2),mean(beta3),var(beta2),var(beta3),cov(beta2,beta3)) 
sta_for_Ftest1$Ur <- coef_groups1 %>% group_by(Ur) %>% summarise(n(),mean(beta2),mean(beta3),var(beta2),var(beta3),cov(beta2,beta3)) 

sta_for_Ftest2$HBsAgB <- coef_groups2 %>% group_by(HBsAgB) %>% summarise(n(),mean(beta2),mean(beta3),var(beta2),var(beta3),cov(beta2,beta3)) 
sta_for_Ftest2$Geo <- coef_groups2 %>% group_by(Geo) %>% summarise(n(),mean(beta2),mean(beta3),var(beta2),var(beta3),cov(beta2,beta3)) 
sta_for_Ftest2$Ill <- coef_groups2 %>% group_by(Ill) %>% summarise(n(),mean(beta2),mean(beta3),var(beta2),var(beta3),cov(beta2,beta3)) 
sta_for_Ftest2$Bed <- coef_groups2 %>% group_by(Bed) %>% summarise(n(),mean(beta2),mean(beta3),var(beta2),var(beta3),cov(beta2,beta3)) 
sta_for_Ftest2$Ur <- coef_groups2 %>% group_by(Ur) %>% summarise(n(),mean(beta2),mean(beta3),var(beta2),var(beta3),cov(beta2,beta3)) 

m=2
p.value <- list()

# Compute the F-statistic and p-value for 2002 NIP
for(i in 1:3){
  num1 <- compare[i,1]
  num2 <- compare[i,2]
  n1 <- sta_for_Ftest1$HBsAgB$`n()`[num1]
  n2 <- sta_for_Ftest1$HBsAgB$`n()`[num2]
  mean1 <- c(sta_for_Ftest1$HBsAgB$`mean(beta2)`[num1],sta_for_Ftest1$HBsAgB$`mean(beta3)`[num1])
  mean2 <- c(sta_for_Ftest1$HBsAgB$`mean(beta2)`[num2],sta_for_Ftest1$HBsAgB$`mean(beta3)`[num2])
  s1 <- matrix(c(sta_for_Ftest1$HBsAgB$`var(beta2)`[num1],sta_for_Ftest1$HBsAgB$`cov(beta2, beta3)`[num1],sta_for_Ftest1$HBsAgB$`cov(beta2, beta3)`[num1],sta_for_Ftest1$HBsAgB$`var(beta3)`[num1]),2)
  s2 <- matrix(c(sta_for_Ftest1$HBsAgB$`var(beta2)`[num2],sta_for_Ftest1$HBsAgB$`cov(beta2, beta3)`[num2],sta_for_Ftest1$HBsAgB$`cov(beta2, beta3)`[num2],sta_for_Ftest1$HBsAgB$`var(beta3)`[num2]),2)
  
  Sc <- ((n1-1)*s1+(n2-1)*s2)/(n1+n2-2)
  T2 <- n1*n2/(n1+n2) * t(mean1-mean2) %*% solve(Sc) %*% (mean1-mean2)
  F_sta <-(n1+n2-m-1)/(n1+n2-2)/m*T2 
  v1 <- m
  v2 <- n1+n2-m-1
  p.value[i] <- 1-pf(F_sta, v1, v2)
}
for(i in 4:18){
  num1 <- compare[i,1]
  num2 <- compare[i,2]
  n1 <- sta_for_Ftest1$Geo$`n()`[num1]
  n2 <- sta_for_Ftest1$Geo$`n()`[num2]
  mean1 <- c(sta_for_Ftest1$Geo$`mean(beta2)`[num1],sta_for_Ftest1$Geo$`mean(beta3)`[num1])
  mean2 <- c(sta_for_Ftest1$Geo$`mean(beta2)`[num2],sta_for_Ftest1$Geo$`mean(beta3)`[num2])
  s1 <- matrix(c(sta_for_Ftest1$Geo$`var(beta2)`[num1],sta_for_Ftest1$Geo$`cov(beta2, beta3)`[num1],sta_for_Ftest1$Geo$`cov(beta2, beta3)`[num1],sta_for_Ftest1$Geo$`var(beta3)`[num1]),2)
  s2 <- matrix(c(sta_for_Ftest1$Geo$`var(beta2)`[num2],sta_for_Ftest1$Geo$`cov(beta2, beta3)`[num2],sta_for_Ftest1$Geo$`cov(beta2, beta3)`[num2],sta_for_Ftest1$Geo$`var(beta3)`[num2]),2)
  
  Sc <- ((n1-1)*s1+(n2-1)*s2)/(n1+n2-2)
  T2 <- n1*n2/(n1+n2) * t(mean1-mean2) %*% solve(Sc) %*% (mean1-mean2)
  F_sta <-(n1+n2-m-1)/(n1+n2-2)/m*T2 
  v1 <- m
  v2 <- n1+n2-m-1
  p.value[i] <- 1-pf(F_sta, v1, v2)
}
for(i in 19){
  num1 <- compare[i,1]
  num2 <- compare[i,2]
  n1 <- sta_for_Ftest1$Ill$`n()`[num1]
  n2 <- sta_for_Ftest1$Ill$`n()`[num2]
  mean1 <- c(sta_for_Ftest1$Ill$`mean(beta2)`[num1],sta_for_Ftest1$Ill$`mean(beta3)`[num1])
  mean2 <- c(sta_for_Ftest1$Ill$`mean(beta2)`[num2],sta_for_Ftest1$Ill$`mean(beta3)`[num2])
  s1 <- matrix(c(sta_for_Ftest1$Ill$`var(beta2)`[num1],sta_for_Ftest1$Ill$`cov(beta2, beta3)`[num1],sta_for_Ftest1$Ill$`cov(beta2, beta3)`[num1],sta_for_Ftest1$Ill$`var(beta3)`[num1]),2)
  s2 <- matrix(c(sta_for_Ftest1$Ill$`var(beta2)`[num2],sta_for_Ftest1$Ill$`cov(beta2, beta3)`[num2],sta_for_Ftest1$Ill$`cov(beta2, beta3)`[num2],sta_for_Ftest1$Ill$`var(beta3)`[num2]),2)
  
  Sc <- ((n1-1)*s1+(n2-1)*s2)/(n1+n2-2)
  T2 <- n1*n2/(n1+n2) * t(mean1-mean2) %*% solve(Sc) %*% (mean1-mean2)
  F_sta <-(n1+n2-m-1)/(n1+n2-2)/m*T2 
  v1 <- m
  v2 <- n1+n2-m-1
  p.value[i] <- 1-pf(F_sta, v1, v2)
}
for(i in 20){
  num1 <- compare[i,1]
  num2 <- compare[i,2]
  n1 <- sta_for_Ftest1$Bed$`n()`[num1]
  n2 <- sta_for_Ftest1$Bed$`n()`[num2]
  mean1 <- c(sta_for_Ftest1$Bed$`mean(beta2)`[num1],sta_for_Ftest1$Bed$`mean(beta3)`[num1])
  mean2 <- c(sta_for_Ftest1$Bed$`mean(beta2)`[num2],sta_for_Ftest1$Bed$`mean(beta3)`[num2])
  s1 <- matrix(c(sta_for_Ftest1$Bed$`var(beta2)`[num1],sta_for_Ftest1$Bed$`cov(beta2, beta3)`[num1],sta_for_Ftest1$Bed$`cov(beta2, beta3)`[num1],sta_for_Ftest1$Bed$`var(beta3)`[num1]),2)
  s2 <- matrix(c(sta_for_Ftest1$Bed$`var(beta2)`[num2],sta_for_Ftest1$Bed$`cov(beta2, beta3)`[num2],sta_for_Ftest1$Bed$`cov(beta2, beta3)`[num2],sta_for_Ftest1$Bed$`var(beta3)`[num2]),2)
  
  Sc <- ((n1-1)*s1+(n2-1)*s2)/(n1+n2-2)
  T2 <- n1*n2/(n1+n2) * t(mean1-mean2) %*% solve(Sc) %*% (mean1-mean2)
  F_sta <-(n1+n2-m-1)/(n1+n2-2)/m*T2 
  v1 <- m
  v2 <- n1+n2-m-1
  p.value[i] <- 1-pf(F_sta, v1, v2)
}
for(i in 21){
  num1 <- compare[i,1]
  num2 <- compare[i,2]
  n1 <- sta_for_Ftest1$Ur$`n()`[num1]
  n2 <- sta_for_Ftest1$Ur$`n()`[num2]
  mean1 <- c(sta_for_Ftest1$Ur$`mean(beta2)`[num1],sta_for_Ftest1$Ur$`mean(beta3)`[num1])
  mean2 <- c(sta_for_Ftest1$Ur$`mean(beta2)`[num2],sta_for_Ftest1$Ur$`mean(beta3)`[num2])
  s1 <- matrix(c(sta_for_Ftest1$Ur$`var(beta2)`[num1],sta_for_Ftest1$Ur$`cov(beta2, beta3)`[num1],sta_for_Ftest1$Ur$`cov(beta2, beta3)`[num1],sta_for_Ftest1$Ur$`var(beta3)`[num1]),2)
  s2 <- matrix(c(sta_for_Ftest1$Ur$`var(beta2)`[num2],sta_for_Ftest1$Ur$`cov(beta2, beta3)`[num2],sta_for_Ftest1$Ur$`cov(beta2, beta3)`[num2],sta_for_Ftest1$Ur$`var(beta3)`[num2]),2)
  
  Sc <- ((n1-1)*s1+(n2-1)*s2)/(n1+n2-2)
  T2 <- n1*n2/(n1+n2) * t(mean1-mean2) %*% solve(Sc) %*% (mean1-mean2)
  F_sta <-(n1+n2-m-1)/(n1+n2-2)/m*T2 
  v1 <- m
  v2 <- n1+n2-m-1
  p.value[i] <- 1-pf(F_sta, v1, v2)
}
coef_compared1 <- cbind(compare, p.value)

# Compute the F-statistic and p-value for 2009 catch-up vaccine
for(i in 1:3){ 
  num1 <- compare[i,1]
  num2 <- compare[i,2]
  n1 <- sta_for_Ftest2$HBsAgB$`n()`[num1]
  n2 <- sta_for_Ftest2$HBsAgB$`n()`[num2]
  mean1 <- c(sta_for_Ftest2$HBsAgB$`mean(beta2)`[num1],sta_for_Ftest2$HBsAgB$`mean(beta3)`[num1])
  mean2 <- c(sta_for_Ftest2$HBsAgB$`mean(beta2)`[num2],sta_for_Ftest2$HBsAgB$`mean(beta3)`[num2])
  s1 <- matrix(c(sta_for_Ftest2$HBsAgB$`var(beta2)`[num1],sta_for_Ftest2$HBsAgB$`cov(beta2, beta3)`[num1],sta_for_Ftest2$HBsAgB$`cov(beta2, beta3)`[num1],sta_for_Ftest2$HBsAgB$`var(beta3)`[num1]),2)
  s2 <- matrix(c(sta_for_Ftest2$HBsAgB$`var(beta2)`[num2],sta_for_Ftest2$HBsAgB$`cov(beta2, beta3)`[num2],sta_for_Ftest2$HBsAgB$`cov(beta2, beta3)`[num2],sta_for_Ftest2$HBsAgB$`var(beta3)`[num2]),2)
  
  Sc <- ((n1-1)*s1+(n2-1)*s2)/(n1+n2-2)
  T2 <- n1*n2/(n1+n2) * t(mean1-mean2) %*% solve(Sc) %*% (mean1-mean2)
  F_sta <-(n1+n2-m-1)/(n1+n2-2)/m*T2 
  v1 <- m
  v2 <- n1+n2-m-1
  p.value[i] <- 1-pf(F_sta, v1, v2)
}
for(i in 4:18){
  num1 <- compare[i,1]
  num2 <- compare[i,2]
  n1 <- sta_for_Ftest2$Geo$`n()`[num1]
  n2 <- sta_for_Ftest2$Geo$`n()`[num2]
  mean1 <- c(sta_for_Ftest2$Geo$`mean(beta2)`[num1],sta_for_Ftest2$Geo$`mean(beta3)`[num1])
  mean2 <- c(sta_for_Ftest2$Geo$`mean(beta2)`[num2],sta_for_Ftest2$Geo$`mean(beta3)`[num2])
  s1 <- matrix(c(sta_for_Ftest2$Geo$`var(beta2)`[num1],sta_for_Ftest2$Geo$`cov(beta2, beta3)`[num1],sta_for_Ftest2$Geo$`cov(beta2, beta3)`[num1],sta_for_Ftest2$Geo$`var(beta3)`[num1]),2)
  s2 <- matrix(c(sta_for_Ftest2$Geo$`var(beta2)`[num2],sta_for_Ftest2$Geo$`cov(beta2, beta3)`[num2],sta_for_Ftest2$Geo$`cov(beta2, beta3)`[num2],sta_for_Ftest2$Geo$`var(beta3)`[num2]),2)
  
  Sc <- ((n1-1)*s1+(n2-1)*s2)/(n1+n2-2)
  T2 <- n1*n2/(n1+n2) * t(mean1-mean2) %*% solve(Sc) %*% (mean1-mean2)
  F_sta <-(n1+n2-m-1)/(n1+n2-2)/m*T2 
  v1 <- m
  v2 <- n1+n2-m-1
  p.value[i] <- 1-pf(F_sta, v1, v2)
}
for(i in 19){
  num1 <- compare[i,1]
  num2 <- compare[i,2]
  n1 <- sta_for_Ftest2$Ill$`n()`[num1]
  n2 <- sta_for_Ftest2$Ill$`n()`[num2]
  mean1 <- c(sta_for_Ftest2$Ill$`mean(beta2)`[num1],sta_for_Ftest2$Ill$`mean(beta3)`[num1])
  mean2 <- c(sta_for_Ftest2$Ill$`mean(beta2)`[num2],sta_for_Ftest2$Ill$`mean(beta3)`[num2])
  s1 <- matrix(c(sta_for_Ftest2$Ill$`var(beta2)`[num1],sta_for_Ftest2$Ill$`cov(beta2, beta3)`[num1],sta_for_Ftest2$Ill$`cov(beta2, beta3)`[num1],sta_for_Ftest2$Ill$`var(beta3)`[num1]),2)
  s2 <- matrix(c(sta_for_Ftest2$Ill$`var(beta2)`[num2],sta_for_Ftest2$Ill$`cov(beta2, beta3)`[num2],sta_for_Ftest2$Ill$`cov(beta2, beta3)`[num2],sta_for_Ftest2$Ill$`var(beta3)`[num2]),2)
  
  Sc <- ((n1-1)*s1+(n2-1)*s2)/(n1+n2-2)
  T2 <- n1*n2/(n1+n2) * t(mean1-mean2) %*% solve(Sc) %*% (mean1-mean2)
  F_sta <-(n1+n2-m-1)/(n1+n2-2)/m*T2 
  v1 <- m
  v2 <- n1+n2-m-1
  p.value[i] <- 1-pf(F_sta, v1, v2)
}
for(i in 20){
  num1 <- compare[i,1]
  num2 <- compare[i,2]
  n1 <- sta_for_Ftest2$Bed$`n()`[num1]
  n2 <- sta_for_Ftest2$Bed$`n()`[num2]
  mean1 <- c(sta_for_Ftest2$Bed$`mean(beta2)`[num1],sta_for_Ftest2$Bed$`mean(beta3)`[num1])
  mean2 <- c(sta_for_Ftest2$Bed$`mean(beta2)`[num2],sta_for_Ftest2$Bed$`mean(beta3)`[num2])
  s1 <- matrix(c(sta_for_Ftest2$Bed$`var(beta2)`[num1],sta_for_Ftest2$Bed$`cov(beta2, beta3)`[num1],sta_for_Ftest2$Bed$`cov(beta2, beta3)`[num1],sta_for_Ftest2$Bed$`var(beta3)`[num1]),2)
  s2 <- matrix(c(sta_for_Ftest2$Bed$`var(beta2)`[num2],sta_for_Ftest2$Bed$`cov(beta2, beta3)`[num2],sta_for_Ftest2$Bed$`cov(beta2, beta3)`[num2],sta_for_Ftest2$Bed$`var(beta3)`[num2]),2)
  
  Sc <- ((n1-1)*s1+(n2-1)*s2)/(n1+n2-2)
  T2 <- n1*n2/(n1+n2) * t(mean1-mean2) %*% solve(Sc) %*% (mean1-mean2)
  F_sta <-(n1+n2-m-1)/(n1+n2-2)/m*T2 
  v1 <- m
  v2 <- n1+n2-m-1
  p.value[i] <- 1-pf(F_sta, v1, v2)
}
for(i in 21){
  num1 <- compare[i,1]
  num2 <- compare[i,2]
  n1 <- sta_for_Ftest2$Ur$`n()`[num1]
  n2 <- sta_for_Ftest2$Ur$`n()`[num2]
  mean1 <- c(sta_for_Ftest2$Ur$`mean(beta2)`[num1],sta_for_Ftest2$Ur$`mean(beta3)`[num1])
  mean2 <- c(sta_for_Ftest2$Ur$`mean(beta2)`[num2],sta_for_Ftest2$Ur$`mean(beta3)`[num2])
  s1 <- matrix(c(sta_for_Ftest2$Ur$`var(beta2)`[num1],sta_for_Ftest2$Ur$`cov(beta2, beta3)`[num1],sta_for_Ftest2$Ur$`cov(beta2, beta3)`[num1],sta_for_Ftest2$Ur$`var(beta3)`[num1]),2)
  s2 <- matrix(c(sta_for_Ftest2$Ur$`var(beta2)`[num2],sta_for_Ftest2$Ur$`cov(beta2, beta3)`[num2],sta_for_Ftest2$Ur$`cov(beta2, beta3)`[num2],sta_for_Ftest2$Ur$`var(beta3)`[num2]),2)
  
  Sc <- ((n1-1)*s1+(n2-1)*s2)/(n1+n2-2)
  T2 <- n1*n2/(n1+n2) * t(mean1-mean2) %*% solve(Sc) %*% (mean1-mean2)
  F_sta <-(n1+n2-m-1)/(n1+n2-2)/m*T2 
  v1 <- m
  v2 <- n1+n2-m-1
  p.value[i] <- 1-pf(F_sta, v1, v2)
}
coef_compared2 <- cbind(compare, p.value)

# write.csv(coef_compared1,"coef_compared.2002.csv", row.names = F) 
# write.csv(coef_compared2,"coef_compared.2009.csv", row.names = F)

# ----------------------
# Compare whether there is a difference in EIR between subgroups ###
# ----------------------
Delta.method<-function(P1,P2,Pop1,Pop2,n){
  MFR.point<-P1/P2
  MFR.log<-log(MFR.point)
  Se.MFR<-sqrt((1-abs(P1))/abs(P1)/Pop1+(1-abs(P2))/abs(P2)/Pop2)
  z <- case_when(
    n == 2 ~ 1.96,
    n == 3 ~ 2.38,
    n == 6 ~ 2.97,
  )
  MFR.upper<-exp(MFR.log+n*Se.MFR)
  MFR.lower<-exp(MFR.log-n*Se.MFR)
  MFR<-data.frame(Point=MFR.point,Lower=MFR.lower,Upper=MFR.upper)
  return(MFR)
}

EIR1 <- read_excel("EIR1.xlsx") # EIRs for subgroup
EIR2 <- read_excel("EIR2.xlsx")

EIR1 <- EIR1 %>% mutate(POP=population11_groups$sum/180)
EIR2 <- EIR2 %>% mutate(POP=population12_groups$sum/180)

compare_result1 <- rbind(
  Delta.method(EIR1$EIR.point[1]/100000,EIR1$EIR.point[2]/100000,EIR1$POP[1],EIR1$POP[2],3),
  Delta.method(EIR1$EIR.point[1]/100000,EIR1$EIR.point[3]/100000,EIR1$POP[1],EIR1$POP[3],3),
  Delta.method(EIR1$EIR.point[2]/100000,EIR1$EIR.point[3]/100000,EIR1$POP[2],EIR1$POP[3],3),
  
  Delta.method(EIR1$EIR.point[4]/100000,EIR1$EIR.point[5]/100000,EIR1$POP[4],EIR1$POP[5],6),
  Delta.method(EIR1$EIR.point[4]/100000,EIR1$EIR.point[6]/100000,EIR1$POP[4],EIR1$POP[6],6),
  Delta.method(EIR1$EIR.point[4]/100000,EIR1$EIR.point[7]/100000,EIR1$POP[4],EIR1$POP[7],6),
  Delta.method(EIR1$EIR.point[4]/100000,EIR1$EIR.point[8]/100000,EIR1$POP[4],EIR1$POP[8],6),
  Delta.method(EIR1$EIR.point[4]/100000,EIR1$EIR.point[9]/100000,EIR1$POP[4],EIR1$POP[9],6),
  Delta.method(EIR1$EIR.point[5]/100000,EIR1$EIR.point[6]/100000,EIR1$POP[5],EIR1$POP[6],6),
  Delta.method(EIR1$EIR.point[5]/100000,EIR1$EIR.point[7]/100000,EIR1$POP[5],EIR1$POP[7],6),
  Delta.method(EIR1$EIR.point[5]/100000,EIR1$EIR.point[8]/100000,EIR1$POP[5],EIR1$POP[8],6),
  Delta.method(EIR1$EIR.point[5]/100000,EIR1$EIR.point[9]/100000,EIR1$POP[5],EIR1$POP[9],6),
  Delta.method(EIR1$EIR.point[6]/100000,EIR1$EIR.point[7]/100000,EIR1$POP[6],EIR1$POP[7],6),
  Delta.method(EIR1$EIR.point[6]/100000,EIR1$EIR.point[8]/100000,EIR1$POP[6],EIR1$POP[8],6),
  Delta.method(EIR1$EIR.point[6]/100000,EIR1$EIR.point[9]/100000,EIR1$POP[6],EIR1$POP[9],6),
  Delta.method(EIR1$EIR.point[7]/100000,EIR1$EIR.point[8]/100000,EIR1$POP[7],EIR1$POP[8],6),
  Delta.method(EIR1$EIR.point[7]/100000,EIR1$EIR.point[9]/100000,EIR1$POP[7],EIR1$POP[9],6),
  Delta.method(EIR1$EIR.point[8]/100000,EIR1$EIR.point[9]/100000,EIR1$POP[8],EIR1$POP[9],6),
  
  Delta.method(EIR1$EIR.point[10]/100000,EIR1$EIR.point[11]/100000,EIR1$POP[10],EIR1$POP[11],2),
  Delta.method(EIR1$EIR.point[12]/100000,EIR1$EIR.point[13]/100000,EIR1$POP[12],EIR1$POP[13],2),
  Delta.method(EIR1$EIR.point[14]/100000,EIR1$EIR.point[15]/100000,EIR1$POP[14],EIR1$POP[15],2)
) %>% cbind(compare)

compare_result2 <- rbind(
  Delta.method(EIR2$EIR.point[1]/100000,EIR2$EIR.point[2]/100000,EIR2$POP[1],EIR2$POP[2],3),
  Delta.method(EIR2$EIR.point[1]/100000,EIR2$EIR.point[3]/100000,EIR2$POP[1],EIR2$POP[3],3),
  Delta.method(EIR2$EIR.point[2]/100000,EIR2$EIR.point[3]/100000,EIR2$POP[2],EIR2$POP[3],3),
  
  Delta.method(EIR2$EIR.point[4]/100000,EIR2$EIR.point[5]/100000,EIR2$POP[4],EIR2$POP[5],6),
  Delta.method(EIR2$EIR.point[4]/100000,EIR2$EIR.point[6]/100000,EIR2$POP[4],EIR2$POP[6],6),
  Delta.method(EIR2$EIR.point[4]/100000,EIR2$EIR.point[7]/100000,EIR2$POP[4],EIR2$POP[7],6),
  Delta.method(EIR2$EIR.point[4]/100000,EIR2$EIR.point[8]/100000,EIR2$POP[4],EIR2$POP[8],6),
  Delta.method(EIR2$EIR.point[4]/100000,EIR2$EIR.point[9]/100000,EIR2$POP[4],EIR2$POP[9],6),
  Delta.method(EIR2$EIR.point[5]/100000,EIR2$EIR.point[6]/100000,EIR2$POP[5],EIR2$POP[6],6),
  Delta.method(EIR2$EIR.point[5]/100000,EIR2$EIR.point[7]/100000,EIR2$POP[5],EIR2$POP[7],6),
  Delta.method(EIR2$EIR.point[5]/100000,EIR2$EIR.point[8]/100000,EIR2$POP[5],EIR2$POP[8],6),
  Delta.method(EIR2$EIR.point[5]/100000,EIR2$EIR.point[9]/100000,EIR2$POP[5],EIR2$POP[9],6),
  Delta.method(EIR2$EIR.point[6]/100000,EIR2$EIR.point[7]/100000,EIR2$POP[6],EIR2$POP[7],6),
  Delta.method(EIR2$EIR.point[6]/100000,EIR2$EIR.point[8]/100000,EIR2$POP[6],EIR2$POP[8],6),
  Delta.method(EIR2$EIR.point[6]/100000,EIR2$EIR.point[9]/100000,EIR2$POP[6],EIR2$POP[9],6),
  Delta.method(EIR2$EIR.point[7]/100000,EIR2$EIR.point[8]/100000,EIR2$POP[7],EIR2$POP[8],6),
  Delta.method(EIR2$EIR.point[7]/100000,EIR2$EIR.point[9]/100000,EIR2$POP[7],EIR2$POP[9],6),
  Delta.method(EIR2$EIR.point[8]/100000,EIR2$EIR.point[9]/100000,EIR2$POP[8],EIR2$POP[9],6),
  
  Delta.method(EIR2$EIR.point[10]/100000,EIR2$EIR.point[11]/100000,EIR2$POP[10],EIR2$POP[11],2),
  Delta.method(EIR2$EIR.point[12]/100000,EIR2$EIR.point[13]/100000,EIR2$POP[12],EIR2$POP[13],2),
  Delta.method(EIR2$EIR.point[14]/100000,EIR2$EIR.point[15]/100000,EIR2$POP[14],EIR2$POP[15],2)
) %>% cbind(compare)

compare_result1 <- compare_result1 %>% mutate(
  CI = sprintf("%.2f (%.2f, %.2f)", Point, Lower, Upper)
)

compare_result2 <- compare_result2 %>% mutate(
  CI = sprintf("%.2f (%.2f, %.2f)", Point, Lower, Upper)
)

