#######1.1 select candidate genes
rm(list=ls())
# 设置环境参数
#work_dir <- "C:/Users/特定突变基因数目VSTMB"
# 设置工作目录
#setwd(work_dir)
clinical <- read.csv("01_Non-Small_Cell_Lung_Cancer_clinical.csv",row.names = 1,check.names=FALSE)
Gene <-read.csv("02_Non-Small_Cell_Lung_Cancer_mutation.csv")

response<-clinical
Da <- data.frame()
for (i in (1:3537)) {
  g <-Gene[i,1]
  r<-which(Gene$Hugo_Symbol==g)
  Data_r<-Gene[r,]
  Data_r<-Data_r[!duplicated(Data_r$Tumor_Sample_Barcode),]
  
  k <- Data_r$Tumor_Sample_Barcode
  row = match(k,response$Tumor_Sample_Barcode)
  row = na.omit(row)
  n=length(row)
  
  response_ =response[row,]
  response_w =response[-row,]
  Resp <- rbind(response_,response_w)
  Resp$Mut <- c(rep("g-Mut",n),rep("g-Wild",350-n))#
  #单个卡方检验
  table(Resp$Mut,Resp$OS_STATUS)
  f=fisher.test(Resp$Mut,Resp$OS_STATUS)
  Result <-data.frame(n,f$p.value,f$estimate)
  Da <- rbind(Da, Result)
}
write.csv(Da,file = "02_NSCLC_mutation_number.csv")

gene<- read.csv("02_NSCLC_mutation_select.csv")
gene_sel<-gene[!duplicated(gene$Hugo_Symbol),]
write.csv(gene_sel,file = "03_gene_select_NSCLC.csv")



#######1.2 Forest
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/14_GMB_H_L"
# 设置工作目录
setwd(work_dir)

library(survival)
library(survminer)

clinical <- read.csv("03_NSCLC_clinical_11.csv",row.names = 1,check.names=FALSE)

####GMB HR
KM <- survfit(Surv(OS_MONTHS, STATUS) ~ Cut_, data =clinical, type = 'kaplan-meier', conf.type = 'log')
KM

###添加HR ,CI ,P
p3 <- ggsurvplot(KM, data = clinical,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 legend.labs = c("Gene_mutation_burden < 11", "Gene_mutation_burden ≥ 11"),
                 palette="jto", ylab = "Overall survival" )
res_cox<-coxph(Surv(OS_MONTHS, STATUS) ~ Cut_, data = clinical)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3

########
###########forest
rm(list=ls())
result<-read.csv("14_forest_UTF.csv",row.names = 1,check.names=FALSE)

result$name<-rownames(result)
#result <-result[order(result[,3]),]

res<-result
#基本图形
#1-2.画森林图的包
#install.packages("forestplot")
library(forestplot)
library(stringr)
fig1<-forestplot(res[,c(8,7,3)],  #1,6,5列显示为变量 HR(CI) p数值形式
                 mean=res[,4],   #表格第2列为HR，要变成森林图的小方块
                 lower=res[,5],  #表格第3列为5%CI，
                 upper=res[,6],  #表格第4列为95%CI，它俩要化作线段，穿过方块
                 zero=1,            #零线或参考线为HR=1即x轴的垂直线
                 boxsize=0.4,       #设置小黑块的大小
                 graph.pos="right") #森林图放在最右侧
fig1

#优化
fig1<-forestplot(res[,c(8,7,3)],  #1,6,5列显示为变量 HR(CI) p数值形式
                 mean=res[,4],   #表格第2列为HR，要变成森林图的小方块
                 lower=res[,5],  #表格第3列为5%CI，
                 upper=res[,6],  #表格第4列为95%CI，它俩要化作线段，穿过方块
                 zero=1,            #零线或参考线为HR=1即x轴的垂直线
                 boxsize=0.1,       #设置小黑块的大小
                 graph.pos="right", #森林图放在最右侧
                 ##fpColors函数设置颜色
                 col=fpColors(box="#1c61b6", lines="#1c61b6", zero = "gray50"),
                 #箱线图中基准线的位置
                 cex=0.9, lineheight = "auto",
                 colgap=unit(8,"mm"),
                 #箱子大小，线的宽度
                 lwd.ci=1,
                 #箱线图两端添加小竖线，高度
                 ci.vertices=TRUE, ci.vertices.height = 0.08,
                 ##定义x轴
                 xlab=" HR (95% CI)",
                 #fpTxtGp函数中的cex参数设置各个组件的大小
                 txt_gp=fpTxtGp(label=gpar(cex=1.25),
                                ticks=gpar(cex=1.1),
                                xlab=gpar(cex = 1.2),
                                title=gpar(cex = 1.2)),) 
fig1



#######1.3 the predictive value of GMB in the discovery cohort
rm(list=ls())
clinical <- read.csv("01_Non-Small_Cell_Lung_Cancer_clinical.csv",row.names = 1,check.names=FALSE)
Gene <- read.csv("02_Non-Small_Cell_Lung_Cancer_mutation.csv")
gene <- read.csv("03_gene_select_NSCLC.csv")
attach(Gene)

D <- data.frame()
for (i in (1:41)) {
  g <-gene[i,2]
  r<-which(Hugo_Symbol==g)
  Gene_r<-Gene[r,]
  D <- rbind(D,Gene_r)
}

Da <- data.frame()
for (i in (1:350)) {
  samp <-clinical[i,6]
  r<-which(D$Tumor_Sample_Barcode==samp)
  Gene_r<-D[r,]
  n=length(r)
  
  Result <-data.frame(n)
  Da <- rbind(Da, Result)
}
write.csv(Da,file = "03_NSCLC_clinical_number.csv")

library(survminer)
res.cut <- surv_cutpoint(clinical,time="OS_MONTHS",event = "STATUS",
                         variables = "n")
summary(res.cut)   
plot(res.cut,"n",palette="npg")
res.cat <- surv_categorize(res.cut)

library("survival")
fit<- survfit(Surv(OS_MONTHS, STATUS) ~ n, data = res.cat)
surv_pvalue(fit)$pval.txt
clinical$cut<-  res.cat$n
write.csv(clinical,file = "03_NSCLC_clinical.csv")

ggsurvplot(fit, data = res.cat,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("High_mutation_number", "Low_mutation_number"),
           palette="lancet", ylab = "Overall survival" )

p3 <- ggsurvplot(fit, data = clinical,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 legend.labs = c("Gene_mutation_burden ≥ 3", "Gene_mutation_burden < 3"),
                 palette="jto", ylab = "Overall survival" )
res_cox<-coxph(Surv(OS_MONTHS, STATUS) ~ n, data = res.cat)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3


#####TMB suvival
library(survival)
attach(clinical)
KM <- survfit(Surv(OS_MONTHS,STATUS) ~ TMB_, data =clinical, type = 'kaplan-meier', conf.type = 'log')
KM
library(survminer)
ggsurvplot(KM, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("TMB < 10", "TMB > 10"),
           palette="lancet", ylab = "Overall survival" )

##TMB VS. GMB
library(survival)
clinical <- read.csv("03_NSCLC_clinical.csv",row.names = 1,check.names=FALSE)
attach(clinical)
KM <- survfit(Surv(OS_MONTHS,STATUS) ~ TMB_Cut, data =clinical, type = 'kaplan-meier', conf.type = 'log')
KM

survdiff(Surv(OS_MONTHS,STATUS) ~ TMB_Mut, data = clinical)
restest <- pairwise_survdiff(Surv(OS_MONTHS,STATUS) ~ TMB_Cut,
                             data = clinical)
restest
ggsurvplot(KM, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           #legend.labs = c("TMB < 10", "TMB > 10"),
           palette="jto", ylab = "Overall survival" )
KM <- survfit(Surv(OS_MONTHS, STATUS) ~ Cut_, data =clinical, type = 'kaplan-meier', conf.type = 'log')
KM
###HR ,CI ,P
p3 <- ggsurvplot(KM, data = clinical,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 legend.labs = c("Gene_mutation_burden ≥ 3", "Gene_mutation_burden < 3"),
                 palette="jto", ylab = "Overall survival" )
res_cox<-coxph(Surv(OS_MONTHS, STATUS) ~ Cut_, data = clinical)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3












