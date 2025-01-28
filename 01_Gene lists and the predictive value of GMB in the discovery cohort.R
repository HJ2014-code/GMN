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



#######2.1 validation of GMB in the Rizvi cohort(2018)
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/Validation"
# 设置工作目录
setwd(work_dir)

###去掉同义突变
counts <- read.csv("data_mutations_mskcc.csv")#,row.names = 1,check.names=FALSE)
attach(counts)
r<-which(Consequence=="synonymous_variant")
expr<- counts[-r,]
write.csv(expr,file = "01_240_data_mutations.csv")

#####临床数据，TMB合并
clinical <- read.csv("data_clinical_patient.csv",row.names = 1,check.names=FALSE)
clin <-clinical[,c(1:3,6:9)]
TMB <-read.csv("data_clinical_sample.csv",row.names = 1,check.names=FALSE)

k <- rownames(clin)
row = match(k,rownames(TMB))
row = na.omit(row)
TMB_ =TMB[row,]

clin$Tumor_Sample_Barcode <- TMB_$SAMPLE_ID
clin$Cancer_Type <- TMB_$CANCER_TYPE
clin$TMB <- TMB_$TMB_NONSYNONYMOUS
write.csv(clin,file = "01_240_all_clinical_TMB.csv")

##与模拟数据集有无重合
clin <- read.csv("03_NSCLC_clinical_test.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("01_240_all_clinical_TMB.csv",row.names = 1,check.names=FALSE)

k <- rownames(clin)
row = match(k,rownames(clinical))
row = na.omit(row)
clinical_ =clinical[row,]
clinical_r =clinical[-row,]
write.csv(clinical_r,file = "01_35_all_clinical_TMB.csv")

###validation
clinical <- read.csv("01_35_all_clinical_TMB.csv",row.names = 1,check.names=FALSE)
Gene <- read.csv("01_240_data_mutations.csv")
gene <- read.csv("03_gene_select_NSCLC.csv")
attach(Gene)

###挑选特定基因所有突变样本
D <- data.frame()
for (i in (1:41)) {
  g <-gene[i,2]
  r<-which(Hugo_Symbol==g)
  Gene_r<-Gene[r,]
  D <- rbind(D,Gene_r)
}

###每个样本含有特定基因突变数目
Da <- data.frame()
for (i in (1:35)) {
  samp <-clinical[i,8]
  r<-which(D$Tumor_Sample_Barcode==samp)
  Gene_r<-D[r,]
  n=length(r)
  
  Result <-data.frame(n)
  Da <- rbind(Da, Result)
}
write.csv(Da,file = "03_NSCLC_clinical_number.csv")

#####GMB OS
library(survival)
#Kaplan-Meier 分析，详情 ?survfit
clinical <- read.csv("01_35_all_clinical_TMB.csv",row.names = 1,check.names=FALSE)
attach(clinical)

KM <- survfit(Surv(PFS_MONTHS,Status) ~ Cut, data =clinical, type = 'kaplan-meier', conf.type = 'log')
KM
#####
#绘制生存曲线，反映了尚在世的患者数量比例和时间的关系
library(survminer)
ggsurvplot(KM, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("Gene_mutation_burden ≥ 3", "Gene_mutation_burden < 3"),
           palette="jto", ylab = "Progress free survival" )
###添加HR ,CI ,P
p3 <- ggsurvplot(KM, data = clinical,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 legend.labs = c("Gene_mutation_burden ≥ 3", "Gene_mutation_burden < 3"),
                 palette="jto", ylab = "Progress free survival" )
res_cox<-coxph(Surv(PFS_MONTHS,Status) ~ Cut_, data = clinical)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3

###GMB OS 240
clinical <- read.csv("01_240_all_clinical_TMB.csv",row.names = 1,check.names=FALSE)
attach(clinical)
KM <- survfit(Surv(PFS_MONTHS,Status) ~ Cut, data =clinical, type = 'kaplan-meier', conf.type = 'log')
KM
library(survminer)
ggsurvplot(KM, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("Gene_mutation_burden ≥ 3", "Gene_mutation_burden < 3"),
           palette="jto", ylab = "Progress free survival" )
###添加HR ,CI ,P
p3 <- ggsurvplot(KM, data = clinical,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 legend.labs = c("Gene_mutation_burden ≥ 3", "Gene_mutation_burden < 3"),
                 palette="jto", ylab = "Progress free survival" )
res_cox<-coxph(Surv(PFS_MONTHS,Status) ~ Cut_, data = clinical)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3



#######2.2 validation of GMB in the Hellmann cohort
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/Validation2"
# 设置工作目录
setwd(work_dir)

###去掉同义突变
counts <- read.csv("data_mutations_mskcc.csv")#,row.names = 1,check.names=FALSE)
attach(counts)
r<-which(Consequence=="synonymous_variant")
expr<- counts[-r,]
write.csv(expr,file = "01_75_data_mutations.csv")

#####临床数据，TMB合并
clinical <- read.csv("data_clinical_patient.csv",row.names = 1,check.names=FALSE)
clin <-clinical[,c(1:9,12)]
TMB <-read.csv("data_clinical_sample.csv",row.names = 1,check.names=FALSE)

k <- rownames(clin)
row = match(k,TMB$PATIENT_ID)
row = na.omit(row)
TMB_ =TMB[row,]

clin$Tumor_Sample_Barcode <- rownames(TMB_)
clin$Cancer_Type <- TMB_$CANCER_TYPE
clin$TMB <- TMB_$TMB_NONSYNONYMOUS
write.csv(clin,file = "01_75_all_clinical_TMB.csv")
###validation
clinical <- read.csv("01_75_all_clinical_TMB.csv",row.names = 1,check.names=FALSE)
Gene <- read.csv("01_75_data_mutations.csv")
gene <- read.csv("03_gene_select_NSCLC.csv")
attach(Gene)

###挑选特定基因所有突变样本
D <- data.frame()
for (i in (1:41)) {
  g <-gene[i,2]
  r<-which(Hugo_Symbol==g)
  Gene_r<-Gene[r,]
  D <- rbind(D,Gene_r)
}

###每个样本含有特定基因突变数目
Da <- data.frame()
for (i in (1:75)) {
  samp <-clinical[i,11]
  r<-which(D$Tumor_Sample_Barcode==samp)
  Gene_r<-D[r,]
  n=length(r)
  
  Result <-data.frame(n)
  Da <- rbind(Da, Result)
}
write.csv(Da,file = "03_NSCLC_clinical_number.csv")

#####mutation OS
library(survival)
#Kaplan-Meier 分析，详情 ?survfit
clinical <- read.csv("01_75_all_clinical_TMB.csv",row.names = 1,check.names=FALSE)
attach(clinical)

KM <- survfit(Surv(PFS_MONTHS,Status) ~ Cut, data =clinical, type = 'kaplan-meier', conf.type = 'log')
KM
#####
#绘制生存曲线，反映了尚在世的患者数量比例和时间的关系
library(survminer)
#生存曲线，详情 ?ggsurvplot
ggsurvplot(KM, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("Gene_mutation_burden ≥ 3", "Gene_mutation_burden < 3"),
           palette="jto", ylab = "Progress free survival" )

###添加HR ,CI ,P
p3 <- ggsurvplot(KM, data = clinical,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 legend.labs = c("Gene_mutation_burden ≥ 3", "Gene_mutation_burden < 3"),
                 palette="jto", ylab = "Progress free survival" )
res_cox<-coxph(Surv(PFS_MONTHS,Status) ~ Cut_, data = clinical)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3



#######2.3 validation of GMB in the Rizvi cohort (2015)
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/Validation3"
# 设置工作目录
setwd(work_dir)

###去掉同义突变
counts <- read.csv("34data_mutations.csv")#,row.names = 1,check.names=FALSE)
attach(counts)
r<-which(Consequence=="synonymous_variant")
expr<- counts[-r,]
write.csv(expr,file = "01_34_data_mutations.csv")

#####临床数据，TMB合并
clinical <- read.csv("34data_clinical_patient.csv",row.names = 1,check.names=FALSE)
clin <-clinical[,c(1:6,8,10:13,16)]
TMB <-read.csv("34data_clinical_sample.csv",row.names = 1,check.names=FALSE)

k <- rownames(clin)
row = match(k,rownames(TMB))
row = na.omit(row)
TMB_ =TMB[row,]

clin$Tumor_Sample_Barcode <- TMB_$SAMPLE_ID
clin$Cancer_Type <- TMB_$CANCER_TYPE
clin$TMB <- TMB_$TMB_NONSYNONYMOUS
write.csv(clin,file = "01_34_all_clinical_TMB.csv")
###validation
clinical <- read.csv("01_34_all_clinical_TMB.csv",row.names = 1,check.names=FALSE)
Gene <- read.csv("01_34_data_mutations.csv")
gene <- read.csv("03_gene_select_NSCLC.csv")
attach(Gene)

###挑选特定基因所有突变样本
D <- data.frame()
for (i in (1:41)) {
  g <-gene[i,2]
  r<-which(Hugo_Symbol==g)
  Gene_r<-Gene[r,]
  D <- rbind(D,Gene_r)
}

###每个样本含有特定基因突变数目
Da <- data.frame()
for (i in (1:35)) {
  samp <-clinical[i,13]
  r<-which(D$Tumor_Sample_Barcode==samp)
  Gene_r<-D[r,]
  n=length(r)
  
  Result <-data.frame(n)
  Da <- rbind(Da, Result)
}
write.csv(Da,file = "03_NSCLC_clinical_number.csv")

#####mutation OS
library(survival)
#Kaplan-Meier 分析，详情 ?survfit
clinical <- read.csv("01_34_all_clinical_TMB.csv",row.names = 1,check.names=FALSE)
attach(clinical)

KM <- survfit(Surv(PFS_MONTHS,Status) ~ Cut, data =clinical, type = 'kaplan-meier', conf.type = 'log')
KM
#####
#绘制生存曲线，反映了尚在世的患者数量比例和时间的关系
library(survminer)
#生存曲线，详情 ?ggsurvplot
ggsurvplot(KM, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("Gene_mutation_burden ≥ 3", "Gene_mutation_burden < 3"),
           palette="jto", ylab = "Progress free survival" )

###添加HR ,CI ,P
p3 <- ggsurvplot(KM, data = clinical,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 legend.labs = c("Gene_mutation_burden ≥ 3", "Gene_mutation_burden < 3"),
                 palette="jto", ylab = "Progress free survival" )
res_cox<-coxph(Surv(PFS_MONTHS,Status) ~ Cut_, data = clinical)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3



#######2.4 vallidation forest
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/11_validation_forest"
# 设置工作目录
setwd(work_dir)

result<-read.csv("Validation_forest.csv",row.names = 1,check.names=FALSE)
result$name<-rownames(result)
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



#######2.5 the predictive value of GMB in the validation cohort
#####mutation OS
library(survival)
library(survminer)
#Kaplan-Meier 分析，详情 ?survfit
clinical <- read.csv("01_75_34_35all_clinical_TMB_.csv",row.names = 1,check.names=FALSE)
attach(clinical)

KM <- survfit(Surv(PFS_MONTHS,Status) ~ Cut, data =clinical, type = 'kaplan-meier', conf.type = 'log')
KM
ggsurvplot(KM, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("High_mutation_number", "Low_mutation_number"),
           palette="lancet", ylab = "Progress free survival" )

###添加HR ,CI ,P
p3 <- ggsurvplot(KM, data = clinical,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 legend.labs = c("Gene_mutation_burden ≥ 3", "Gene_mutation_burden < 3"),
                 palette="jto", ylab = "Progress free survival" )
res_cox<-coxph(Surv(PFS_MONTHS,Status) ~ Cut_, data = clinical)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3

#####TMB  生存分析
library(survival)
#Kaplan-Meier 分析，详情 ?survfit
clinical$TMB_ <- ''
clinical[clinical$TMB >= 10,]$TMB_ <- "TMB>10"
clinical[clinical$TMB < 10,]$TMB_ <- "TMB<10"
attach(clinical)
write.csv(clinical,file = "03_75_34_35all_clinical_TMB.csv")
KM <- survfit(Surv(PFS_MONTHS,Status) ~ TMB_, data =clinical, type = 'kaplan-meier', conf.type = 'log')
KM
#绘制生存曲线，反映了尚在世的患者数量比例和时间的关系
library(survminer)
#生存曲线，详情 ?ggsurvplot
ggsurvplot(KM, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("TMB < 10", "TMB > 10"),
           palette="lancet", ylab = "Progress free survival" )

###添加HR ,CI ,P
p3 <- ggsurvplot(KM, data = clinical,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 legend.labs = c("TMB < 10", "TMB > 10"),
                 palette="lancet", ylab = "Progress free survival" )
res_cox<-coxph(Surv(PFS_MONTHS,Status) ~ TMB_, data = clinical)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3

##TMB及GMB
###
library(survival)
#Kaplan-Meier 分析，详情 ?survfit
clinical <- read.csv("03_75_34_35all_clinical_TMB_.csv",row.names = 1,check.names=FALSE)
attach(clinical)

KM <- survfit(Surv(PFS_MONTHS,Status) ~ TMB_Cut, data =clinical, type = 'kaplan-meier', conf.type = 'log')
KM

restest <- pairwise_survdiff(Surv(PFS_MONTHS,Status) ~ TMB_Cut,
                             data = clinical)
restest
#绘制生存曲线，反映了尚在世的患者数量比例和时间的关系
library(survminer)
#生存曲线，详情 ?ggsurvplot
ggsurvplot(KM, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           palette="lancet", ylab = "Overall survival" )



#######2.6 the C index of GMB,TMB and other single mutatin gene
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/10_Contrast"
# 设置工作目录
setwd(work_dir)

###点突变基因
gene <-read.csv("02_Non-Small_Cell_Lung_Cancer_mutation.csv")
clinical<-read.csv("03_NSCLC_clinical.csv",row.names = 1,check.names=FALSE)

###350个样本中的突变基因
D <- data.frame()
for (i in (1:350)) {
  k <- clinical[i,6]
  r<-which(gene$Tumor_Sample_Barcode==k)
  Data_r<-gene[r,]
  D <- rbind(D, Data_r)
}

genes <- c(
  "B2M",
  "ZFHX3",
  "KRAS", 
  "TP53",
  "EPHA5",
  "STK11"
)

response<-clinical
Da <- data.frame()
library(survminer)
library(rms)

for (i in 1:length(genes)) {
  g <-genes[i]
  r<-which(D$Hugo_Symbol==g)
  Data_r<-D[r,]
  Data_r<-Data_r[!duplicated(Data_r$Tumor_Sample_Barcode),]
  
  k <- Data_r$Tumor_Sample_Barcode
  row = match(k,response$Tumor_Sample_Barcode)
  row = na.omit(row)
  n=length(row)
  
  response_ =response[row,]
  response_w =response[-row,]
  Resp <- rbind(response_,response_w)
  Resp$Mut <- c(rep("g-Mut",n),rep("g-Wild",350-n))
  
  attach(Resp)
  cox <- coxph(Surv(OS_MONTHS,STATUS)~Mut,data = Resp)
  #BiocManager::install("survcomp")
  library(survcomp)
  ci4 <- concordance.index(
    predict(cox), surv.time=OS_MONTHS, surv.event=STATUS,
    method="noether"
  )
  #ci4$c.index
  #ci4$lower
  #ci4$upper
  Result <-data.frame(ci4$c.index,ci4$lower,ci4$upper)
  Da <- rbind(Da, Result)
}
write.csv(Da,file = "01_single_gene_index.csv")

###GMB Cindex
cox <- coxph(Surv(OS_MONTHS,STATUS)~cut,data = clinical)
attach(clinical)
ci4 <- concordance.index(
  predict(cox), surv.time=OS_MONTHS, surv.event=STATUS,
  method="noether"
)
ci4$c.index
ci4$lower
ci4$upper
###TMB Cindex
cox <- coxph(Surv(OS_MONTHS,STATUS)~TMB_,data = clinical)
ci4 <- concordance.index(
  predict(cox), surv.time=OS_MONTHS, surv.event=STATUS,
  method="noether"
)
ci4$c.index
ci4$lower
ci4$upper


###validation cohort
clinical<-read.csv("03_75_34_35all_clinical_TMB.csv",row.names = 1,check.names=FALSE)
###Cut Cindex
library(survcomp)
attach(clinical)
cox <- coxph(Surv(PFS_MONTHS,Status)~Cut,data = clinical)
ci4 <- concordance.index(
  predict(cox), surv.time=PFS_MONTHS, surv.event=Status,
  method="noether"
)
ci4$c.index
ci4$lower
ci4$upper

###TMB
cox <- coxph(Surv(PFS_MONTHS,Status)~TMB_,data = clinical)
ci4 <- concordance.index(
  predict(cox), surv.time=PFS_MONTHS, surv.event=Status,
  method="noether"
)
ci4$c.index
ci4$lower
ci4$upper

###条形图可视化Cindex
Cindex <-read.csv("01_single_gene_index_.csv")
library(ggplot2)
Cindex <-Cindex[c(1:3,5:8),]

###添加条形图的文本
p1 <- ggplot(data=Cindex, aes(x=Gene_name,y=ci4.c.index,fill= Gene_name ))+
  geom_bar(stat="identity",width=0.5)+
  theme_minimal()

library(dplyr)
Cindex %>%
  mutate(Gene_name = factor(Gene_name, levels = c("KRAS", "TP53", "STK11","EPHA5","B2M","ZFHX3","GMB"))) %>%
  ggplot(aes(fill = Gene_name, y = ci4.c.index, x = Gene_name)) +
  geom_bar(position = "dodge", stat = "identity",width=0.3)+
  theme_minimal()

###
####TMB与GMB
dataxx <-read.csv("01_GMB_TMB.csv",row.names = 1,check.names=FALSE)
datax <- as.matrix(dataxx)
barplot(datax,beside = TRUE,horiz =T)
barplot(datax,beside = TRUE,angle = 15+10*1:5, density = 10,col = rainbow(2),horiz =T)
barplot(datax,beside = TRUE,angle = 15+10*1:5, density = 10)


attach(dataxx)
x <- barplot(dataxx, xlim = Cohort, offset = 0, axis.lty = 1, names.arg = Cohort,
             col = Marker, beside = TRUE, horiz = T)
box()
legend("bottomright", legend = legs, fill = cols, box.col = "transparent")
title(xlab = "Value", ylab = "Sample")
sd <- dataxx * 0.1
for (i in 1:3) plot.error(x[i, ], dataxx[i, ], sd = sd[i, ], horiz = TRUE)
data <- cbind(a = 1:4, b = 1:4)






















