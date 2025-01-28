#######3.1 GMB was not a prognostic factor in the Zehir Cohort
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/prognosis"
# 设置工作目录
setwd(work_dir)
###去掉同义突变
counts <- read.csv("data_mutations.csv")#,row.names = 1,check.names=FALSE)
attach(counts)
r<-which(Consequence=="synonymous_variant")
expr<- counts[-r,]
#write.csv(expr,file = "01_240_data_mutations.csv")

#####临床数据，TMB合并
clinical <- read.csv("data_clinical_patient.csv",row.names = 1,check.names=FALSE)

TMB <-read.csv("data_clinical_sample.csv")
r<-which(TMB$CANCER_TYPE=="Non-Small Cell Lung Cancer")
TMB_NSCLC <- TMB[r,]
write.csv(TMB_NSCLC,file = "01_1668_TMB.csv")
TMB_NSCLC_ <- TMB_NSCLC[!duplicated(TMB_NSCLC$PATIENT_ID),]

k <- TMB_NSCLC_$PATIENT_ID
row = match(k,rownames(clinical))
row = na.omit(row)
clin =clinical[row,]

clin$Tumor_Sample_Barcode <- TMB_NSCLC_$SAMPLE_ID
clin$Cancer_Type <- TMB_NSCLC_$CANCER_TYPE
clin$TMB <- TMB_NSCLC_$TMB_NONSYNONYMOUS
write.csv(clin,file = "01_1567_all_clinical_TMB.csv")

###预后判断
clinical <- read.csv("01_1443_all_clinical_TMB.csv",row.names = 1,check.names=FALSE)
Gene <- read.csv("data_mutations.csv")
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
for (i in (1:1443)) {
  samp <-clinical[i,6]
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
clinical <- read.csv("01_1443_all_clinical_TMB.csv",row.names = 1,check.names=FALSE)
attach(clinical)

KM <- survfit(Surv(OS_MONTHS,Status) ~ Cut, data =clinical, type = 'kaplan-meier', conf.type = 'log')
KM
library(survminer)
ggsurvplot(KM, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("High_mutation_number", "Low_mutation_number"),
           palette="jto", ylab = "overall survival" )

#计算log rank p
sdiff <- survdiff(Surv(OS_MONTHS,Status)~Cut, data=clinical)
p.val = 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
#计算HR
library(survival)
data.survdiff <- survdiff(Surv(OS_MONTHS,Status) ~ Cut)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 

#######3.2 GMB was not a prognostic factor in the TCGA cohort
##########################LUSC#################################
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/TCGA"
# 设置工作目录
setwd(work_dir)

###去掉同义突变
counts <- read.csv("data_mutations.csv")#,row.names = 1,check.names=FALSE)
attach(counts)
r<-which(Consequence=="synonymous_variant")
expr<- counts[-r,]
write.csv(expr,file = "01_LUSC_data_mutations.csv")

#####临床数据，TMB合并
clinical <- read.csv("data_clinical_patient.csv",row.names = 1,check.names=FALSE)
TMB <-read.csv("data_clinical_sample.csv")

k <- TMB$PATIENT_ID
row = match(k,rownames(clinical))
row = na.omit(row)
clin =clinical[row,]

clin$Tumor_Sample_Barcode <- TMB$SAMPLE_ID
clin$Cancer_Type <- TMB$CANCER_TYPE
clin$TMB <- TMB$TMB_NONSYNONYMOUS
write.csv(clin,file = "01_487_all_clinical_TMB.csv")

###预后判断
clinical <- read.csv("01_487_all_clinical_TMB.csv",row.names = 1,check.names=FALSE)
Gene <- read.csv("01_LUSC_data_mutations.csv")
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
for (i in (1:487)) {
  samp <-clinical[i,17]
  r<-which(D$Tumor_Sample_Barcode==samp)
  Gene_r<-D[r,]
  n=length(r)
  
  Result <-data.frame(n)
  Da <- rbind(Da, Result)
}
write.csv(Da,file = "03_NSCLC_clinical_number.csv")


##########################LUAD#################################
###去掉同义突变
counts <- read.csv("data_mutations.csv")#,row.names = 1,check.names=FALSE)
attach(counts)
r<-which(Consequence=="synonymous_variant")
expr<- counts[-r,]
write.csv(expr,file = "01_LUAD_data_mutations.csv")

#####临床数据，TMB合并
clinical <- read.csv("data_clinical_patient.csv",row.names = 1,check.names=FALSE)
clin <- clinical[,c(-5,-6,-7)]

TMB <-read.csv("data_clinical_sample.csv")

k <- TMB$PATIENT_ID
row = match(k,rownames(clinical))
row = na.omit(row)
clin =clinical[row,]

clin$Tumor_Sample_Barcode <- TMB$SAMPLE_ID
clin$Cancer_Type <- TMB$CANCER_TYPE
clin$TMB <- TMB$TMB_NONSYNONYMOUS
write.csv(clin,file = "01_566_all_clinical_TMB.csv")

###预后判断
clinical <- read.csv("01_566_all_clinical_TMB.csv",row.names = 1,check.names=FALSE)
Gene <- read.csv("01_LUAD_data_mutations.csv")
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
for (i in (1:566)) {
  samp <-clinical[i,20]
  r<-which(D$Tumor_Sample_Barcode==samp)
  Gene_r<-D[r,]
  n=length(r)
  
  Result <-data.frame(n)
  Da <- rbind(Da, Result)
}
write.csv(Da,file = "03_NSCLC_clinical_number.csv")


##########################NSCLC#################################
library(survival)
#Kaplan-Meier 分析，详情 ?survfit
clinical <- read.csv("01_566_487_all_clinical_TMB.csv",row.names = 1,check.names=FALSE)
attach(clinical)

KM <- survfit(Surv(OS_MONTHS,Status) ~ Cut, data =clinical, type = 'kaplan-meier', conf.type = 'log')
KM

#绘制生存曲线，反映了尚在世的患者数量比例和时间的关系
library(survminer)
ggsurvplot(KM, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("High_mutation_number", "Low_mutation_number"),
           palette="jto", ylab = "Progress free survival" )
#计算log rank p
sdiff <- survdiff(Surv(OS_MONTHS,Status)~Cut, data=clinical)
p.val = 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
#计算HR
library(survival)
data.survdiff <- survdiff(Surv(OS_MONTHS,Status) ~ Cut)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 



