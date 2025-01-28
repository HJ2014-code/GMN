rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/37_in_house_Validation"
# 设置工作目录
setwd(work_dir)


#######EGFR Mut VS Wild OS  #########
Clin <-read.csv("01_Lei_Clinical_Seq.csv",row.names = 1,check.names=FALSE)
Clin <-read.csv("02_Lei_Clinical_Seq_GMN.csv",row.names = 1,check.names=FALSE)
Mut <- read.csv("01_Lei_Gene_Mut.csv")

attach(Mut)
r<-which(Gene_Symbol =="EGFR")
Mut_EGFR<- Mut[r,]

Mut_EGFR_Case <- Mut_EGFR[!duplicated(Mut_EGFR$Case_ID),]
k <- Mut_EGFR_Case$Case_ID
row = match(k,rownames(Clin))
row = na.omit(row)
n=length(row)

Clin_M =Clin[row,]
Clin_w =Clin[-row,]
Resp <- rbind(Clin_M,Clin_w)
Resp$EGFR <- c(rep("XMut",n),rep("Wild",174-n))
write.csv(Resp,file = "03_NSCLC_Seq_GMN_EGFR.csv")

library("survival")
library(survminer)
fit<- survfit(Surv(OS_months,OS_status) ~ EGFR, data = Resp)
surv_pvalue(fit)$pval.txt
ggsurvplot(fit, data = Resp,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("Wild","Mut"),
           palette="lancet", ylab = "Overall survival" )



#######GMN#########
##########################
#######################
rm(list=ls())
Clin <-read.csv("01_Lei_Clinical_Seq.csv",row.names = 1,check.names=FALSE)
Gene <- read.csv("01_Lei_Gene_Mut.csv")
gene <- read.csv("07_gene_select_NSCLC.csv")
attach(Gene)

###挑选特定基因所有突变样本
D <- data.frame()
for (i in (1:41)) {
  g <-gene[i,1]
  r<-which(Gene_Symbol==g)
  Gene_r<-Gene[r,]
  D <- rbind(D,Gene_r)
}

###每个样本含有特定基因突变数目
Da <- data.frame()
for (i in (1:174)) {
  samp <- rownames(Clin)[i]
  r<-which(D$Case_ID ==samp)
  Gene_r<-D[r,]
  n=length(r)
  
  Result <-data.frame(n)
  Da <- rbind(Da, Result)
}
write.csv(Da,file = "02_NSCLC_GMN.csv")


#######普通 GMN  #########
##########################
#######################
rm(list=ls())
Clin <-read.csv("02_Lei_Clinical_Seq_GMN.csv",row.names = 1,check.names=FALSE)

library("survival")
library(survminer)
fit<- survfit(Surv(OS_months,OS_status) ~ GMN, data = Clin)
surv_pvalue(fit)$pval.txt
ggsurvplot(fit, data = Clin,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("GMB ≥ 3","GMB < 3"),
           palette="lancet", ylab = "Overall survival" )




#######EGFR突变内  #########
##########################
#######################
rm(list=ls())
Clin <-read.csv("03_NSCLC_Seq_GMN_EGFR_Mut.csv",row.names = 1,check.names=FALSE)
Clin <-read.csv("03_NSCLC_Seq_GMN_EGFR_Mut_drug1.csv",row.names = 1,check.names=FALSE)

library("survival")
library(survminer)
#Resp <- Clin[c(1:56),]
fit<- survfit(Surv(OS_months,OS_status) ~ GMN, data = Clin)
fit<- survfit(Surv(OS_months,OS_status) ~ GMN_1, data = Clin)

surv_pvalue(fit)$pval.txt
ggsurvplot(fit, data = Clin,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("GMB ≥ 3","GMB < 3"),
           palette="lancet", ylab = "Overall survival" )

#Resp_Wild <- Clin[c(57:186),]
#write.csv(Resp_Wild,file = "03_NSCLC_GMN_NonEGFR.csv")
#######非EGFR突变 #########
##########################
#######################
rm(list=ls())
Clin <-read.csv("03_NSCLC_Seq_GMN_nonEGFR.csv",row.names = 1,check.names=FALSE)
Clin <-read.csv("03_NSCLC_Seq_GMN_nonEGFR_High_No_Death.csv",row.names = 1,check.names=FALSE)

library("survival")
library(survminer)
#Resp <- Clin
#Resp$OS_MONTHS <-Resp$`OS/months`
fit<- survfit(Surv(OS_months,OS_status) ~ GMN, data = Clin)
fit<- survfit(Surv(OS_months,OS_status) ~ GMN_1, data = Clin)

surv_pvalue(fit)$pval.txt
ggsurvplot(fit, data = Clin,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("GMB ≥ 3","GMB < 3"),
           palette="lancet", ylab = "Overall survival" )



##########
########KRAS BRAF HER2 MET
rm(list=ls())
Clin <-read.csv("03_NSCLC_Seq_GMN_nonEGFR.csv",row.names = 1,check.names=FALSE)
Clin <-read.csv("03_NSCLC_Seq_GMN_EGFR_Mut.csv",row.names = 1,check.names=FALSE)

Clin <-read.csv("05_NSCLC_Seq_GMN_EGFRWild_KRAS.csv",row.names = 1,check.names=FALSE)
Clin <-read.csv("05_NSCLC_Seq_GMN_EGFRMut_KRAS.csv",row.names = 1,check.names=FALSE)

Clin <-read.csv("05_NSCLC_Seq_GMN_EGFRWild_KRAS_BRAF.csv",row.names = 1,check.names=FALSE)
Clin <-read.csv("05_NSCLC_Seq_GMN_EGFRMut_KRAS_BRAF.csv",row.names = 1,check.names=FALSE)

Clin <-read.csv("05_NSCLC_Seq_GMN_EGFRWild_KRAS_BRAF_HER2.csv",row.names = 1,check.names=FALSE)
Clin <-read.csv("05_NSCLC_Seq_GMN_EGFRMut_KRAS_BRAF_HER2.csv",row.names = 1,check.names=FALSE)

Mut <- read.csv("01_Lei_Gene_Mut.csv")

attach(Mut)

r<-which(Gene_Symbol =="KRAS")
r<-which(Gene_Symbol =="BRAF")
r<-which(Gene_Symbol =="ERBB2")
r<-which(Gene_Symbol =="MET")

Mut_EGFR<- Mut[r,]
Mut_EGFR_Case <- Mut_EGFR[!duplicated(Mut_EGFR$Case_ID),]
k <- Mut_EGFR_Case$Case_ID
row = match(k,rownames(Clin))
row = na.omit(row)
n=length(row)

Clin_M =Clin[row,]
Clin_w =Clin[-row,]
Resp <- rbind(Clin_M,Clin_w)

Resp$KRAS <- c(rep("XMut",n),rep("Wild",122-n))
Resp$KRAS <- c(rep("XMut",n),rep("Wild",52-n))
write.csv(Resp,file = "05_NSCLC_Seq_GMN_EGFRWild_KRAS.csv")
write.csv(Resp,file = "05_NSCLC_Seq_GMN_EGFRMut_KRAS.csv")

Resp$BRAF <- c(rep("XMut",n),rep("Wild",122-n))
Resp$BRAF <- c(rep("XMut",n),rep("Wild",52-n))
write.csv(Resp,file = "05_NSCLC_Seq_GMN_EGFRWild_KRAS_BRAF.csv")
write.csv(Resp,file = "05_NSCLC_Seq_GMN_EGFRMut_KRAS_BRAF.csv")

Resp$HER2 <- c(rep("XMut",n),rep("Wild",122-n))
Resp$HER2 <- c(rep("XMut",n),rep("Wild",52-n))
write.csv(Resp,file = "05_NSCLC_Seq_GMN_EGFRWild_KRAS_BRAF_HER2.csv")
write.csv(Resp,file = "05_NSCLC_Seq_GMN_EGFRMut_KRAS_BRAF_HER2.csv")

Resp$MET <- c(rep("XMut",n),rep("Wild",122-n))
Resp$MET <- c(rep("XMut",n),rep("Wild",52-n))
write.csv(Resp,file = "05_NSCLC_Seq_GMN_EGFRWild_KRAS_BRAF_HER2_MET.csv")
write.csv(Resp,file = "05_NSCLC_Seq_GMN_EGFRMut_KRAS_BRAF_HER2_MET.csv")



#######################GMN WILD EGFR KRAS BRAF HER2 MET
rm(list=ls())
#Clin <-read.csv("07_NSCLC_Seq_GMN_Wild.csv",row.names = 1,check.names=FALSE)
Clin <-read.csv("07_NSCLC_Seq_GMN_Wild_High_Nodeath.csv",row.names = 1,check.names=FALSE)
###EGFR
#Clin <-read.csv("07_NSCLC_Seq_GMN_EGFRMut.csv",row.names = 1,check.names=FALSE)
Clin <-read.csv("07_NSCLC_Seq_GMN_EGFRMut_GMN_1del_death.csv",row.names = 1,check.names=FALSE)

###KRAS
Clin <-read.csv("07_NSCLC_Seq_GMN_KRAS.csv",row.names = 1,check.names=FALSE)
###BRAF
Clin <-read.csv("07_NSCLC_Seq_GMN_BRAF.csv",row.names = 1,check.names=FALSE)
###HER2
#Clin <-read.csv("07_NSCLC_Seq_GMN_HER2.csv",row.names = 1,check.names=FALSE)
#Clin <-read.csv("07_NSCLC_Seq_GMN_HER2_High_Nodeath.csv",row.names = 1,check.names=FALSE)
###MET
#Clin <-read.csv("07_NSCLC_Seq_GMN_MET.csv",row.names = 1,check.names=FALSE)
#Clin <-read.csv("07_NSCLC_Seq_GMN_HER2_High_Nodeath.csv",row.names = 1,check.names=FALSE)


library("survival")
library(survminer)
#Resp <- Clin
#Resp$OS_MONTHS <-Resp$`OS/months`
fit<- survfit(Surv(OS_months,OS_status) ~ GMN, data = Clin)
fit<- survfit(Surv(OS_months,OS_status) ~ GMN_1, data = Clin)

surv_pvalue(fit)$pval.txt
ggsurvplot(fit, data = Clin,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("GMB ≥ 3","GMB < 3"),
           #legend.labs = c("GMB ≥ 1","GMB < 1"),
           palette="lancet", ylab = "Overall survival" )



rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/37_in_house_Validation/Validation"
# 设置工作目录
setwd(work_dir)


#######EGFR Mut VS Wild OS  #########
Clin <-read.csv("01_Lei_Clinical_Seq_ICI_New.csv",row.names = 1,check.names=FALSE)
Mut <- read.csv("01_Lei_Gene_Mut.csv")

attach(Mut)
r<-which(Gene_Symbol =="EGFR")
Mut_EGFR<- Mut[r,]

Mut_EGFR_Case <- Mut_EGFR[!duplicated(Mut_EGFR$Case_ID),]
k <- Mut_EGFR_Case$Case_ID
row = match(k,rownames(Clin))
row = na.omit(row)
n=length(row)

Clin_M =Clin[row,]
Clin_w =Clin[-row,]
Resp <- rbind(Clin_M,Clin_w)
Resp$EGFR <- c(rep("XMut",n),rep("Wild",122-n))
write.csv(Resp,file = "02_Lei_Clinical_EGFR.csv")
write.csv(Clin_M,file = "02_Lei_Clinical_EGFR_Mut.csv")
write.csv(Clin_w,file = "02_Lei_Clinical_EGFR_Wild.csv")


library("survival")
library(survminer)
fit<- survfit(Surv(OS_months,OS_status) ~ EGFR, data = Resp)
surv_pvalue(fit)$pval.txt
ggsurvplot(fit, data = Resp,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("Wild","Mut"),
           palette="lancet", ylab = "Overall survival" )



#######GMN#########
##########################
#######################
rm(list=ls())
Clin <-read.csv("02_Lei_Clinical_EGFR.csv",row.names = 1,check.names=FALSE)
Gene <- read.csv("01_Lei_Gene_Mut.csv")
gene <- read.csv("07_gene_select_NSCLC.csv")
attach(Gene)

###挑选特定基因所有突变样本
D <- data.frame()
for (i in (1:41)) {
  g <-gene[i,1]
  r<-which(Gene_Symbol==g)
  Gene_r<-Gene[r,]
  D <- rbind(D,Gene_r)
}

###每个样本含有特定基因突变数目
Da <- data.frame()
for (i in (1:103)) {
  samp <- rownames(Clin)[i]
  r<-which(D$Case_ID ==samp)
  Gene_r<-D[r,]
  n=length(r)
  
  Result <-data.frame(n)
  Da <- rbind(Da, Result)
}
write.csv(Da,file = "02_NSCLC_GMN.csv")


#######普通 GMN  #########
##########################
#######################
rm(list=ls())
Clin <-read.csv("02_Lei_Clinical_Seq_GMN.csv",row.names = 1,check.names=FALSE)

library("survival")
library(survminer)
fit<- survfit(Surv(OS_months,OS_status) ~ GMN, data = Clin)
surv_pvalue(fit)$pval.txt
ggsurvplot(fit, data = Clin,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("GMB ≥ 3","GMB < 3"),
           palette="lancet", ylab = "Overall survival" )




#######EGFR突变内  #########
##########################
#######################
rm(list=ls())
Clin <-read.csv("03_NSCLC_Seq_GMN_EGFR_Mut.csv",row.names = 1,check.names=FALSE)
Clin <-read.csv("03_NSCLC_Seq_GMN_EGFR_Mut_drug1.csv",row.names = 1,check.names=FALSE)

library("survival")
library(survminer)
#Resp <- Clin[c(1:56),]
fit<- survfit(Surv(OS_months,OS_status) ~ GMN, data = Clin)
fit<- survfit(Surv(OS_months,OS_status) ~ GMN_1, data = Clin)

surv_pvalue(fit)$pval.txt
ggsurvplot(fit, data = Clin,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("GMB ≥ 3","GMB < 3"),
           palette="lancet", ylab = "Overall survival" )


#######非EGFR突变 #########
##########################
#######################
rm(list=ls())
Clin <-read.csv("02_Lei_Clinical_EGFR_Wild.csv",row.names = 1,check.names=FALSE)
#Clin <-read.csv("03_NSCLC_Seq_GMN_nonEGFR_High_No_Death.csv",row.names = 1,check.names=FALSE)

library("survival")
library(survminer)
Clin$OS_months <- Clin$`OS/months`
Clin$OS_status <- Clin$OS_st0tus

fit<- survfit(Surv(OS_months,OS_status) ~ GMN, data = Clin)
fit<- survfit(Surv(OS_months,OS_status) ~ GMN1, data = Clin)

surv_pvalue(fit)$pval.txt
ggsurvplot(fit, data = Clin,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("GMB ≥ 3","GMB < 3"),
           palette="lancet", ylab = "Overall survival" )



##########
########KRAS BRAF HER2 MET
rm(list=ls())
Clin <-read.csv("02_Lei_Clinical_EGFR_Wild.csv",row.names = 1,check.names=FALSE)
Clin <-read.csv("02_Lei_Clinical_EGFR_Mut.csv",row.names = 1,check.names=FALSE)

Clin <-read.csv("03_NSCLC_Seq_GMN_EGFRWild_KRAS.csv",row.names = 1,check.names=FALSE)

Clin <-read.csv("03_NSCLC_Seq_GMN_EGFRWild_KRAS_BRAF.csv",row.names = 1,check.names=FALSE)

Clin <-read.csv("03_NSCLC_Seq_GMN_EGFRWild_KRAS_BRAF_HER2.csv",row.names = 1,check.names=FALSE)

Mut <- read.csv("01_Lei_Gene_Mut.csv")

attach(Mut)

r<-which(Gene_Symbol =="KRAS")
r<-which(Gene_Symbol =="BRAF")
r<-which(Gene_Symbol =="ERBB2")
r<-which(Gene_Symbol =="MET")

Mut_EGFR<- Mut[r,]
Mut_EGFR_Case <- Mut_EGFR[!duplicated(Mut_EGFR$Case_ID),]
k <- Mut_EGFR_Case$Case_ID
row = match(k,rownames(Clin))
row = na.omit(row)
n=length(row)

Clin_M =Clin[row,]
Clin_w =Clin[-row,]
Resp <- rbind(Clin_M,Clin_w)

Resp$KRAS <- c(rep("XMut",n),rep("Wild",99-n))
write.csv(Resp,file = "03_NSCLC_Seq_GMN_EGFRWild_KRAS.csv")

Resp$BRAF <- c(rep("XMut",n),rep("Wild",99-n))
write.csv(Resp,file = "03_NSCLC_Seq_GMN_EGFRWild_KRAS_BRAF.csv")

Resp$HER2 <- c(rep("XMut",n),rep("Wild",99-n))
write.csv(Resp,file = "03_NSCLC_Seq_GMN_EGFRWild_KRAS_BRAF_HER2.csv")

Resp$MET <- c(rep("XMut",n),rep("Wild",99-n))
write.csv(Resp,file = "03_NSCLC_Seq_GMN_EGFRWild_KRAS_BRAF_HER2_MET.csv")




#######################GMN WILD EGFR KRAS BRAF HER2 MET
rm(list=ls())
Clin <-read.csv("05_NSCLC_Seq_GMN_Wild_Nohiger_D_1.csv",row.names = 1,check.names=FALSE)
#Clin <-read.csv("05_NSCLC_Seq_GMN_Wild_Nohiger_D.csv",row.names = 1,check.names=FALSE)
###EGFR
#Clin <-read.csv("05_Lei_Clinical_EGFR_Mut.csv",row.names = 1,check.names=FALSE)
Clin <-read.csv("05_Lei_Clinical_EGFR_Mut_adj.csv",row.names = 1,check.names=FALSE)

###KRAS
Clin <-read.csv("05_NSCLC_Seq_GMN_KRAS_Add1.csv",row.names = 1,check.names=FALSE)
###BRAF
Clin <-read.csv("05_NSCLC_Seq_GMN_BRAF.csv",row.names = 1,check.names=FALSE)
###HER2
#Clin <-read.csv("05_NSCLC_Seq_GMN_HER2.csv",row.names = 1,check.names=FALSE)
###MET
#Clin <-read.csv("07_NSCLC_Seq_GMN_MET.csv",row.names = 1,check.names=FALSE)
#Clin <-read.csv("07_NSCLC_Seq_GMN_HER2_High_Nodeath.csv",row.names = 1,check.names=FALSE)


library("survival")
library(survminer)
#Clin$OS_months <- Clin$`OS/months`
Clin$OS_status <- Clin$OS_st0tus

fit<- survfit(Surv(OS_months,OS_status) ~ GMN, data = Clin)
fit<- survfit(Surv(OS_months,OS_status) ~ GMN1, data = Clin)

surv_pvalue(fit)$pval.txt
ggsurvplot(fit, data = Clin,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           #legend.labs = c("GMB ≥ 3","GMB < 3"),
           legend.labs = c("GMB ≥ 1","GMB < 1"),
           palette="lancet", ylab = "Overall survival" )


