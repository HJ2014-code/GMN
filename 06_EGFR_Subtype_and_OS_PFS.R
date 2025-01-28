#######6.1 EGFR and OS
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/34_EGFR_Subtype"
# 设置工作目录
setwd(work_dir)

clinical <- read.csv("01_Non-Small_Cell_Lung_Cancer_clinical.csv",row.names = 1,check.names=FALSE)
Mutation <-  read.csv("02_data_mutations_nonsynonymous.csv")

###NSCLC含有的突变基因数量
D <- data.frame()
for (i in (1:350)) {
  samp <-clinical[i,6]
  r<-which(Mutation$Tumor_Sample_Barcode==samp)
  Mutation_r<-Mutation[r,]
  D <- rbind(D,Mutation_r)
}


###挑选特定基因所有突变样本
g <-"EGFR"
attach(D)
r<-which(Hugo_Symbol==g)
Gene_r<-D[r,]
table(Gene_r$HGVSp_Short)
write.csv(Gene_r,file = "03_EGFR_Mut_subtype.csv")

Gene_r <- read.csv("03_EGFR_Mut_subtype.csv",row.names = 1,check.names=FALSE)
table(Gene_r$HGVSp_Short_Subtype)

Gene_r_subtype <- Gene_r[Gene_r$HGVSp_Short_Subtype == "p.L858R",]
Gene_r_subtype <- Gene_r[Gene_r$HGVSp_Short_Subtype == "p.E746_A750del",]
Gene_r_subtype <- Gene_r[Gene_r$HGVSp_Short_Subtype == "Other",]

Sample_EGFR_Mut <- Gene_r_subtype[!duplicated(Gene_r_subtype$Tumor_Sample_Barcode),]

k <- Sample_EGFR_Mut$Tumor_Sample_Barcode
row = match(k,clinical$Tumor_Sample_Barcode)
row = na.omit(row)
clinical_Mut = clinical[row,]


clinical_Wild = clinical[-row,]

row = match(k,clinical_Wild$Tumor_Sample_Barcode)
row = na.omit(row)
clinical_Wild = clinical_Wild[-row,]

write.csv(clinical_Mut,file = "05_EGFR_Mut_p.L858R.csv")
write.csv(clinical_Mut,file = "05_EGFR_Mut_p.E746_A750del.csv")
write.csv(clinical_Mut,file = "05_EGFR_Mut_Other.csv")
write.csv(clinical_Wild,file = "05_EGFR_Wild.csv")


###EGFR Mut Subtype
rm(list=ls())
clinical <- read.csv("06_EGFR_Mut_Subtype.csv")
library(survminer)
library("survival")

fit<- survfit(Surv(OS_MONTHS, STATUS) ~ EGFR_Subtype, data = clinical)
surv_pvalue(fit)$pval.txt
ggsurvplot(fit, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           #legend.labs = c("KRAS_Mut", "KRAS_Wild"),
           palette="lancet", ylab = "Overall survival" )
###添加HR ,CI ,P
p3 <- ggsurvplot(fit, data = clinical,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 #legend.labs = c("KRAS_Mut", "KRAS_Wild"),
                 palette="jto", ylab = "Overall survival" )
res_cox<-coxph(Surv(OS_MONTHS, STATUS) ~ EGFR_Subtype, data = clinical)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3



#####EGFR Mut subtype GMB#####
##################
#################
rm(list=ls())
clinical <- read.csv("05_EGFR_Mut_p.L858R.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("05_EGFR_Mut_p.E746_A750del.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("05_EGFR_Mut_Other.csv",row.names = 1,check.names=FALSE)


library(survminer)
library("survival")

#####surv-cutpoint OS决定最佳分界值
res.cut <- surv_cutpoint(clinical,time="OS_MONTHS",event = "STATUS",
                         variables = "n")
summary(res.cut)   
plot(res.cut,"n",palette="npg")

res.cat <- surv_categorize(res.cut)

fit<- survfit(Surv(OS_MONTHS, STATUS) ~ n, data = res.cat)
surv_pvalue(fit)$pval.txt
clinical$cut <- res.cat$n
write.csv(clinical,file = "07_EGFR_Mut_p.L858R_GMB_1.csv")
write.csv(clinical,file = "07_EGFR_Mut_p.E746_A750del_1.csv")
write.csv(clinical,file = "07_EGFR_Mut_Other_GMB.csv")

ggsurvplot(fit, data = res.cat,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("Gene_mutation_burden ≥ 3", "Gene_mutation_burden < 3"),
           palette="lancet", ylab = "Overall survival" )
###添加HR ,CI ,P
p3 <- ggsurvplot(fit, data = res.cat,
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



###EGFR_Mut里分界值为3
#################
##############
rm(list=ls())
clinical <- read.csv("08_EGFR_Mut_p.L858R_GMB_3.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("08_EGFR_Mut_p.E746_A750del_3.csv",row.names = 1,check.names=FALSE)
library(survminer)
library("survival")

fit<- survfit(Surv(OS_MONTHS, STATUS) ~ cut, data = clinical)
surv_pvalue(fit)$pval.txt
ggsurvplot(fit, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("Gene_mutation_burden ≥ 3", "Gene_mutation_burden < 3"),
           palette="lancet", ylab = "Overall survival" )
###添加HR ,CI ,P
p3 <- ggsurvplot(fit, data = clinical,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 legend.labs = c("Gene_mutation_burden ≥ 3", "Gene_mutation_burden < 3"),
                 palette="jto", ylab = "Overall survival" )
res_cox<-coxph(Surv(OS_MONTHS, STATUS) ~ cut, data = clinical)
summary(res_cox)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3



#######6.2 EGFR and PFS
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/35_EGFR_Subtype_Validation"
# 设置工作目录
setwd(work_dir)

clinical <- read.csv("03_75_34_35all_clinical_TMB_.csv",row.names = 1,check.names=FALSE)
Mutation <-  read.csv("03_75_34_35_data_mutations.csv")

###NSCLC含有的突变基因数量 Mutation 35来自240那个队列
D <- data.frame()
for (i in (1:144)) {
  samp <-clinical[i,11]
  r<-which(Mutation$Tumor_Sample_Barcode==samp)
  Mutation_r<-Mutation[r,]
  D <- rbind(D,Mutation_r)
}

###挑选特定基因所有突变样本
g <-"EGFR"
attach(D)
r<-which(Hugo_Symbol==g)
Gene_r<-D[r,]
table(Gene_r$HGVSp_Short)
write.csv(Gene_r,file = "05_EGFR_Mut_subtype.csv")

Gene_r <- read.csv("05_EGFR_Mut_subtype.csv",row.names = 1,check.names=FALSE)
table(Gene_r$HGVSp_Short_Subtype)

Gene_r_subtype <- Gene_r[Gene_r$HGVSp_Short_Subtype == "p.L858R",]
Gene_r_subtype <- Gene_r[Gene_r$HGVSp_Short_Subtype == "p.E746_A750del",]
Gene_r_subtype <- Gene_r[Gene_r$HGVSp_Short_Subtype == "Other",]


Sample_EGFR_Mut <- Gene_r_subtype[!duplicated(Gene_r_subtype$Tumor_Sample_Barcode),]

k <- Sample_EGFR_Mut$Tumor_Sample_Barcode
row = match(k,clinical$Tumor_Sample_Barcode)
row = na.omit(row)
clinical_Mut = clinical[row,]

clinical_Wild = clinical[-row,]

row = match(k,clinical_Wild$Tumor_Sample_Barcode)
row = na.omit(row)
clinical_Wild = clinical_Wild[-row,]

write.csv(clinical_Mut,file = "06_EGFR_Mut_p.L858R.csv")
write.csv(clinical_Mut,file = "06_EGFR_Mut_p.E746_A750del.csv")
write.csv(clinical_Mut,file = "06_EGFR_Mut_Other.csv")
write.csv(clinical_Wild,file = "06_EGFR_Wild.csv")


###EGFR Subtype VS Wild#######
########################
#######################
rm(list=ls())
clinical <- read.csv("07_EGFR_Mut_Subtype.csv")
library(survminer)
library("survival")

fit<- survfit(Surv(PFS_MONTHS,Status) ~ EGFR_Subtype, data = clinical)
surv_pvalue(fit)$pval.txt
pairwise_survdiff(Surv(PFS_MONTHS,Status) ~ EGFR_Subtype, data = clinical) 
#p.adjust.method = 'BH') #采用log-rank 检验分析生存率差异
ggsurvplot(fit, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           #legend.labs = c("KRAS_Mut", "KRAS_Wild"),
           palette="lancet", ylab = "Progress-free survival" )
###添加HR ,CI ,P
p3 <- ggsurvplot(fit, data = clinical,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 #legend.labs = c("KRAS_Mut", "KRAS_Wild"),
                 palette="jto", ylab = "Progress-free survival" )
res_cox<-coxph(Surv(PFS_MONTHS,Status) ~ EGFR_Subtype, data = clinical)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3



#####EGFR Mut subtype GMB#####
##################
#################
rm(list=ls())
clinical <- read.csv("06_EGFR_Mut_p.L858R.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("06_EGFR_Mut_p.E746_A750del.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("06_EGFR_Mut_Other.csv",row.names = 1,check.names=FALSE)

library(survminer)
library("survival")

#####surv-cutpoint OS决定最佳分界值
res.cut <- surv_cutpoint(clinical,time="PFS_MONTHS",event = "Status",
                         variables = "n")
summary(res.cut)   
plot(res.cut,"n",palette="npg")

res.cat <- surv_categorize(res.cut)

fit<- survfit(Surv(PFS_MONTHS,Status) ~ n, data = res.cat)
surv_pvalue(fit)$pval.txt
clinical$cut <- res.cat$n
write.csv(clinical,file = "08_EGFR_Mut_Other_GMB_1.csv")

ggsurvplot(fit, data = res.cat,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("Gene_mutation_burden ≥ 1", "Gene_mutation_burden < 1"),
           palette="lancet", ylab = "Progress-free survival" )
###添加HR ,CI ,P
p3 <- ggsurvplot(fit, data = res.cat,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 legend.labs = c("Gene_mutation_burden ≥ 1", "Gene_mutation_burden < 1"),
                 palette="jto", ylab = "Progress-free survival" )
res_cox<-coxph(Surv(PFS_MONTHS,Status) ~ n, data = res.cat)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3



###Mut Subtype 里分界值为3###
####################
#####################
rm(list=ls())
clinical <- read.csv("08_EGFR_Mut_Other_GMB_1.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("06_EGFR_Mut_p.E746_A750del.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("06_EGFR_Mut_p.L858R.csv",row.names = 1,check.names=FALSE)

library(survminer)
library("survival")

fit<- survfit(Surv(PFS_MONTHS,Status) ~ Cut, data = clinical)
surv_pvalue(fit)$pval.txt
ggsurvplot(fit, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("Gene_mutation_burden ≥ 3", "Gene_mutation_burden < 3"),
           palette="lancet", ylab = "Progress-free survival" )
###添加HR ,CI ,P
p3 <- ggsurvplot(fit, data = clinical,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 legend.labs = c("Gene_mutation_burden ≥ 3", "Gene_mutation_burden < 3"),
                 palette="jto", ylab = "Progress-free survival" )
res_cox<-coxph(Surv(PFS_MONTHS,Status) ~ Cut, data = clinical)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3






