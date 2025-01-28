######4.1 KRAS and OS
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/21_KRAS_Immune_GMB"
# 设置工作目录
setwd(work_dir)

clinical <- read.csv("01_Non-Small_Cell_Lung_Cancer_clinical.csv",row.names = 1,check.names=FALSE)
#Gene <- read.csv("02_Non-Small_Cell_Lung_Cancer_mutation.csv")
Mutation <-  read.csv("02_data_mutations_nonsynonymous.csv")
gene <- read.csv("03_gene_select_NSCLC.csv")

###NSCLC含有的突变基因数量
D <- data.frame()
for (i in (1:350)) {
  samp <-clinical[i,6]
  r<-which(Mutation$Tumor_Sample_Barcode==samp)
  Mutation_r<-Mutation[r,]
  D <- rbind(D,Mutation_r)
}
#data <- D[!duplicated(D$Hugo_Symbol),]

###挑选特定基因所有突变样本
g <-"KRAS"
attach(D)
r<-which(Hugo_Symbol==g)
Gene_r<-D[r,]
Sample_KRAS_Mut <- Gene_r[!duplicated(Gene_r$Tumor_Sample_Barcode),]

k <- Sample_KRAS_Mut$Tumor_Sample_Barcode
row = match(k,clinical$Tumor_Sample_Barcode)
row = na.omit(row)
clinical_Mut = clinical[row,]
clinical_Wild = clinical[-row,]
write.csv(clinical_Mut,file = "05_KRAS_Mut.csv")
write.csv(clinical_Wild,file = "05_KRAS_Wild.csv")


###KRAS Mut VS Wild
rm(list=ls())
clinical <- read.csv("05_KRAS_Mut_Wild.csv",row.names = 1,check.names=FALSE)
library(survminer)
library("survival")

fit<- survfit(Surv(OS_MONTHS, STATUS) ~ KRAS_Status, data = clinical)
surv_pvalue(fit)$pval.txt
ggsurvplot(fit, data = clinical,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("KRAS_Mut", "KRAS_Wild"),
           palette="lancet", ylab = "Overall survival" )
###添加HR ,CI ,P
p3 <- ggsurvplot(fit, data = clinical,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 legend.labs = c("KRAS_Mut", "KRAS_Wild"),
                 palette="jto", ylab = "Overall survival" )
res_cox<-coxph(Surv(OS_MONTHS, STATUS) ~ KRAS_Status, data = clinical)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3



#####KRAS Mut GMB
rm(list=ls())
clinical <- read.csv("05_KRAS_Mut.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("05_KRAS_Wild.csv",row.names = 1,check.names=FALSE)
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
write.csv(clinical,file = "06_KRAS_Mut_GMB.csv")
write.csv(clinical,file = "06_KRAS_Wild_GMB.csv")

ggsurvplot(fit, data = res.cat,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("Gene_mutation_burden ≥ 2", "Gene_mutation_burden < 2"),
           palette="lancet", ylab = "Overall survival" )
###添加HR ,CI ,P
p3 <- ggsurvplot(fit, data = res.cat,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 legend.labs = c("Gene_mutation_burden ≥ 2", "Gene_mutation_burden < 2"),
                 palette="jto", ylab = "Overall survival" )
res_cox<-coxph(Surv(OS_MONTHS, STATUS) ~ n, data = res.cat)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3



###Wild 里分界值为3
###KRAS Mut VS Wild
rm(list=ls())
clinical <- read.csv("07_KRAS_Wild_GMB_3.csv",row.names = 1,check.names=FALSE)
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
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3



######4.2 KRAS and PFS
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/21_KRAS_Immune_GMB"
# 设置工作目录
setwd(work_dir)

clinical <- read.csv("05_KRAS_Mut_Wild.csv",row.names = 1,check.names=FALSE)
PFS <- read.csv("09_205_PFS.csv")

k <- PFS$X
row = match(k,rownames(clinical))
row = na.omit(row)
clinical_Res = clinical[row,]

clinical_Res$PFS_Status = PFS$Status
clinical_Res$PFS_MONTHS = PFS$PFS_MONTHS

library(survminer)
library("survival")

fit<- survfit(Surv(PFS_MONTHS, PFS_Status) ~ KRAS_Status, data = clinical_Res)
surv_pvalue(fit)$pval.txt
ggsurvplot(fit, data = clinical_Res,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("KRAS_Mut", "KRAS_Wild"),
           palette="lancet", ylab = "Progress-free survival" )
###添加HR ,CI ,P
p3 <- ggsurvplot(fit, data = clinical_Res,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 legend.labs = c("KRAS_Mut", "KRAS_Wild"),
                 palette="jto", ylab = "Progress-free survival" )
res_cox<-coxph(Surv(PFS_MONTHS, PFS_Status) ~ KRAS_Status, data = clinical_Res)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3



#####KRAS Mut GMB########
#####################
####################
rm(list=ls())
clinical <- read.csv("05_KRAS_Mut.csv",row.names = 1,check.names=FALSE)
PFS <- read.csv("09_205_PFS.csv")

k <- rownames(clinical) 
row = match(k,PFS$X)
row = na.omit(row)
PFS_final = PFS[row,]

k <- PFS_final$X
row = match(k,rownames(clinical))
row = na.omit(row)
clinical_Res = clinical[row,]


rm(list=ls())
clinical <- read.csv("05_KRAS_Wild.csv",row.names = 1,check.names=FALSE)
PFS <- read.csv("09_205_PFS.csv")

k <- PFS$X
row = match(k,rownames(clinical))
row = na.omit(row)
clinical_Res = clinical[row,]

k <- rownames(clinical_Res) 
row = match(k,PFS$X)
row = na.omit(row)
PFS_final = PFS[row,]

clinical_Res$PFS_Status = PFS_final$Status
clinical_Res$PFS_MONTHS = PFS_final$PFS_MONTHS

library(survminer)
library("survival")

#####surv-cutpoint OS决定最佳分界值
res.cut <- surv_cutpoint(clinical_Res,time="PFS_MONTHS",event = "PFS_Status",
                         variables = "n")
summary(res.cut)   
plot(res.cut,"n",palette="npg")

res.cat <- surv_categorize(res.cut)

fit<- survfit(Surv(PFS_MONTHS, PFS_Status) ~ n, data = res.cat)
surv_pvalue(fit)$pval.txt
clinical_Res$cut <- res.cat$n
write.csv(clinical_Res,file = "09_KRAS_Mut_GMB_PFS.csv")
write.csv(clinical_Res,file = "09_KRAS_Wild_GMB_PFS_1.csv")

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
res_cox<-coxph(Surv(PFS_MONTHS, PFS_Status) ~ n, data = res.cat)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3



###Wild 里分界值为3
####################
######################
rm(list=ls())
clinical <- read.csv("10_KRAS_Wild_GMB_PFS_3.csv",row.names = 1,check.names=FALSE)
library(survminer)
library("survival")

fit<- survfit(Surv(PFS_MONTHS, PFS_Status) ~ cut, data = clinical)
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
res_cox<-coxph(Surv(PFS_MONTHS, PFS_Status) ~ cut, data = clinical)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3



######4.3 KRAS and Response
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/21_KRAS_Immune_GMB"
# 设置工作目录
setwd(work_dir)

clinical <- read.csv("05_KRAS_Mut_Wild.csv",row.names = 1,check.names=FALSE)
response <- read.csv("08_MSKCC_NTRK3_response.csv")

k <- response$PATIENT_ID
row = match(k,rownames(clinical))
row = na.omit(row)
clinical_Res = clinical[row,]

clinical_Res$Best_treatment_response = response$Best_treatment_response
####免疫应答
table(clinical_Res$Best_treatment_response,clinical_Res$KRAS_Status)
fisher.test(clinical_Res$Best_treatment_response,clinical_Res$KRAS_Status)

# bar_plot
df<-clinical_Res[,c(13,14)]
df$count <-c(1)
ggplot(df, aes(x=KRAS_Status, y=count))+
  geom_bar(stat="identity", position="fill", aes(fill=Best_treatment_response)) +
  scale_y_continuous(labels = scales::percent)


###Mut GMB Response#######
#############################
##########################
rm(list=ls())
clinical <- read.csv("06_KRAS_Mut_GMB.csv",row.names = 1,check.names=FALSE)
response <- read.csv("08_MSKCC_NTRK3_response.csv")

k <- rownames(clinical) 
row = match(k,response$PATIENT_ID)
row = na.omit(row)
response_final = response[row,]

k <- response_final$PATIENT_ID
row = match(k,rownames(clinical))
row = na.omit(row)
clinical_Res = clinical[row,]

clinical_Res$Best_treatment_response = response_final$Best_treatment_response
####免疫应答
table(clinical_Res$Best_treatment_response,clinical_Res$cut)
fisher.test(clinical_Res$Best_treatment_response,clinical_Res$cut)

# bar_plot
df<-clinical_Res[,c(13,14)]
df$count <-c(1)
ggplot(df, aes(x=cut, y=count))+
  geom_bar(stat="identity", position="fill", aes(fill=Best_treatment_response)) +
  scale_y_continuous(labels = scales::percent)

###Waterfall plots
#安装R包
install.packages("waterfalls")
library(waterfalls)

#读取数据
data<-clinical_Res[,c(12,14)]
write.csv(data,"08_KRAS_Mut_GMB_ORR_Waterfall.csv")
data <- read.csv("08_KRAS_Mut_GMB_ORR_Waterfall.csv",row.names = 1,check.names=FALSE)
View(data)

#基本作图
barplot(data$number,col="blue", border="blue", space=0.5, ylim=c(-5,15),
        #图标题
        #main = "Waterfall plot for gene mutation number",
        ylab="Gene mutation number",cex.axis=1.2, cex.lab=1.4)

# 美化
col <- ifelse(data$Best_treatment_response == "CR/PR", "#BC5A42", "#009296")
barplot(data$number, col=col, border=col, space=0.5, ylim=c(-5,15),
        #main = "Waterfall plot for changes in QoL scores", 
        ylab="Gene mutation number",
        #xlab="ID",
        cex.axis=1.2, cex.lab=1.4, legend.text=c('CR/PR','PD/SD'),
        #图例显示及设定
        args.legend=list(title="Best_treatment_response", fill=c("#BC5A42", "#009296"), border=NA, cex=0.9,bty = "n"))



###Wild GMB Response#####
######################
########################
rm(list=ls())
clinical <- read.csv("06_KRAS_Wild_GMB.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("07_KRAS_Wild_GMB_3.csv",row.names = 1,check.names=FALSE)
response <- read.csv("08_MSKCC_NTRK3_response.csv")

k <- response$PATIENT_ID
row = match(k,rownames(clinical))
row = na.omit(row)
clinical_Res = clinical[row,]

k <- rownames(clinical_Res) 
row = match(k,response$PATIENT_ID)
row = na.omit(row)
response_final = response[row,]

clinical_Res$Best_treatment_response = response_final$Best_treatment_response
####免疫应答
table(clinical_Res$Best_treatment_response,clinical_Res$cut)
fisher.test(clinical_Res$Best_treatment_response,clinical_Res$cut)
# bar_plot
df<-clinical_Res[,c(13,14)]
df$count <-c(1)
ggplot(df, aes(x=cut, y=count))+
  geom_bar(stat="identity", position="fill", aes(fill=Best_treatment_response)) +
  scale_y_continuous(labels = scales::percent)










