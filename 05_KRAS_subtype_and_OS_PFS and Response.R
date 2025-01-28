######5.1 KRAS subtype and OS
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/25_KRAS_subtype_GMB"
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
g <-"KRAS"
attach(D)
r<-which(Hugo_Symbol==g)
Gene_r<-D[r,]
table(Gene_r$HGVSp_Short)
write.csv(Gene_r,file = "03_KRAS_Mut_subtype.csv")

Gene_r <- read.csv("03_KRAS_Mut_subtype.csv",row.names = 1,check.names=FALSE)
table(Gene_r$HGVSp_Short_Subtype)

Gene_r_subtype <- Gene_r[Gene_r$HGVSp_Short_Subtype == "p.G12C",]
Gene_r_subtype <- Gene_r[Gene_r$HGVSp_Short_Subtype == "p.G12V",]
Gene_r_subtype <- Gene_r[Gene_r$HGVSp_Short_Subtype == "p.G12D",]
Gene_r_subtype <- Gene_r[Gene_r$HGVSp_Short_Subtype == "Other",]

Sample_KRAS_Mut <- Gene_r_subtype[!duplicated(Gene_r_subtype$Tumor_Sample_Barcode),]

k <- Sample_KRAS_Mut$Tumor_Sample_Barcode
row = match(k,clinical$Tumor_Sample_Barcode)
row = na.omit(row)
clinical_Mut = clinical[row,]


clinical_Wild = clinical[-row,]

row = match(k,clinical_Wild$Tumor_Sample_Barcode)
row = na.omit(row)
clinical_Wild = clinical_Wild[-row,]

write.csv(clinical_Mut,file = "05_KRAS_Mut_p.G12C.csv")
write.csv(clinical_Mut,file = "05_KRAS_Mut_p.G12V.csv")
write.csv(clinical_Mut,file = "05_KRAS_Mut_p.G12D.csv")
write.csv(clinical_Mut,file = "05_KRAS_Mut_Other.csv")
write.csv(clinical_Wild,file = "05_KRAS_Wild.csv")


###KRAS Mut Subtype
rm(list=ls())
clinical <- read.csv("06_KRAS_Mut_Subtype.csv")
library(survminer)
library("survival")

fit<- survfit(Surv(OS_MONTHS, STATUS) ~ KRAS_Subtype, data = clinical)
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
res_cox<-coxph(Surv(OS_MONTHS, STATUS) ~ KRAS_Subtype, data = clinical)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3



#####KRAS Mut subtype GMB#####
##################
#################
rm(list=ls())
clinical <- read.csv("05_KRAS_Mut_p.G12C.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("05_KRAS_Mut_p.G12D.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("05_KRAS_Mut_p.G12V.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("05_KRAS_Mut_Other.csv",row.names = 1,check.names=FALSE)


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
write.csv(clinical,file = "07_KRAS_Mut_p.G12C_GMB.csv")
write.csv(clinical,file = "07_KRAS_Mut_p.G12D_GMB.csv")
write.csv(clinical,file = "07_KRAS_Mut_p.G12V_GMB.csv")
write.csv(clinical,file = "07_KRAS_Mut_Other_GMB_2.csv")

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



###KRAS_Mut_Other 里分界值为3
#################
##############
rm(list=ls())
clinical <- read.csv("07_KRAS_Mut_Other_GMB_3.csv",row.names = 1,check.names=FALSE)
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



###Mut Subtype Response#######
#############################
##########################
rm(list=ls())
clinical <- read.csv("05_KRAS_Mut_p.G12C.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("05_KRAS_Mut_p.G12D.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("05_KRAS_Mut_p.G12V.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("05_KRAS_Mut_Other.csv",row.names = 1,check.names=FALSE)

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


######## 可视化###
##########
#########
clinical_Res=clinical_Res[order(clinical_Res$Best_treatment_response),]
table(clinical_Res$Best_treatment_response)

##云雨图
######
{
  library(ggplot2)
  library(grid)
  library(RColorBrewer)
  library(dplyr)
  library(SuppDists) 
}

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x,
                     xmax = x + width / 2)
            
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, xminv = x,
                              xmaxv = x + violinwidth * (xmax - x)) #利用transform函数为数据框mydata增加数据
            
            newdata <- rbind(plyr::arrange(transform(data, x = xmaxv), -y),plyr::arrange(transform(data, x = xminv), y))
            newdata_Polygon <- rbind(newdata, newdata[1,])
            newdata_Polygon$colour<-NA
            
            newdata_Path <- plyr::arrange(transform(data, x = xmaxv), -y)
            
            ggplot2:::ggname("geom_flat_violin", grobTree(
              GeomPolygon$draw_panel(newdata_Polygon, panel_scales, coord),
              GeomPath$draw_panel(newdata_Path, panel_scales, coord))
            )
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),
          
          required_aes = c("x", "y")
  )

library(dplyr)
library(ggpubr)
attach(clinical_Res)
clinical_Res<- clinical_Res[!is.na(clinical_Res$n),]
clinical_Res$Best_treatment_response

d <- group_by(clinical_Res, Best_treatment_response) %>%  summarize(mean = mean(n),sd = sd(n))

compaired <- list(c("CR/PR", "PD/SD"))
ggplot(clinical_Res, aes(Best_treatment_response,n, fill=Best_treatment_response))  +  
  geom_flat_violin(position=position_nudge(x=0.2)) + 
  geom_jitter(aes(color=Best_treatment_response), width=0.1) +  
  geom_pointrange(aes(y=mean, ymin=mean-sd, ymax=mean+sd),data=d, size=1, position=position_nudge(x=.2)) +
  stat_compare_means(comparisons = compaired,method = "wilcox.test",label = "p.value")+#,label.y =90
  #coord_flip() +  
  theme_bw() +
  ylab("Gene Mutation Number")+xlab(" ")+  
  theme(axis.text = element_text(size=13), axis.title =  element_text(size=15), legend.position="none") 



######5.2 KRAS subtype and PFS
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/27_KRAS_subtype_GMB_PFS"
# 设置工作目录
setwd(work_dir)

clinical <- read.csv("06_KRAS_Mut_Subtype.csv")
PFS <- read.csv("09_205_PFS.csv")

k <- PFS$X
row = match(k,clinical$X)
row = na.omit(row)
clinical_Res = clinical[row,]

clinical_Res$PFS_Status = PFS$Status
clinical_Res$PFS_MONTHS = PFS$PFS_MONTHS

library(survminer)
library("survival")

fit<- survfit(Surv(PFS_MONTHS, PFS_Status) ~ KRAS_Subtype, data = clinical_Res)
surv_pvalue(fit)$pval.txt
ggsurvplot(fit, data = clinical_Res,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           #legend.labs = c("KRAS_Mut", "KRAS_Wild"),
           palette="lancet", ylab = "Progress-free survival" )
###添加HR ,CI ,P
p3 <- ggsurvplot(fit, data = clinical_Res,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 #legend.labs = c("KRAS_Mut", "KRAS_Wild"),
                 palette="jto", ylab = "Progress-free survival" )
res_cox<-coxph(Surv(PFS_MONTHS, PFS_Status) ~ KRAS_Subtype, data = clinical_Res)
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
clinical <- read.csv("07_KRAS_Mut_p.G12C_GMB.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("07_KRAS_Mut_p.G12D_GMB.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("07_KRAS_Mut_p.G12V_GMB.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("07_KRAS_Mut_Other_GMB_3.csv",row.names = 1,check.names=FALSE)

PFS <- read.csv("09_205_PFS.csv")

k <- rownames(clinical) 
row = match(k,PFS$X)
row = na.omit(row)
PFS_final = PFS[row,]

k <- PFS_final$X
row = match(k,rownames(clinical))
row = na.omit(row)
clinical_Res = clinical[row,]

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
write.csv(clinical_Res,file = "10_KRAS_p.G12C_GMB_PFS.csv")
write.csv(clinical_Res,file = "10_KRAS_p.G12D_GMB_PFS_2.csv")
write.csv(clinical_Res,file = "10_KRAS_p.G12V_GMB_PFS.csv")
write.csv(clinical_Res,file = "10_KRAS_Other_GMB_PFS_2.csv")

ggsurvplot(fit, data = res.cat,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("Gene_mutation_burden ≥ 2", "Gene_mutation_burden < 2"),
           palette="lancet", ylab = "Progress-free survival" )
###添加HR ,CI ,P
p3 <- ggsurvplot(fit, data = res.cat,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 legend.labs = c("Gene_mutation_burden ≥ 2", "Gene_mutation_burden < 2"),
                 palette="jto", ylab = "Progress-free survival" )
res_cox<-coxph(Surv(PFS_MONTHS, PFS_Status) ~ n, data = res.cat)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3



###Subtype 里分界值为3
####################
######################
rm(list=ls())
clinical <- read.csv("10_KRAS_p.G12D_GMB_PFS_3.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("10_KRAS_Other_GMB_PFS_3.csv",row.names = 1,check.names=FALSE)
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



######5.3 Validation of KRAS subtype and PFS
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/26_KRAS_subtype_GMB_validation"
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
g <-"KRAS"
attach(D)
r<-which(Hugo_Symbol==g)
Gene_r<-D[r,]
table(Gene_r$HGVSp_Short)
write.csv(Gene_r,file = "05_KRAS_Mut_subtype.csv")

Gene_r <- read.csv("05_KRAS_Mut_subtype.csv",row.names = 1,check.names=FALSE)
table(Gene_r$HGVSp_Short)

Gene_r_subtype <- Gene_r[Gene_r$HGVSp_Short == "p.G12C",]
Gene_r_subtype <- Gene_r[Gene_r$HGVSp_Short == "p.G12V",]
Gene_r_subtype <- Gene_r[Gene_r$HGVSp_Short == "p.G12D",]
Gene_r_subtype <- Gene_r[Gene_r$HGVSp_Short == "Other",]


Sample_KRAS_Mut <- Gene_r_subtype[!duplicated(Gene_r_subtype$Tumor_Sample_Barcode),]

k <- Sample_KRAS_Mut$Tumor_Sample_Barcode
row = match(k,clinical$Tumor_Sample_Barcode)
row = na.omit(row)
clinical_Mut = clinical[row,]

clinical_Wild = clinical[-row,]

row = match(k,clinical_Wild$Tumor_Sample_Barcode)
row = na.omit(row)
clinical_Wild = clinical_Wild[-row,]

write.csv(clinical_Mut,file = "06_KRAS_Mut_p.G12C.csv")
write.csv(clinical_Mut,file = "06_KRAS_Mut_p.G12V.csv")
write.csv(clinical_Mut,file = "06_KRAS_Mut_p.G12D.csv")
write.csv(clinical_Mut,file = "06_KRAS_Mut_Other.csv")
write.csv(clinical_Wild,file = "06_KRAS_Wild.csv")


###KRAS Subtype VS Wild#######
########################
#######################
rm(list=ls())
clinical <- read.csv("07_KRAS_Mut_Subtype.csv")
library(survminer)
library("survival")

fit<- survfit(Surv(PFS_MONTHS,Status) ~ KRAS_Subtype, data = clinical)
surv_pvalue(fit)$pval.txt
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
res_cox<-coxph(Surv(PFS_MONTHS,Status) ~ KRAS_Subtype, data = clinical)
p3$plot = p3$plot + ggplot2::annotate("text",x = 5, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 5, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 5, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3



#####KRAS Mut subtype GMB#####
##################
#################
rm(list=ls())
clinical <- read.csv("06_KRAS_Mut_p.G12C.csv",row.names = 1,check.names=FALSE)
#clinical <- read.csv("06_KRAS_Mut_p.G12D.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("06_KRAS_Mut_p.G12V.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("06_KRAS_Mut_Other.csv",row.names = 1,check.names=FALSE)

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
write.csv(clinical,file = "08_KRAS_Mut_p.G12C_GMB_2.csv")
write.csv(clinical,file = "08_KRAS_Mut_p.G12D_GMB.csv")
write.csv(clinical,file = "08_KRAS_Mut_p.G12V_GMB_5.csv")
write.csv(clinical,file = "08_KRAS_Mut_Other_GMB.csv")

ggsurvplot(fit, data = res.cat,
           pval = TRUE,
           pval.coord = c(0.5,0.05),
           risk.table = TRUE,
           xlab = "Follow up time (months)",
           legend = c(0.8,0.75),
           legend.title = "",
           legend.labs = c("Gene_mutation_burden ≥ 3", "Gene_mutation_burden < 3"),
           palette="lancet", ylab = "Progress-free survival" )
###添加HR ,CI ,P
p3 <- ggsurvplot(fit, data = res.cat,
                 pval = TRUE,
                 pval.coord = c(0.5,0.05),
                 risk.table = TRUE,
                 xlab = "Follow up time (months)",
                 legend = c(0.8,0.75),
                 legend.title = "",
                 legend.labs = c("Gene_mutation_burden ≥ 3", "Gene_mutation_burden < 3"),
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
clinical <- read.csv("08_KRAS_Mut_p.G12C_GMB_2.csv",row.names = 1,check.names=FALSE)
clinical <- read.csv("08_KRAS_Mut_p.G12V_GMB_5.csv",row.names = 1,check.names=FALSE)

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



