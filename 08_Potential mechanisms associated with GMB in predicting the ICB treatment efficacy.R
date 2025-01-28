#######6.1 Neoantigen in GMB-H
##########################Rizvi 2015#################################
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/12_mechanism"
# 设置工作目录
setwd(work_dir)

# bar_plot 2021/01/19
# 导入所需的包
library(ggplot2)
#install.packages("devtools")
library(devtools)
#install_github('cttobin/ggthemr')
library(ggthemr)
library(ggsignif)
library(tidyverse)
library(dplyr)
library(ggpubr)
#install.packages("devEMF")
library(devEMF)

# 导入并处理数据,需要两张表，一个是原始汇总表格，另外一些是每一个测量数据的均值和标准差
data1 <- read.csv("01_34_all_clinical_TMB.csv")
data <- data1[,c(13,18)]

Score_mean <- data %>% 
  dplyr::group_by(Cut) %>% 
  dplyr::summarize(
    count=n(),
    mean = mean(NEOANTIGEN_BURDEN),
    sd = sd(NEOANTIGEN_BURDEN)
  )
plot_data1 <- Score_mean
plot_data2 <- data

plot_data1$Cut <- as.factor(plot_data1$Cut)
plot_data2$Cut <- as.factor(plot_data2$Cut)
p4 <- ggplot()+
  geom_errorbar(data=plot_data1,mapping=aes(x = Cut,ymin = mean-sd, ymax = mean+sd), # 误差线添加
                width = 0.1, #误差线的宽度
                color = 'black', #颜色
                size=0.8)+ #粗细
  geom_bar(data=plot_data1,mapping=aes(x=Cut,y=mean,fill=Cut), # fill填充
           position="dodge", # 柱状图格式
           stat="identity", # 数据格式
           width = 0.7)+  # 柱状图尺寸
  scale_fill_manual(values = c("#4E4E56","#DA635D"))+ # 柱状图颜色,, "#DA635D","#B1938B"
  geom_signif(data=plot_data2,mapping=aes(x=Cut,y=NEOANTIGEN_BURDEN), # 不同组别的显著性
              comparisons = list(c("1", "2")), # 哪些组进行比较
              annotation=c("ns"), # 显著性差异做标记
              map_signif_level=T, # T为显著性，F为p value
              tip_length=c(0.04,0.04,0.05,0.05), # 修改显著性那个线的长短
              #y_position = c(4100,3000), # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 10, # 修改*标记的大小
              test = "t.test")+ # 检验的类型
  #geom_errorbar(data=plot_data1,mapping=aes(x = group,ymin = mean-sd, ymax = mean+sd), # 误差线添加
  #              width = 0.1, #误差线的宽度
  #              color = 'black', #颜色
  #              size=0.8)+ #粗细
  scale_y_continuous(limits =c(0,600) ,expand = c(0,0))+ # y轴的范围
  theme_classic(  # 主题设置，这个是无线条主题
    base_line_size = 1 # 坐标轴的粗细
  )+
  labs(title="",x="",y="")+ # 添加标题，x轴，y轴内容
  theme(plot.title = element_text(size = 20,
                                  colour = "red",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   # family = "myFont", # 类型
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 13,  # 修改y轴上字体大小，
                                   # family = "myFont", # 类型
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0) #角度
  ) 
p4
emf(file = "SOD.emf") # 打开一个矢量图画布，这种格式的图片放在word里不会失真
print(p4) # 打印图片
dev.off() #关闭画布


##########################Hellmann 2018#################################
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/12_mechanism"
# 设置工作目录
setwd(work_dir)

# bar_plot 2021/01/19
# 导入所需的包
library(ggplot2)
#install.packages("devtools")
library(devtools)
#install_github('cttobin/ggthemr')
library(ggthemr)
library(ggsignif)
library(tidyverse)
library(dplyr)
library(ggpubr)
#install.packages("devEMF")
library(devEMF)

# 导入并处理数据,需要两张表，一个是原始汇总表格，另外一些是每一个测量数据的均值和标准差
data1 <- read.csv("01_75_all_clinical_TMB.csv")
data <- data1[,c(11,16)]

Score_mean <- data %>% 
  dplyr::group_by(Cut) %>% 
  dplyr::summarize(
    count=n(),
    mean = mean(PREDICTED_NEOANTIGEN_BURDEN),
    sd = sd(PREDICTED_NEOANTIGEN_BURDEN)
  )

plot_data1 <- Score_mean
plot_data2 <- data

plot_data1$Cut <- as.factor(plot_data1$Cut)
plot_data2$Cut <- as.factor(plot_data2$Cut)
p4 <- ggplot()+
  geom_errorbar(data=plot_data1,mapping=aes(x = Cut,ymin = mean-sd, ymax = mean+sd), # 误差线添加
                width = 0.1, #误差线的宽度
                color = 'black', #颜色
                size=0.8)+ #粗细
  geom_bar(data=plot_data1,mapping=aes(x=Cut,y=mean,fill=Cut), # fill填充
           position="dodge", # 柱状图格式
           stat="identity", # 数据格式
           width = 0.7)+  # 柱状图尺寸
  scale_fill_manual(values = c("#4E4E56","#DA635D"))+ # 柱状图颜色,, "#DA635D","#B1938B"
  geom_signif(data=plot_data2,mapping=aes(x=Cut,y=PREDICTED_NEOANTIGEN_BURDEN), # 不同组别的显著性
              comparisons = list(c("1", "2")), # 哪些组进行比较
              annotation=c("ns"), # 显著性差异做标记
              map_signif_level=T, # T为显著性，F为p value
              tip_length=c(0.04,0.04,0.05,0.05), # 修改显著性那个线的长短
              #y_position = c(4100,3000), # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 10, # 修改*标记的大小
              test = "t.test")+ # 检验的类型
  #geom_errorbar(data=plot_data1,mapping=aes(x = group,ymin = mean-sd, ymax = mean+sd), # 误差线添加
  #              width = 0.1, #误差线的宽度
  #              color = 'black', #颜色
  #              size=0.8)+ #粗细
  scale_y_continuous(limits =c(0,2200) ,expand = c(0,0))+ # y轴的范围
  theme_classic(  # 主题设置，这个是无线条主题
    base_line_size = 1 # 坐标轴的粗细
  )+
  labs(title="",x="",y="")+ # 添加标题，x轴，y轴内容
  theme(plot.title = element_text(size = 20,
                                  colour = "red",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   # family = "myFont", # 类型
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 13,  # 修改y轴上字体大小，
                                   # family = "myFont", # 类型
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0) #角度
  ) 
p4
emf(file = "SOD.emf") # 打开一个矢量图画布，这种格式的图片放在word里不会失真
print(p4) # 打印图片
dev.off() #关闭画布


##########################Rizvi + Hellmann#################################
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/12_mechanism"
# 设置工作目录
setwd(work_dir)

# bar_plot 2021/01/19
# 导入所需的包
library(ggplot2)
#install.packages("devtools")
library(devtools)
#install_github('cttobin/ggthemr')
library(ggthemr)
library(ggsignif)
library(tidyverse)
library(dplyr)
library(ggpubr)
#install.packages("devEMF")
library(devEMF)

# 导入并处理数据,需要两张表，一个是原始汇总表格，另外一些是每一个测量数据的均值和标准差
data1 <- read.csv("01_34_75_neoantigen.csv")
data <- data1[,c(2,3)]

Score_mean <- data %>% 
  dplyr::group_by(Cut) %>% 
  dplyr::summarize(
    count=n(),
    mean = mean(NEOANTIGEN_BURDEN),
    sd = sd(NEOANTIGEN_BURDEN)
  )
plot_data1 <- Score_mean
plot_data2 <- data

plot_data1$Cut <- as.factor(plot_data1$Cut)
plot_data2$Cut <- as.factor(plot_data2$Cut)
p4 <- ggplot()+
  geom_errorbar(data=plot_data1,mapping=aes(x = Cut,ymin = mean-sd, ymax = mean+sd), # 误差线添加
                width = 0.1, #误差线的宽度
                color = 'black', #颜色
                size=0.8)+ #粗细
  geom_bar(data=plot_data1,mapping=aes(x=Cut,y=mean,fill=Cut), # fill填充
           position="dodge", # 柱状图格式
           stat="identity", # 数据格式
           width = 0.7)+  # 柱状图尺寸
  scale_fill_manual(values = c("#4E4E56","#DA635D"))+ # 柱状图颜色,, "#DA635D","#B1938B"
  geom_signif(data=plot_data2,mapping=aes(x=Cut,y=NEOANTIGEN_BURDEN), # 不同组别的显著性
              comparisons = list(c("1", "2")), # 哪些组进行比较
              annotation=c("ns"), # 显著性差异做标记
              map_signif_level=T, # T为显著性，F为p value
              tip_length=c(0.04,0.04,0.05,0.05), # 修改显著性那个线的长短
              #y_position = c(4100,3000), # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 10, # 修改*标记的大小
              test = "t.test")+ # 检验的类型
  #geom_errorbar(data=plot_data1,mapping=aes(x = group,ymin = mean-sd, ymax = mean+sd), # 误差线添加
  #              width = 0.1, #误差线的宽度
  #              color = 'black', #颜色
  #              size=0.8)+ #粗细
  scale_y_continuous(limits =c(0,2000) ,expand = c(0,0))+ # y轴的范围
  theme_classic(  # 主题设置，这个是无线条主题
    base_line_size = 1 # 坐标轴的粗细
  )+
  labs(title="",x="",y="")+ # 添加标题，x轴，y轴内容
  theme(plot.title = element_text(size = 20,
                                  colour = "red",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   # family = "myFont", # 类型
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 13,  # 修改y轴上字体大小，
                                   # family = "myFont", # 类型
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0) #角度
  ) 
p4
emf(file = "SOD.emf") # 打开一个矢量图画布，这种格式的图片放在word里不会失真
print(p4) # 打印图片
dev.off() #关闭画布



#######6.2 Pathway analysis by GSEA
##########################LUAD#################################
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/12_mechanism"
# 设置工作目录
setwd(work_dir)

Symbole_trans <- read.csv("02_gencode.v22.annotation.LUAD.probeMap.csv",row.names = 1,check.names=FALSE)
LUAD <-read.csv("02_TCGA_LUAD.htseq_fpkm.tsv.csv",row.names = 1,check.names=FALSE)

head(LUAD)
head(Symbole_trans)
k <- rownames(Symbole_trans)
row = match(k,rownames(LUAD))
row = na.omit(row)
LUAD_ =LUAD[row,]

LUAD_$gene_name <- Symbole_trans$gene
write.csv(LUAD_,file = "02_LUAD_gene.csv")

count_matrix_ <-read.csv("02_LUAD_gene.csv")
Cut <- read.csv("02_LUAD_Cut.csv",row.names = 1,check.names=FALSE)
count_matrix<-count_matrix_
##重复基因 平均值
exprSet_symbol1 <- aggregate(x = count_matrix[,2:ncol(count_matrix)],
                             by = list(count_matrix$gene_name),
                             FUN = mean)
head(exprSet_symbol1)
rownames(exprSet_symbol1) <- exprSet_symbol1$Group.1
exprSet_symbol <- exprSet_symbol1[,-1]

names(exprSet_symbol) <- sapply(strsplit(names(exprSet_symbol),"[.]"),function(x) paste0(x[1:3],collapse="-"))

k = rownames(Cut)
row = match(k,colnames(exprSet_symbol))
row = na.omit(row)
exprSet =exprSet_symbol[,row]

k = colnames(exprSet)
row = match(k,rownames(Cut))
row = na.omit(row)
Cut_ =Cut[row,]
table(Cut_$Cut)
write.csv(exprSet,file = "12_GSEA_LUAD.csv")

##limma
LUAD <-read.csv("12_GSEA_LUAD.csv",row.names = 1,check.names=FALSE)
library(limma)
group <- factor(c(rep("Lower_GMB",361),rep("Higher_GMB",149)))

design <- model.matrix(~0+group)
rownames(design)<-colnames(LUAD)
colnames(design)<-levels(group)

library("edgeR")
dge <- DGEList(counts=LUAD)
dge <- calcNormFactors(dge)
logCPM <-cpm(dge,log=TRUE, prior.count = 3)

v<-voom(dge,design,normalize="quantile")
fit <-lmFit(v,design)

contrasts <-paste(rev(levels(group)),collapse="-")
cont.matrix <-makeContrasts(contrasts = contrasts,levels= design)
fit2 =contrasts.fit(fit,cont.matrix)
fit2 <-eBayes(fit2)

DEG = topTable(fit2,coef =contrasts, n=Inf )
DEG =na.omit(DEG)

###GSEA
library(clusterProfiler)
library(dplyr)

options(connectionObserver = NULL)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)

data <-DEG
data$SYMBOL <-rownames(DEG)
#使用MsigDb数据库分析
KEGG <- read.gmt("c2.cp.kegg.v7.4.symbols.gmt")
geneList <-data$logFC
names(geneList) <-data$SYMBOL
geneList <-sort(geneList,decreasing = T)
gsea <- GSEA (geneList,TERM2GENE = KEGG)
write.csv(gsea,file = "12_GSEA_LUAD_Results.csv")


##########################LUSC#################################
Symbole_trans <- read.csv("02_gencode.v22.annotation.LUSC.probeMap.csv",row.names = 1,check.names=FALSE)
LUSC <-read.csv("02_TCGA_LUSC.htseq_fpkm.tsv.csv",row.names = 1,check.names=FALSE)

head(LUSC)
head(Symbole_trans)
k <- rownames(Symbole_trans)
row = match(k,rownames(LUSC))
row = na.omit(row)
LUSC_ =LUSC[row,]

LUSC_$gene_name <- Symbole_trans$gene
write.csv(LUSC_,file = "02_LUSC_gene.csv")

count_matrix_ <-read.csv("02_LUSC_gene.csv")
Cut <- read.csv("02_LUSC_Cut.csv",row.names = 1,check.names=FALSE)
count_matrix<-count_matrix_
##重复基因 平均值
exprSet_symbol1 <- aggregate(x = count_matrix[,2:ncol(count_matrix)],
                             by = list(count_matrix$gene_name),
                             FUN = mean)
head(exprSet_symbol1)
rownames(exprSet_symbol1) <- exprSet_symbol1$Group.1
exprSet_symbol <- exprSet_symbol1[,-1]

names(exprSet_symbol) <- sapply(strsplit(names(exprSet_symbol),"[.]"),function(x) paste0(x[1:3],collapse="-"))

k = rownames(Cut)
row = match(k,colnames(exprSet_symbol))
row = na.omit(row)
exprSet =exprSet_symbol[,row]

k = colnames(exprSet)
row = match(k,rownames(Cut))
row = na.omit(row)
Cut_ =Cut[row,]
table(Cut_$Cut)
write.csv(exprSet,file = "12_GSEA_LUSC.csv")

rm(list=ls())
LUSC <-read.csv("12_GSEA_LUSC.csv",row.names = 1,check.names=FALSE)
library(limma)
group <- factor(c(rep("Lower_GMB",369), rep("Higher_GMB", 115)))

design <- model.matrix(~0+group)
rownames(design)<-colnames(LUSC)
colnames(design)<-levels(group)

library("edgeR")
dge <- DGEList(counts= LUSC)
dge <- calcNormFactors(dge)
logCPM <-cpm(dge,log=TRUE, prior.count = 3)

v<-voom(dge,design,normalize="quantile")
fit <-lmFit(v,design)

contrasts <-paste(rev(levels(group)),collapse="-")
cont.matrix <-makeContrasts(contrasts = contrasts,levels= design)
fit2 =contrasts.fit(fit,cont.matrix)
fit2 <-eBayes(fit2)

DEG = topTable(fit2,coef =contrasts, n=Inf )
DEG =na.omit(DEG)

###GSEA
library(clusterProfiler)
library(dplyr)

options(connectionObserver = NULL)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
data <-DEG
data$SYMBOL <-rownames(DEG)

#使用MsigDb数据库分析
KEGG <- read.gmt("c2.cp.kegg.v7.4.symbols.gmt")
geneList <-data$logFC
names(geneList) <-data$SYMBOL
geneList <-sort(geneList,decreasing = T)
gsea <- GSEA (geneList,TERM2GENE = KEGG)
write.csv(gsea,file = "12_GSEA_LUSC_Results_.csv")


##########################NSCLC#################################
rm(list=ls())
NSCLC <-read.csv("12_GSEA_NSCLC.csv",row.names = 1,check.names=FALSE)
library(limma)
group <- factor(c(rep("Lower_GMB",361),rep("Higher_GMB",149),
                  rep("Lower_GMB",369), rep("Higher_GMB", 115)))

design <- model.matrix(~0+group)
rownames(design)<-colnames(NSCLC)
colnames(design)<-levels(group)

library("edgeR")
dge <- DGEList(counts= NSCLC)
dge <- calcNormFactors(dge)
logCPM <-cpm(dge,log=TRUE, prior.count = 3)

v<-voom(dge,design,normalize="quantile")
fit <-lmFit(v,design)

contrasts <-paste(rev(levels(group)),collapse="-")
cont.matrix <-makeContrasts(contrasts = contrasts,levels= design)
fit2 =contrasts.fit(fit,cont.matrix)
fit2 <-eBayes(fit2)

DEG = topTable(fit2,coef =contrasts, n=Inf )
DEG =na.omit(DEG)

###GSEA
library(clusterProfiler)
library(dplyr)

options(connectionObserver = NULL)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)

data <-DEG
data$SYMBOL <-rownames(DEG)

#使用MsigDb数据库分析
KEGG <- read.gmt("c2.cp.kegg.v7.4.symbols.gmt")

geneList <-data$logFC
names(geneList) <-data$SYMBOL
geneList <-sort(geneList,decreasing = T)
gsea <- GSEA (geneList,TERM2GENE = KEGG)
write.csv(gsea,file = "12_GSEA_NSCLC_Results_.csv")

####可视化
dataxx <-read.csv("13_GSEA_pathway.csv",row.names = 1,check.names=FALSE)
datax <- as.matrix(dataxx)
barplot(datax,beside = TRUE,horiz =T)
barplot(datax,beside = TRUE,angle = 15+10*1:5, density = 10,col = rainbow(2),horiz =T)
barplot(datax,beside = TRUE,angle = 15+10*1:5, density = 10,,horiz =T)



#######6.3 tumor immune microenvironment in GMB-H vs GMB-L
##########################LUAD#################################
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/13_CD8_cell"
# 设置工作目录
setwd(work_dir)
Immunecell <- read.csv("LUAD_CIBERSORTx_Job1_Results.csv",row.names = 1,check.names=FALSE)
LUAD <-read.csv("02_LUAD_Cut.csv",row.names = 1,check.names=FALSE)

k = rownames(LUAD)
row = match(k,rownames(Immunecell))
row = na.omit(row)
Immunecell_ =Immunecell[row,]

k = rownames(Immunecell_)
row = match(k,rownames(LUAD))
row = na.omit(row)
LUAD_ =LUAD[row,]
table(LUAD_$Cut)

Immunecell_$Cut <- LUAD_$Cut
write.csv(Immunecell_,file = "13_Tcell_LUAD.csv")

##误差条图+箱线图+散点图
library(ggplot2)
library(RColorBrewer)

#使用ggplot2包生成箱线图
##CD8+ T cell
compaired <- list(c("High", "Low")) 
P1 <- ggplot(Immunecell_,aes(x=Cut,y=Immunecell_[,4]))+ #”fill=“设置填充颜色
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ #由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
  geom_boxplot(size=0.5,outlier.fill="white",outlier.color="white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
  geom_jitter(aes(fill=Cut),width =0.2,shape = 21,size=2.5)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
  scale_fill_manual(values = c("#E69F00", "#0072B2"))+  #设置填充的颜色"#F0E442"
  scale_color_manual(values=c("black","black"))+ #设置散点图的圆圈的颜色为黑色"black"
  ggtitle(" ")+ #设置总的标题
  theme_bw()+ #背景变为白色
  theme(legend.position="none", #不需要图例
        axis.text.x=element_text(colour="black",family="Times",size=14), #设置x轴刻度标签的字体属性
        axis.text.y=element_text(family="Times",size=14,face="plain"), #设置x轴刻度标签的字体属性
        axis.title.y=element_text(family="Times",size = 14,face="plain"), #设置y轴的标题的字体属性
        axis.title.x=element_text(family="Times",size = 14,face="plain"), #设置x轴的标题的字体属性
        plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5), #设置总标题的字体属性
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank())+
  ylab("CD8+ T Cell")+xlab(" ") #设置x轴和y轴的标题
P1+stat_compare_means(comparisons = compaired, method = "wilcox.test",label = "p.value",label.y = 0.5) 

## Treg cell
compaired <- list(c("High", "Low")) 
P1 <- ggplot(Immunecell_,aes(x=Cut,y=Immunecell_[,9]))+ #”fill=“设置填充颜色
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ #由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
  geom_boxplot(size=0.5,outlier.fill="white",outlier.color="white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
  geom_jitter(aes(fill=Cut),width =0.2,shape = 21,size=2.5)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
  scale_fill_manual(values = c("#E69F00", "#0072B2"))+  #设置填充的颜色"#F0E442"
  scale_color_manual(values=c("black","black"))+ #设置散点图的圆圈的颜色为黑色"black"
  ggtitle(" ")+ #设置总的标题
  theme_bw()+ #背景变为白色
  theme(legend.position="none", #不需要图例
        axis.text.x=element_text(colour="black",family="Times",size=14), #设置x轴刻度标签的字体属性
        axis.text.y=element_text(family="Times",size=14,face="plain"), #设置x轴刻度标签的字体属性
        axis.title.y=element_text(family="Times",size = 14,face="plain"), #设置y轴的标题的字体属性
        axis.title.x=element_text(family="Times",size = 14,face="plain"), #设置x轴的标题的字体属性
        plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5), #设置总标题的字体属性
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank())+
  ylab("Treg Cell")+xlab(" ") #设置x轴和y轴的标题
P1+stat_compare_means(comparisons = compaired, method = "wilcox.test",label = "p.value",label.y = 0.2) 


##########################LUSC#################################
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/13_CD8_cell"
# 设置工作目录
setwd(work_dir)

##LUSC 
Immunecell <- read.csv("LUSC_CIBERSORTx_Job2_Results.csv",row.names = 1,check.names=FALSE)
LUSC <-read.csv("02_LUSC_Cut.csv",row.names = 1,check.names=FALSE)

k = rownames(LUSC)
row = match(k,rownames(Immunecell))
row = na.omit(row)
Immunecell_ =Immunecell[row,]

k = rownames(Immunecell_)
row = match(k,rownames(LUSC))
row = na.omit(row)
LUSC_ =LUSC[row,]
table(LUSC_$Cut)

Immunecell_$Cut <- LUSC_$Cut
write.csv(Immunecell_,file = "13_Tcell_LUSC.csv")

##误差条图+箱线图+散点图
#rm(list = ls()) #清除工作区
library(ggplot2)
library(RColorBrewer)

#使用ggplot2包生成箱线图
##CD8+ T cell
compaired <- list(c("High", "Low")) 
P1 <- ggplot(Immunecell_,aes(x=Cut,y=Immunecell_[,4]))+ #”fill=“设置填充颜色
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ #由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
  geom_boxplot(size=0.5,outlier.fill="white",outlier.color="white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
  geom_jitter(aes(fill=Cut),width =0.2,shape = 21,size=2.5)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
  scale_fill_manual(values = c("#E69F00", "#0072B2"))+  #设置填充的颜色"#F0E442"
  scale_color_manual(values=c("black","black"))+ #设置散点图的圆圈的颜色为黑色"black"
  ggtitle(" ")+ #设置总的标题
  theme_bw()+ #背景变为白色
  theme(legend.position="none", #不需要图例
        axis.text.x=element_text(colour="black",family="Times",size=14), #设置x轴刻度标签的字体属性
        axis.text.y=element_text(family="Times",size=14,face="plain"), #设置x轴刻度标签的字体属性
        axis.title.y=element_text(family="Times",size = 14,face="plain"), #设置y轴的标题的字体属性
        axis.title.x=element_text(family="Times",size = 14,face="plain"), #设置x轴的标题的字体属性
        plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5), #设置总标题的字体属性
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank())+
  ylab("CD8+ T Cell")+xlab(" ") #设置x轴和y轴的标题
P1+stat_compare_means(comparisons = compaired, method = "wilcox.test",label = "p.value",label.y = 0.5) 

## Treg cell
compaired <- list(c("High", "Low")) 
P1 <- ggplot(Immunecell_,aes(x=Cut,y=Immunecell_[,9]))+ #”fill=“设置填充颜色
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ #由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
  geom_boxplot(size=0.5,outlier.fill="white",outlier.color="white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
  geom_jitter(aes(fill=Cut),width =0.2,shape = 21,size=2.5)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
  scale_fill_manual(values = c("#E69F00", "#0072B2"))+  #设置填充的颜色"#F0E442"
  scale_color_manual(values=c("black","black"))+ #设置散点图的圆圈的颜色为黑色"black"
  ggtitle(" ")+ #设置总的标题
  theme_bw()+ #背景变为白色
  theme(legend.position="none", #不需要图例
        axis.text.x=element_text(colour="black",family="Times",size=14), #设置x轴刻度标签的字体属性
        axis.text.y=element_text(family="Times",size=14,face="plain"), #设置x轴刻度标签的字体属性
        axis.title.y=element_text(family="Times",size = 14,face="plain"), #设置y轴的标题的字体属性
        axis.title.x=element_text(family="Times",size = 14,face="plain"), #设置x轴的标题的字体属性
        plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5), #设置总标题的字体属性
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank())+
  ylab("Treg Cell")+xlab(" ") #设置x轴和y轴的标题
P1+stat_compare_means(comparisons = compaired, method = "wilcox.test",label = "p.value",label.y = 0.18) 


##########################NSCLC#################################
rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/特定突变基因数目VSTMB/13_CD8_cell"
# 设置工作目录
setwd(work_dir)

##LUSC 
Immunecell_ <- read.csv("13_Tcell_NSCLC.csv",row.names = 1,check.names=FALSE)
#使用ggplot2包生成箱线图
##CD8+ T cell
compaired <- list(c("High", "Low")) 
P1 <- ggplot(Immunecell_,aes(x=Cut,y=Immunecell_[,4]))+ #”fill=“设置填充颜色
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ #由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
  geom_boxplot(size=0.5,outlier.fill="white",outlier.color="white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
  geom_jitter(aes(fill=Cut),width =0.2,shape = 21,size=2.5)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
  scale_fill_manual(values = c("#E69F00", "#0072B2"))+  #设置填充的颜色"#F0E442"
  scale_color_manual(values=c("black","black"))+ #设置散点图的圆圈的颜色为黑色"black"
  ggtitle(" ")+ #设置总的标题
  theme_bw()+ #背景变为白色
  theme(legend.position="none", #不需要图例
        axis.text.x=element_text(colour="black",family="Times",size=14), #设置x轴刻度标签的字体属性
        axis.text.y=element_text(family="Times",size=14,face="plain"), #设置x轴刻度标签的字体属性
        axis.title.y=element_text(family="Times",size = 14,face="plain"), #设置y轴的标题的字体属性
        axis.title.x=element_text(family="Times",size = 14,face="plain"), #设置x轴的标题的字体属性
        plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5), #设置总标题的字体属性
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank())+
  ylab("CD8+ T Cell")+xlab(" ") #设置x轴和y轴的标题
P1+stat_compare_means(comparisons = compaired, method = "wilcox.test",label = "p.value",label.y = 0.5) 

## Treg cell
compaired <- list(c("High", "Low")) 
P1 <- ggplot(Immunecell_,aes(x=Cut,y=Immunecell_[,9]))+ #”fill=“设置填充颜色
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ #由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
  geom_boxplot(size=0.5,outlier.fill="white",outlier.color="white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
  geom_jitter(aes(fill=Cut),width =0.2,shape = 21,size=2.5)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
  scale_fill_manual(values = c("#E69F00", "#0072B2"))+  #设置填充的颜色"#F0E442"
  scale_color_manual(values=c("black","black"))+ #设置散点图的圆圈的颜色为黑色"black"
  ggtitle(" ")+ #设置总的标题
  theme_bw()+ #背景变为白色
  theme(legend.position="none", #不需要图例
        axis.text.x=element_text(colour="black",family="Times",size=14), #设置x轴刻度标签的字体属性
        axis.text.y=element_text(family="Times",size=14,face="plain"), #设置x轴刻度标签的字体属性
        axis.title.y=element_text(family="Times",size = 14,face="plain"), #设置y轴的标题的字体属性
        axis.title.x=element_text(family="Times",size = 14,face="plain"), #设置x轴的标题的字体属性
        plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5), #设置总标题的字体属性
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank())+
  ylab("Treg Cell")+xlab(" ") #设置x轴和y轴的标题
P1+stat_compare_means(comparisons = compaired, method = "wilcox.test",label = "p.value",label.y = 0.18) 







