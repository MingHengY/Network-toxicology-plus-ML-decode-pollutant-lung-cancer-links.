#提取风险评分并保存
load("TCGA_LUAD.RData")
data <- res_luad[["riskscore"]][["Lasso + RSF"]][["TCGA_LUAD"]]

save(data,file = "TCGA_LUAD_riskscore.RData")
load("TCGA_LUAD_riskscore.RData")

#根据风险评分分组
data$Group <- ifelse(data$RS > median(data$RS),"High Risk","Low Risk")

highrisk <- data[data$Group == "High Risk",]
lowrisk <- data[data$Group == "Low Risk",]


#合并表达数据
library(tidyverse)
library(limma)
data=read.table("TCGA-LUAD.star_tpm.tsv",header = T,sep = '\t',check.names = F,row.names = 1)
platform_file=read.delim("gencode.v36.annotation.gtf.gene.probemap", header = TRUE, sep = "\t")
plat<-platform_file[,c(1,2)]
data$id<-rownames(data)
luad<-merge(plat,data,by="id")
#对重复基因合并取平均值
luad_gset = as.data.frame(avereps(luad[,-1],ID = luad$gene) )
rownames(luad_gset) <- luad_gset$gene
luad_gset <- luad_gset[,-1]
expr_data <- luad_gset



data_highrisk<-expr_data[,highrisk$ID]
data_lowrisk<-expr_data[,lowrisk$ID]
data_highrisk<-t(data_highrisk)
data_lowrisk<-t(data_lowrisk)

# 确保数据是字符型
expr_data <- as.data.frame(lapply(expr_data, as.character))

# 将字符型数据转换为数值型
expr_data <- as.data.frame(lapply(expr_data, as.numeric))
rownames(expr_data) <- rownames(luad_gset)
colnames(expr_data) <- colnames(luad_gset)
# 检查转换后的数据类型
sapply(expr_data, class)

keyset=expr_data[rowMeans(expr_data)>1,]

high=keyset[,highrisk$ID]
low=keyset[,lowrisk$ID]
high<-log2(high+1)
low<-log2(low+1)
data2<-cbind(high,low)
write.table(data2,"TPM_with_log2.txt",sep = "\t",quote = F, row.names = T,col.names = T)

#导出
write.table(data_highrisk,"Highrisk_all.txt",sep = "\t",quote = F, row.names = T,col.names = T)
write.table(data_lowrisk,"Lowrisk_all.txt",sep = "\t",quote = F, row.names = T,col.names = T)
write.table(highrisk,"Highriskgroup.txt",sep = "\t",quote = F, row.names = T,col.names = T)
write.table(lowrisk,"Lowriskgroup.txt",sep = "\t",quote = F, row.names = T,col.names = T)
write.table(high,"High_keySet.txt",sep = "\t",quote = F, row.names = T,col.names = T)
write.table(low,"Low_keySet.txt",sep = "\t",quote = F, row.names = T,col.names = T)


##免疫浸润----
#CIBERSORT
#生成表达矩阵文件
highRisk<-read.table("Highrisk_all.txt",
                     header = T, sep = '\t', check.names = F, row.names = 1)
lowRisk<-read.table("Lowrisk_all.txt",
                    header = T, sep = '\t', check.names = F, row.names = 1)


library(tidyverse)
library(limma)
data=read.table("TCGA-LUAD.star_tpm.tsv",header = T,sep = '\t',check.names = F,row.names = 1)
platform_file=read.delim("gencode.v36.annotation.gtf.gene.probemap", header = TRUE, sep = "\t")
plat<-platform_file[,c(1,2)]
data$id<-rownames(data)
luad<-merge(plat,data,by="id")
#对重复基因合并取平均值
luad_gset = as.data.frame(avereps(luad[,-1],ID = luad$gene) )
rownames(luad_gset) <- luad_gset$gene
luad_gset <- luad_gset[,-1]
expr_data <- luad_gset

# 将字符型数据转换为数值型
expr_data <- as.data.frame(lapply(expr_data, as.numeric))
rownames(expr_data) <- rownames(luad_gset)
colnames(expr_data) <- colnames(luad_gset)

expr_data <- 2^expr_data

highRiskGroup<-expr_data[,colnames(highRisk)]
lowRiskGroup<-expr_data[,colnames(lowRisk)]

write.table(highRiskGroup,"high_forCIBERSORT.txt",sep = "\t",quote = F, row.names = T,col.names = T)
write.table(lowRiskGroup,"low_forCIBERSORT.txt",sep = "\t",quote = F, row.names = T,col.names = T)



source("D:/Documents/Documents/CIBERSORT/sourcecibersort.R")
high_cibersort<-CIBERSORT(sig_matrix = "LM22.txt",mixture_file = "high_forCIBERSORT.txt",perm = 1000,QN=F)
low_cibersort<-CIBERSORT(sig_matrix = "LM22.txt",mixture_file = "low_forCIBERSORT.txt",perm = 1000,QN=F)

write.table(high_cibersort,"high_cibersort.txt",sep = "\t",quote = F, row.names = T,col.names = T)
write.table(low_cibersort,"low_cibersort.txt",sep = "\t",quote = F, row.names = T,col.names = T)

library(reshape2)
library(dplyr)
#高风险组数据
HC_data<-as.data.frame(high_cibersort[,1:22])
HC_data$sample<-rownames(HC_data)
HC_data$group<-"High Risk"
HC_melt<-melt(HC_data)
colnames(HC_melt)=c("Sample","Group","Celltype","Composition")
head(HC_melt)
#低风险组数据
LC_data<-as.data.frame(low_cibersort[,1:22])
LC_data$sample<-rownames(LC_data)
LC_data$group<-"Low Risk"
LC_melt<-melt(LC_data)
colnames(LC_melt)=c("Sample","Group","Celltype","Composition")
head(LC_melt)

#合并
IMM<-rbind(HC_melt,LC_melt)

#作图
library(ggplot2)
library(ggpubr)
if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"),
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12))
}

#箱线图
box_IMM <- ggplot(IMM, aes(x = Celltype, y = Composition))+
  labs(y="Cell composition",x= NULL,title = "High & Low Risk Group Cell composition From CIBERSORT")+
  geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+
  scale_fill_manual(values = c("#EB7369", "#1CB4B8"))+
  theme_classic() + mytheme +
  stat_compare_means(aes(group =  Group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)
box_IMM;ggsave("Immume_CIBERSORT_LUAD.pdf",box_IMM,height=15,width=25,unit="cm")


##QuanTiseq----

library(reshape2)
library(ggpubr)
library(ggplot2)
library(tidyverse)
library(ggsci)
library(quantiseqr)
library(dplyr)
library(tidyr)
library(tibble)
library(GEOquery)
library(SummarizedExperiment)

highRisk<-read.table("Highrisk_all.txt",
                     header = T, sep = '\t', check.names = F, row.names = 1)
lowRisk<-read.table("Lowrisk_all.txt",
                    header = T, sep = '\t', check.names = F, row.names = 1)

H_ti_racle <- quantiseqr::run_quantiseq(
  expression_data = highRisk,
  signature_matrix = "TIL10",
  is_arraydata = FALSE,
  is_tumordata = TRUE,
  scale_mRNA = TRUE
)

L_ti_racle <- quantiseqr::run_quantiseq(
  expression_data = lowRisk,
  signature_matrix = "TIL10",
  is_arraydata = FALSE,
  is_tumordata = TRUE,
  scale_mRNA = TRUE
)

quantiplot(H_ti_racle)
quantiplot(L_ti_racle)

library(reshape2)
library(dplyr)
#高风险组数据
HC_data<-as.data.frame(H_ti_racle[,2:12])
HC_data$sample<-rownames(HC_data)
HC_data$group<-"High Risk"
HC_melt<-melt(HC_data)
colnames(HC_melt)=c("Sample","Group","Celltype","Composition")
head(HC_melt)
#低风险组数据
LC_data<-as.data.frame(L_ti_racle[,2:12])
LC_data$sample<-rownames(LC_data)
LC_data$group<-"Low Risk"
LC_melt<-melt(LC_data)
colnames(LC_melt)=c("Sample","Group","Celltype","Composition")
head(LC_melt)

#合并
IMM<-rbind(HC_melt,LC_melt)

#作图
library(ggplot2)
library(ggpubr)
if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"),
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12))
}

#箱线图
box_IMM <- ggplot(IMM, aes(x = Celltype, y = Composition))+
  labs(y="Cell composition",x= NULL,title = "High & Low Risk Group Cell composition From CIBERSORT")+
  geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+
  scale_fill_manual(values = c("#EB7369", "#3399FF"))+
  theme_classic() + mytheme +
  stat_compare_means(aes(group =  Group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)
box_IMM;ggsave("Immume_QuanTiseq_LUAD.pdf",box_IMM,height=15,width=25,unit="cm")


##TIMER----
library(IOBR)

highRisk<-read.table("Highrisk_all.txt",
                     header = T, sep = '\t', check.names = F, row.names = 1)
lowRisk<-read.table("Lowrisk_all.txt",
                    header = T, sep = '\t', check.names = F, row.names = 1)

T_high<-deconvo_tme(eset = highRisk, method = "timer", group_list = rep("luad",dim(highRisk)[2]))
T_low<-deconvo_tme(eset = lowRisk, method = "timer", group_list = rep("luad",dim(lowRisk)[2]))
library(reshape2)
library(dplyr)
#高风险组数据
T_high<-as.data.frame(T_high)
colnames(T_high)<-gsub("_TIMER", "", colnames(T_high))
T_high$Group<-"High Risk"
T_high1<-melt(T_high)
colnames(T_high1)<-c("ID","Group","Celltype","Composition")

#低风险组数据
T_low<-as.data.frame(T_low)
colnames(T_low)<-gsub("_TIMER", "", colnames(T_low))
T_low$Group<-"Low Risk"
T_low1<-melt(T_low)
colnames(T_low1)<-c("ID","Group","Celltype","Composition")

Timer<-rbind(T_high1,T_low1)
head(Timer)

if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"),
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12))
}

box_TIMER <- ggplot(Timer, aes(x = Celltype, y = Composition))+
  labs(y="Cell composition",x= NULL,title = "High & Low Risk Group Cell composition From Timer")+
  geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+
  scale_fill_manual(values = c("#CC6633", "#3399FF"))+
  theme_classic() + mytheme +
  stat_compare_means(aes(group =  Group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)
ggsave("Immume_Timer.pdf",box_TIMER,height=15,width=25,unit="cm");box_TIMER


##ssGSEA----
library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(GSVA)
library(GSEABase)
library(pheatmap)
library(limma)
library(stats)
library(BiocParallel)

#表达数据
highRisk<-read.table("Highrisk_all.txt",
                     header = T, sep = '\t', check.names = F, row.names = 1)
lowRisk<-read.table("Lowrisk_all.txt",
                    header = T, sep = '\t', check.names = F, row.names = 1)

sample<-cbind(highRisk,lowRisk)
sample<-sample[rowMeans(sample)>1,]

hg<-data.frame(Group=rep("High Risk",288))
rownames(hg)<-colnames(highRisk)

lg<-data.frame(Group=rep("Low Risk",288))
rownames(lg)<-colnames(lowRisk)

Group<-rbind(hg,lg)
Group$sample<-rownames(Group)
#基因集
load("TISIDB肿瘤浸润淋巴细胞基因集.rdata")


# 创建参数对象
gsvaPar <- ssgseaParam(exprData = as.matrix(sample), geneSets = tisidb_cell, normalize = T)

# 调用 gsva 函数
gsva_mat <- gsva(param = gsvaPar,
                 verbose = TRUE,
                 BPPARAM = BiocParallel::SnowParam(parallel::detectCores()))

library(dplyr)
library(reshape2)
library(rstatix)

gsva_mat <- t(gsva_mat) %>% as.data.frame()
colnames(gsva_mat)

gsva_mat$sample <- rownames(gsva_mat)
longer_data <- gsva_mat %>%
  pivot_longer(
    cols = -c("sample"),
    names_to = "ImmuneCell",
    values_to = "Score"
  )
ssgsea_data <- merge(longer_data, Group, by = "sample")
write.csv(ssgsea_data, file = 'ssgsea_results.csv')

wilcox_res <- ssgsea_data %>%
  group_by(ImmuneCell) %>%
  wilcox_test(Score ~ Group) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p")

write.csv(wilcox_res, file = 'wilcox_res_sample.csv')

library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)

p_values <- ssgsea_data %>%
  group_by(ImmuneCell) %>%
  wilcox_test(Score ~ Group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  select(ImmuneCell, p.adj, p.adj.signif)

# 计算各组的均值和中位数
summary_data <- ssgsea_data %>%
  group_by(ImmuneCell, Group) %>%
  summarise(
    mean_score = mean(Score),
    sd_score = sd(Score),
    .groups = "drop"
  )

# 合并统计结果
plot_data <- left_join(summary_data, p_values, by = "ImmuneCell")

# 绘制柱状图
p<-ggplot(plot_data, aes(x = ImmuneCell, y = mean_score, fill = Group)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(0.8),
    width = 0.7
  ) +
  geom_errorbar(
    aes(ymin = mean_score - sd_score, ymax = mean_score + sd_score),
    position = position_dodge(0.8),
    width = 0.2
  ) +
  geom_text(
    aes(label = p.adj.signif, y = mean_score + sd_score + 0.05),
    position = position_dodge(0.8),
    vjust = 0
  ) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  labs(
    x = "Immune Cell Type",
    y = "Enrichment Score",
    fill = "Risk Group"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.margin = unit(c(1,1,1,2), "cm")
  )
ggsave("ssgsea.pdf",plot = p,height = 8,width = 12,units = "in");p


##TIDE----
#设置分组
highRisk<-read.table("Highrisk_all.txt",
                     header = T, sep = '\t', check.names = F, row.names = 1)
lowRisk<-read.table("Lowrisk_all.txt",
                    header = T, sep = '\t', check.names = F, row.names = 1)

pat <- cbind(highRisk,lowRisk)
write.table(pat,"tide_array.txt",sep = "\t",quote = F, row.names = T,col.names = T)

hg<-data.frame(Group=rep("High Risk",288))
rownames(hg)<-colnames(highRisk)

lg<-data.frame(Group=rep("Low Risk",288))
rownames(lg)<-colnames(lowRisk)

Group<-rbind(hg,lg)
Group$Patient<-rownames(Group)

#分析TIDE计算结果
res_tide<-read.csv("TIDE0726.csv",header = T)
tide<-merge(res_tide,Group,by="Patient")

#绘图（云雨图）
library(ggplot2)
library(ggrain)
library(gghalves)

#TIDE评分
tidy<-tide[,c("TIDE","Group")]
pt<-ggplot(tidy,aes(Group,TIDE,fill = Group))+
  geom_rain(alpha=1,
            cov = "TIDE",
            boxplot.args.pos = list(width=0.05,position=position_nudge(x=0.13)),
            violin.args.pos = list(side="r",width=0.7,position=position_nudge(x=0.2)))+
  theme_classic()+
  scale_color_brewer(palette = "Set1")+
  guides(fill="none",color="none")+
  scale_color_viridis_c(option = "A",direction = -1)+
  ggsignif::geom_signif(comparisons = list(c("High Risk","Low Risk")),
                        test = "wilcox.test",
                        map_signif_level = c("***"=0.001,"**"=0.01,"*"=0.05),
                        size = 1,
                        textsize = 5)
ggsave(filename = "TIDEScore.pdf",plot = pt,width = 10,height = 8,units = "in")


#Dysfunction
tidy1<-tide[,c("Dysfunction","Group")]
pt1<-ggplot(tidy1,aes(Group,Dysfunction,fill = Group))+
  geom_rain(alpha=1,
            cov = "Dysfunction",
            boxplot.args.pos = list(width=0.05,position=position_nudge(x=0.13)),
            violin.args.pos = list(side="r",width=0.7,position=position_nudge(x=0.2)))+
  theme_classic()+
  scale_color_brewer(palette = "Set1")+
  guides(fill="none",color="none")+
  scale_color_viridis_c(option = "A",direction = -1)+
  ggsignif::geom_signif(comparisons = list(c("High Risk","Low Risk")),
                        test = "wilcox.test",
                        map_signif_level = c("***"=0.001,"**"=0.01,"*"=0.05),
                        size = 1,
                        textsize = 5)
ggsave(filename = "TIDE_Dysfunction.pdf",plot = pt1,width = 10,height = 8,units = "in")

#Exclusion
tidy2<-tide[,c("Exclusion","Group")]
pt2<-ggplot(tidy2,aes(Group,Exclusion,fill = Group))+
  geom_rain(alpha=1,
            cov = "Exclusion",
            boxplot.args.pos = list(width=0.05,position=position_nudge(x=0.13)),
            violin.args.pos = list(side="r",width=0.7,position=position_nudge(x=0.2)))+
  theme_classic()+
  scale_color_brewer(palette = "Set1")+
  guides(fill="none",color="none")+
  scale_color_viridis_c(option = "A",direction = -1)+
  ggsignif::geom_signif(comparisons = list(c("High Risk","Low Risk")),
                        test = "wilcox.test",
                        map_signif_level = c("***"=0.001,"**"=0.01,"*"=0.05),
                        size = 1,
                        textsize = 5)
ggsave(filename = "TIDE_Exclusion.pdf",plot = pt2,width = 10,height = 8,units = "in")

#TAM
tidy3<-tide[,c("TAM.M2","Group")]
pt3<-ggplot(tidy3,aes(Group,TAM.M2,fill = Group))+
  geom_rain(alpha=1,
            cov = "TAM.M2",
            boxplot.args.pos = list(width=0.05,position=position_nudge(x=0.13)),
            violin.args.pos = list(side="r",width=0.7,position=position_nudge(x=0.2)))+
  theme_classic()+
  scale_color_brewer(palette = "Set1")+
  guides(fill="none",color="none")+
  scale_color_viridis_c(option = "A",direction = -1)+
  ggsignif::geom_signif(comparisons = list(c("High Risk","Low Risk")),
                        test = "wilcox.test",
                        map_signif_level = c("***"=0.001,"**"=0.01,"*"=0.05),
                        size = 1,
                        textsize = 5)
ggsave(filename = "TIDE_M2Macrophage.pdf",plot = pt3,width = 10,height = 8,units = "in")


#MDSC
tidy4<-tide[,c("MDSC","Group")]
pt4<-ggplot(tidy4,aes(Group,MDSC,fill = Group))+
  geom_rain(alpha=1,
            cov = "MDSC",
            boxplot.args.pos = list(width=0.05,position=position_nudge(x=0.13)),
            violin.args.pos = list(side="r",width=0.7,position=position_nudge(x=0.2)))+
  theme_classic()+
  scale_color_brewer(palette = "Set1")+
  guides(fill="none",color="none")+
  scale_color_viridis_c(option = "A",direction = -1)+
  ggsignif::geom_signif(comparisons = list(c("High Risk","Low Risk")),
                        test = "wilcox.test",
                        map_signif_level = c("***"=0.001,"**"=0.01,"*"=0.05),
                        size = 1,
                        textsize = 5)
ggsave(filename = "TIDE_MDSC.pdf",plot = pt4,width = 10,height = 8,units = "in")

#CD8
tidy5<-tide[,c("CD8","Group")]
pt5<-ggplot(tidy5,aes(Group,CD8,fill = Group))+
  geom_rain(alpha=1,
            cov = "CD8",
            boxplot.args.pos = list(width=0.05,position=position_nudge(x=0.13)),
            violin.args.pos = list(side="r",width=0.7,position=position_nudge(x=0.2)))+
  theme_classic()+
  scale_color_brewer(palette = "Set1")+
  guides(fill="none",color="none")+
  scale_color_viridis_c(option = "A",direction = -1)+
  ggsignif::geom_signif(comparisons = list(c("High Risk","Low Risk")),
                        test = "wilcox.test",
                        map_signif_level = c("***"=0.001,"**"=0.01,"*"=0.05),
                        size = 1,
                        textsize = 5)
ggsave(filename = "TIDE_CD8.pdf",plot = pt5,width = 10,height = 8,units = "in")


##药物敏感性分析----
# install.packages("oncoPredict")
library(oncoPredict)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)

highRisk<-read.table("high_keySet.txt",
                     header = T, sep = '\t', check.names = F, row.names = 1)
lowRisk<-read.table("low_keySet.txt",
                    header = T, sep = '\t', check.names = F, row.names = 1)

highexprSet<-as.matrix(t(highRisk))
lowexprSet<-as.matrix(t(lowRisk))
exprSet<-rbind(highexprSet,lowexprSet)
exprSet<-t(exprSet)


GDSC2_Expr <- readRDS(file = "D:/Documents/Documents/DataFiles/DataFiles/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
GDSC2_Res <- readRDS(file = "D:/Documents/Documents/DataFiles/DataFiles/Training Data/GDSC2_Res.rds")
GDSC2_Res<-exp(GDSC2_Res)

#计算
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = exprSet,#需要matrix
              batchCorrect = 'eb',
              #IC50是对数转换的，所以表达矩阵也用对数转换过的
              powerTransformPhenotype = F,
              minNumSamples = 20,
              printOutput = T,
              removeLowVaryingGenes = 0.2,
              removeLowVaringGenesFrom = "homogenizeData"
)


##GDSC2----
#读取数据并处理
res <- read.csv("DrugPredictions.csv")
rownames(res)<-colnames(exprSet)
res1<-res[,c(2:ncol(res))]
colnames(res1)<-gsub("(.*)\\_(\\d+)","\\1",colnames(res1))
res1[,dim(res1)[2]+1]=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3\\",
                           rownames(res1))
res2=res1[!duplicated(res1$V199),]
rownames(res2)=res1[,dim(res1)[2]]
res3<-as.data.frame(res2[,-199])
rownames(res3)<-rownames(res2)
#分组
sample_group<-rep(c("High Risk","Low Risk"),c(288,288))
group_df <- data.frame(group = sample_group)
rownames(group_df)<-colnames(exprSet)
rt=cbind(group_df,res3)
#比较组，用于差异分析
rt$group=factor(rt$group,levels = c("High Risk","Low Risk"))
type=levels(factor(rt[,"group"]))
comp=combn(type,2)
comparisons=list()
for(i in 1:ncol(comp)){comparisons[[i]]<-comp[,i]}

#寻找具有显著差异的药物
sigGene=c()
for(i in colnames(rt)[2:(ncol(rt))]){
  if(sd(rt[,i])<0.05){next}
  wilcoxTest=wilcox.test(rt[,i]~rt[,"group"])
  pvalue=wilcoxTest$p.value
  if(wilcoxTest$p.value<0.01){
    sigGene=c(sigGene,i)
  }
}
sigGene=c(sigGene,"group")
rt1=rt[,sigGene]

#数据转化为ggplot2以用于绘图
library(tidyverse)
library(reshape2)
rt2=melt(rt1,id.vars = c("group"))
colnames(rt2)=c("risk","Gene","Expression")

#再次分组
group=levels(factor(rt2$risk))
rt2$risk=factor(rt2$risk,levels = c("High Risk","Low Risk"))
comp=combn(group,2)
comparisons=list()
for(i in 1:ncol(comp)){comparisons[[i]]<-comp[,i]}
rt2$risk<-factor(rt2$risk,levels = rev(levels(rt2$risk)))
write.csv(rt2,"Drugrt2.csv",row.names = F)
#箱线图
boxplot=ggboxplot(rt2,x="Gene",y="Expression",fill = "risk",
                  xlab = "",
                  ylab = "Drug Sensitivity",
                  legend.title = "Risk",
                  width = 0.8,
                  palette = c("DodgerBlue1","Firebrick2"))+
  rotate_x_text(50)+
  stat_compare_means(aes(group=risk),method = "wilcox.test",
                     symnum.args =list(cutpoints=c(0,0.001,0.01,0.05,1),
                                       symbols=c("***","**","*","ns")),label = "p.signif" )+
  theme(axis.text = element_text(face = "bold.italic",colour = "#441718",size = 16),
        axis.title = element_text(face = "bold.italic",colour = "#441718",size = 16),
        axis.line = element_blank(),
        plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16),
        legend.text = element_text(face = "bold.italic"),
        panel.border = element_rect(fill = NA,color = "#35A79D",size = 1.5,linetype = "solid"),
        panel.background = element_rect(fill = "#F1F6FC"),
        panel.grid.major = element_line(color = "#CFD3D6",size = 0.5,linetype = "dotdash"),
        legend.title = element_text(face = "bold.italic",size = 13)
  )

pdf(file = "drugSensitivity_GDSC2.pdf",width = 60,height = 10)
print(boxplot)
dev.off()

#再次挑选
# 计算每个显著组合的中位数
median_values <- rt2 %>%
  group_by(Gene, risk) %>%
  summarise(median_expr = median(Expression, na.rm = TRUE))

#选择每个组中的最大值
max_values <- median_values %>%
  group_by(Gene) %>%
  summarise(max_expr = max(median_expr,na.rm = TRUE))

#排序后选择IC50值最小的前9个基因
top_values<-max_values%>%arrange(max_expr)
top_Gene<-as.character(top_values$Gene)[1:12]

#在原列表中获取
rt3<-rt2%>%filter(Gene %in% c(top_Gene))
rt3$Expression<-exp(rt3$Expression)

# # 四分位数范围（IQR）法去除异常值
# Q1 <- quantile(rt3$Expression, 0.25, na.rm = TRUE)
# Q3 <- quantile(rt3$Expression, 0.75, na.rm = TRUE)
# IQR <- Q3 - Q1
#
# # 定义异常值的界限
# lower_bound <- Q1 - 1.5 * IQR
# upper_bound <- Q3 + 1.5 * IQR
#
# # 过滤异常值
# filtered_rt3 <- rt3[rt3$Expression >= lower_bound & rt3$Expression <= upper_bound, ]

#绘图
#挑取每种药物的数据

rt3_Rapamycin<-rt3[rt3$Gene=="Rapamycin",]

rt3_Dactinomycin.1<-rt3[rt3$Gene=="Dactinomycin.1",]

rt3_Sepantronium.bromide<-rt3[rt3$Gene=="Sepantronium.bromide",]

rt3_MG.132<-rt3[rt3$Gene=="MG.132",]

rt3_Paclitaxel<-rt3[rt3$Gene=="Paclitaxel",]

rt3_Sabutoclax<-rt3[rt3$Gene=="Sabutoclax",]

rt3_Vinorelbine<-rt3[rt3$Gene=="Vinorelbine",]

rt3_Dactolisib<-rt3[rt3$Gene=="Dactolisib",]

rt3_AZD8055<-rt3[rt3$Gene=="AZD8055",]

rt3_Sabutoclax<-rt3[rt3$Gene=="Buparlisib",]

rt3_Luminespib<-rt3[rt3$Gene=="Luminespib",]

rt3_Staurosporine<-rt3[rt3$Gene=="AZD7762",]

#组合为列表
MedFile<-list(rt3_Rapamycin, rt3_Dactinomycin.1, rt3_Sepantronium.bromide,
              rt3_MG.132, rt3_Paclitaxel, rt3_Sabutoclax,
              rt3_Vinorelbine,rt3_Dactolisib,rt3_AZD8055,rt3_Sabutoclax,
              rt3_Luminespib,rt3_Staurosporine)

names(MedFile)<-c("Rapamycin", "Dactinomycin.1", "Sepantronium.bromide",
                  "MG.132", "Paclitaxel", "Sabutoclax",
                  "Vinorelbine","Dactolisib","AZD8055","Buparlisib",
                  "Luminespib","AZD7762")

#过滤异常值
for (i in names(MedFile)) {
  # 四分位数范围（IQR）法去除异常值
  Q1 <- quantile(MedFile[[i]]$Expression, 0.25, na.rm = TRUE)
  Q3 <- quantile(MedFile[[i]]$Expression, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1

  # 定义异常值的界限
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR

  # 过滤异常值
  MedFile[[i]] <- MedFile[[i]][MedFile[[i]]$Expression >= lower_bound & MedFile[[i]]$Expression <= upper_bound, ]

}


for(i in names(MedFile)){
  boxplot=ggboxplot(MedFile[[i]],x="Gene",y="Expression",fill = "risk",
                    xlab = "",
                    ylab = "IC50",
                    legend.title = "Risk",
                    width = 0.4,
                    palette = c("DodgerBlue1","Firebrick2"))+
    stat_compare_means(aes(group=risk),method = "wilcox.test",
                       symnum.args =list(cutpoints=c(0,0.001,0.01,0.05,1),
                                         symbols=c("***","**","*","ns")),label = "p.signif" )+
    theme(axis.text = element_text(face = "bold.italic",colour = "#441718",size = 12),
          axis.title = element_text(face = "bold.italic",colour = "#441718",size = 12),
          axis.line = element_blank(),
          plot.title = element_text(face = "bold.italic",colour = "#441718",size = 12),
          legend.text = element_text(face = "bold.italic"),
          panel.border = element_rect(fill = NA,color = "#35A79D",size = 1.5,linetype = "solid"),
          panel.background = element_rect(fill = "#F1F6FC"),
          panel.grid.major = element_line(color = "#CFD3D6",size = 0.5,linetype = "dotdash"),
          legend.title = element_text(face = "bold.italic",size = 10)
    )
  print(boxplot)
  ggsave(paste0(i,".pdf"),plot = boxplot,width = 10,height = 8,units = "in")
}





