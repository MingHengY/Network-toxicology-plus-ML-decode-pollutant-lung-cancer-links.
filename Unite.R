##导入并合并所有数据----

geoluad2_gse50081<-read.table("gse50081_geneset.txt",
                              header = T, sep = '\t', check.names = F, row.names = 1)

tcga_luad<-read.table("LUAD_geneset.txt",
                      header = T, sep = '\t', check.names = F, row.names = 1)
tcga_lusc<-read.table("lusc_geneset.txt",
                      header = T, sep = '\t', check.names = F, row.names = 1)

library(dplyr)
str1<-intersect(colnames(geoluad1_gse68465),colnames(geoluad2_gse50081))
str2<-intersect(colnames(geolusc1_gse30219),colnames(geolusc2_gse37745))

ppitable<-read.table("AP_LC Target.txt",
                     header = F, sep = '\t', check.names = F, row.names = 1)
tar<-rownames(ppitable)



library(VennDiagram)

venn.diagram(x=list('GSE68465'=colnames(geoluad1_gse68465)[4:130],
                    'GSE50081'=colnames(geoluad2_gse50081)[4:131],
                    'ppi'=tar),
             scaled = F, # 根据比例显示大小
             alpha= 0.5, #透明度
             lwd=1,lty=1,col=c("dodgerblue","darkorange1","#FF99FF"), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             cex = 2, # 数字大小
             fontface = "bold",  # 字体粗细；加粗bold
             fill=c("dodgerblue","darkorange1","#FF99FF"), # 填充色 配色https://www.58pic.com/
             category.names = c("GSE68465","GSE50081", "PPI") , #标签名
             cat.dist = 0.02, # 标签距离圆圈的远近
             cat.pos = -180, # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             cat.cex = 2, #标签字体大小
             cat.fontface = "bold",  # 标签字体加粗
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             output=TRUE,
             filename='Venn.png',# 文件保存
             imagetype="png",  # 类型（tiff png svg）
             resolution = 400,  # 分辨率
             compression = "lzw"# 压缩算法
)

#luad
#统一基因名
geoluad1_gse68465_Valid<-geoluad1_gse68465
geoluad2_gse50081_Valid<-geoluad2_gse50081
luadset<-tcga_luad
#统一生存时间单位
geoluad2_gse50081_Valid$OS.time<-geoluad2_gse50081_Valid$OS.time*12
luadset$OS.time<-luadset$OS.time/30

#抽样
set.seed(1234)#随机种子，确保结果可以复现
df1<-sample(nrow(luadset),103,replace = F)

TCGA_Vali<-luadset[df1,]
TCGA_Train<-luadset[-df1,]
write.table(TCGA_Train,file = "LUAD_geneset_80.txt",sep = "\t",quote = F, row.names = T)
#lusc
#统一基因名
geolusc1_gse30219_Valid<-geolusc1_gse30219[,str2]
geolusc2_gse37745_Valid<-geolusc2_gse37745[,str2]
luscset<-tcga_lusc[,str2]
#统一生存时间单位
geolusc2_gse37745_Valid$OS.time<-geolusc2_gse37745_Valid$OS.time/30
luscset$OS.time<-luscset$OS.time/30


library(Mime1)
##luad分析----
vali_Data_LUAD<-list(TCGA_Train,TCGA_Vali,geoluad2_gse50081_Valid)

names(vali_Data_LUAD)<-c("TCGA_LUAD_TRAIN","TCGA_LUAD_VALI","GSE50081")

ppitable<-read.table("gene.txt",
                     header = T, sep = '\t', check.names = F, row.names = 1)
tar<-ppitable$x

gene1<-tar


res_luad <- ML.Dev.Prog.Sig(train_data = vali_Data_LUAD[["TCGA_LUAD_TRAIN"]],
                        list_train_vali_Data = vali_Data_LUAD,
                        unicox.filter.for.candi = T,
                        unicox_p_cutoff = 0.05,
                        candidate_genes = gene1,
                        mode = 'all',nodesize =5,seed = 1234567 )


##C-index热图----
# 加载必要的包
library(tidyr)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(ggplotify)

#luad
# 1. 数据准备与清洗
# 假设数据已读取为data数据框
load("TCGA_LUAD_0409.RData")
Cindex<-res_luad[["Cindex.res"]]
data<-Cindex

# 统一ID格式（处理下划线和短横线）
data$ID <- gsub("_", "-", data$ID)

# 2. 转换数据为宽格式
wide_data <- data %>%
  pivot_wider(
    id_cols = Model,
    names_from = ID,
    values_from = Cindex
  )

# 3. 按TCGA-LIHC的Cindex值降序排列
sorted_models <- wide_data %>%
  arrange(desc(`TCGA-LUAD-TRAIN`))

# 4. 调整列顺序（TCGA-LUAD第一列，其他按原始顺序）
id_order <- c("TCGA-LUAD-TRAIN", setdiff(colnames(sorted_models)[-1], "TCGA-LUAD-TRAIN"))
sorted_data <- sorted_models[, c("Model", id_order)]

# #添加平均值
numeric_cols<-sapply(sorted_data,is.numeric)
Average<-rowMeans(sorted_data[,numeric_cols])
Cindex_mean<-data.frame(Model=sorted_data$Model,Mean=Average)
Cindex_mean$Mean<-round(Cindex_mean$Mean,digits = 3)

# 5. 创建热图矩阵
mat <- as.matrix(sorted_data[, -1])
rownames(mat) <- sorted_data$Model

# 6. 创建列注释（ID颜色块）
annotation_col <- data.frame(ID = factor(colnames(mat)))
rownames(annotation_col) <- colnames(mat)

#创建行注释
row_anno<-data.frame(Mean=rowMeans(mat,na.rm = T))

# 7. 定义颜色方案
id_colors <- list(
  ID = setNames(
    colorRampPalette(brewer.pal(8, "Set2"))(ncol(mat)),
    colnames(mat)
  )
)

# 8. 绘制热图

library(ggplot2)

# 设置条带长度的范围
min_value <- 0.55  # 最小值
max_value <- 0.75  # 最大值

# 创建条形图
p3 <- ggplot(data = Cindex_mean, aes(x = Model, y = Mean, fill = Mean)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#0084A7", high = "#E05D00") +  # 使用颜色渐变
  geom_text(aes(label = round(Mean, 3)), vjust = 0.5, hjust = 1.2, size = 2) +
  theme(
    axis.title = element_text(size = 8),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white")
  ) +
  labs(y = "") +
  coord_flip() +
  coord_cartesian(ylim = c(min_value, max_value))  # 设置 y 轴范围

# 打印图表
print(p3)

p_heat<-pheatmap(mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.3f",
         annotation_col = annotation_col,
         #annotation_row = row_anno,
         annotation_colors = id_colors,
         color = colorRampPalette(c("#457fff", "#ffffbf", "#d73027"))(100),
         main = "C-index Heatmap",
         fontsize_row = 8,
         fontsize_col = 8,
         fontsize_number = 7,
         number_color = "black",
         gaps_col = 1)

p_heat<-as.ggplot(p_heat)
p_combined<-p_heat %>% insert_left(p3,width = 0.3)


pdf("Cindex_of_all_models_LUAD_0409_2.pdf",width = 8,height = 13.5)


dev.off()


#lusc
# 1. 数据准备与清洗
# 假设数据已读取为data数据框
Cindex<-res_lusc[["Cindex.res"]]
data<-Cindex

# 统一ID格式（处理下划线和短横线）
data$ID <- gsub("_", "-", data$ID)

# 2. 转换数据为宽格式
wide_data <- data %>%
  pivot_wider(
    id_cols = Model,
    names_from = ID,
    values_from = Cindex
  )

# 3. 按TCGA-LUSC的Cindex值降序排列
sorted_models <- wide_data %>%
  arrange(desc(`TCGA-LUSC`))

# 4. 调整列顺序（TCGA-LUSC第一列，其他按原始顺序）
id_order <- c("TCGA-LUSC", setdiff(colnames(sorted_models)[-1], "TCGA-LUSC"))
sorted_data <- sorted_models[, c("Model", id_order)]

# 5. 创建热图矩阵
mat <- as.matrix(sorted_data[, -1])
rownames(mat) <- sorted_data$Model

# 6. 创建列注释（ID颜色块）
annotation_col <- data.frame(ID = factor(colnames(mat)))
rownames(annotation_col) <- colnames(mat)

#创建行注释
row_anno<-data.frame(Mean=rowMeans(mat,na.rm = T))

# 7. 定义颜色方案
id_colors <- list(
  ID = setNames(
    colorRampPalette(brewer.pal(8, "Set2"))(ncol(mat)),
    colnames(mat)
  )
)

# 8. 绘制热图
pdf("Cindex_of_all_models_LUSC.pdf",width = 8,height = 13.5)
pheatmap(mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.3f",
         annotation_col = annotation_col,
         #annotation_row = row_anno,
         annotation_colors = id_colors,
         color = colorRampPalette(c("#457f99", "#ffffbf", "#d73027"))(100),
         main = "C-index Heatmap",
         fontsize_row = 8,
         fontsize_col = 8,
         fontsize_number = 7,
         number_color = "black",
         gaps_col = 1)
dev.off()






#逐步回归
##SetpCox_LUAD----

load("TCGA_LUAD_0407.RData")


library(My.stepwise)
library(survival)
library(survminer)
rt<-read.table("LUAD_geneset_80.txt",header=T,sep="\t",row.names=1,check.names=F)
rt1<-rt
rownames(rt)<-rt$ID
rt<-rt[,-1]
# 单因素COX回归分析

#设置p值的阈值
pfilter1 <- 0.05
#新建空白数据框
uniresult <- data.frame()
#使用for循环对输入数据中的100个基因依次进行单因素COX分析
#单因素COX回归分析中p值＜0.05的基因，其分析结果输入到之前新建的空白数据框uniresult中
for(i in colnames(rt[,3:ncol(rt)])){
  unicox <- coxph(Surv(time = OS.time, event = OS) ~ rt[,i], data = rt)
  unisum<- summary(unicox)
  pvalue <- round(unisum$coefficients[,5],3)
  if(pvalue<pfilter1){
    uniresult <- rbind(uniresult,
                       cbind(gene=i,
                             HR=unisum$coefficients[,2],
                             L95CI=unisum$conf.int[,3],
                             H95CI=unisum$conf.int[,4],
                             pvalue=unisum$coefficients[,5]
                       ))
  }
}
#保存单因素COX回归分析结果
library(dplyr)
unids<-dplyr::select(rt1,ID,OS.time,OS,uniresult$gene)

write.table(uniresult,file="UniCoxSurvival_LUAD_1.txt",sep="\t",row.names=F,quote=F)
write.table(unids,file="UniCoxData_LUAD_1.txt",sep="\t",row.names=F,quote=F)

##单因素COX回归森林图----
tducs <- read.table("UniCoxSurvival_LUAD_1.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(tducs)
#将HR及其置信区间数据组合成一个因子
hr <- sprintf("%.3f",tducs$"HR")
hrLow  <- sprintf("%.3f",tducs$"L95CI")
hrHigh <- sprintf("%.3f",tducs$"H95CI")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pValue <- ifelse(tducs$pvalue<0.001, "<0.001", sprintf("%.3f", tducs$pvalue))
#绘图
pdf(file="UniCoxSurForestPlot_LUAD_1.pdf", width = 6,height = 3)
n <- nrow(tducs)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(2.5,2))

#开始画图的左边部分，即图1，森林图的文字及数值标注部分
xlim = c(0,2.5)
par(mar=c(4,2.5,2,1))
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.7
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.1,n:1,pValue,adj=1,cex=text.cex);text(1.1,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(2.5,n:1,Hazard.ratio,adj=1,cex=text.cex);text(2.5,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

#开始画图的右边部分，即图2，森林图的图形
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 1, col = boxcolor, cex=1.3)
axis(1)

dev.off()

##LASSO分析----
library(glmnet)
rt<-read.table("UniCoxData_LUAD_1.txt",header=T,sep="\t",row.names=1,check.names=F)
#随机种子
set.seed(123456)
#构建LASSO回归模型
library(survival)
library(survminer)
x=as.matrix(rt[,c(3:ncol(rt))])
scaleX<-scale(x)
y=data.matrix(Surv(rt$OS.time,rt$OS))
fit=glmnet(x,y,family = "cox",nfolds=10)

#c-index
cvfit<-cv.glmnet(x,y,family= "cox",type.measure = "C",nfolds = 10)
pdf("lasso_c-index_1.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()
#deviance图形
cvfit<-cv.glmnet(x,y,family= "cox",type.measure = "deviance",nfolds = 10)
pdf("lasso_cvfit_1.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()
#Coefficients图形
pdf("lasso_lambda_1.pdf")
plot(fit,xvar = "lambda",label = T)
abline(v=log(cvfit$lambda.min),lty="dashed")
dev.off()
#输出lasso显著基因表达量
coef=coef(fit, s=cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
row.names(coef)[index]
actCoef
index
lassoSigExp=rt[,c("OS.time","OS",lassoGene)]
lassoSigExpOut=cbind(id=row.names(lassoSigExp),lassoSigExp)
write.csv(lassoSigExpOut,file = "LassoResults_1.csv",row.names = F)
LASSOsigGene<-data.frame(gene=row.names(coef)[index],Coef=actCoef)
write.table(LASSOsigGene,file = "LASSO_Significant_Gene_1.txt",sep = "\t",quote = F, row.names = F)




##随机森林----
library(dplyr)
lasso=read.csv("LassoResults_1.csv")
rownames(lasso)

rt<-lasso
rownames(rt)<-rt$id
rt$id<-NULL

library(randomForestSRC)
library(rfPermute)
library(survival)
#构建RSF生存模型
fit <- rfsrc(Surv(OS.time,OS)~.,data = rt,
             ntree = 1000,
             nodesize = 10,
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = 1234)
#筛选重要基因
library(dplyr)
library(tibble)
importance_gene <- data.frame(fit$importance) %>%
  rownames_to_column("gene") %>%
  arrange(- fit.importance) %>%
  head(20)

importance_gene

#可视化
pdf(file = "trees_250409_lasso.pdf",width = 10,height = 6)
plot(fit,10)
dev.off()

#使用ggplot2绘制柱形图，使用reorder函数进行排序
library(ggplot2)
pr<-ggplot(data=importance_gene, aes(x = reorder(gene,  fit.importance),
                                     y=fit.importance,fill=gene)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(legend.position = 'none') +
  coord_flip()
pr
ggsave(filename = "RFS_gene_lasso_0409.pdf", plot = pr, width = 10, height = 8, units = "in")

RFS_gene<-importance_gene
RFS_res<-dplyr::select(rt, OS.time,OS,RFS_gene$gene)

write.table(RFS_res,"RSFresults_lasso_250409.txt",sep = "\t",quote = F, row.names = T)

# 计算变量重要性的阈值
var_importance<-vimp(fit)$importance
# 计算变量重要性的阈值（保留前20%的变量）
median<-median(var_importance)
# 选择重要性得分高于阈值的变量
selected_vars <- names(var_importance)[var_importance > median]

# 查看选中的变量
print(selected_vars)
s<-importance_gene[importance_gene$fit.importance>median,]
write.table(s,"Lasso_RSF_selected_0409.txt",sep = "\t",quote = F, row.names = T)




