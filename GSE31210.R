library(limma)
library(sva)
library(GEOquery)
##整理下载数据----
id5<-"GSE31210"
gse5<-getGEO(id5,getGPL = F)
gse5 <- gse5[[1]]
summary(exprs(gse5))
exprs5 <- exprs(gse5)
#提取注释和表达信息
annotation5 <- featureData(gse5)
# 读取芯片平台文件txt
platform_file5 <- read.delim("GPL570-55999.txt", header = TRUE, sep = "\t", comment.char = "#")
#查看平台文件列名
colnames(platform_file5)

##数据处理----
# 假设芯片平台文件中有两列，一列是探针ID，一列是基因名
probe_names <- platform_file5$ID
gene_symbol <- platform_file5$Gene.Symbol

platform_file_set=platform_file5[,c(1,11)]

#一个探针对应多个基因名，保留第一个基因名
ids = platform_file_set
library(tidyverse)
test_function <- apply(ids,
                       1,
                       function(x){
                         paste(x[1],
                               str_split(x[2],'///',simplify=T),
                               sep = "...")
                       })
x = tibble(unlist(test_function))

colnames(x) <- "ttt"
ids <- separate(x,ttt,c("ID","Gene.Symbol"),sep = "\\...")
dim(ids)


#将Matrix格式表达矩阵转换为data.frame格式
exprSet <- data.frame(exprs5)
dim(exprSet)

#给表达矩阵新增加一列ID
exprSet$ID <- rownames(exprSet) # 得到表达矩阵，行名为ID，需要转换，新增一列
dim(exprSet)
#矩阵表达文件和平台文件有相同列‘ID’，使用merge函数合并
express <- merge(x = exprSet, y = ids, by.x = "ID")

#删除探针ID列
express$ID =NULL


dim(express)

matrix = express
dim(matrix)
#查看多少个基因重复了
table(duplicated(matrix$Gene.Symbol))


#把重复的Symbol取平均值
matrix <- aggregate(.~Gene.Symbol, matrix, mean)  ##把重复的Symbol取平均值
row.names(matrix) <- matrix$Gene.Symbol  #把行名命名为SYMBOL

dim(matrix)

matrix_na = na.omit(matrix)   #去掉缺失值
dim(matrix_na)

matrix_final = matrix_na[matrix_na$Gene.Symbol != "",]
dim(matrix_final)

matrix_final <- subset(matrix_final, select = -1)  #删除Symbol列（一般是第一列）
dim(matrix_final)
#+  经过注释、探针名基因名处理、删除基因名为空值、删除缺失值 得到最终 matrix_final

#数据离群处理
#处理极端值
#定义向量极端值处理函数
#用于处理异常值，将超出一定范围的值替换为中位数，以减少异常值对后续分析的影响。
dljdz=function(x) {
  DOWNB=quantile(x,0.25)-1.5*(quantile(x,0.75)-quantile(x,0.25))
  UPB=quantile(x,0.75)+1.5*(quantile(x,0.75)-quantile(x,0.25))
  x[which(x<DOWNB)]=quantile(x,0.5)
  x[which(x>UPB)]=quantile(x,0.5)
  return(x)
}

#第一列设置为行名
matrix_leave=matrix_final

boxplot(matrix_leave,outline=FALSE, notch=T, las=2)  ##出箱线图
dim(matrix_leave)

#处理离群值
matrix_leave_res=apply(matrix_leave,2,dljdz)

boxplot(matrix_leave_res,outline=FALSE, notch=T, las=2)  ##出箱线图
dim(matrix_leave_res)


### limma的函数归一化，矫正差异  ，表达矩阵自动log2化

#1.归一化不是绝对必要的，但是推荐进行归一化。
#有重复的样本中，应该不具备生物学意义的外部因素会影响单个样品的表达，
#例如中第一批制备的样品会总体上表达高于第二批制备的样品，假设所有样品表达值的范围和分布都应当相似，
#需要进行归一化来确保整个实验中每个样本的表达分布都相似。
#2.归一化要在log2标准化之前做


library(limma)

exprSet=normalizeBetweenArrays(matrix_leave_res)

boxplot(exprSet,outline=FALSE, notch=T, las=2)  ##出箱线图

## 这步把矩阵转换为数据框很重要
class(exprSet)   ##注释：此时数据的格式是矩阵（Matrix）
exprSet <- as.data.frame(exprSet)


#标准化 表达矩阵自动log2化
qx <- as.numeric(quantile(exprSet, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)


## 开始判断
if (LogC) {
  # 假设exprSet是一个多列的数据框，我们需要对每一列进行操作
  exprSet_clean <- exprSet  # 创建一个新的数据框来存储清理后的数据

  # 对exprSet中的每一列进行操作
  for (i in seq_along(exprSet)) {
    # 使用ifelse函数来替换小于或等于0的值为NaN
    exprSet_clean[[i]] <- ifelse(exprSet[[i]] <= 0, NaN, exprSet[[i]])

    # 对非NaN值取log2
    exprSet_clean[[i]][!is.nan(exprSet_clean[[i]])] <- log2(exprSet_clean[[i]][!is.nan(exprSet_clean[[i]])])
  }

  print("log2 transform finished")
  exprSet <- exprSet_clean  # 将清理后的数据框赋值回exprSet
} else {
  print("log2 transform not needed")
}

boxplot(exprSet_clean,outline=FALSE, notch=T, las=2)  ##出箱线图
boxplot(exprSet,outline=FALSE, notch=T, las=2)  ##如果数据已经log2化了
### limma的函数归一化，矫正差异  ，表达矩阵自动log2化 得到exprSet_clean

# saving all data to the path
save(exprSet_clean, file ="exprSet_gse31210.RData")


# #临床信息
# patient2<-data.frame(ID=gse2@phenoData@data[["geo_accession"]],OS.time=gse2@phenoData@data[["characteristics_ch1.8"]],
#                      OS=gse2@phenoData@data[["status:ch1"]],Type=gse2@phenoData@data[["characteristics_ch1.1"]])
#
# library(stringr)
# library(dplyr)
#
# patient2<-patient2[patient2$Type=="histology: adenocarcinoma",]
# patient2$Type<-NULL
# patient2$OS.time <- as.numeric(str_extract(patient2$OS.time, "\\d+\\.?\\d*"))
#
# patient2$OS<-ifelse(patient2$OS == "alive",0,1)
# write.table(patient2,"GSE50081Patient.txt",sep = "\t",quote = F, row.names = T)
#

geneset<-read.table("Lasso_RSF_selected_0409.txt",header=T,sep="\t",row.names=1,check.names=F)
tar<-geneset$gene

#分组
load("exprSet_gse31210.RData")

info<-data.frame(ID=gse5@phenoData@data[["geo_accession"]],
                 sample=gse5@phenoData@data[["tissue:ch1"]])

tumor_group<-info[info$sample=="primary lung tumor",]
normal_group<-info[info$sample=="normal lung",]

exprSet1<-exprSet_clean[tar,]
tumor<-exprSet1[,tumor_group$ID]
normal<-exprSet1[,normal_group$ID]

conNum1<-ncol(normal)
treatNum1<-ncol(tumor)

expr1<-cbind(normal,tumor)
expr1<-as.matrix(expr1)

##几个基因的表达热图
library(limma)
library(pheatmap)
library(stringr)
library(edgeR)
library(ComplexHeatmap)

annotation_col1 = data.frame(
  Type=c(rep("Normal",conNum1),rep("Tumor",treatNum1))
)
rownames(annotation_col1) = colnames(expr1)
#构建列的注释信息

#绘制热图
pdf(file = "candidate_genes_GSE31210.pdf",height = 5,width = 8)
ComplexHeatmap::pheatmap(expr1,
                         annotation_col = annotation_col1,
                         annotation_colors = list(Type=c( "Tumor"= "#E31A1C", "Normal"= "#5A8FCA")),
                         color = colorRampPalette(c(rep('blue',5),'white',rep('red',5)))(50),
                         cluster_rows = T,
                         cluster_cols = F,
                         show_colnames = F,
                         scale = "row",
                         fontsize = 5,
                         fontsize_row = 8,
                         fontsize_col = 8)
dev.off()

