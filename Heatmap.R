##绘制热图----
#TCGA
library(DESeq2)
library(limma)
library(pheatmap)
library(stringr)
library(edgeR)

data=read.table("TCGA-LUAD.gene_expression_TPM.tsv",header = T,sep = '\t',check.names = F,row.names = 1)

dimnames=list(rownames(data),colnames(data))
data=matrix(as.numeric(as.matrix(data)),nrow = nrow(data),dimnames = dimnames)

#正常组和肿瘤组数目,第14、15字符，01-09：癌症，10-19：正常；20-29：癌旁
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
#将2变为1
group=gsub("2","1",group)
#获取正常与肿瘤组数目
conNum=length(group[group==1])
treatNum=length(group[group==0])
Type=c(rep(1,conNum),rep(2,treatNum))
#根据正常组和肿瘤组排序
data1=data[,group==1]
data2=data[,group==0]
data=cbind(data1,data2)

data<-as.data.frame(data)
geneset<-read.table("Lasso_RSF_selected_0409.txt",header=T,sep="\t",row.names=1,check.names=F)
expr<-data[geneset$gene,]
expr<-log2(expr+1)
##几个基因的表达热图
library(limma)
library(pheatmap)
library(stringr)
library(edgeR)
library(ComplexHeatmap)

annotation_col = data.frame(
  Type=c(rep("Normal",conNum),rep("Tumor",treatNum))
)
rownames(annotation_col) = colnames(expr)
#构建列的注释信息

#绘制热图
pdf(file = "candidate_genes_TCGALUAD_2.pdf",height = 5,width = 8)
ComplexHeatmap::pheatmap(expr,
                         annotation_col = annotation_col,
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


#GEO

geneset<-read.table("Lasso_RSF_selected_0409.txt",header=T,sep="\t",row.names=1,check.names=F)
tar<-geneset$gene

#GSE68465
load("exprSet_gse68465.RData")

library(limma)
library(sva)
library(GEOquery)
##整理下载数据，获得样本分类信息
id1<-"GSE68465"
gse1<-getGEO(id1,getGPL = F)
gse1 <- gse1[[1]]
summary(exprs(gse1))

info<-data.frame(ID=gse1@phenoData@data[["geo_accession"]],
                 sample=gse1@phenoData@data[["characteristics_ch1"]])

tumor_group<-info[info$sample=="disease_state: Lung Adenocarcinoma",]
normal_group<-info[info$sample=="disease_state: Normal",]

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
pdf(file = "candidate_genes_GSE68465_1.pdf",height = 5,width = 8)
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



#GSE50081
# load("exprSet_gse50081.RData")
# ##整理下载数据，获得样本分类信息
# library(limma)
# library(sva)
# library(GEOquery)
#
# id2<-"GSE50081"
# gse2<-getGEO(id2,getGPL = F)
# gse2 <- gse2[[1]]
# summary(exprs(gse2))

data=read.table("TCGA-LUAD.star_tpm.tsv",header = T,sep = '\t',check.names = F,row.names = 1)
platform_file=read.delim("gencode.v36.annotation.gtf.gene.probemap", header = TRUE, sep = "\t")
plat<-platform_file[,c(1,2)]
data$id<-rownames(data)
luad<-merge(plat,data,by="id")
library(limma)
library(dplyr)

luad<-luad[,-1]

#对重复基因合并取平均值
luad_gset = as.data.frame(avereps(luad[,-1],ID = luad$gene) )
luad_gset1<-luad_gset[tar,]
luad_gset1<-as.data.frame(t(luad_gset1))

dimnames=list(rownames(luad_gset1),colnames(luad_gset1))
luad_gset2=matrix(as.numeric(as.matrix(luad_gset1)),nrow = nrow(luad_gset1),dimnames = dimnames)
luad_gset2<-as.data.frame(t(luad_gset2))

#正常组和肿瘤组数目,第14、15字符，01-09：癌症，10-19：正常；20-29：癌旁
group=sapply(strsplit(colnames(luad_gset2),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
#将2变为1
group=gsub("2","1",group)
#获取正常与肿瘤组数目
conNum=length(group[group==1])
treatNum=length(group[group==0])
Type=c(rep(1,conNum),rep(2,treatNum))
#根据正常组和肿瘤组排序
con<-luad_gset2[,group==1]
treat<-luad_gset2[,group==0]
data<-cbind(con,treat)
expr2<-as.matrix(data)

##几个基因的表达热图
library(limma)
library(pheatmap)
library(stringr)
library(edgeR)
library(ComplexHeatmap)

annotation_col2 = data.frame(
  Type=c(rep("Normal",conNum),rep("Tumor",treatNum))
)
rownames(annotation_col2) = colnames(expr2)
#构建列的注释信息

#绘制热图
pdf(file = "candidate_genes_TCGALUAD_1.pdf",height = 5,width = 8)
ComplexHeatmap::pheatmap(expr2,
                         annotation_col = annotation_col2,
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

#合并的数据集
load("F:/TCGA_GTEx_pancancer_mrna_pheno.rdata")
LUAD_proj<-tcga_gtex_mrna_pheno[tcga_gtex_mrna_pheno$project=="LUAD",]
rownames(LUAD_proj)<-LUAD_proj$sample_id

LUAD_tumor<-LUAD_proj[LUAD_proj$sample_type=="TCGA_tumor",]
LUAD_normal<-LUAD_proj[LUAD_proj$sample_type!="TCGA_tumor",]

LUAD_tumor1<-LUAD_tumor[,c(5:ncol(LUAD_tumor))]
LUAD_normal1<-LUAD_normal[,c(5:ncol(LUAD_normal))]

LUAD_expr<-rbind(LUAD_normal1,LUAD_tumor1)

geneset<-read.table("Lasso_RSF_selected_0409.txt",header=T,sep="\t",row.names=1,check.names=F)
tar<-geneset$gene

LUAD_expr1<-LUAD_expr[,tar]
LUAD_expr2<-as.matrix(t(LUAD_expr1))


##几个基因的表达热图
library(limma)
library(pheatmap)
library(stringr)
library(edgeR)
library(ComplexHeatmap)

annotation_col2 = data.frame(
  Type=c(rep("Normal",nrow(LUAD_normal1)),rep("Tumor",nrow(LUAD_tumor1)))
)
rownames(annotation_col2) = colnames(LUAD_expr2)
#构建列的注释信息

#绘制热图
pdf(file = "candidate_genes_LUAD_TCGA_withGTE_1.pdf",height = 5,width = 8)
ComplexHeatmap::pheatmap(LUAD_expr2,
                         annotation_col = annotation_col2,
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



KEGG=read.table("enrichment.KEGG.tsv",header = F,sep = '\t',check.names = F,row.names = 1)
colnames(KEGG)<-c("description","count in network","total genes in pathway","strength","signal","fdr","Ensemble ID","Gene Symbol")
write.csv(KEGG,file = "KEGG.csv")
