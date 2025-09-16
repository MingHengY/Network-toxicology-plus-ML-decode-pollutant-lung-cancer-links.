##差异分析----
library(DESeq2)
library(limma)
library(pheatmap)
library(stringr)
library(edgeR)
library(dplyr)
data=read.table("TCGA-LUAD.star_counts.tsv",header = T,sep = '\t',check.names = F,row.names = 1)
platform_file=read.delim("gencode.v36.annotation.gtf.gene.probemap", header = TRUE, sep = "\t")
plat<-platform_file[,c(1,2)]
data$id<-rownames(data)
data<-merge(plat,data,by="id")

data<-data[,-1]
#对重复基因合并取平均值
luad_gset = as.data.frame(avereps(data[,-1],ID = data$gene) )
luad_gset1<-as.data.frame(t(luad_gset))

data<-luad_gset1
data<-2^data-1

dimnames=list(rownames(data),colnames(data))
data=matrix(as.numeric(as.matrix(data)),nrow = nrow(data),dimnames = dimnames)
data<-t(data)

#去除低表达的基因
data=data[rowMeans(data)>1,]

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
data=cbind(data2,data1)

#DESeq
#分组信息

group <- c(rep('tumor', treatNum), rep('normal', conNum))
group <- factor(group, levels = c("normal", "tumor"))

colData <- data.frame(row.names = colnames(data),
                      group = group)
colData$group <- factor(colData$group, levels = c("normal", "tumor"))

#整数化
exp_int <- apply(data, 2, as.integer)
rownames(exp_int) <- rownames(data)

# 构建DESeqDataSet对象，也就是dds矩阵，将基因计数数据、样本分组信息和设计矩阵关联起来
dds <- DESeqDataSetFromMatrix(countData = exp_int, # 表达矩阵
                              colData = colData,        # 表达矩阵列名和分组信息的对应关系
                              design = ~ group)         # group为colData中的group，也就是分组信息
# 查看一下构建好的矩阵
head(dds)

# 进行差异表达分析
dds <- DESeq(dds)

# 查看结果的名称
resultsNames(dds)

# 提取差异表达结果，进行对比，这里contrast参数指定了对比的组别
# contrast参数必须写成下面三个元素的向量格式，且顺序不能反
res <- results(dds, contrast = c("group", rev(levels(group))))
# res <- results(dds, contrast = c("group", levels(group)[2], levels(group)[1]))

# 按照padj（调整后的p值）的大小对差异结果进行排序（只有DESeq2需要，limma和edgeR会自动排好）
resOrdered <- res[order(res$padj), ]

# 将差异表达结果转换为数据框
DEG <- as.data.frame(resOrdered)

# 去除缺失值，如果没有这一步，一些表达量很低的基因计算后会出现NA，给后续分析和绘图带来麻烦，远离麻烦！
DEG_deseq2 <- na.omit(DEG)

# 将处理后的差异表达结果保存为R数据文件
save(DEG_deseq2, file = 'DEG_deseq2.Rdata')

load("./data/DEG_deseq2.Rdata")

# 加change列,标记上下调基因，可根据需求设定阈值
logFC = 1.2
P.Value = 0.05
k1 <- (DEG_deseq2$pvalue < P.Value) & (DEG_deseq2$log2FoldChange < -logFC)
k2 <- (DEG_deseq2$pvalue < P.Value) & (DEG_deseq2$log2FoldChange > logFC)
DEG_deseq2 <- mutate(DEG_deseq2, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))

table(DEG_deseq2$change)

# 火山图
library(ggplot2)
p <- ggplot(data = DEG_deseq2,
            aes(x = log2FoldChange,
                y = -log10(pvalue))) +
  geom_point(alpha = 0.4, size = 3.5,
             aes(color = change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values = c("blue4", "grey", "red3"))+
  geom_vline(xintercept = c(-logFC, logFC), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(P.Value), lty = 4, col = "black", lwd = 0.8) +
  theme_bw()
p

ggsave(filename = "volcano_plot_deseq2.pdf", plot = p, device = "pdf", width = 6, height = 5)
dev.off()

# 差异基因热图
data<-as.data.frame(data)
deg_opt <- DEG_deseq2 %>% filter(DEG_deseq2$change != "stable")
exp_brca_heatmap <- data %>% filter(rownames(data) %in% rownames(deg_opt))
annotation_col <- data.frame(group = group)
rownames(annotation_col) <- colnames(exp_brca_heatmap)

p1 <- pheatmap(exp_brca_heatmap, show_colnames = F, show_rownames = F,
               scale = "row",
               cluster_cols = F,
               annotation_col = annotation_col,
               breaks = seq(-3, 3, length.out = 100))
p1

ggsave(filename = "heatmap_plot_deseq2.pdf", plot = p1, device = "pdf", width = 5, height = 6)
dev.off()

write.table(deg_opt,file = "TCGA_LUAD_diff_count_DESeq.txt",sep = "\t",quote = F,row.names = T)


#取交集
ppitable<-read.table("AP_LC Target.txt",
                     header = F, sep = '\t', check.names = F, row.names = 1)
tar<-rownames(ppitable)

library(VennDiagram)

diffdata<-read.table("TCGA_LUAD_diff_count_DESeq.txt",header = T,sep = '\t',check.names = F,row.names = 1)
diff<-rownames(diffdata)

venn.diagram(x=list('DEGs'=diff,
                    'ppi'=tar),
             scaled = F, # 根据比例显示大小
             alpha= 0.5, #透明度
             lwd=1,lty=1,col=c("dodgerblue","darkorange1"), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             cex = 2, # 数字大小
             fontface = "bold",  # 字体粗细；加粗bold
             fill=c("dodgerblue","darkorange1"), # 填充色 配色https://www.58pic.com/
             category.names = c("DEGs", "PPI") , #标签名
             cat.dist = 0.02, # 标签距离圆圈的远近
             cat.pos = -180, # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             cat.cex = 2, #标签字体大小
             cat.fontface = "bold",  # 标签字体加粗
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             output=TRUE,
             filename='Venn_diffvsppi_0408.png',# 文件保存
             imagetype="png",  # 类型（tiff png svg）
             resolution = 400,  # 分辨率
             compression = "lzw"# 压缩算法
)

selected_gene<-intersect(diff,tar)
write.table(selected_gene,file = "gene.txt",sep = "\t",quote = F, row.names = T)



