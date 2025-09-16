
ppitable<-read.table("gene.txt",
                     header = T, sep = '\t', check.names = F, row.names = 1)
tar<-ppitable$x


##TCGA-LUAD----
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
luad_gset3=luad_gset2[,group==0]
luad_gset3<-as.data.frame(t(luad_gset3))
luad_gset3$ID<-rownames(luad_gset3)

#临床数据
luad_cli=read.table("TCGA-LUAD.survival.tsv",header = T,sep = '\t',check.names = F,row.names = 1)
luad_cli<-luad_cli[,c(1,2)]
luad_gset3$ID<-rownames(luad_gset3)
luad_cli$ID<-rownames(luad_cli)
luad_cli<-luad_cli[,c("ID","OS.time","OS")]

luad<-merge(luad_cli,luad_gset3,by="ID")
write.table(luad,file = "LUAD_geneset.txt",sep = "\t",quote = F, row.names = T)


##LUSC----
data=read.table("TCGA-LUSC.star_tpm.tsv",header = T,sep = '\t',check.names = F,row.names = 1)
platform_file=read.delim("gencode.v36.annotation.gtf.gene.probemap", header = TRUE, sep = "\t")
plat<-platform_file[,c(1,2)]
data$id<-rownames(data)
chol<-merge(plat,data,by="id")
library(limma)
library(dplyr)

chol<-chol[,-1]
#对重复基因合并取平均值
chol_gset = as.data.frame(avereps(chol[,-1],ID = chol$gene) )
chol_gset1<-chol_gset[tar,]
chol_gset1<-as.data.frame(t(chol_gset1))

#临床数据
chol_cli=read.table("TCGA-LUSC.survival.tsv",header = T,sep = '\t',check.names = F,row.names = 1)
chol_cli<-chol_cli[,c(1,2)]
chol_gset1$ID<-rownames(chol_gset1)
chol_cli$ID<-rownames(chol_cli)
chol_cli<-chol_cli[,c("ID","OS.time","OS")]

chol<-merge(chol_cli,chol_gset1,by="ID")
write.table(chol,file = "LUSC_geneset.txt",sep = "\t",quote = F, row.names = T)



##差异分析----
library(DESeq2)
library(limma)
library(pheatmap)
library(stringr)
library(edgeR)
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
data=cbind(data1,data2)
#设置分组信息
Type=factor(Type)
design<-model.matrix(~0+Type)
rownames(design)=colnames(data)
colnames(design)<-levels(Type)
#创建对象
DGElist<-DGEList(counts = data,group = Type)
#过滤CPM小于等于1的基因
keep_gene<-rowSums(cpm(DGElist)>1)>=2
table(keep_gene)
#生成DGElist对象
DGElist<-DGElist[keep_gene,,keep.lib.sizes = FALSE]
#计算比例因子，以便将原始库大小转化为有效库大小
DGElist<-calcNormFactors(DGElist)
#转换logCPM
v<-voom(DGElist,design,plot = TRUE,normalize.method ="quantile")
#非线性最小二乘法
fit<-lmFit(v,design)
colnames(design)=c("normal","tumor")
#构建比较矩阵
cont.matrix<-makeContrasts(contrasts = c('tumor-normal'),levels = design)
#线性拟合模型构建
fit2<-contrasts.fit(fit,cont.matrix)
#用经验贝叶斯调整t-test中方差的部分
fit2<-eBayes(fit2)
nrDEG_limma_voom=topTable(fit2, coef = 'tumor-normal',n=Inf)
nrDEG_limma_voom=na.omit(nrDEG_limma_voom)
head(nrDEG_limma_voom)
#筛选标准
padj=0.05
logFC=1
#筛选
outDiff=nrDEG_limma_voom[(nrDEG_limma_voom$adj.P.Val<padj&
                            (nrDEG_limma_voom$logFC>logFC|nrDEG_limma_voom$logFC<(-logFC))),]
outDiff=outDiff[order(outDiff$logFC),]
write.table(data.frame(ID=rownames(outDiff),outDiff),file = "TCGA_LUAD_diff_count.txt",sep = "\t",quote = F,row.names = F)
#火山图
logFCfilter=1
fdrFilter=0.05
pdf(file = "vol.pdf",width = 5,height = 5)
xMax=6
yMax=max(-log10(nrDEG_limma_voom$adj.P.Val))+1
logFC=nrDEG_limma_voom$logFC
fdr=nrDEG_limma_voom$adj.P.Val
plot(as.numeric(as.vector(nrDEG_limma_voom$logFC)),-log10(nrDEG_limma_voom$adj.P.Val),
     xlab = "logFC",ylab = "-log10(adj.P.Val)",main = "Volcano",
     ylim = c(0,yMax),xlim = c(-xMax,xMax),yaxs = "i",pch=20,cex=1.2)
diffsub=subset(nrDEG_limma_voom,fdr<fdrFilter&as.numeric(as.vector(logFC))>logFCfilter)
points(as.numeric(as.vector(diffsub$logFC)),-log10(diffsub$adj.P.Val),pch=20,col='red',cex=1.5)
diffsub=subset(nrDEG_limma_voom,fdr<fdrFilter&as.numeric(as.vector(logFC))<(-logFCfilter))
points(as.numeric(as.vector(diffsub$logFC)),-log10(diffsub$adj.P.Val),pch=20,col='green',cex=1.5)
abline(v=0,lty=2,lwd=3)
dev.off()

##GO与KEGG----

library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(ggplot2)
library(stringr)
library(openxlsx)
library(KEGGREST)

# sig<-read.table("TCGA_LIHC_diff_count.txt",
#                 header = T, sep = '\t', check.names = F, row.names = 1)
# key=rownames(sig)

ppitable<-read.table("AP_LC Target.txt",
                     header = F, sep = '\t', check.names = F, row.names = 1)
tar<-rownames(ppitable)
key=tar

library(org.Hs.eg.db)
#基因ID转换#
keytypes(org.Hs.eg.db) #查看所有可转化类型
entrezid_all = mapIds(x = org.Hs.eg.db,  #id转换的比对基因组（背景基因）的物种，以人为例
                      keys = key, #将输入的gene_name列进行数据转换
                      keytype = "SYMBOL", #输入数据的类型
                      column = "ENTREZID")#输出数据的类型
entrezid_all  = na.omit(entrezid_all)  #na省略entrezid_all中不是一一对应的数据情况
entrezid_all = data.frame(entrezid_all) #将entrezid_all变成数据框格式
head(entrezid_all)
#GO分析
GO_enrich = enrichGO(gene = entrezid_all[,1], #表示前景基因，即待富集的基因列表;[,1]表示对entrezid_all数据集的第1列进行处理
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID", #输入数据的类型
                     ont = "ALL", #可以输入CC/MF/BP/ALL
                     #universe = 背景数据集 # 表示背景基因，无参的物种选择组装出来的全部unigenes作为背景基因；有参背景基因则不需要。
                     pvalueCutoff = 1,qvalueCutoff = 1, #表示筛选的阈值，阈值设置太严格可导致筛选不到基因。可指定 1 以输出全部
                     readable = T) #是否将基因ID映射到基因名称。
GO_enrich  = data.frame(GO_enrich) #将GO_enrich导成数据框格式
write.csv(GO_enrich,'GO_enrich_aplcTar.csv')


#KEGG富集分析
KEGG_enrich = enrichKEGG(gene = entrezid_all[,1], #即待富集的基因列表
                         keyType = "kegg",
                         pAdjustMethod = 'fdr',  #指定p值校正方法
                         organism= "human",  #hsa，可根据你自己要研究的物种更改，可在https://www.kegg.jp/brite/br08611中寻找
                         qvalueCutoff = 1, #指定 p 值阈值（可指定 1 以输出全部）
                         pvalueCutoff=1) #指定 q 值阈值（可指定 1 以输出全部）
KEGG_enrich  = data.frame(KEGG_enrich)
write.csv(KEGG_enrich,'KEGG_enrich_aplcTar.csv')



GO_enrich$term <- GO_enrich$Description
# GO_enrich$term <- sub("^(.*?\\:.*?\\:.*?).*", "\\1", GO_enrich$term)

for(i in 1:nrow(GO_enrich)){
  description_splite=strsplit(GO_enrich$term[i],split = " ")
  #选择前五个单词GO term名字
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ")
  GO_enrich$term[i]=description_collapse
  #gsub()可以用于字段的删减、增补、替换和切割，
  #可以处理一个字段也可以处理由字段组成的向量。gsub(“目标字符”, “替换字符”, 对象)
  GO_enrich$term=gsub(pattern = "NA","",GO_enrich$term)
}
unique_terms <- unique(GO_enrich$term[GO_enrich$term != ""])
GO_enrich$term <- factor(GO_enrich$term, levels = unique_terms, ordered = TRUE)


#每一种选择几个进行绘图
library(dplyr)
# 筛选BP的前10个结果
bp_top10 <- GO_enrich %>%
  filter(ONTOLOGY == "BP") %>%
  slice_head(n = 10)

# 筛选CC的前10个结果
cc_top10 <- GO_enrich %>%
  filter(ONTOLOGY == "CC") %>%
  slice_head(n = 10)

# 筛选MF的前10个结果
mf_top10 <- GO_enrich %>%
  filter(ONTOLOGY == "MF") %>%
  slice_head(n = 10)

combined_top10 <- bind_rows(bp_top10, cc_top10, mf_top10)

#纵向柱状图#
p<-ggplot(combined_top10,
          aes(x=term,y=Count, fill=ONTOLOGY)) + #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666") ) +  #柱状图填充颜色
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  coord_flip() +  #让柱状图变为纵向
  xlab("GO term") +  #x轴标签
  ylab("Gene_Number") +  #y轴标签
  labs(title = "GO Terms Enrich")+  #设置标题
  theme_bw()
p
ggsave("GO_barplot_aplc.pdf",p,width= 10 , height= 8 , units="in")

#横向柱状图#
ggplot(combined_top10,
       aes(x=term,y=Count, fill=ONTOLOGY)) +  #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666") ) + #柱状图填充颜色
  facet_grid(.~ONTOLOGY, scale = 'free_x', space = 'free_x')+
  xlab("GO term") + #x轴标签
  ylab("Gene_Number") +  #y轴标签
  labs(title = "GO Terms Enrich")+ #设置标题
  theme_bw() +
  theme(axis.text.x=element_text(family="sans",face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 )) #对字体样式、颜色、还有横坐标角度（）

#气泡图#
p2<-ggplot(combined_top10,
           aes(y=term,x=Count))+
  geom_point(aes(size=Count,color=p.adjust))+
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"),
       x="Gene Ratio",y="GO term",title="GO Enrichment")+
  theme_bw()
p2
ggsave("GO_bubble_plot_aplc.pdf",p2,width= 10 , height= 6 , units="in")

#KEGG可视化
###柱状图
hh <- as.data.frame(KEGG_enrich)#保存结果
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
keg<-hh[(hh$pvalue<0.05),]

kp1<-ggplot(keg,aes(y=order,x=Count,fill=p.adjust))+
  geom_bar(stat = "identity",width=0.7)+####柱子宽度
  #coord_flip()+##颠倒横纵轴
  scale_fill_gradient(low = "red",high ="blue" )+#颜色自己可以换
  labs(title = "KEGG Pathways Enrichment",
       x = "Gene numbers",
       y = "Pathways")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()
kp1
ggsave("KEGG_bar_plot_aplc.pdf",kp1,width= 10 , height= 18 , units="in")

###气泡图
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
keg<-hh[(hh$pvalue<0.05),]
kp2<-ggplot(keg,aes(y=order,x=Count))+
  geom_point(aes(size=Count,color=-1*p.adjust))+# 修改点的大小
  scale_color_gradient(low="green",high = "red")+
  labs(color=expression(p.adjust,size="Count"),
       x="Gene Number",y="Pathways",title="KEGG Pathway Enrichment")+
  theme_bw()
kp2
ggsave("KEGG_bubble_plot_aplc.pdf",kp2,width= 10 , height= 18 , units="in")







##取交集----
ppitable<-read.table("AP_LC Target.txt",
                     header = F, sep = '\t', check.names = F, row.names = 1)
tar<-rownames(ppitable)

library(VennDiagram)

diffdata<-read.table("TCGA_LUAD_diff_count.txt",header = T,sep = '\t',check.names = F,row.names = 1)
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
             filename='Venn_diffvsppi.png',# 文件保存
             imagetype="png",  # 类型（tiff png svg）
             resolution = 400,  # 分辨率
             compression = "lzw"# 压缩算法
)

selected_gene<-intersect(diff,tar)
write.table(selected_gene,file = "gene.txt",sep = "\t",quote = F, row.names = T)

