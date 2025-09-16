library(TCGAbiolinks)
library(maftools)
library(SummarizedExperiment)
library(ggplot2)
library(ggpubr)
library(TCGAmutations)

#突变----
#导入突变数据
maf <- TCGAmutations::tcga_load(study = "LUAD")
genes<-c("CYP17A1", "PTGES", "CDK1", "KIF11", "F2", "HNF4A")

oncoplot(
  maf = maf,
  genes = genes,
  top = 20,  # 显示前10个高频突变基因
  showTumorSampleBarcodes = FALSE
)

# 打开 PDF 设备，指定文件名和尺寸
pdf("LUAD_Oncoplot.pdf", width = 10, height = 7)  # 可调整宽高适应图形

# 绘制 oncoplot
oncoplot(
  maf = maf,
  genes = genes,
  top = 20,
  showTumorSampleBarcodes = FALSE
)

# 关闭图形设备（保存文件）
dev.off()


# 统计突变频率
gene_mut_freq <- subsetMaf(
  maf,
  genes = genes
)@gene.summary[, .(Hugo_Symbol, MutatedSamples)]

# 可视化突变频率
p <- ggplot(gene_mut_freq, aes(x = Hugo_Symbol, y = MutatedSamples)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Mutation Frequency", x = "Gene", y = "Mutated Samples") +
  theme_minimal()


##拷贝数变异(CNV)----
library(TCGAbiolinks)
f = "masked.Rdata"
if(!file.exists(f)){
  query <- GDCquery(project = "TCGA-LUAD",
                    data.category = "Copy Number Variation",
                    data.type = "Masked Copy Number Segment")

  # 设置分块下载，每块100个文件
  GDCdownload(query, files.per.chunk = 100)

  seg <- GDCprepare(query)
  colnames(seg)
  seg = seg[,c(7,2:6)]
  save(seg,file = f)
  write.table(seg,file = "seg.txt",row.names = F,quote = F,sep = "\t")
}

load(f)
head(seg)
length(unique(seg$Sample)) #人数

marker = data.table::fread("snp6.na35.remap.hg38.subset.txt.gz",data.table = F)
head(marker)
table(marker$freqcnv)

marker = marker[!marker$freqcnv,] #只保留FALSE
marker =marker[,1:3]
colnames(marker) = c("Marker Name","Chromosome","Marker Position")
head(marker)

write.table(marker,file = "marker.txt",row.names = F,quote = F,sep = "\t")

#从GISTIC2.0获取结果后继续分析
library(maftools)
library(tidyverse)

laml.gistic = readGistic(gisticAllLesionsFile = "648367/all_lesions.conf_99.txt",
                         gisticAmpGenesFile = "648367/amp_genes.conf_99.txt",
                         gisticDelGenesFile = "648367/del_genes.conf_99.txt",
                         gisticScoresFile = "648367/scores.gistic",
                         isTCGA = TRUE)

head(laml.gistic@data)

gisticChromPlot(gistic = laml.gistic, ref.build = 'hg38',markBands = "all")
head(laml.gistic@gis.scores)
gisticBubblePlot(gistic = laml.gistic)

#找出感兴趣的区域里面的基因
tmp = as.data.frame(laml.gistic@data[,c(1,5)])
library(tidyverse)
tmp = separate(tmp,col = Cytoband,into = c("cnv","cyto"),sep = ":")
head(tmp)
table(tmp$cnv)

dir("648367/",pattern = "genes.txt$")
all_data = read.delim("648367/all_data_by_genes.txt",row.names = 1,check.names = F)
all_data[1:4,1:5]

all_thr = read.delim("648367/all_thresholded.by_genes.txt",row.names = 1,check.names = F)
all_thr[1:4,1:5]

#差异分析
if(!require(GeoTcgaData))BiocManager::install("GeoTcgaData")
library(GeoTcgaData)
cnv = all_thr[,-(1:2)]
rownames(cnv) = 1:nrow(cnv)
cnv[cnv >= 1] = 1
cnv[cnv <= -1] = -1
library(tinyarray)
group = make_tcga_group(cnv)
diffcnv = differential_CNV(cnv,group)

diffcnv$gene = rownames(all_thr)
rownames(cnv) = rownames(all_thr)
head(diffcnv)

#针对目标基因
my_gene = c("CYP17A1", "PTGES", "CDK1", "KIF11", "F2", "HNF4A")
cnv = cnv[my_gene,group=="tumor"]

library(tidyverse)
cnv_df = as.data.frame(t(cnv)) %>% rownames_to_column("Sample")
cnv_long <- gather(cnv_df, key = "Gene", value = "CNV", 2:7)

# 统计每个基因的拷贝数变化样本比例
cnv_summary <- cnv_long %>%
  filter(CNV != 0) %>%  # 排除正常的样本
  group_by(Gene, CNV) %>%
  summarise(Freq = n()/ncol(cnv) *100) %>%
  ungroup() %>%
  mutate(CNV = ifelse(CNV>0,"Amp","Del"))
# 控制基因顺序
tmp = cnv_summary %>%
  group_by(Gene) %>%
  summarise(max = max(Freq)) %>%
  arrange(desc(max))
cnv_summary$Gene = factor(cnv_summary$Gene,levels = tmp$Gene)
# 使用 ggplot2 绘制点图
library(ggplot2)
p <- ggplot(cnv_summary, aes(x = Gene, y = Freq)) +
  geom_segment( aes(x=Gene, xend=Gene, y=0, yend=Freq),color = "grey",linewidth = 3) +
  geom_point(aes(color = CNV),size = 5) +
  scale_color_manual(values = c("Del" = "#2fa1dd", "Amp" = "#f87669")) +
  theme_minimal()
ggsave("CNVpoint.pdf", plot = p, width = 11.69, height = 8.27)



##甲基化1----
library(wateRmelon)
library(minfi)
library(methylationArrayAnalysis)
library(ChAMP)
library(data.table)
library(missMethyl)
# 读取 tsv.gz 文件

pd.all <- data.table::fread("TCGA-LUAD.methylation450.tsv.gz", sep = "\t")
dim(pd.all)
#读取样本信息
sample <- read.table("TCGA.LUAD.sampleMap_LUAD_clinicalMatrix",
                    header = T, sep = '\t', check.names = F, row.names = 1)

clinical <- read.table("TCGA-LUAD.survival.tsv",
                       header = T, sep = '\t', check.names = F, row.names = 1)

#格式化并合并
pd <- sample[,c("bcr_sample_barcode","sample_type")]
pd[pd ==""] <- NA
pd1 <- na.omit(pd)

pd.all <- as.data.frame(pd.all)
rownames(pd.all) <- pd.all$`Composite Element REF`
pd.matrix <- pd.all[,c(2:ncol(pd.all))]

pd.matrix <- t(pd.matrix)
pd.matrix <-as.data.frame(pd.matrix)

pd.matrix$ID <- rownames(pd.matrix)
colnames(pd1)[1] <- "ID"

matrix <- merge(pd.matrix, pd1, by= "ID")
rownames(matrix) <- matrix$ID
matrix1 <- matrix
matrix1[is.na(matrix1)] <- 0

#保留Normal-Tumor信息
library(dplyr)
keep_samples <- pd1[pd1$sample_type %in% c("Primary Tumor", "Solid Tissue Normal"),]
matrix2 <- matrix1[keep_samples$ID,]
matrix2 <- na.omit(matrix2)
#注释文件
colnames(keep_samples) <- c("SampleName","SampleGroup")
pd2 <- keep_samples[keep_samples$SampleName %in% rownames(matrix2),]
#甲基化文件
matrix2$sample_type <- NULL
matrix2$ID <- NULL
matrix3 <- t(matrix2)

write.table(matrix3, file = "matrix3.txt", sep = "\t", row.names = TRUE)
write.csv(pd2, file = "pd2.csv")

matrix3 <- read.table("matrix3.txt",
           header = T, sep = '\t', check.names = F, row.names = 1)
matrix3 <- as.matrix(matrix3)
pd2 <- read.csv("pd2.csv",header = T)
rownames(pd2) <- pd2$X
pd2$X <- NULL
colnames(pd2) <- c("Sample_Name","Sample_Group")


# 载入为ChAMP对象
myLoad=champ.filter(beta = matrix3 ,pd = pd2)
champ.QC(beta = myLoad$beta,pheno = myLoad$pd$Sample_Group)
#归一化
myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=10)

#主成分分析
library(FactoMineR)
dat <- t(myNorm)
pd=myLoad$pd
group_list=pd$Sample_Group
table(group_list)

dat.pca <- PCA(dat, graph = FALSE,ncp = 10)

# 加载必要的包
if (!require(factoextra)) install.packages("factoextra")
library(factoextra)

# 绘制 PCA 图
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", # 注释掉此行代码以查看离群点
                         col.ind = group_list,
                         addEllipses = TRUE,
                         legend.title = "Groups",
                         # ellipse.type = "confidence" # 置信椭圆
)

# 保存图形
ggsave("pca.pdf", plot = pca_plot, width = 11.69, height = 8.27)

#热图聚类
cg=names(tail(sort(apply(myNorm,1,sd)),1000))

ac=data.frame(group=group_list)
rownames(ac)=colnames(myNorm)

library(pheatmap)

# 打开 PDF 图形设备
pdf("pca_heatmap.pdf", width = 11.69, height = 8.27)

# 绘制热图
pheatmap(myNorm[cg, ],
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = ac)

# 关闭图形设备
dev.off()


##甲基化2----
library(wateRmelon)
library(minfi)
library(methylationArrayAnalysis)
library(ChAMP)
library(data.table)
library(missMethyl)


# 加载注释数据
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# 定义目标基因列表
target_genes <- c("CYP17A1", "PTGES", "CDK1", "KIF11", "F2", "HNF4A")

# 提取每个基因相关的探针
gene_probes <- lapply(target_genes, function(gene) {
  # 获取基因相关的所有探针
  probes <- rownames(anno)[grep(gene, anno$UCSC_RefGene_Name)]

  # 添加位置信息
  probe_info <- data.frame(
    Probe = probes,
    Gene = gene,
    Relation = anno[probes, "UCSC_RefGene_Group"],
    Chromosome = anno[probes, "chr"],
    Position = anno[probes, "pos"],
    stringsAsFactors = FALSE
  )

  # 按位置关系排序（TSS区域优先）
  relation_order <- c("TSS200", "TSS1500", "5'UTR", "1stExon", "Body", "3'UTR")
  probe_info$Relation <- factor(probe_info$Relation, levels = relation_order)
  probe_info <- probe_info[order(probe_info$Relation), ]

  return(probe_info)
})

# 合并所有基因的探针信息
all_probes <- do.call(rbind, gene_probes)

# 去除重复探针（某些探针可能对应多个基因）
all_probes <- all_probes[!duplicated(all_probes$Probe), ]



# 提取目标探针的甲基化数据
# 过滤掉不在 myNorm 行名中的探针
valid_probes <- all_probes$Probe %in% rownames(myNorm)
target_beta <- myNorm[all_probes$Probe[valid_probes], ]

# 转置数据并添加探针信息
plot_data <- as.data.frame(t(target_beta))
plot_data$Sample <- rownames(plot_data)

# 添加分组信息
plot_data$Group <- pd$Sample_Group[match(plot_data$Sample, pd$Sample_Name)]

# 转换为长格式
library(tidyr)
plot_long <- pivot_longer(
  plot_data,
  cols = -c(Sample, Group),
  names_to = "Probe",
  values_to = "Beta"
)

# 添加基因和位置信息
plot_long <- dplyr::left_join(plot_long, all_probes, by = "Probe")

# 计算每个样本每个基因的平均甲基化水平
gene_avg <- plot_long %>%
  dplyr::group_by(Sample, Gene, Group) %>%
  dplyr::summarise(Avg_Beta = mean(Beta, na.rm = TRUE)) %>%
  dplyr::ungroup()


# 对每个基因进行统计检验
stat_results <- gene_avg %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    p_value = t.test(Avg_Beta ~ Group)$p.value,
    delta_beta = mean(Avg_Beta[Group == "Tumor"]) - mean(Avg_Beta[Group == "Normal"]),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    adj_p_value = p.adjust(p_value, method = "BH"),
    Significance = dplyr::case_when(
      adj_p_value < 0.001 ~ "***",
      adj_p_value < 0.01 ~ "**",
      adj_p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# 查看统计结果
print(stat_results)

# 对每个探针进行统计检验
probe_stat <- plot_long %>%
  dplyr::group_by(Probe, Gene, Relation) %>%
  dplyr::summarise(
    # 检查每个分组中的数据是否满足条件
    valid_tumor = var(Beta[Group == "Tumor"], na.rm = TRUE) > 0 & length(Beta[Group == "Tumor"]) > 1,
    valid_normal = var(Beta[Group == "Normal"], na.rm = TRUE) > 0 & length(Beta[Group == "Normal"]) > 1,
    p_value = if (valid_tumor & valid_normal) {
      t.test(Beta ~ Group)$p.value
    } else {
      NA  # 如果不满足条件，返回 NA
    },
    delta_beta = mean(Beta[Group == "Tumor"], na.rm = TRUE) - mean(Beta[Group == "Normal"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    adj_p_value = p.adjust(p_value, method = "BH"),
    Significance = dplyr::case_when(
      adj_p_value < 0.001 ~ "***",
      adj_p_value < 0.01 ~ "**",
      adj_p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

library(ggplot2)
library(ggpubr)

# 合并统计结果到绘图数据
gene_plot_data <- dplyr::left_join(gene_avg, stat_results, by = "Gene")


# 绘制箱线图并添加显著性标记
pl<-ggplot(gene_plot_data, aes(x = Gene, y = Avg_Beta)) +
  geom_boxplot(aes(fill = Group), position = position_dodge(width = 0.8), alpha = 0.8, outlier.shape = NA) +
  geom_point(aes(color = Group), position = position_jitterdodge(jitter.width = 0.2),
             size = 1.5, alpha = 0.6) +
  geom_text(data = dplyr::distinct(gene_plot_data, Gene, Significance),
            aes(x = Gene, y = max(gene_plot_data$Avg_Beta) + 0.1, label = Significance),  # 使用固定值调整 y 坐标
            vjust = -1,  # 调整标签的垂直位置
            size = 5) +
  scale_fill_manual(values = c("Solid Tissue Normal" = "#1b9e77", "Primary Tumor" = "#d95f02")) +
  scale_color_manual(values = c("Solid Tissue Normal" = "#1b9e77", "Primary Tumor" = "#d95f02")) +
  labs(title = "Gene-level Methylation Comparison",
       subtitle = "LUAD Tumor vs Normal Tissue",
       y = "Average Beta Value", x = "Gene") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        legend.position = "top") +
  ylim(0, max(gene_plot_data$Avg_Beta) + 0.1)

ggsave("Gene-level Methylation Comparison.pdf",plot = pl,width = 11.69, height = 8.27)





