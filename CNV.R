

# 加载必要的R包
library(tidyverse)
library(ggplot2)
library(ggpubr)

# 设置数据目录（根据实际文件夹路径修改）
data_dir <- "648367"  # GISTIC2输出数据所在文件夹

# 1. 读取CNV数据 (GISTIC2阈值处理后的数据)
# 文件名为 all_thresholded.by_genes.txt
cnv_file <- file.path(data_dir, "all_thresholded.by_genes.txt")
cnv_data <- read.delim(cnv_file, row.names = 1, check.names = FALSE)

# 提取目标基因的CNV数据
target_genes <- c("CYP17A1", "PTGES", "CDK1", "KIF11", "F2", "HNF4A")
cnv_target <- cnv_data[rownames(cnv_data) %in% target_genes, ]

# 转置数据：行为样本，列为基因
cnv_df <- as.data.frame(t(cnv_target)) %>%
  rownames_to_column(var = "sample_id") %>%
  # 简化样本ID为前15字符（TCGA标准）
  mutate(sample_id = substr(sample_id, 1, 15))

# 2. 读取表达数据
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
# 提取目标基因的表达数据并转置
expr_target <- expr_data[rownames(expr_data) %in% target_genes, ]
expr_df <- as.data.frame(t(expr_target)) %>%
  rownames_to_column(var = "sample_id") %>%
  mutate(sample_id = substr(sample_id, 1, 15))

# 3. 合并CNV和表达数据
combined_data <- inner_join(cnv_df, expr_df, by = "sample_id", suffix = c("_cnv", "_expr"))



# 4. 计算每个基因的CNV与表达的相关性
results <- list()
for (gene in target_genes) {
  # 提取当前基因的数据
  df_gene <- combined_data %>%
    select(sample_id, cnv = all_of(paste0(gene, "_cnv")), expr = all_of(paste0(gene, "_expr")))

  # 将 cnv 和 expr 转换为数值型
  df_gene$cnv <- as.numeric(as.character(df_gene$cnv))
  df_gene$expr <- as.numeric(as.character(df_gene$expr))

  # 计算Spearman相关系数
  cor_test <- cor.test(df_gene$cnv, df_gene$expr, method = "spearman")

  # 存储结果
  results[[gene]] <- data.frame(
    Gene = gene,
    Rho = cor_test$estimate,
    P_value = cor_test$p.value
  )

  # 绘制散点图
  p <- ggplot(df_gene, aes(x = cnv, y = log2(expr + 1))) +
    geom_point(alpha = 0.6, color = "blue") +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(
      title = paste("CNV vs Expression:", gene),
      x = "CNV Status (GISTIC2)",
      y = "log2(Expression + 1)"
    ) +
    theme_bw() +
    annotate("text", x = -1, y = max(log2(df_gene$expr + 1)),
             label = paste0("Rho = ", round(cor_test$estimate, 3), "\nP = ", format.pval(cor_test$p.value, digits = 2)))

  print(p)  # 显示图形
  ggsave(paste0(gene, "_cnv_vs_expr.pdf"), plot = p, width = 9, height = 8)
}

# 5. 汇总相关性结果
results_df <- do.call(rbind, results)
print(results_df)

# 可选：保存结果到CSV
write.csv(results_df, "CNV_Expression_Correlation_Results.csv", row.names = FALSE)
