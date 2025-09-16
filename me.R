##甲基化----
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


