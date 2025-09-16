library(ggplot2)
library(pheatmap)
library(gridExtra)
library(reshape2)
library(scales)

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
bar_data<-data.frame(Model=sorted_data$Model,Mean=Average)
bar_data$Mean<-round(Cindex_mean$Mean,digits = 3)

# 5. 创建热图矩阵
heat_data <- as.matrix(sorted_data[, -1])
rownames(heat_data) <- sorted_data$Model

# 确保顺序一致
heat_data <- heat_data[match(bar_data$Model, rownames(heat_data)), ]

# 创建列注释（关键修改部分）
annotation_col <- data.frame(
  Dataset = factor(colnames(heat_data),
                   levels = c("TCGA-LUAD-TRAIN", "TCGA-LUAD-VALI", "GSE50081"))
)
rownames(annotation_col) <- colnames(heat_data)

# 定义注释颜色
annotation_colors <- list(
  Dataset = c("TCGA-LUAD-TRAIN" = "#1B9E77",
              "TCGA-LUAD-VALI" = "#D95F02",
              "GSE50081" = "#7570B3")
)

# 绘制条形图（保持原有）
bar_plot <- ggplot(bar_data, aes(x = Mean, y = reorder(Model, Mean))) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_text(aes(label = format(Mean, nsmall = 3)),
            hjust = -0.1, size = 3, color = "black") +
  scale_x_continuous(limits = c(0, 0.8), expand = c(0, 0)) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(r = 0)) +
  labs(x = "Mean")

# 绘制带注释的热图（关键修改部分）
heatmap_plot <- pheatmap(heat_data,
                         color = colorRampPalette(c("#ffffbf","darkorange1","firebrick3"))(100),
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         display_numbers = TRUE,
                         number_format = "%.3f",
                         number_color = "black",
                         fontsize_row = 8,
                         fontsize_col = 8,
                         angle_col = 45,
                         annotation_col = annotation_col,
                         annotation_colors = annotation_colors,
                         annotation_names_col = TRUE,
                         silent = TRUE)



# 组合图形
pdf("Cindex_of_all_models_LUAD_0409_2.pdf",width = 8,height = 14.5)

grid.arrange(
  bar_plot,
  heatmap_plot[[4]],
  ncol = 2,
  widths = c(1, 2.5),
  top = "Model C-index Performance Comparison"
)

dev.off()


load("TCGA_LUAD_0409.RData")


##ROC曲线----
library(timeROC)
library(survival)
library(ggplot2)
riskgroup<-res_luad[["riskscore"]][["Lasso + RSF"]][["TCGA_LUAD_TRAIN"]]
riskgroup$OS.time<-riskgroup$OS.time/30
#计算
time_ROC<-timeROC(T=riskgroup$OS.time,
                  delta=riskgroup$OS,
                  marker=riskgroup$RS,
                  cause=1,
                  weighting = "marginal",
                  times = c(12,36,60),
                  ROC=TRUE,
                  iid=TRUE
)

time_ROC

summary(time_ROC)
time_ROC$AUC
#列表以便作图
time_ROC.res<-data.frame(TP_1year=time_ROC$TP[,1],
                         FP_1year=time_ROC$FP[,1],
                         TP_3year=time_ROC$TP[,2],
                         FP_3year=time_ROC$FP[,2],
                         TP_5year=time_ROC$TP[,3],
                         FP_5year=time_ROC$FP[,3])
#作图
p<-ggplot()+
  geom_line(data=time_ROC.res,aes(x=FP_1year,y=TP_1year),size=.5,color="red")+
  geom_line(data=time_ROC.res,aes(x=FP_3year,y=TP_3year),size=.5,color="blue")+
  geom_line(data=time_ROC.res,aes(x=FP_5year,y=TP_5year),size=.5,color="purple")+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey",size=1, linetype = 2
  )+
  theme_bw()+
  annotate("text",x=0.75,y=0.30,size=5,
           label=paste0("AUC at 1 years = ",round(time_ROC$AUC[[1]],3)),color="red")+
  annotate("text",x=0.75,y=0.20,size=5,
           label=paste0("AUC at 3 years = ",round(time_ROC$AUC[[2]],3)),color="blue")+
  annotate("text",x=0.75,y=0.10,size=5,
           label=paste0("AUC at 5 years = ",round(time_ROC$AUC[[3]],3)),color="purple")+
  labs(x="False positive rate",y="True positive rate")+
  theme(axis.text=element_text( size=12,  color="black"),
        axis.title=element_text( size=12, color="black"))
p
ggsave(filename = "ROC_LUAD80.pdf", plot = p, width = 10, height = 8, units = "in")


#LUAD-Vali
riskgroup1<-res_luad[["riskscore"]][["Lasso + RSF"]][["TCGA_LUAD_VALI"]]
riskgroup1$OS.time<-riskgroup1$OS.time/30
#计算
time_ROC<-timeROC(T=riskgroup1$OS.time,
                  delta=riskgroup1$OS,
                  marker=riskgroup1$RS,
                  cause=1,
                  weighting = "marginal",
                  times = c(12,36,60),
                  ROC=TRUE,
                  iid=TRUE
)

time_ROC

summary(time_ROC)
time_ROC$AUC
#列表以便作图
time_ROC.res<-data.frame(TP_1year=time_ROC$TP[,1],
                         FP_1year=time_ROC$FP[,1],
                         TP_3year=time_ROC$TP[,2],
                         FP_3year=time_ROC$FP[,2],
                         TP_5year=time_ROC$TP[,3],
                         FP_5year=time_ROC$FP[,3])
#作图
p<-ggplot()+
  geom_line(data=time_ROC.res,aes(x=FP_1year,y=TP_1year),size=.5,color="red")+
  geom_line(data=time_ROC.res,aes(x=FP_3year,y=TP_3year),size=.5,color="blue")+
  geom_line(data=time_ROC.res,aes(x=FP_5year,y=TP_5year),size=.5,color="purple")+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey",size=1, linetype = 2
  )+
  theme_bw()+
  annotate("text",x=0.75,y=0.30,size=5,
           label=paste0("AUC at 1 years = ",round(time_ROC$AUC[[1]],3)),color="red")+
  annotate("text",x=0.75,y=0.20,size=5,
           label=paste0("AUC at 3 years = ",round(time_ROC$AUC[[2]],3)),color="blue")+
  annotate("text",x=0.75,y=0.10,size=5,
           label=paste0("AUC at 5 years = ",round(time_ROC$AUC[[3]],3)),color="purple")+
  labs(x="False positive rate",y="True positive rate")+
  theme(axis.text=element_text( size=12,  color="black"),
        axis.title=element_text( size=12, color="black"))
p
ggsave(filename = "ROC_LUAD20_1.pdf", plot = p, width = 10, height = 8, units = "in")


#GSE50081
riskgroup2<-res_luad[["riskscore"]][["Lasso + RSF"]][["GSE50081"]]

#计算
time_ROC<-timeROC(T=riskgroup2$OS.time,
                  delta=riskgroup2$OS,
                  marker=riskgroup2$RS,
                  cause=1,
                  weighting = "marginal",
                  times = c(12,36,60),
                  ROC=TRUE,
                  iid=TRUE
)

time_ROC

summary(time_ROC)
time_ROC$AUC
#列表以便作图
time_ROC.res<-data.frame(TP_1year=time_ROC$TP[,1],
                         FP_1year=time_ROC$FP[,1],
                         TP_3year=time_ROC$TP[,2],
                         FP_3year=time_ROC$FP[,2],
                         TP_5year=time_ROC$TP[,3],
                         FP_5year=time_ROC$FP[,3])
#作图
p<-ggplot()+
  geom_line(data=time_ROC.res,aes(x=FP_1year,y=TP_1year),size=.5,color="red")+
  geom_line(data=time_ROC.res,aes(x=FP_3year,y=TP_3year),size=.5,color="blue")+
  geom_line(data=time_ROC.res,aes(x=FP_5year,y=TP_5year),size=.5,color="purple")+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey",size=1, linetype = 2
  )+
  theme_bw()+
  annotate("text",x=0.75,y=0.30,size=5,
           label=paste0("AUC at 1 years = ",round(time_ROC$AUC[[1]],3)),color="red")+
  annotate("text",x=0.75,y=0.20,size=5,
           label=paste0("AUC at 3 years = ",round(time_ROC$AUC[[2]],3)),color="blue")+
  annotate("text",x=0.75,y=0.10,size=5,
           label=paste0("AUC at 5 years = ",round(time_ROC$AUC[[3]],3)),color="purple")+
  labs(x="False positive rate",y="True positive rate")+
  theme(axis.text=element_text( size=12,  color="black"),
        axis.title=element_text( size=12, color="black"))
p
ggsave(filename = "ROC_GSE50081.pdf", plot = p, width = 10, height = 8, units = "in")


