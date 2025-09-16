library(tidyverse)
library(sf)
library(ggspatial)
library(RColorBrewer)
library(sjmisc)
library(ggplot2)
library(ggpubr)

china <- st_read("china_country.shp")
china_prov <- st_read("china.shp")

china_line = st_read("china_nine_dotted_line.shp")

region_code<- data.frame(
  FCNAME = c("北京市", "天津市", "河北省", "山西省", "内蒙古自治区", "黑龙江省", "吉林省", "辽宁省",
           "上海市", "江苏省", "浙江省", "山东省", "安徽省","福建省","江西省","河南省","湖北省","湖南省",
           "广东省","广西壮族自治区","海南省","重庆市","四川省","贵州省","云南省","西藏自治区","陕西省",
           "甘肃省","青海省","宁夏回族自治区","新疆维吾尔自治区","台湾省","香港特别行政区","澳门特别行政区"),
  LUCA_rate = c(35.0, 35.0, 35.0, 35.0, 35.0, 41.6, 41.6, 41.6, 35.3, 35.3, 35.3, 35.3, 35.3, 35.3, 35.3, 36.5, 36.5, 36.5,
                34.9, 34.9, 34.9, 38.6, 38.6, 38.6, 38.6, 38.6, 28.9, 28.9, 28.9, 28.9, 28.9, 35.3, 34.9, 34.9),
  region = c("华北", "华北", "华北", "华北", "华北", "东北", "东北", "东北",
             "华东", "华东", "华东", "华东", "华东", "华东", "华东", "华中", "华中", "华中",
             "华南", "华南", "华南", "西南", "西南", "西南", "西南", "西南", "西北", "西北",
             "西北", "西北", "西北", "华东", "华南", "华南"),
  region_en = c("North China","North China","North China","North China","North China","Northeast China","Northeast China","Northeast China",
                "East China","East China","East China","East China","East China","East China","East China","Central China","Central China","Central China",
                "South China","South China","South China","Southwest China","Southwest China","Southwest China","Southwest China","Southwest China",
                "Northwest China","Northwest China","Northwest China","Northwest China","Northwest China","East China",
                "South China","South China")
)

# 合并GDP数据
map_data1 <- left_join(china_prov, region_code, by = "FCNAME")


# 合并并验证拓扑
map_data1 <- left_join(china_prov, region_code, by = "FCNAME") %>%
  st_make_valid() %>%
  filter(!is.na(region))  # 移除未定义区域

# 计算区域指标（无几何操作）
region_avg <- map_data1 %>%
  st_drop_geometry() %>%
  group_by(region_en) %>%
  summarise(LUCA_rate = mean(LUCA_rate))

# 单独处理区域几何
region_geom <- map_data1 %>%
  group_by(region_en) %>%
  summarise(geometry = st_union(geometry)) %>%
  st_make_valid()  # 关键修复

# 计算中心点
region_centers <- region_geom %>%
  mutate(centroid = st_centroid(geometry),
         label_pos = st_coordinates(centroid)) %>%
  st_drop_geometry() %>%
  select(region_en, label_pos)

# 创建最终数据集
region_data <- region_geom %>%
  left_join(region_avg, by = "region_en") %>%
  left_join(region_centers, by = "region_en")



# 创建南海小地图
south_china_sea <- ggplot() +
  geom_sf(data = china, fill = "white", color = "gray60", size = 0.3) +
  geom_sf(data = china_line, color = "red", size = 0.5) +
  coord_sf(
    xlim = c(105, 125),   # 东经105-125度
    ylim = c(3, 25),      # 北纬3-25度
    expand = FALSE,
    datum = sf::st_crs(4326)  # 使用WGS84坐标系
  ) +
  theme_void() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    panel.background = element_rect(fill = "white")
  )

# 绘制主地图并添加小地图
main_plot <- ggplot() +
  # 主地图 - 使用自定义连续色标
  geom_sf(data = region_data, aes(fill = LUCA_rate)) +
  geom_text(data = region_data,
            aes(x = label_pos[,1], y = label_pos[,2],
                label = paste0(region_en, "\n", round(LUCA_rate,1), "%")),
            size = 3) +
  # 自定义连续色标，在关键值附近创建明显断点
  scale_fill_gradientn(
    name = "Prevalence(%)",
    colors = c("#F7F7F7", "#D1E5F0", "#92C5DE", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"),
    values = scales::rescale(c(28.9, 34.9, 35.0, 35.3, 36.5, 38.6, 41.6)),
    breaks = c(28.9, 35.0, 36.5, 38.6, 41.6),
    guide = guide_colorbar(
      barheight = unit(3, "cm"),
      barwidth = unit(0.5, "cm"))
  ) +

  # 添加国界线
  geom_sf(data = china, fill = NA, color = "black", size = 0.5) +

  # 添加小地图
  annotation_custom(
    grob = ggplotGrob(south_china_sea),
    xmin = 125,   # 小地图左下角经度
    xmax = 137.5,   # 小地图右上角经度
    ymin = 15,     # 小地图左下角纬度
    ymax = 5     # 小地图右上角纬度
  ) +

  # 添加比例尺和指北针
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering
  ) +

  # 设置主题和标题
  labs(title = "Distribution Map of Lung Cancer Prevalence in China") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = c(0.95, 0.3),  # 将图例移到右上角
    legend.justification = c(1, 0)
  )

# 保存为PDF
ggsave(
  filename = "China_Lung_Cancer_Prevalence0624_2.pdf",  # 文件名
  plot = main_plot,                             # 要保存的图形对象
  device = "pdf",                               # 输出格式
  width = 10.5,                                   # 宽度（英寸）
  height = 10.5,                                   # 高度（英寸）
  units = "in",                                 # 尺寸单位
  dpi = 300                                     # 分辨率（适用于嵌入的栅格元素）
)

