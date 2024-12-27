library(sf)
library(ggplot2)
library(dplyr)
heilongj_sf <- st_read("./佳木斯市/佳木斯市.shp")
data <- read.csv("佳市数据.csv", fileEncoding = "UTF-8")

heilongj_sf <- heilongj_sf %>%
  left_join(data, by = c("name" = "name"))

# 绘制热图并添加地级市名称
ggplot(data = heilongj_sf) +
  geom_sf(aes(fill = 数据), color = "black") +  # 根据产量填充颜色
  scale_fill_gradient(low = "blue", high = "red") +  # 蓝色到红色的渐变
  geom_sf_text(aes(label = name), size = 2.6, color = "black", fontface = "bold") +  # 添加地级市名称
  theme_minimal() +
  ggtitle("佳木斯数据展示图") +
  theme(legend.position = "right")
