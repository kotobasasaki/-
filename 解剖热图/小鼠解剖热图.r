library(gganatogram)
library(dplyr)
library(viridis)
library(gridExtra)
library(patchwork)
library(openxlsx)
library(ggplot2)


data = readxl::read_excel("./mmFemale_key.xlsx",sheet = "value")
lung_data <- data %>% filter(organ %in% c("lung", "trachea"))

gganatogram(data = lung_data, 
            outline = TRUE, 
            fillOutline = 'white', 
            organism = 'mouse', 
            sex = 'female', 
            fill = "value") +
  theme_void() +  # 使用空白主题
  scale_fill_gradient(low = "red", high = "pink") +  # 使用viridis配色方案
  labs(title = "小鼠呼吸系统器官重量") +  # 添加标题
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # 居中标题
  ) +
  annotate("text", 
           x = 1, y = -160,  # 设置文本位置，您可以根据需要调整这些值
           label = "单位：克（g）", 
           size = 5, 
           hjust = -5, vjust = -1, 
           fontface = "bold", 
           color = "black")  # 字体设置


