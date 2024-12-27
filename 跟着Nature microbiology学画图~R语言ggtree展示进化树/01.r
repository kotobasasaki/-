library(ggtree)
library(treeio)
library(ggh4x)
tree<-read.tree("figure1b.txt")
tree
tree$root.edge<-0.1
p = ggtree(tree)+
  geom_tiplab(fontface="italic")+
  theme_tree2()+
  geom_rootedge()
p

p+
  geom_nodelab(aes(label=node))
p1 = ggtree::rotate(p,node = 15)

p1+
  scale_x_continuous(breaks =c (seq(0,3.2,by=0.45),3.5),
                     limits = c(-0.3,5.5),
                     labels = c(8:0))+
  guides(x=guide_axis_truncated(trunc_lower = 0,
                                trunc_upper = 3.5))+
  labs(x="Age (million years)")+
  theme(axis.title.x = element_text(hjust=0.3))+
  #geom_nodelab(aes(label=node))+
  geom_cladelab(node=18,label="I",textcolor="red",
                offset=1)+
  geom_cladelab(node=20,label="II",textcolor="red",
                offset=1)+
  geom_cladelab(node=23,label="III",textcolor="red",
                offset=1)+
  geom_cladelab(node=24,label="IV",textcolor="red",
                offset=1.7) -> pp

ggplot_build(pp)$data[[5]]
a = pp+
  annotate(geom = "segment",x=1.8,y=1.5,xend = 2.6,yend=3,
           arrow=arrow(angle = 20,type="closed",length = unit(3,'mm')))+
  annotate(geom="label",x=1.8,y=1.5,label="Mean wild-cultivated\ntomato divergence\n0.71 Ma",hjust=0.5)+
  annotate(geom = "segment",x=1,y=4,xend = 1.7,yend=5.2,
           arrow=arrow(angle = 20,type="closed",length = unit(3,'mm')))+
  annotate(geom="label",x=1,y=4,label="Mean divergence\n1.73 Ma",hjust=0.5)+
  annotate(geom = "segment",x=-0.1,y=10,xend = -0.05,yend=12.2,
           arrow=arrow(angle = 20,type="closed",length = unit(3,'mm')))+
  annotate(geom="label",x=-0.1,y=10,label="Fix at\n7.5 Ma",hjust=0.5)


library(readxl)
library(tidyverse)
fig1a.df<-read_excel("41588_2023_1340_MOESM5_ESM.xlsx",
                     sheet = "Fig1",
                     skip = 1)

# 解决 Species 列重复问题，确保唯一性
fig1a.df <- fig1a.df %>%
  mutate(Species = make.unique(as.character(Species)))  # 添加唯一标识防止重复

# 调整 Species 为因子类型，并逆序排列
fig1a.df <- fig1a.df %>%
  mutate(Species = factor(Species, levels = unique(rev(Species))))

# 数据长格式转换
fig1a.df <- fig1a.df %>%
  pivot_longer(
    cols = -Species, 
    names_to = "name", 
    values_to = "value"
  ) %>%
  mutate(
    name = factor(name, levels = c(
      "Non repetitive", "Unknown repeats", "Other repeats", 
      "DNA transpons", "Other retranspons", "ClassI/LTR", 
      "ClassI/LTR/Gypsy", "ClassI/LTR/Copia"
    ))
  )

p.stacked = ggplot(fig1a.df, aes(x = value, y = Species)) +
  geom_bar(stat = "identity", aes(fill = name)) +  # 堆叠条形图
  theme_bw() +  # 设置主题
  theme(
    legend.position = c(0.8, 0.2),             # 设置图例位置
    panel.border = element_blank(),           # 移除面板边框
    axis.text.x = element_text(),             # x轴文本格式
    axis.line.x = element_line(),             # x轴线条
    axis.text.y = element_blank(),            # 隐藏y轴文本
    axis.ticks.y = element_blank(),           # 隐藏y轴刻度
    panel.grid = element_blank(),             # 移除网格
    legend.title = element_blank()            # 隐藏图例标题
  ) +
  scale_x_continuous(
    limits = c(0, 2000000000),                # x轴范围
    labels = function(x) x / 1000000,         # 将值缩放为百万单位
    breaks = seq(0, 1500000000, by = 300000000)  # 设置刻度
  ) +
  guides(
    x = guide_axis_truncated(trunc_lower = 0, trunc_upper = 1500000000)  # 截断x轴
  ) +
  labs(x = NULL, y = NULL) +                  # 去除轴标签
  scale_fill_manual(
    values = c(
      "#34a56c", "#26c2e2", "#37519e", "#3fbc9c",
      "#eedf42", "#907135", "#f05435", "#67696b"
    ),
    limits = rev(c(
      "Non repetitive", "Unknown repeats", "Other repeats", 
      "DNA transpons", "Other retranspons", "ClassI/LTR", 
      "ClassI/LTR/Gypsy", "ClassI/LTR/Copia"
    )),
    labels = rev(c(
      "Non repetitive", "Unknown repeats", "Other repeats", 
      "DNA transpons", "Other retranspons", "LTR/Others", 
      "LTR/Gypsy", "LTR/Copia"
    ))
  )

