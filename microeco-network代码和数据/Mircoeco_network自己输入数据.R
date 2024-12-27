# 清理工作环境中的所有对象
rm(list = ls())

#安装microeco包
if(!require("microeco")){install.packages("microeco")}
library(microeco)  
#安装并加载ape包
if(!require("ape")){install.packages("ape")}
library(ape)

##################数据基本操作#######################
# 加载示例数据；16S rRNA基因扩增子测序数据
#读取样本分组信息表
sample_info_16S <- read.csv("sample_info_16S.csv")
#读取特征表
otu_table_16S <- read.csv("otu_table_16S.csv",row.names =1)
#读取系统发育树
phylo_tree_16S <- read.tree("phylo_tree_16S.tre", tree.names = NULL)

###############创建微表microtable
# 如果你只提供丰度表，该类可以帮助你创建一个样本信息表
mt <- microtable$new(otu_table = otu_table_16S)
# 通常情况下添加元数据
mt <- microtable$new(otu_table = otu_table_16S, sample_table = sample_info_16S)
mt
# 让我们创建一个包含更多信息的 microtable 对象
mt <- microtable$new(sample_table = sample_info_16S, otu_table = otu_table_16S, tax_table = taxonomy_table_16S, phylo_tree = phylo_tree_16S)
mt

##################microeco网络分析###################

#计算相关性矩阵
# 参数 cor_method 用于选择相关性计算方法，默认使用 Pearson 或 Spearman 相关性，

# 调用 R 基础的 cor.test 函数，这个函数执行速度较慢(可选)
#t1 <- trans_network$new(
#  dataset = mt,         # 使用的输入数据集
#  cor_method = "spearman",  # 选择 Spearman 相关性方法
#  filter_thres = 0.001     # 设置相关性过滤阈值
#)

#使用WGCNA包进行计算，该包提升了计算效率，适合大规模数据
#加载并安装WGCNA包
if(!require("WGCNA")) install.packages("WGCNA", repos = BiocManager::repositories())
#使用WGCNA进行相关性矩阵计算
t1 <- trans_network$new(
  dataset = mt,                      # 输入数据集 mt
  cor_method = "spearman",            # 选择 Spearman 相关性计算方法
  use_WGCNA_pearson_spearman = TRUE,  # 使用 WGCNA 包来计算相关性
  filter_thres = 0.0001               # 设置过滤阈值为 0.0001
)

# 构建相关性网络
#安装并加载igraph包
if(!require("igraph")){install.packages("igraph")}
# 使用 p 值阈值和网络优化进行网络构建
#t1$cal_network(
#  COR_p_thres = 0.01,      # 设定相关性 p 值的阈值，只有 p 值小于 0.01 的相关性才会用于构建网络
#  COR_optimization = TRUE  # 启用网络优化，使用算法构建更稀疏的网络，提高网络可解释性
#)

# 使用自定义的相关系数阈值来构建网络(可选)
t1$cal_network(
  COR_p_thres = 0.01,  # 设定相关性 p 值的阈值，过滤掉不显著的相关性
  COR_cut = 0.7        # 使用相关系数阈值 0.7，只有相关系数大于 0.7 的变量对会被保留在网络中
)



# 对无向网络应用 igraph 的 cluster_fast_greedy 函数进行模块划分
t1$cal_module(method = "cluster_fast_greedy")  # 调用模块划分方法

#将网络保存为gexf文件，用于Gephi软件进行可视化
#安装并加载rgexf包,xfun包兼容问题，建议先安装xfun
if(!require(xfun)){install.packages("xfun")}
if(!require(rgexf)){install.packages('rgexf')}
#保存网络
t1$save_network(filepath = "network.gexf")

# 计算网络属性
t1$cal_network_attr()
t1$res_network_attr

#函数get_node_table、get_edge_table和get_adjacency_matrix分别用于从网络获取节点属性表、边属性表和邻接矩阵。
# 获得节点属性
t1$get_node_table(node_roles = TRUE)
# 获得边属性表
t1$get_edge_table()
# 获得邻接矩阵 
t1$get_adjacency_matrix()

#输出节点和边数据文件，可用于cytoscope软件绘制网络图
write.csv(t1$res_edge_table,"res_edge_table.csv")
write.csv(t1$res_node_table,"res_node_table.csv")

#绘制网络图
t1$plot_network(method = "ggraph", node_color = "module",node_label = NA)

#可选其他的绘图方式
#t1$plot_network(method = "networkD3", node_color = "module")
#t1$plot_network(method = "igraph")

#模块内连接和模块间连接绘制节点分类
#add_label = TRUE 可以用于直接为点添加文本标签
t1$plot_taxa_roles(use_type = 1)

# plot node roles with phylum information
t1$plot_taxa_roles(use_type = 2)

#提取子网络
#提取特定模块的子网络
#t1$subset_network(node = t1$res_node_table %>% base::subset(module == "M1") %>% rownames, rm_single = TRUE)
#提取特定类型边的子网络(正相关)
#t1$subset_network(edge = "+")

#实例
# 提取样本 'S1' 的子网络
sub1 <- t1$subset_network(node = mt$otu_table %>% .[.[, "S1"] != 0, ] %>% rownames, rm_single = TRUE)
# 通过直接赋值创建 t2 对象
t2 <- t1
t2$res_network <- sub1  # 将提取的 'S1' 子网络赋值给 t2 的网络属性

# 现在 t2 包含 'S1' 的网络，可以用于进一步分析
t2$cal_module()  # 计算 t2 网络的模块
t2$save_network("S1.gexf")  # 保存 'S1' 的网络为 .gexf 格式文件

#绘制子网络图
t2$plot_network(method = "ggraph", node_color = "module",node_label = NA)

#对于微生物网络图的绘制多种方法，可以综合其他方法的选择进行选取
