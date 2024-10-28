library(haven)
library(tidyverse)
violence_data<-read_spss("C:\\Users\\16124\\OneDrive\\桌面\\network analysis\\CHN for wenjie.sav")
colnames(violence_data)

# 选择需要的变量，并过滤无效数据
violence_data_clean <- violence_data %>%
  select(actage,BMI_class,q8_unemployed,q014_educlass,q00901_alone,
         q175_relation,q174_stress,aspd_3_imp,q022_physical,child_mal,druguse,
         hadsanx1plus,hadsdep1plus,psychosis,drugscore,AUDIT,hads_anx,
         hads_dep,PSQ,condpdscore,nopsymbd) %>%
  filter(complete.cases(.)) # 移除缺失值

cor_matrix <- cor(violence_data_clean, 
                  method = "spearman", 
                  use = "pairwise.complete.obs")
# 设置相关性矩阵对角线为 0
diag(cor_matrix) <- 0

#相关性可视化
library(corrplot)
corrplot(cor_matrix, method = "color", 
         type = "upper",           # 只显示上三角
         tl.col = "black",        # 标签颜色
         tl.srt = 45,             # 标签旋转角度
         addCoef.col = "black",   # 在热图上添加相关系数
         number.cex = 0.3,       # 数字的大小
         tl.cex = 0.4)           # 调整标签大小

# 设置阈值并将低于阈值的相关性设为 0
threshold <- 0.1
cor_matrix[abs(cor_matrix) < threshold] <- 0

library(igraph)
# 创建网络图
network <- graph_from_adjacency_matrix(cor_matrix, 
                                       mode = "undirected", weighted = TRUE, 
                                       diag = FALSE)
#将可能存在的边的负权重转化为其绝对值
E(network)$weight <- abs(E(network)$weight)


library(ggraph)
library(tidygraph)
#将igraph产生的结果转化为ggraph能够识别的结果
tidy_network <- as_tbl_graph(network)

ggraph(tidy_network, layout = "fr") +  # 使用 Fruchterman-Reingold 布局
  geom_edge_link(aes(edge_width = weight), alpha = 0.5) +  # 根据边的权重调整宽度
  geom_node_point(color = "lightblue", size = 5) +  # 节点的颜色和大小
  geom_node_label(aes(label = name), size = 3) +  # 节点标签
  theme_void() +  # 去掉背景网格和坐标轴
  ggtitle("Network Analysis of Spearman Correlation")  # 添加标题

library(qgraph)

# 假设你的相关性矩阵是 cor_matrix
# 这里直接使用qgraph进行绘图
qgraph(cor_matrix,
       layout = "spring",  # 布局类型
       edge.width = abs(cor_matrix),  # 边的宽度根据相关性绝对值调整
       labels = colnames(cor_matrix),  # 使用列名作为节点标签
       label.cex = 1,  # 调整标签大小
       node.color = "lightblue",  # 节点颜色
       vsize = 8,  # 节点大小
       title = "Network Analysis of Spearman Correlation",  # 添加标题
       legend = TRUE  # 显示图例
)