#加载数据包
library(reshape2)
library(haven)
library(ggplot2)
library(bootnet)
library(qgraph)
library(BiocGenerics)
library(graph)
library(ggm)
library(dplyr)
library(tidyr)
library(psych)
library(igraph)
library(scales)
library(networktools)
library(ggpubr)
library(nodeIdentifyR)
##########################################################################################################################################################
#数据准备
#加载数据
violence_data <- read_spss("C:\\Users\\16124\\OneDrive\\桌面\\network analysis\\case for Ising model.sav")
names(violence_data)

# 选择需要的变量，并过滤无效数据
violence_data_clean <- violence_data %>%
  select(q8_unemployed,q00901_alone,q175_relation,q174_stress,aspd_3_imp,
         q022_physical,child_mal,AUDIT_class20,druguse,
         hadsanx1plus,hadsdep1plus,psychosis,stalk_c,
         gambpathgy,pornaddic,standard1,standard2) %>%
  filter(complete.cases(.)) # 移除缺失值
node_names <- names(violence_data_clean)

# 检查缺失值
any(is.na(violence_data_clean))
# 检查无穷值
any(sapply(violence_data_clean, is.infinite))
#检查是否数据属于二进制数据
data_check <- apply(violence_data_clean, 2, function(x) all(unique(x) %in% c(0, 1)))
print(data_check)
#############################################################################################################################################################
#网络模型估计
# 构建 Ising 网络
isingnet <- bootnet::estimateNetwork( violence_data_clean,
                                      rule = "OR",
                                      default = "IsingFit",
                                      tuning = 0.25
                                     )


#############################################################################################################################################################
isingnet_matrix <- isingnet$graph
nrow(isingnet_matrix)  # 检查矩阵的行数

#这里可以创建副本
isingnet_matrix_filtered <- isingnet_matrix

print(isingnet_matrix_filtered)

# 假设邻接矩阵是 symmetric，取上三角部分进行检查
conflicts <- isingnet_matrix_filtered * t(isingnet_matrix_filtered) < 0  # 正负相乘小于 0 表示冲突
if (any(conflicts[upper.tri(conflicts)])) {
  print("存在冲突的正负连接")
} else {
  print("无冲突连接")
}

# 定义节点名称
node_names <- c("q8_unemployed", "q00901_alone", "q175_relation", "q174_stress", 
                "aspd_3_imp", "q022_physical", "child_mal", "AUDIT_class20", 
                "druguse", "hadsanx1plus", "hadsdep1plus", "psychosis", 
                "stalk_c", "gambpathgy","pornaddic", "standard1", "standard2")

# 定义每个节点的组
node_groups <- c("group2", "group2", "group2", "group2", "group2", 
                 "group2", "group2", "group2", "group2", "group2",
                 "group2", "group2", "group1", "group1", 
                 "group1", "group1", "group1")  
# group1: gambpathgy, stalk_c, pornaddic, standard1, standard2

length(node_names)  # 检查 node_names 的长度
nrow(isingnet_matrix_filtered)  # 检查矩阵的行数

# 绘图
C_undirected_network <- qgraph(
  isingnet_matrix_filtered,
  layout = "spring",
  vsize = 9,
  groups = node_groups,
  palette = "colorblind",
  diag = FALSE,
  negDashed = TRUE, 
  legend = FALSE,
  fade = TRUE,         # 开启渐变效果
  maximum = 0.9,       # 最高权重边缘透明度设为完全不透明
  cut = 0.2,           # 忽略权重小于 0.2 的边
  labels = node_names,
  label.cex = 0.9,
  repulsion = 1,
  posCol = "#006400",
  negCol = "#FF9E9E"
)

C_undirected_network_centrality <- centrality(
  C_undirected_network,
  all.shortest.paths = TRUE  # 设置为 TRUE 计算所有节点间的最短路径
  )
print(C_undirected_network_centrality)

################################################################################################################################################################
#在ising模型中不会使用closeness和betweenness
centralityPlot(C_undirected_network, 
               include = c( "Closeness", "Betweenness"), 
               scale = "z-scores")
################################################################################################################################################################

# 定义节点的分组（即所属维度或类别）
Impact_factor_membership <- c(2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1)

# 计算桥中心性
bridge_centrality <- bridge(C_undirected_network, Impact_factor_membership)

# 创建数据框，合并桥强度和中心性指标
strength_centrality_merged <- data.frame(
  "Impact_factor" = names(bridge_centrality$`Bridge Strength`),
  "Strength" = C_undirected_network_centrality$InExpectedInfluence,
  "BridgeStrength" = bridge_centrality$`Bridge Strength`
)

# 按 Strength 排序
strength_centrality_merged <- strength_centrality_merged[
  order(strength_centrality_merged$Strength),]

# 确保 Impact_factor 列是因子，且按排序后的顺序绘制
strength_centrality_merged$Impact_factor <- factor(strength_centrality_merged$Impact_factor, 
                                              levels = strength_centrality_merged$Impact_factor)
str(strength_centrality_merged)

# 绘制 Strength 图
Strength_g <- ggplot(strength_centrality_merged, aes(y = Impact_factor, x = Strength, group = 1)) +geom_point() +geom_line() +theme_bw(base_size = 16)
# 显示图形对象
print(Strength_g)

# 绘制 Bridge Strength 图
BridgeStrength_g <- ggplot(strength_centrality_merged, aes(y = Impact_factor, x = BridgeStrength, group = 1)) +
  geom_point() +
  geom_line() +
  theme_bw(base_size = 16) +
  xlab("Bridge Strength")
print(BridgeStrength_g)

strength_g_merged <- ggarrange(
  Strength_g, BridgeStrength_g,
  nrow = 1, ncol = 2,
  widths = c(7, 7),  # 扩大宽度
  heights = c(3)         # 设置高度
)

print(strength_g_merged)

set.seed(123)
bootnet_nonpar <- bootnet(isingnet, 
                          nBoots = 1000,
                          nCores = 7)
edge_stability_g <- plot(bootnet_nonpar,labels = FALSE)
print(edge_stability_g)

set.seed(123)
bootnet_case_dropping <- bootnet(isingnet, 
                                 nBoots = 1000,
                                 type = "case",
                                 nCores = 6,
                                 statistics = c('edge', 'strength', 'bridgeStrength'),
                                 communities = Impact_factor_membership)
centrality_stability_g <- plot(bootnet_case_dropping, 'all')

print(centrality_stability_g)

# 提取上三角矩阵（不包含对角线）并取绝对值
gs_total <- sum(abs(edgeWeightMatrix[upper.tri(edgeWeightMatrix)]))

# 输出 GS 总分
print(gs_total)
##########################################################################################################################################################################
library(usethis)
library(devtools)

# 从 GitHub 安装 nodeIdentifyR
devtools::install_github("JasperNaberman/nodeIdentifyR", force = TRUE)


adj_matrix <- isingnet$graph  # 提取邻接矩阵
print(adj_matrix)

#计算GS总分
# 提取上三角矩阵（不包含对角线）并取绝对值
gs_total <- sum(abs(adj_matrix[upper.tri(adj_matrix)]))

# 输出 GS 总分
print(gs_total)

thresholdVector <- isingnet$intercepts  # 提取截距作为阈值

#对此模型进行模拟干预
gs_IsingSamples <- simulateResponses(adj_matrix,
                                     thresholdVector,
                                     "aggravating", 
                                     2)
gs_sumIsingSamples <- calculateSumScores(gs_IsingSamples)

#转为长格式方便绘图
gs_sumIsingSamplesLong <- prepareDFforPlottingAndANOVA(gs_sumIsingSamples)

plotSumScores(sum_scores_long = gs_sumIsingSamplesLong,
              perturbation_type = "aggravating",
              x_label_size = 12, y_label_size = 12)

?plotSumScores

# 生成图形并存储在变量中
result <- plotSumScores(sum_scores_long = gs_sumIsingSamplesLong, 
                        perturbation_type = "aggravating")

# 获取ggplot对象
sumScoresPlot <- result$sumScoresPlot  # 假设返回的对象中包含图形

# 调整标签和主题
sumScoresPlot <- sumScoresPlot + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # X轴标签
        axis.text.y = element_text(size = 12),                       # Y轴标签
        plot.title = element_text(size = 14, hjust = 0.5),          # 标题大小
        legend.text = element_text(size = 12))                       # 图例文本大小

# 打印修改后的图形
print(sumScoresPlot)


