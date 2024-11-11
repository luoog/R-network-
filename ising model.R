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
print(isingnet)

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
gs_IsingSamples_aggravating  <- simulateResponses(adj_matrix,
                                     thresholdVector,
                                     "aggravating", 
                                     2)
gs_sumIsingSamples_aggravating  <- calculateSumScores(gs_IsingSamples_aggravating )
print(gs_sumIsingSamples_aggravating)
print(gs_sumIsingSamplesLong_aggravating)

#转为长格式方便绘图
gs_sumIsingSamplesLong_aggravating <- prepareDFforPlottingAndANOVA(gs_sumIsingSamples_aggravating)

# 调用 plotSumScores 绘制增强干预的图形
result_aggravating  <- plotSumScores(sum_scores_long = gs_sumIsingSamplesLong_aggravating ,
                        perturbation_type = "aggravating")

# 获取ggplot对象
sumScoresPlot_aggravating  <- result_aggravating $sumScoresPlot  # 假设返回的对象中包含图形
print(sumScoresPlot_aggravating)

# 调整标签和主题
sumScoresPlot_aggravating  <- sumScoresPlot_aggravating  + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # X轴标签
        axis.text.y = element_text(size = 12),                       # Y轴标签
        plot.title = element_text(size = 14, hjust = 0.5),          # 标题大小
        legend.text = element_text(size = 12))                       # 图例文本大小

# 打印修改后的图形
print(sumScoresPlot_aggravating )

# 模拟减弱（缓解）干预
gs_IsingSamples_alleviating <- simulateResponses(adj_matrix,
                                                 thresholdVector,
                                                 "alleviating",  # 缓解
                                                 2)

# 计算减弱干预的总分
gs_sumIsingSamples_alleviating <- calculateSumScores(gs_IsingSamples_alleviating)

# 转换为长格式以便绘图
gs_sumIsingSamplesLong_alleviating <- prepareDFforPlottingAndANOVA(gs_sumIsingSamples_alleviating)

# 调用 plotSumScores 绘制减弱干预的图形
result_alleviating <- plotSumScores(sum_scores_long = gs_sumIsingSamplesLong_alleviating,
                        perturbation_type = "alleviating")
print(result_alleviating)

# 获取ggplot对象
sumScoresPlot_alleviating <- result_alleviating$sumScoresPlot  

# 调整标签和主题
sumScoresPlot_alleviating <- sumScoresPlot_alleviating + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # X轴标签
        axis.text.y = element_text(size = 12),                       # Y轴标签
        plot.title = element_text(size = 14, hjust = 0.5),          # 标题大小
        legend.text = element_text(size = 12))                       # 图例文本大小

# 打印修改后的图形
print(sumScoresPlot_alleviating)

# 从增强和减弱干预结果中提取数据
print(result_aggravating)
data_aggravating <- result_aggravating$plottedInformation
data_alleviating <- result_alleviating$plottedInformation
print(data_aggravating)

# 添加干预类型列
data_aggravating$perturbation_type <- "Aggravating"
data_alleviating$perturbation_type <- "Alleviating"

# 合并数据框
combined_data <- rbind(data_aggravating, data_alleviating)
print(combined_data)

# 设置 X 轴顺序
combined_data$thresholdIteration <- factor(combined_data$thresholdIteration, 
                                           levels = c("original", "child_mal", "aspd_3_imp", 
                                                      "standard2", "standard1", "hadsanx1plus",
                                                      "q174_stress", "psychosis", "hadsdep1plus",
                                                      "q8_unemployed", "q00901_alone", "q022_physical",
                                                      "AUDIT_class20", "q175_relation", "gambpathgy",
                                                      "stalk_c", "druguse", "pornaddic"))
print(combined_data)

ggplot(combined_data, aes(x = thresholdIteration, y = meanSumscore, group = perturbation_type, linetype = perturbation_type)) +
  geom_line(aes(color = perturbation_type), size = 1) +       # 绘制线条并设置颜色
  geom_point(aes(color = perturbation_type), size = 2) +      # 添加数据点
  labs(x = "Symptom of which the threshold is altered", y = "Sum score", title = "Effect of Aggravating vs Alleviating Interventions") +
  scale_color_manual(values = c("Aggravating" = "black", "Alleviating" = "black")) + # 设置颜色为黑色
  scale_linetype_manual(values = c("Aggravating" = "solid", "Alleviating" = "dashed")) + # 设置线条类型
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  # X轴标签的字体大小和角度
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5)
  )





