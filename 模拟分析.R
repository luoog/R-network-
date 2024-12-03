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
library(Rmisc)
library(grid)
library(RColorBrewer)
library(nodeIdentifyR)
library(gridExtra)


#数据准备
#加载数据
violence_data <- read_spss("C:\\Users\\16124\\OneDrive\\桌面\\network analysis\\case for Ising model.sav")
names(violence_data)

# 选择需要的变量，并过滤无效数据
violence_data_clean <- violence_data %>%
  select(q8_unemployed,q00901_alone,q175_relation,q174_stress,aspd_3_imp,
         q022_physical,child_mal,AUDIT_class20,druguse,
         hadsanx1plus,hadsdep1plus,psychosis,stalk_c,
         gambpathgy,pornprob,vio_others,vio_self) %>%
  filter(complete.cases(.)) # 移除缺失值
node_names <- names(violence_data_clean)
print(violence_data_clean)

#网络模型估计
# 构建 Ising 网络
isingnet <- bootnet::estimateNetwork( violence_data_clean,
                                      rule = "OR",
                                      default = "IsingFit",
                                      tuning = 0.25
)
print(isingnet)

# 提取邻接矩阵
adj_matrix <- as.matrix(isingnet$graph)

# 将邻接矩阵转换为 igraph 图对象
g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)

# 运行 Spinglass 算法
spinglass_result <- cluster_spinglass(g, implementation = "neg")

# 查看聚类结果
print(spinglass_result)  # 输出社区分配
print(membership(spinglass_result))  # 查看每个节点的社区编号

# 获取社区分配
community_membership <- membership(spinglass_result)  # 获取每个节点的社区编号

# 为四个社区选择颜色，使用RColorBrewer生成调色板
community_colors <- brewer.pal(4, "Set3")  # 选择四种颜色

# 或者你也可以直接指定颜色
# community_colors <- c("#FF5733", "#33FF57", "#3357FF", "#F8FF33")

# 根据社区编号分配颜色
V(g)$color <- community_colors[community_membership]

# 继续绘制图形
C_undirected_network <- qgraph(
  adj_matrix,
  layout = "spring",       # 使用 Spring Layout
  vsize = 9,               # 节点大小
  cut = 0.2,               # 忽略权重小于 0.2 的边
  labels = node_names,     # 节点标签
  diag = FALSE,            # 不显示对角线
  negDashed = TRUE,        # 负边为虚线
  legend = FALSE,          # 不显示图例
  fade = TRUE,             # 渐变效果
  maximum = 0.9,           # 最大权重边完全不透明
  label.cex = 1.1,           # 节点标签大小
  posCol = "#006400",      # 正权重边颜色
  negCol = "#FF9E9E",      # 负权重边颜色
  groups = community_membership,  # 节点分组
  color = V(g)$color,      # 使用指定的颜色
  repulsion = 1            # 节点间排斥力
)

###########################################################################################################################################################################################################################################################################################################
adj_matrix <- isingnet$graph  # 提取邻接矩阵
thresholdVector <- isingnet$intercepts  # 提取截距作为阈值
print(thresholdVector)

set.seed(123)
gs_IsingSamples_aggravating2  <- simulateResponses(adj_matrix,
                                                  thresholdVector,
                                                  "aggravating", 
                                                  2)
str(gs_IsingSamples_aggravating2)
print(gs_IsingSamples_aggravating2)

# 去除 vio_others 和 vio_self
remaining_data_2 <- gs_IsingSamples_aggravating2[-c(17, 18)]
str(remaining_data_2)

# 提取这16个数据框中的第16列和第17列
column_16_2 <- lapply(remaining_data_2, function(df) df[, 16])
column_17_2 <- lapply(remaining_data_2, function(df) df[, 17])
str(column_16_2)

# 计算每个变量的激活率（即每个向量的1的比例）
activation_rates_others2 <- sapply(column_16_2, mean)
str(activation_rates_others2)
activation_rates_self2 <- sapply(column_17_2,mean)
str(activation_rates_self2)

# 计算标准误差（SE）
SE_self_2 <- sapply(column_17_2, function(x) {
  p <- mean(x)  # 激活率（均值）
  n <- length(x)  # 样本大小
  sqrt((p * (1 - p)) / n)  # 计算标准误差
})

SE_others_2 <- sapply(column_16_2, function(x) {
  p <- mean(x)  # 激活率（均值）
  n <- length(x)  # 样本大小
  sqrt((p * (1 - p)) / n)  # 计算标准误差
})

# 计算95%置信区间
CI_self_lower_2 <- activation_rates_self2 - 1.96 * SE_self_2
CI_self_upper_2 <- activation_rates_self2 + 1.96 * SE_self_2

CI_others_lower_2 <- activation_rates_others2 - 1.96 * SE_others_2
CI_others_upper_2 <- activation_rates_others2 + 1.96 * SE_others_2

# 创建一个数据框，包括激活率及其置信区间
df_self_2 <- data.frame(
  Variable = names(activation_rates_self2),
  Mean = activation_rates_self2,
  LowerCI = CI_self_lower_2,
  UpperCI = CI_self_upper_2
)

df_others_2 <- data.frame(
  Variable = names(activation_rates_others2),
  Mean = activation_rates_others2,
  LowerCI = CI_others_lower_2,
  UpperCI = CI_others_upper_2
)

set.seed(123)
gs_IsingSamples_aggravating3  <- simulateResponses(adj_matrix,
                                                   thresholdVector,
                                                   "aggravating", 
                                                   3)
# 去除 vio_others 和 vio_self
remaining_data_3 <- gs_IsingSamples_aggravating3[-c(17, 18)]

# 提取这16个数据框中的第16列和第17列
column_16_3 <- lapply(remaining_data_3, function(df) df[, 16])
column_17_3 <- lapply(remaining_data_3, function(df) df[, 17])

# 计算每个变量的激活率（即每个向量的1的比例）
activation_rates_others3 <- sapply(column_16_3, mean)
str(activation_rates_others3)
str(activation_rates_others2)
activation_rates_self3 <- sapply(column_17_3,mean)
str(activation_rates_self3)
str(activation_rates_self2)

# 计算标准误差（SE）
SE_self_3 <- sapply(column_17_3, function(x) {
  p <- mean(x)  # 激活率（均值）
  n <- length(x)  # 样本大小
  sqrt((p * (1 - p)) / n)  # 计算标准误差
})

SE_others_3 <- sapply(column_16_3, function(x) {
  p <- mean(x)  # 激活率（均值）
  n <- length(x)  # 样本大小
  sqrt((p * (1 - p)) / n)  # 计算标准误差
})

# 计算95%置信区间
CI_self_lower_3 <- activation_rates_self3 - 1.96 * SE_self_3
CI_self_upper_3 <- activation_rates_self3 + 1.96 * SE_self_3

CI_others_lower_3 <- activation_rates_others3 - 1.96 * SE_others_3
CI_others_upper_3 <- activation_rates_others3 + 1.96 * SE_others_3

# 创建一个数据框，包括激活率及其置信区间
df_self_3 <- data.frame(
  Variable = names(activation_rates_self3),
  Mean = activation_rates_self3,
  LowerCI = CI_self_lower_3,
  UpperCI = CI_self_upper_3
)

df_others_3 <- data.frame(
  Variable = names(activation_rates_others3),
  Mean = activation_rates_others3,
  LowerCI = CI_others_lower_3,
  UpperCI = CI_others_upper_3
)

# 合并 self 数据框
df_self <- rbind(
  transform(df_self_2, Intervention = 2),  # 添加干预强度列为 2
  transform(df_self_3, Intervention = 3)   # 添加干预强度列为 3
)

# 合并 others 数据框
df_others <- rbind(
  transform(df_others_2, Intervention = 2),  # 添加干预强度列为 2
  transform(df_others_3, Intervention = 3)   # 添加干预强度列为 3
)

# 确保 'original' 在第一个位置，不修改其他变量的排序
df_self$Variable <- factor(df_self$Variable, levels = c("original", setdiff(unique(df_self$Variable), "original")))

# 绘制 self 激活率
p_self <- ggplot(df_self, aes(x = Variable, y = Mean, color = as.factor(Intervention))) +
  geom_point(position = position_dodge(width = 0.5)) +  # 增加 width，变宽
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, position = position_dodge(width = 0.7)) +
  labs(x = "Variable", y = "Activation Rate (Self)", color = "Intervention Intensity") +
  scale_color_manual(values = c("2" = "#6fa3e3",        # 干预强度 2 - 柔和蓝色
                                "3" = "#f46d6d")) +     # 干预强度 3 - 柔和红色
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # 倾斜 45 度
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

# 绘制 others 激活率，确保 'original' 排在第一
df_others$Variable <- factor(df_others$Variable, levels = c("original", setdiff(unique(df_others$Variable), "original")))

p_others <- ggplot(df_others, aes(x = Variable, y = Mean, color = as.factor(Intervention))) +
  geom_point(position = position_dodge(width = 0.5)) +  # 稍微错开位置
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, position = position_dodge(width = 0.5)) +
  labs(x = "Variable", y = "Activation Rate (Others)", color = "Intervention Intensity") +
  scale_color_manual(values = c("2" = "#6fa3e3",        # 干预强度 2 - 柔和蓝色
                                "3" = "#f46d6d")) +     # 干预强度 3 - 柔和红色
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # 倾斜 45 度
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

# 显示图表
print(p_self)
print(p_others)

grid.arrange(
  arrangeGrob(p_self, top = textGrob("A", gp = gpar(fontsize = 20, fontface = "bold"))),
  arrangeGrob(p_others, top = textGrob("B", gp = gpar(fontsize = 20, fontface = "bold"))),
  ncol = 2
)

# 模拟减弱（缓解）干预
set.seed(123)
gs_IsingSamples_alleviating2 <- simulateResponses(adj_matrix,
                                                 thresholdVector,
                                                 "alleviating",  # 缓解
                                                 2)
set.seed(123)
gs_IsingSamples_alleviating3 <- simulateResponses(adj_matrix,
                                                  thresholdVector,
                                                  "alleviating",  # 缓解
                                                  3)

str(gs_IsingSamples_alleviating2)
str(gs_IsingSamples_alleviating3)

# 去除 vio_others 和 vio_self
remaining_data_alleviating_2 <- gs_IsingSamples_alleviating2[-c(17, 18)]
remaining_data_alleviating_3 <- gs_IsingSamples_alleviating3[-c(17, 18)]
str(remaining_data_alleviating_2)
str(remaining_data_alleviating_3)

# 提取这16个数据框中的第16列和第17列
column_16_alleviating_2 <- lapply(remaining_data_alleviating_2, function(df) df[, 16])
column_17_alleviating_2 <- lapply(remaining_data_alleviating_2, function(df) df[, 17])

column_16_alleviating_3 <- lapply(remaining_data_alleviating_3, function(df) df[, 16])
column_17_alleviating_3 <- lapply(remaining_data_alleviating_3, function(df) df[, 17])

activation_rates_alleviating_others2 <- sapply(column_16_alleviating_2, mean)
activation_rates_alleviating_self2 <- sapply(column_17_alleviating_2,mean)
activation_rates_alleviating_others3 <- sapply(column_16_alleviating_3, mean)
activation_rates_alleviating_self3 <- sapply(column_17_alleviating_3,mean)

# 计算标准误差（SE）
SE_self_alleviating_2 <- sapply(column_17_alleviating_2, function(x) {
  p <- mean(x)  # 激活率（均值）
  n <- length(x)  # 样本大小
  sqrt((p * (1 - p)) / n)  # 计算标准误差
})
SE_self_alleviating_3 <- sapply(column_17_alleviating_3, function(x) {
  p <- mean(x)  # 激活率（均值）
  n <- length(x)  # 样本大小
  sqrt((p * (1 - p)) / n)  # 计算标准误差
})

SE_others_alleviating_2 <- sapply(column_16_alleviating_2, function(x) {
  p <- mean(x)  # 激活率（均值）
  n <- length(x)  # 样本大小
  sqrt((p * (1 - p)) / n)  # 计算标准误差
})
SE_others_alleviating_3 <- sapply(column_16_alleviating_3, function(x) {
  p <- mean(x)  # 激活率（均值）
  n <- length(x)  # 样本大小
  sqrt((p * (1 - p)) / n)  # 计算标准误差
})

# 计算95%置信区间
CI_self_lower_alleviating_3 <- activation_rates_alleviating_self3 - 1.96 * SE_self_alleviating_3
CI_self_upper_alleviating_3 <- activation_rates_alleviating_self3 + 1.96 * SE_self_alleviating_3

CI_others_lower_alleviating_3 <- activation_rates_alleviating_others3 - 1.96 * SE_others_alleviating_3
CI_others_upper_alleviating_3 <- activation_rates_alleviating_others3 + 1.96 * SE_others_alleviating_3

CI_self_lower_alleviating_2 <- activation_rates_alleviating_self2 - 1.96 * SE_self_alleviating_2
CI_self_upper_alleviating_2 <- activation_rates_alleviating_self2 + 1.96 * SE_self_alleviating_2

CI_others_lower_alleviating_2 <- activation_rates_alleviating_others2 - 1.96 * SE_others_alleviating_2
CI_others_upper_alleviating_2 <- activation_rates_alleviating_others2 + 1.96 * SE_others_alleviating_2

# 创建一个数据框，包括激活率及其置信区间
df_self_alleviating_3 <- data.frame(
  Variable = names(activation_rates_alleviating_self3),
  Mean = activation_rates_alleviating_self3,
  LowerCI = CI_self_lower_alleviating_3,
  UpperCI = CI_self_upper_alleviating_3
)

df_self_alleviating_2 <- data.frame(
  Variable = names(activation_rates_alleviating_self2),
  Mean = activation_rates_alleviating_self2,
  LowerCI = CI_self_lower_alleviating_2,
  UpperCI = CI_self_upper_alleviating_2
)

df_others_alleviating_2 <- data.frame(
  Variable = names(activation_rates_alleviating_others2),
  Mean = activation_rates_alleviating_others2,
  LowerCI = CI_others_lower_alleviating_2,
  UpperCI = CI_others_upper_alleviating_2
)

df_others_alleviating_3 <- data.frame(
  Variable = names(activation_rates_alleviating_others3),
  Mean = activation_rates_alleviating_others3,
  LowerCI = CI_others_lower_alleviating_3,
  UpperCI = CI_others_upper_alleviating_3
)

df_alleviating_self <- rbind(
  transform(df_self_alleviating_2, Intervention = 2),  # 添加干预强度列为 2
  transform(df_self_alleviating_3, Intervention = 3)   # 添加干预强度列为 3
)

df_alleviating_others <- rbind(
  transform(df_others_alleviating_2, Intervention = 2),  # 添加干预强度列为 2
  transform(df_others_alleviating_3, Intervention = 3)   # 添加干预强度列为 3
)

# 确保 'original' 在第一个位置，不修改其他变量的排序
df_alleviating_self$Variable <- factor(df_alleviating_self$Variable, levels = c("original", setdiff(unique(df_alleviating_self$Variable), "original")))

# 绘制 self 激活率
p_alleviating_self <- ggplot(df_alleviating_self, aes(x = Variable, y = Mean, color = as.factor(Intervention))) +
  geom_point(position = position_dodge(width = 0.5)) +  # 增加 width，变宽
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, position = position_dodge(width = 0.7)) +
  labs(x = "Variable", y = "Activation Rate (Self)", color = "Intervention Intensity") +
  scale_color_manual(values = c("2" = "#6fa3e3",        # 干预强度 2 - 柔和蓝色
                                "3" = "#f46d6d")) +     # 干预强度 3 - 柔和红色
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # 倾斜 45 度
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

# 确保 'original' 排在第一
df_alleviating_others$Variable <- factor(df_alleviating_others$Variable, levels = c("original", setdiff(unique(df_alleviating_others$Variable), "original")))

# 绘制 others 激活率
p_alleviating_others <- ggplot(df_alleviating_others, aes(x = Variable, y = Mean, color = as.factor(Intervention))) +
  geom_point(position = position_dodge(width = 0.5)) +  # 稍微错开位置
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, position = position_dodge(width = 0.5)) +
  labs(x = "Variable", y = "Activation Rate (Others)", color = "Intervention Intensity") +
  scale_color_manual(values = c("2" = "#6fa3e3",        # 干预强度 2 - 柔和蓝色
                                "3" = "#f46d6d")) +     # 干预强度 3 - 柔和红色
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # 倾斜 45 度
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

# 显示图表
print(p_alleviating_self)
print(p_alleviating_others)

# 将两个图表合并在一起
grid.arrange(
  arrangeGrob(p_alleviating_self, top = textGrob("A", gp = gpar(fontsize = 20, fontface = "bold"))),
  arrangeGrob(p_alleviating_others, top = textGrob("B", gp = gpar(fontsize = 20, fontface = "bold"))),
  ncol = 2
)

###########################################################################################################################################################################################################################################################################################################




