#required packages (install them first = install.packages(name of the package))
library(qgraph)
library(stats)
library(readr)
library(ggplot2)
library(bootnet)
library(igraph)
library(dplyr)
library(mgm)
library(reshape2)
library(data.table)
library(haven)
library(pbapply)  # 用于显示进度条
library(parallel)
##################################################################################################################################################################

violence_data<-read_spss("C:\\Users\\16124\\OneDrive\\桌面\\network analysis\\CHN for wenjie.sav")
colnames(violence_data)

# 选择需要的变量，并过滤无效数据
violence_data_clean <- violence_data %>%
  select(standard3_total,audit_1,aspd_tot,stalk,gambling,pornadscore,
         drugscore,q126_total,hads_anx,hads_dep,PSQ,condpdscore) %>%
  filter(complete.cases(.)) # 移除缺失值

# 检查缺失值
any(is.na(violence_data_clean))
# 检查无穷值
any(sapply(violence_data_clean, is.infinite))

# 定义多个 tuning 值
tuning_values <- seq(0.01, 1, by = 0.01)

# 使用 pblapply 显示进度条来计算每个 tuning 值下的网络
models <- pblapply(tuning_values, function(t) {
  bootnet::estimateNetwork(
                          data = violence_data_clean,
                          default = "EBICglasso",
                          corMethod = "cor",
                          corArgs = list(method = "spearman",
                                         use = "pairwise.complete.obs"),
                          tuning = t
                           )
                    })

# 假设评估标准是网络边数，可以选边数适中的模型
target_edges <- 10  # 假设期望的边数
model_edges <- sapply(models, function(m) sum(m$graph != 0))
best_tuning <- tuning_values[which.min(abs(model_edges - target_edges))]

print(best_tuning)
print(data.frame(tuning = tuning_values, edges = model_edges))


network1 <- bootnet::estimateNetwork(
                                     violence_data_clean,
                                     default = "EBICglasso",
                                     corMethod = "cor",
                                     corArgs = list(method = "spearman",
                                                    use = "pairwise.complete.obs"),
                                     tuning = best_tuning
                                    )

#violence_data_clean格式为数据框，mgm输入需要矩阵
data_matrix <- as.matrix(violence_data_clean)
# 确定节点数量
num_nodes <- ncol(data_matrix)

# 定义节点类型（假设所有节点为连续型）
type <- rep("g", num_nodes)

# 定义每个节点的水平数（连续型节点设为 1,离散变量根据分类数来判断）
level <- rep(1, num_nodes)
# 将数据框转换为矩阵
data_matrix <- as.matrix(violence_data_clean)

# 拟合模型
fit1 <- mgm(data = data_matrix,
            type = type, 
            level = rep(1, num_nodes))  # 每个节点的最大邻居数


# 预测
pred1 <- predict(fit1, violence_data_clean)

# 查看预测误差
pred1$error

# 打印每个节点的可预测性值
print(pred1$error)

# 提取R2信息
R2_values <- pred1$errors$R2

# 打印R2信息
print(R2_values)
warning()

# 从 network1 提取加权邻接矩阵
Network1 <- getWmat(network1)

class(Network1)

colnames(violence_data_clean)
items <- colnames(violence_data_clean)

# 假设 pred1$error 是包含预测误差的数据框，并且 pred1$error[,2] 是饼图数据
plot_network <- qgraph(Network1,  
                       layout = "spring",
                       labels = items,
                       label.cex = 0.9,
                       label.color = "black",
                       negDashed = TRUE,              # 负相关的边使用虚线
                       pie = pred1$error[, 3],
                       pieColor = "skyblue",
                       color = "#fdae61",
                       fade = FALSE,
                       vsize = 10
                       )

# 计算中心性
centrality_results <- centrality_auto(plot_network)

centralityPlot(plot_network, 
               include = c("Strength", "Closeness", "Betweenness"), 
               scale = "z-scores")

##############################################################################################################################################
#复杂化计算中心性
#  network1 是由 estimateNetwork 生成的对象,转化为对象邻接矩阵
adj_matrix <- as.matrix(network1$graph)

network_qgraph <- qgraph(adj_matrix)

# 计算中心性
centrality_results <- centrality_auto(network_qgraph)
str(centrality_results)

print(centrality_results$node.centrality)

# 提取需要的中心性指标
node_centrality <- centrality_results$node.centrality[, c("Betweenness",
                                                          "Closeness", 
                                                          "Strength")]

# 将节点名称添加为列，以便在绘图时标识
node_centrality$Node <- rownames(node_centrality)

# 转换为长格式，方便绘图
library(tidyr)
centrality_data_long <- pivot_longer(node_centrality, 
                                     cols = c("Betweenness", "Closeness", "Strength"), 
                                     names_to = "Centrality", 
                                     values_to = "Value")

# 绘制 Betweenness 中介中心性折线图
ggplot(node_centrality, aes(x = factor(Node), y = Betweenness, group = 1)) +
  geom_line(color = "skyblue", size = 1) +
  geom_point(color = "skyblue", size = 2) +
  labs(x = "Node", y = "Betweenness Centrality") +
  ggtitle("Betweenness Centrality of Nodes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 绘制 Closeness 紧密中心性折线图
ggplot(node_centrality, aes(x = factor(Node), y = Closeness, group = 1)) +
  geom_line(color = "salmon", size = 1) +
  geom_point(color = "salmon", size = 2) +
  labs(x = "Node", y = "Closeness Centrality") +
  ggtitle("Closeness Centrality of Nodes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 绘制 Strength 强度中心性折线图
ggplot(node_centrality, aes(x = factor(Node), y = Strength, group = 1)) +
  geom_line(color = "lightgreen", size = 1) +
  geom_point(color = "lightgreen", size = 2) +
  labs(x = "Node", y = "Strength Centrality") +
  ggtitle("Strength Centrality of Nodes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
##############################################################################################################################################

#计算边线的Bootstrap置信区间
class(network1)  # 查看对象类型
str(plot_network)
is.na(plot_network$Arguments$input)

# 提取邻接矩阵
adj_matrix1 <- plot_network$Arguments$input
# 确保它是数值矩阵
adj_matrix_real <- as.matrix(adj_matrix1)

#确定电脑核心数目
num_cores <- detectCores()
print(num_cores)

# 使用 bootnet 进行引导估计
Result1 <- bootnet::bootnet(network1, 
                           statistics = c("Strength", "Closeness", "Betweenness", "edge"), 
                           nBoots = 1000,
                           nCores = num_cores-1
                           )
pdf("C:/Users/16124/OneDrive/桌面/置信区间1.pdf", width = 5, height = 10)
plot(Result1,  order = "sample")
dev.off()

plot(Result1, "strength", plot = "difference")    # 分析节点强度上的差异
plot(Result1, "closeness", plot = "difference")   # 分析节点紧密性上的差异
plot(Result1, "betweenness", plot = "difference") # 分析节点中介性上的差异

pdf("C:/Users/16124/OneDrive/桌面/网络分析结果/置信区间差异.pdf", width = 50, height = 50)
plot(Result1, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")
dev.off()

#节点中心指标的稳定性
Result2 <- bootnet(
                   network1, 
                   statistics = c("Strength", "Closeness", "Betweenness"), 
                   nBoots = 1000, 
                   nCores = num_cores, 
                   type = "case"
                    )

plot(Result2, statistics = c("Strength", "Closeness", "Betweenness"))
corStability(Result2)

# 获取 CS 分数
cs_scores <- corStability(Result2)
# 比如 cs_scores$strength、cs_scores$closeness、cs_scores$betweenness
# 假设 cs_scores 是一个包含 3 个值的向量
cs_data <- data.frame(
                      Metric = c("Strength", "Closeness", "Betweenness"),
                      CS_Score = cs_scores  # cs_scores 是一个向量，直接使用
                      )

# 绘制 CS 分数的柱状图
ggplot(cs_data, aes(x = Metric, y = CS_Score)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = round(CS_Score, 2)), vjust = -0.5) +
  ylim(0, 1) +
  labs(title = "CS 分数（中心性稳定性系数）", x = "中心性指标", y = "CS 分数") +
  theme_minimal()


