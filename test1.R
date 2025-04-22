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
library(png)

############################################################################################################################
############################################################################################################################
#build network structure
violence_data <- read_spss("C:\\Users\\16124\\OneDrive\\桌面\\network analysis\\case for Ising model.sav")
names(violence_data)

violence_data_clean <- violence_data %>%
  select(q8_unemployed, q00901_alone, q175_relation, q174_stress, aspd_3_imp,
         q022_physical, child_mal, AUDIT_class20, druguse,
         hadsanx1plus, hadsdep1plus, psychosis, stalk_c,
         gambpathgy, pornprob, vio_others, vio_self) %>%
  filter(complete.cases(.)) 

node_names <- names(violence_data_clean)
print(violence_data_clean)

isingnet <- bootnet::estimateNetwork(violence_data_clean,
                                     rule = "OR",           
                                     default = "IsingFit",   
                                     tuning = 0.25)          
print(isingnet)

Names <-c("Unemployment","Living alone","No close relationship","Adulthood stress event","Impulsivity",
          "Chronic illness/disability","Child maltreatment","Alcohol dependence","Drug abuse","Anxiety","Depression",
          "Psychosis","Stalking","Pathological Gambling","Problem Pornography","Interpersonal violence","Self-directied violence")

adj_matrix_named <- isingnet$graph
rownames(adj_matrix_named) <- Names
colnames(adj_matrix_named) <- Names
print(adj_matrix_named)
g <- graph_from_adjacency_matrix(adj_matrix_named, mode = "undirected", weighted = TRUE)
spinglass_result <- cluster_spinglass(g, implementation = "neg")

print(spinglass_result)  
print(membership(spinglass_result))  

community_membership <- membership(spinglass_result) 
community_colors <- brewer.pal(4, "Set3")  

Labels<-c("A1","A2","A3","B1","B2",
          "C1","B3","D1","D2","C2",
          "C3","C4","B4","D3","D4",
          "B5","C5")


clusrer_group <-c("Community A: Social isolation","Community A: Social isolation","Community A: Social isolation",
                  "Community B: Externalizing violence","Community B: Externalizing violence",
                  "Community C: Internalizing violence","Community B: Externalizing violence","Community D: Addictive behavior",
                  "Community D: Addictive behavior","Community C: Internalizing violence","Community C: Internalizing violence",
                  "Community C: Internalizing violence","Community B: Externalizing violence","Community D: Addictive behavior",
                  "Community D: Addictive behavior","Community B: Externalizing violence","Community C: Internalizing violence")

graph1 <- plot(isingnet, layout = 'spring',groups = clusrer_group,maximum=0.9,cut =0.1,vsize=6,
               legend=TRUE,borders=TRUE,palette="pastel",posCol = "#008000",negCol = "#FF4500",
               fade = T,
               legend.cex=0.4,label.cex=1.2,height=10,width=10,
               labels=Labels, nodeNames=Names,mar=c(4,1,2,2))

set.seed(123)
bootnet_nonpar <- bootnet::bootnet(isingnet, 
                                   nBoots = 1000,
                                   nCores = 7)
edge_stability_g <- plot(bootnet_nonpar,labels = FALSE)
print(edge_stability_g)
print(community_membership)
set.seed(123)
bootnet_case_dropping <- bootnet::bootnet(isingnet, 
                                          nBoots = 1000,
                                          type = "case",
                                          nCores = 7,
                                          statistics = c("edge","bridgeStrength","expectedInfluence"),
                                          communities = community_membership)
centrality_stability_g <- plot(bootnet_case_dropping, 'all')

print(centrality_stability_g)
print(bootnet_case_dropping)

cs_result <- corStability(bootnet_case_dropping)
print(cs_result)

############################################################################################################################
#ABC community network analysis
abc_data <- violence_data_clean %>%
  select(q8_unemployed, q00901_alone, q175_relation,q174_stress, aspd_3_imp, child_mal, stalk_c, vio_others,
         q022_physical, hadsanx1plus, hadsdep1plus, psychosis, vio_self)
names(abc_data)

abc_net <- bootnet::estimateNetwork(abc_data,
                                   rule = "OR",
                                   default = "IsingFit",
                                   tuning = 0.25)

abc_names <- c("Unemployment","Living alone","No close relationship",  # Community A
               "Adulthood stress event","Impulsivity","Child maltreatment","Stalking","Interpersonal violence",  # Community B
               "Chronic illness/disability","Anxiety","Depression","Psychosis","Self-Directed violence")  # Community C

abc_labels <- c("A1","A2","A3",
                "B1","B2","B3","B4","B5",
                "C1","C2","C3","C4","C5")

abc_group <- c("Community A: Social isolation","Community A: Social isolation","Community A: Social isolation",
               "Community B: Externalizing violence","Community B: Externalizing violence","Community B: Externalizing violence",
               "Community B: Externalizing violence","Community B: Externalizing violence",
               "Community C: Internalizing violence","Community C: Internalizing violence",
               "Community C: Internalizing violence","Community C: Internalizing violence",
               "Community C: Internalizing violence")


png("abc_network_wide.png", width = 1800, height = 1000, res = 150)
graph_abc <- plot(abc_net, layout = 'spring',groups = abc_group,maximum=0.9,cut =0.1,vsize=6,
                  legend=TRUE,borders=TRUE,palette="pastel",posCol = "#008000",negCol = "#FF4500",
                  fade = T,
                  legend.cex=0.4,label.cex=1.2,height=20,width=20,
                  labels=abc_labels, nodeNames=abc_names,mar=c(4,1,2,4))
dev.off()
getwd()

set.seed(123)
bootnet_nonpar_abc <- bootnet::bootnet(abc_net, 
                                      nBoots = 1000,
                                      nCores = 7)
edge_stability_g_abc <- plot(bootnet_nonpar_abc,labels = FALSE)
print(edge_stability_g_abc)
set.seed(123)
abc_boot <- bootnet(abc_net, 
                   nBoots = 1000, 
                   nCores = 7,
                   type = "case", 
                   statistics = c("edge","bridgeStrength","expectedInfluence"),
                   communities = abc_group
)
centrality_stability_abc <- plot(abc_boot, 'all')
print(centrality_stability_abc)
cs_result_abc <- corStability(abc_boot)
print(cs_result_abc)

#BCD community network analysis
bcd_data <- violence_data_clean %>%
  select(q174_stress, aspd_3_imp, child_mal, stalk_c, vio_others,
         q022_physical, hadsanx1plus, hadsdep1plus, psychosis, vio_self,
         AUDIT_class20, druguse, gambpathgy, pornprob)

bcd_net <- bootnet::estimateNetwork(bcd_data,
                                    rule = "OR",
                                    default = "IsingFit",
                                    tuning = 0.25)

bcd_names <- c("Adulthood stress event","Impulsivity","Child maltreatment","Stalking","Interpersonal violence",   # B
               "Chronic illness/disability","Anxiety","Depression","Psychosis","Self-directed violence",         # C
               "Alcohol dependence","Drug abuse","Pathological Gambling","Problem Pornography")                  # D

bcd_labels <- c("B1","B2","B3","B4","B5",
                "C1","C2","C3","C4","C5",
                "D1","D2","D3","D4")

bcd_group <- c("Community B: Externalizing violence","Community B: Externalizing violence","Community B: Externalizing violence",
               "Community B: Externalizing violence","Community B: Externalizing violence",
               "Community C: Internalizing violence","Community C: Internalizing violence",
               "Community C: Internalizing violence","Community C: Internalizing violence","Community C: Internalizing violence",
               "Community D: Addictive behavior","Community D: Addictive behavior",
               "Community D: Addictive behavior","Community D: Addictive behavior")

png("bcd_network_wide.png", width = 1800, height = 1000, res = 150)
graph_bcd <- plot(bcd_net, layout = 'spring',groups = bcd_group,maximum=0.9,cut =0.1,vsize=6,
                  legend=TRUE,borders=TRUE,palette="pastel",posCol = "#008000",negCol = "#FF4500",
                  fade = T,
                  legend.cex=0.4,label.cex=1.2,height=10,width=10,
                  labels=bcd_labels, nodeNames=bcd_names,mar=c(4,1,2,4))
dev.off()

set.seed(123)
bootnet_nonpar_bcd <- bootnet::bootnet(bcd_net, 
                                       nBoots = 1000,
                                       nCores = 7)
edge_stability_g_bcd <- plot(bootnet_nonpar_bcd,labels = FALSE)
print(edge_stability_g_bcd)
set.seed(123)
bcd_boot <- bootnet(bcd_net, 
                    nBoots = 1000, 
                    nCores = 7,
                    type = "case", 
                    statistics = c("edge","bridgeStrength","expectedInfluence"),
                    communities = bcd_group
)
centrality_stability_bcd <- plot(bcd_boot, 'all')
print(centrality_stability_bcd)
cs_result_bcd <- corStability(bcd_boot)
print(cs_result_bcd)

############################################################################################################################
out1 <- expectedInf(graph_abc, step = c("both", 1, 2), directed = FALSE)

out1$step1 <- out1$step1[!names(out1$step1) %in% c("Interpersonal violence", "Self-directed violence")]
out1$step2 <- out1$step2[!names(out1$step2) %in% c("Interpersonal violence", "Self-directed violence")]
plot(out1,order = "value",zscore = TRUE)

############################################################################################################################
bridge_centrality_abc <- bridge(graph_abc, communities=abc_group)
bridge_centrality_abc
str(bridge_centrality_abc)
print(bridge_centrality_abc[["Bridge Strength"]])
bridge_strength <- bridge_centrality_abc[["Bridge Strength"]]
bridge_strength_filtered <- bridge_strength[!names(bridge_strength) %in% 
                                              c("Interpersonal violence", "Self-Directed violence")]
print(bridge_strength_filtered)
plot(bridge_centrality_abc, 
     include=c("Bridge Strength"),
     order="value", 
     zscore=TRUE) 

############################################################################################################################
adj_matrix <- abc_net$graph
thresholdVector <- abc_net$intercepts  
str(adj_matrix)

new_labels <-c("Original","Unemployment","Living alone","No close relationship",  # Community A
          "Adulthood stress event","Impulsivity","Child maltreatment","Stalking", # Community B
          "Chronic illness/disability","Anxiety","Depression","Psychosis")

############################################################################################################################
set.seed(123)
gs_IsingSamples_aggravating2  <- simulateResponses(adj_matrix,
                                                   thresholdVector,
                                                   "aggravating", 
                                                   2)
set.seed(123)
gs_IsingSamples_aggravating3  <- simulateResponses(adj_matrix,
                                                   thresholdVector,
                                                   "aggravating", 
                                                   3)

names(gs_IsingSamples_aggravating2)
str(gs_IsingSamples_aggravating2)

remaining_data_2 <- gs_IsingSamples_aggravating2[!(names(gs_IsingSamples_aggravating2) %in% 
                                                     c("vio_self", "vio_others"))]
remaining_data_3 <- gs_IsingSamples_aggravating3[!(names(gs_IsingSamples_aggravating3) %in% 
                                                     c("vio_self", "vio_others"))]
str(remaining_data_3 )

column_interpersonal_2 <- lapply(remaining_data_2, function(df) df[, 8])
column_self_2 <- lapply(remaining_data_2, function(df) df[, 13])
column_interpersonal_3 <- lapply(remaining_data_3, function(df) df[, 8])
column_self_3 <- lapply(remaining_data_3, function(df) df[, 13])

activation_rates_others2 <- sapply(column_interpersonal_2, mean)
activation_rates_self2 <- sapply(column_self_2,mean)
activation_rates_others3 <- sapply(column_interpersonal_3, mean)
activation_rates_self3 <- sapply(column_self_3,mean)

SE_self_2 <- sapply(column_self_2, function(x) {
  p <- mean(x)  
  n <- length(x)  
  sqrt((p * (1 - p)) / n)  
})
SE_others_2 <- sapply(column_interpersonal_2, function(x) {
  p <- mean(x) 
  n <- length(x)  
  sqrt((p * (1 - p)) / n) 
})
SE_self_3 <- sapply(column_self_3, function(x) {
  p <- mean(x)  
  n <- length(x) 
  sqrt((p * (1 - p)) / n)  
})

SE_others_3 <- sapply(column_interpersonal_3, function(x) {
  p <- mean(x)  
  n <- length(x)  
  sqrt((p * (1 - p)) / n) 
})

CI_self_lower_2 <- activation_rates_self2 - 1.96 * SE_self_2
CI_self_upper_2 <- activation_rates_self2 + 1.96 * SE_self_2
CI_others_lower_2 <- activation_rates_others2 - 1.96 * SE_others_2
CI_others_upper_2 <- activation_rates_others2 + 1.96 * SE_others_2

CI_self_lower_3 <- activation_rates_self3 - 1.96 * SE_self_3
CI_self_upper_3 <- activation_rates_self3 + 1.96 * SE_self_3
CI_others_lower_3 <- activation_rates_others3 - 1.96 * SE_others_3
CI_others_upper_3 <- activation_rates_others3 + 1.96 * SE_others_3

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

df_self <- rbind(
  transform(df_self_2, Intervention = 2),  
  transform(df_self_3, Intervention = 3)   
)
df_others <- rbind(
  transform(df_others_2, Intervention = 2),  
  transform(df_others_3, Intervention = 3)   
)

df_self$Variable <- factor(df_self$Variable, levels = c("original", setdiff(unique(df_self$Variable), "original")))
df_others$Variable <- factor(df_others$Variable, levels = c("original", setdiff(unique(df_others$Variable), "original")))
df_self$Variable <- factor(df_self$Variable,levels = levels(df_self$Variable),labels = new_labels)
df_others$Variable <- factor(df_others$Variable,levels = levels(df_others$Variable),labels = new_labels)
print(df_self$Variable)

p_self <- ggplot(df_self, aes(x = Variable, y = Mean, color = as.factor(Intervention))) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, position = position_dodge(width = 0.7)) +
  labs(x = "Variables", 
       y = "Activation Rate ", 
       title = "Simulation of Activation Rate under Aggravating Intervention(Self)",
       color = "Intervention Intensity") +
  scale_color_manual(values = c("2" = "#6fa3e3", "3" = "#f46d6d")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30,hjust = 1),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12)
  )
print(p_self)

p_others <- ggplot(df_others, aes(x = Variable, y = Mean, color = as.factor(Intervention))) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, position = position_dodge(width = 0.5)) +
  labs(x = "Variables", 
       y = "Activation Rate",
       title = "Simulation of Activation Rate under Aggravating Intervention (Others)",
       color = "Intervention Intensity") +
  scale_color_manual(values = c("2" = "#6fa3e3", "3" = "#f46d6d")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30,hjust = 1),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12)
  )
print(p_others)

############################################################################################################################
set.seed(123)
gs_IsingSamples_alleviating2 <- simulateResponses(adj_matrix,
                                                  thresholdVector,
                                                  "alleviating",
                                                  2)
set.seed(123)
gs_IsingSamples_alleviating3 <- simulateResponses(adj_matrix,
                                                  thresholdVector,
                                                  "alleviating", 
                                                  3)

remaining_data_alleviating_2 <- gs_IsingSamples_alleviating2[!(names(gs_IsingSamples_alleviating2) %in% 
                                                     c("vio_self", "vio_others"))]
remaining_data_alleviating_3 <- gs_IsingSamples_alleviating3[!(names(gs_IsingSamples_alleviating3) %in% 
                                                     c("vio_self", "vio_others"))]

column_interpersonal_alleviating_2 <- lapply(remaining_data_alleviating_2, function(df) df[, 8])
column_self_alleviating_2 <- lapply(remaining_data_alleviating_2, function(df) df[, 13])
column_interpersonal_alleviating_3 <- lapply(remaining_data_alleviating_3, function(df) df[, 8])
column_self_alleviating_3 <- lapply(remaining_data_alleviating_3, function(df) df[, 13])

activation_rates_alleviating_others2 <- sapply(column_interpersonal_alleviating_2, mean)
activation_rates_alleviating_self2 <- sapply(column_self_alleviating_2,mean)
activation_rates_alleviating_others3 <- sapply(column_interpersonal_alleviating_3, mean)
activation_rates_alleviating_self3 <- sapply(column_self_alleviating_3 ,mean)

SE_self_alleviating_2 <- sapply(column_self_alleviating_2, function(x) {
  p <- mean(x)  
  n <- length(x)  
  sqrt((p * (1 - p)) / n)  
})
SE_self_alleviating_3 <- sapply(column_self_alleviating_3, function(x) {
  p <- mean(x) 
  n <- length(x)  
  sqrt((p * (1 - p)) / n)
})

SE_others_alleviating_2 <- sapply(column_interpersonal_alleviating_2, function(x) {
  p <- mean(x)  
  n <- length(x) 
  sqrt((p * (1 - p)) / n) 
})
SE_others_alleviating_3 <- sapply(column_interpersonal_alleviating_3, function(x) {
  p <- mean(x)  
  n <- length(x)  
  sqrt((p * (1 - p)) / n)  
})

CI_self_lower_alleviating_3 <- activation_rates_alleviating_self3 - 1.96 * SE_self_alleviating_3
CI_self_upper_alleviating_3 <- activation_rates_alleviating_self3 + 1.96 * SE_self_alleviating_3
CI_self_lower_alleviating_2 <- activation_rates_alleviating_self2 - 1.96 * SE_self_alleviating_2
CI_self_upper_alleviating_2 <- activation_rates_alleviating_self2 + 1.96 * SE_self_alleviating_2

CI_others_lower_alleviating_3 <- activation_rates_alleviating_others3 - 1.96 * SE_others_alleviating_3
CI_others_upper_alleviating_3 <- activation_rates_alleviating_others3 + 1.96 * SE_others_alleviating_3
CI_others_lower_alleviating_2 <- activation_rates_alleviating_others2 - 1.96 * SE_others_alleviating_2
CI_others_upper_alleviating_2 <- activation_rates_alleviating_others2 + 1.96 * SE_others_alleviating_2

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
  transform(df_self_alleviating_2, Intervention = 2), 
  transform(df_self_alleviating_3, Intervention = 3)   
)
df_alleviating_others <- rbind(
  transform(df_others_alleviating_2, Intervention = 2),  
  transform(df_others_alleviating_3, Intervention = 3)   
)

df_alleviating_self$Variable <- factor(df_alleviating_self$Variable, 
                                       levels = c("original", 
                                                  setdiff(unique(df_alleviating_self$Variable), "original")))
df_alleviating_others$Variable <- factor(df_alleviating_others$Variable, 
                                         levels = c("original", 
                                                    setdiff(unique(df_alleviating_others$Variable), "original")))
df_alleviating_self$Variable <- factor(df_alleviating_self$Variable,
                                       levels = levels(df_alleviating_self$Variable),labels = new_labels)
df_alleviating_others$Variable <- factor(df_alleviating_others$Variable,
                                         levels = levels(df_alleviating_others$Variable),labels = new_labels)


p_alleviating_self <- ggplot(df_alleviating_self, aes(x = Variable, y = Mean, color = as.factor(Intervention))) +
  geom_point(position = position_dodge(width = 0.5)) +  
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, position = position_dodge(width = 0.7)) +
  labs(x = "Variable", 
       y = "Activation Rate ", 
       title = "Simulation of Activation Rate under Alleviating Intervention(Self)",
       color = "Intervention Intensity") +
  scale_color_manual(values = c("2" = "#6fa3e3",       
                                "3" = "#f46d6d")) +     
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), 
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 12))

p_alleviating_others <- ggplot(df_alleviating_others, aes(x = Variable, y = Mean, color = as.factor(Intervention))) +
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, position = position_dodge(width = 0.5)) +
  labs(x = "Variable", 
       y = "Activation Rate ", 
       title = "Simulation of Activation Rate under Alleviating Intervention(Others)",
       color = "Intervention Intensity") +
  scale_color_manual(values = c("2" = "#6fa3e3",       
                                "3" = "#f46d6d")) +     
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),  
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 12))

print(p_alleviating_self)
print(p_alleviating_others)


