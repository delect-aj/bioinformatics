library(vegan)
library(tidyverse)

metadata <- read.table("data/metadata.txt", header=T, sep="\t", 
                       comment.char="", stringsAsFactors=F)

# 主坐标轴分析，可选距离矩阵bray_curtis、unifrac、unifrac_binary、jaccard、manhatten、euclidean
# 设置距离矩阵类似，常用bray_curtis或unifrac
distance_type <- "bray_curtis"
# Data reading
distance_mat <- read.table(paste0("data/", distance_type,".txt"), 
                           header=T, sep="\t", 
                           comment.char="")
# 重命名第一列
colnames(distance_mat)[1] <- "SampleID"

meta_distance <- inner_join(distance_mat, metadata, by="SampleID")

all_dist <- meta_distance %>% 
  select(all_of(.[["SampleID"]])) %>%
  as.dist()

set.seed(1234)
all_test <- adonis2(all_dist ~ Group, 
                   data=meta_distance,
                   permutations = 1000)
all_test$`Pr(>F)`[1]

# 两组之间的比较
comparison_adonis <- function(meta_distance, group_list){
    pariwise_p <- numeric()
    comp_list <- t(combn(group_list, 2))
    for (i in 1 : nrow(comp_list)){
        group1 <- comp_list[i, 1]
        group2 <- comp_list[i, 2]
        meta_twogroup <- meta_distance %>%
          filter(Group == paste0(group1) | Group == paste(group2))
        dis_twogroup <- meta_twogroup %>%
          select(all_of(.[["SampleID"]])) %>%
          as.dist()
        
        twogroup_test <- adonis2(dis_twogroup ~ Group, 
                           data=meta_twogroup,
                           permutations = 1000)
        pariwise_p[paste0(group1, "_", group2)] <- twogroup_test$`Pr(>F)`[1]
    }
    pariwise_p <- p.adjust(pariwise_p, method = "BH")
    return(pariwise_p)
}
comparison_adonis(meta_distance, group_list)
