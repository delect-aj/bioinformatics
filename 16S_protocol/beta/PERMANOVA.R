library(vegan)
metadata <- read.table("data/metadata.txt", header=T, row.names=1, sep="\t", 
                       comment.char="", stringsAsFactors=F)

# 主坐标轴分析，可选距离矩阵bray_curtis、unifrac、unifrac_binary、jaccard、manhatten、euclidean
# 设置距离矩阵类似，常用bray_curtis或unifrac
distance_type = "bray_curtis"
# Data reading
distance_mat = read.table(paste0("data/", distance_type,".txt"), header=T, row.names=1, sep="\t", comment.char="")

# 仅查看分组对群落的解释
adonis_var <- adonis (as.dist(distance_mat) ~ Group, data = metadata, by=NULL, parallel=4)
adonis_var$aov.tab


write.table(adonis_var$aov.tab, file="beta/adonis_var_group.txt", append = F, sep="\t", quote=F, row.names=T, col.names=T)

# 计算Group, Date, Site三类变量及交互对群落结构差异的解释
adonis_var <- adonis (as.dist(distance_mat) ~ Group*Date*Site, data = metadata, by=NULL, parallel=4)
adonis_var$aov.tab

# 变量不能重合，否则无效。如本次Date与Site一致，所以Site无结果
# 计算Group和Site的解释率
adonis_var <- adonis(as.dist(distance_mat) ~ Group*Site, data = metadata, by=NULL, parallel=4)
adonis_var$aov.tab

write.table(adonis_var$aov.tab, file="beta/adonis_var.txt", append = F, sep="\t", quote=F, row.names=T, col.names=T)