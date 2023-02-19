heaatmap_diff <- function(
    otu,
    tax,
    group_df, ## mwtadata
    group # 分组信息
    ){
    if (!require("BiocManager"))
      install.packages('BiocManager') 
    if (!require("ComplexHeatmap"))
      BiocManager::install('ComplexHeatmap') 
    # 加载包
    library(ComplexHeatmap)
    library(RColorBrewer)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    #处理后有73种差异还比较明显的颜色，基本够用
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
    
    colnames(otu)[1] <- "OTU"
    #修改物种注释信息
    colnames(tax)[1] <- "OTU"
    tax <- tax %>% mutate(Species = paste0(Species, "(", OTU, ")", sep = ""))
    df <- merge(tax, otu)
    
    group_df <- group_df[group]
    colnames(group_df) <- "Group"
    #选取列
    df <- df[, c("Class", "Family", "Species", rownames(group_df))]
    #按Class和Family列对df的所有行进行排序
    df<-df[order(df$Class,df$Family),]
    #获得相对丰度矩阵
    mat<-df[-c(1,2,3)]
    mat<-as.matrix(mat) #将数据框转化为矩阵
    row.names(mat)<-df[,3]
    #获得class-family-species对应表
    df_annotation<-df[c(1,2,3)]
    
    #生成列分组注释
    top_col <- sample(col_vector, length(unique(group_df[, "Group"])))
    names(top_col) <- unique(group_df[, "Group"])
    top_anno<-HeatmapAnnotation(State = group_df[, "Group"], 
                                col = list(State = top_col), #设置颜色
                                annotation_legend_param = list(State = list(at=unique(group_df[, "Group"])))) #设置图例顺序
    
    
    #生成行分组注释
    left_col <- sample(col_vector, length(unique(df[, "Family"]))) #生成行分组颜色的具名向量
    names(left_col) <- unique(df[, "Family"])
    left_anno<-rowAnnotation(Family = df_annotation[,2],
                             show_annotation_name = F, #不展示行分组注释名
                             col = list(Family =left_col),
                             annotation_legend_param = list(Family = list(at=unique(df[, "Family"]))))
    
    
    ht<-Heatmap(mat,
                cluster_rows = F,  #不按行聚类
                show_column_names = F, #不展示列名
                heatmap_legend_param = list(title = "Log2 relative abundance"), #设置热图图例名称
                col = c("#FFFFFF","#D32F2F"), #设置热图颜色
                top_annotation = top_anno, #添加列注释
                left_annotation = left_anno) #添加行注释
    
    # class<-unique(df_annotation$Class)
    
    
    #绘制class单独的图例
    # lgd<-Legend(labels = class, title = "Class",
    #             legend_gp = gpar(fill = c("#1976D2","#F8BBD0","#D32F2F","#4CAF50","#FFEB3B","#673AB7")))
    
    #将单独图例与其他部分合并
    # pdf("Figure 1B.pdf",width = 11, height = 10)
    # draw(ht, ht_gap = unit(7, "mm"),annotation_legend_list = lgd)
    # dev.off()
    return(ht)
}

#读取相对丰度矩阵及class-family-species对应表
otu<-read.table("data/otutab.txt", sep = "\t", header = T, comment.char="", stringsAsFactors=F)
otu <- otu[sample(nrow(otu), 30), ]

tax <- read.table("data/taxonomy.txt", sep = "\t", header = T)

#读取分组数据表
group_df <- read.table("data/metadata.txt", header = T, sep = "\t", row.names = 1)

res <- heaatmap_diff(otu, tax, group_df, "Group")
