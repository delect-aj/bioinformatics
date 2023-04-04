ternary_plot <- function(
    otu=NULL,
    map=NULL,
    tax=NULL,
    list_group=NULL, # 指定分组名称，vector
    Group="Group",
    j="Phylum", # 使用门水平绘制丰度图表
    Top=10){
    # phyloseq导出特征表函数
    vegan_otu = function(physeq){
      OTU= otu_table(physeq)
      if(taxa_are_rows(OTU)){
        OTU= t(OTU)
      }
      return(as(OTU,"matrix"))
    }
    
    p_list = c("BiocManager", "tidyverse", "ggtern", "phyloseq", "reshape2")
    for(p in p_list){if (!requireNamespace(p)){install.packages(p)}
      library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}
    
    # phyloseq导出物种注释函数
    vegan_tax <-  function(physeq){
      tax <-  tax_table(physeq)
      return(as(tax,"matrix"))
    }
    
    # 数据交叉筛选
    # Extract only those ID in common between the two tables
    idx=rownames(otu) %in% rownames(tax)
    otu=otu[idx,]
    tax=tax[rownames(otu),]
    
    # 分组列重命名为Group
    map = map[Group]
    colnames(map)="Group"
    map$ID=row.names(map)
    
    # 数据导入phyloseq
    ps=phyloseq(sample_data(map),otu_table(as.matrix(otu), taxa_are_rows=TRUE), tax_table(as.matrix(tax)))
    # phyloseq(ps)对象标准化
    ps1_rela=phyloseq::transform_sample_counts(ps, function(x) x / sum(x) )
    
    # 导出OTU表
    otu=as.data.frame(t(vegan_otu(ps1_rela)))
    
    # 按照分类学门(j)合并
    psdata=ps1_rela %>% tax_glom(taxrank=j)
    
    # 转化丰度值
    # if (tran == TRUE) {
    #   psdata=psdata%>% transform_sample_counts(function(x) {x/sum(x)} )
    # }
    
    #--提取otu和物种注释表格
    otu=otu_table(psdata)
    tax=tax_table(psdata)
    
    #--按照指定的Top数量进行筛选与合并
    for (i in 1:dim(tax)[1]) {
      if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing=TRUE)[1:Top])) {
        tax[i,j] =tax[i,j]
      } else {
        tax[i,j]= "Other"
      }
    }
    tax_table(psdata) = tax
    
    ##转化为表格
    Taxonomies <- psdata %>% psmelt() %>% rename('tax' = sprintf(j))
    # head(Taxonomies)
    Taxonomies$Abundance=Taxonomies$Abundance * 100
    
    Taxonomies <- Taxonomies %>% group_by(tax, Group) %>% summarise(Abundance = mean(Abundance)) %>%
      dcast(Group ~ tax) %>% column_to_rownames("Group") %>% apply(1, function(x) x / sum(x) * 100) %>%
      as.data.frame() %>% rownames_to_column("tax") 
    
    Taxonomies <- Taxonomies[, c("tax", list_group)]
    colnames(Taxonomies) <- c("tax", "x", "y", "z")
    p <- ggtern(data=Taxonomies,
                aes(x=x, y=y, z=z)) + 
      labs(x = list_group[1],y = list_group[2],z = list_group[3], xarrow = list_group[1],
           yarrow = list_group[2], zarrow = list_group[3], color = "Taxonomy")+
      geom_mask() + # 可将超出边界的点正常显示出来
      geom_point(aes(color=tax),alpha=0.8, size = 5) +
      theme_bvbw()
    return(p)
}
