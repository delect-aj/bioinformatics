tax_stack_clust <- function(
    otu=NULL,
    map=NULL,
    tax=NULL,
    dist="bray",
    Group="Group",
    j="Phylum", # 使用门水平绘制丰度图表
    rep=6 ,# 重复数量是6个
    Top=10, # 提取丰度前十的物种注释
    tran=TRUE, # 转化为相对丰度值
    hcluter_method="complete",
    cuttree=3){
  
  # 加载默认参数调试函数
  # otu=otutab
  # map=metadata
  # tax=taxonomy
  # dist="bray"
  # Group="Group"
  # j="Phylum" # 使用门水平绘制丰度图表
  # rep=6 # 重复数量是6个
  # Top=10 # 提取丰度前十的物种注释
  # tran=TRUE # 转化为相对丰度值
  # hcluter_method="complete"
  # cuttree=3
  
  # 加载包
  p_list = c("ggplot2", "BiocManager", "tidyverse", "ggdendro", "ggstance")
  for(p in p_list){if (!requireNamespace(p)){install.packages(p)}
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}
  p_list = c("ggtree", "phyloseq")
  for(p in p_list){if (!requireNamespace(p)){install.packages(p)}
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}
  
  # phyloseq导出特征表函数
  vegan_otu = function(physeq){
    OTU= otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU= t(OTU)
    }
    return(as(OTU,"matrix"))
  }
  
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
  
  idx = rownames(map) %in% colnames(otu)
  map = map[idx, , drop = F]
  otu = otu[, rownames(map)]
  
  # 分组列重命名为Group
  map <- map[Group]
  colnames(map)="Group"
  map$ID=row.names(map)
  
  # 数据导入phyloseq
  ps=phyloseq(sample_data(map),otu_table(as.matrix(otu), taxa_are_rows=TRUE), tax_table(as.matrix(tax)))
  # phyloseq(ps)对象标准化
  ps1_rela=phyloseq::transform_sample_counts(ps, function(x) x / sum(x) )
  
  # 导出OTU表
  otu=as.data.frame(t(vegan_otu(ps1_rela)))
  
  #计算距离矩阵
  unif=phyloseq::distance(ps1_rela , method=dist)
  # 聚类树，method默认为complete
  hc <- stats::hclust(unif, method=hcluter_method )
  
  #  take grouping with hcluster tree
  clus <- cutree(hc, cuttree)
  # 提取树中分组的标签和分组编号
  d=data.frame(label=names(clus),
               member=factor(clus))
  # Extract mapping file
  map=as.data.frame(sample_data(ps))
  # 合并树信息到样本元数据
  dd=merge(d,map,by="row.names",all=F)
  row.names(dd)=dd$Row.names
  dd$Row.names=NULL
  
  # ggtree绘图 #----
  p =ggtree(hc) %<+% dd +
    geom_tippoint(size=5, shape=21, aes(fill= Group, x=x)) +
    geom_tiplab(aes(color=Group,x=x * 1.2), hjust=1)
  
  # 按照分类学门(j)合并
  psdata=ps1_rela %>% tax_glom(taxrank=j)
  
  # 转化丰度值
  if (tran == TRUE) {
    psdata=psdata%>% transform_sample_counts(function(x) {x/sum(x)} )
  }
  
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
  tax_table(psdata)= tax
  
  ##转化为表格
  Taxonomies <- psdata %>% psmelt()
  # head(Taxonomies)
  Taxonomies$Abundance=Taxonomies$Abundance * 100
  
  Taxonomies$OTU=NULL
  colnames(Taxonomies)[1]="id"
  # head(Taxonomies)
  
  p <- p + ggnewscale::new_scale_fill()
  p1 <- facet_plot(p, panel='Stacked Barplot', data=Taxonomies, geom=geom_barh,mapping=aes(x=Abundance, fill=Phylum),color="black",stat='identity' )
  
  grotax <- Taxonomies %>%
    group_by(Group,Phylum) %>%
    summarise(Abundance=mean(Abundance))
  
  #--绘制分组的聚类结果
  ps1_rela=phyloseq::transform_sample_counts(ps, function(x) x / sum(x) )
  #按照分组合并OTU表格
  hc=phyloseq::merge_samples(ps1_rela, "Group",fun=mean) %>%
    phyloseq::distance(method=dist) %>%
    stats::hclust( method=hcluter_method )
  
  #  take grouping with hcluster tree
  clus <- cutree(hc, cuttree)
  # 提取树中分组的标签和分组编号
  d=data.frame(label=names(clus), member=factor(clus))
  # Extract mapping file
  map=as.data.frame(sample_data(phyloseq::merge_samples(ps1_rela, "Group",fun=mean)))
  # 合并树信息到样本元数据
  dd=merge(d,map,by="row.names",all=F)
  row.names(dd)=dd$Row.names
  dd$Row.names=NULL
  
  # ggtree绘图 #----
  p3 =ggtree(hc) %<+% dd +
    geom_tippoint(size=5, shape=21, aes(fill= member, x=x)) +
    geom_tiplab(aes(color=member,x=x * 1.2), hjust=1)
  p3 <- p3 + ggnewscale::new_scale_fill()
  p4 <- facet_plot(p3, panel='Stacked Barplot', data=grotax, geom=geom_barh,mapping=aes(x=Abundance, fill=Phylum),color="black",stat='identity' )
  
  return(list(p,p1,p3,p4,Taxonomies))
}
