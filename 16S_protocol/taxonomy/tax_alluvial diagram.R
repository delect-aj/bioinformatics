# 加载包
p_list = c("BiocManager", "tidyverse", "ggalluvial", "phyloseq", "reshape2")
for(p in p_list){if (!requireNamespace(p)){install.packages(p)}
  library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}

tax_alluvial <- function(
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
    tax_table(psdata)= tax
    
    ##转化为表格
    Taxonomies <- psdata %>% psmelt() %>% rename('tax' = sprintf(j))
    # head(Taxonomies)
    Taxonomies$Abundance=Taxonomies$Abundance * 100
    
    Taxonomies <- Taxonomies %>% group_by(tax, Group) %>% summarise(Abundance = mean(Abundance)) %>%
      dcast(Group ~ tax) %>% column_to_rownames("Group") %>% apply(1, function(x) x / sum(x) * 100) %>%
      as.data.frame() %>% rownames_to_column("tax") %>% melt(id.vars = "tax", variable.name = "Group", value.name = "Abundance")
    # 筛选需要的组别，并指定设定的顺序
    Taxonomies <- Taxonomies %>% filter(Group %in% list_group) %>%
                  mutate(Group, factor(Group, levels = list_group))
    
    p <- ggplot(Taxonomies, aes(x = Group, y = Abundance, fill = tax,
                                stratum = tax, alluvium = tax)) +
         geom_stratum() +  #代替 geom_col() 绘制堆叠柱形图
         geom_flow(alpha = 0.5) +
         theme_bw() +
         scale_fill_brewer(palette="Set3") +  #填充颜色赋值
         labs(x = '', y = 'Relative Abundance(%)') +
         theme(panel.grid = element_blank(), strip.text = element_text(size = 12),
               panel.background = element_rect(color = 'black', fill = 'transparent')) +
         theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13),
               legend.title = element_blank(), legend.text = element_text(size = 11))
    
    return(p)
}

otutab <- read.table("data/otutab.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors=F)
metadata <- read.table("data/metadata.txt", header=T, row.names=1, sep="\t", 
                       comment.char="", stringsAsFactors=F)
taxonomy=read.table("data/taxonomy.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors=F)

p <- tax_alluvial(otu = otutab, map = metadata, tax = taxonomy, list_group = c("OE", "KO"),
                  j = "Phylum", Top = 10)

zoom=1.5 # 控制图片缩放比例
ggsave(paste0("tax_alluvial.jpg"), p, width=89*zoom, height=56*zoom, units="mm")