### 距离统计检验
MicroTest = function(otu = NULL, map = NULL, ps = ps,
                     group = "Group", Micromet = "MRPP", dist = "bray"){
  library(phyloseq)
  # otu = otutab_rare
  # map = metadata
  # ps = ps
  # group = "Group"
  # Micromet = "MRPP"
  # dist = "bray"
  
  # dist_methods = unlist(phyloseq::distanceMethodList)
  # dist = dist_methods[dist]
  
  # # 数据导入 otutab 和 map文件
  if (is.null(otu)&is.null(map)) {
    ps = ps
  }else {
    # 数据导入PhyloSeq格式
    if (is.null(otu)&is.null(map)) {
      ps = ps
    }else{
      otu = as.matrix(otu)
      map$Group = as.factor(map[, group])
      ps = phyloseq(otu_table(otu, taxa_are_rows=TRUE),sample_data(map))
    }
    # 只有使用树相关的距离算法时，读取树
    if (dist %in% c("unifrac" , "wunifrac",  "dpcoa")) {
      phy_tree(ps) = tree
    }
  }
  if (dist %in% c("unifrac" , "wunifrac",  "dpcoa")) {
    phy_tree(ps) = tree
  }
  
  # 求取相对丰度
  ps1_rela  = transform_sample_counts(ps, function(x) x / sum(x) );ps1_rela
  
  library(vegan)
  # -准备矩阵和分组文件
  map = as.data.frame(sample_data(ps1_rela))
  # ?distance
  # unif= distance(ps1_rela , method=dist)
  # unif = vegdist(t(otu), methodC="bray")
  # unif = distance(otu, method=dist)
  
  
  unif = phyloseq::distance(ps, method=dist)
  
  # adonis#----
  if (Micromet == "adonis") {
    ado =  vegan::adonis(unif ~ map$Group,permutations = 999)
    a = round(as.data.frame(ado$aov.tab[5])[1,1],3)
    R2 = paste("adonis:R ",a, sep = "")
    b = as.data.frame(ado$aov.tab[6])[1,1]
    p_v = paste("p: ",b, sep = "")
    title1 = paste(R2," ",p_v, sep = "")
    title1
  }
  
  # MRPP#----
  if (Micromet == "MRPP") {
    dat.mrpp = mrpp(unif, map$Group)
    a = round(dat.mrpp$delta,3)
    R2 = paste("MRPP.delta ",a, sep = "")
    p_v = paste("p: ",round(dat.mrpp$Pvalue,3), sep = "")
    title1 = paste(R2," ",p_v, sep = "")
    title1
  }
  
  # anosim#----
  if (Micromet == "anosim") {
    dat.ano = anosim(unif, map$Group)
    a = round(dat.ano$statistic,3)
    R2 = paste("ANOSIM.r ",a, sep = "")
    p_v = paste("p: ",round(dat.ano$signif,3), sep = "")
    title1 = paste(R2," ",p_v, sep = "")
    title1
  }
  return(title1)
}

### 两两比较的统计检验
pairMicroTest=function(ps=ps, Micromet="anosim", dist="bray"){
  
  if (!requireNamespace("vegan", quietly=TRUE))
    install.packages("vegan")
  library(vegan)
  # 安装Bioconductor的R包phyloseq
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
  suppressWarnings(suppressMessages(library(BiocManager)))
  if (!requireNamespace("phyloseq", quietly=TRUE))
    BiocManager::install("phyloseq")
  library(phyloseq)
  
  # 生成phyloseq对象
  # ps=phyloseq(otu_table(otutab_rare, taxa_are_rows=TRUE),
  #               sample_data(metadata))
  ps1_rela =transform_sample_counts(ps, function(x) x / sum(x) );ps
  
  #-准备矩阵和分组文件
  map=as.data.frame(sample_data(ps1_rela))# ,stringsAsFactors=T
  # 转换为后，levels消失了，改用unique添加levels
  map$Group=factor(map$Group, levels=unique(map$Group))
  (aa=levels(map$Group))
  (aaa=combn(aa,2))
  dim(aaa)[2]
  
  # 构建三个空列
  ID=rep("a",dim(aaa)[2])
  R=rep("a",dim(aaa)[2])
  P=rep("a",dim(aaa)[2])
  # i=1
  for (i in 1:dim(aaa)[2]) {
    # print(i)
    Desep_group=aaa[,i]
    map=as.data.frame(sample_data(ps1_rela))
    # head(map)
    map$ID=row.names(map)
    # maps=dplyr::filter(map, Group %in% Desep_group)
    # 取子集 #----
    maps=subset(map, Group %in% Desep_group)
    row.names(maps)=maps$ID
    ps_sub=ps1_rela
    sample_data(ps_sub)=maps
    ps_sub=phyloseq::filter_taxa(ps_sub, function(x) sum(x ) > 0 , TRUE);ps_sub
    map=as.data.frame(sample_data(ps_sub))
    unif <- phyloseq::distance(ps_sub, method=dist)
    
    if (Micromet == "MRPP") {
      mrpp=vegan::mrpp(unif, map$Group)
      as1=round(mrpp$delta,3)
      R2 <- paste("MRPP.delta ",as1, sep="")
      # R[i]=R2
      R2
      p_v=paste("p: ",round(mrpp$Pvalue,3), sep="")
      # p_v
      # P[i]=p_v
    }
    
    if (Micromet == "anosim") {
      dat.ano=anosim(unif, map$Group)
      a=round(dat.ano$statistic,3)
      R2 <- paste("ANOSIM.r ",a, sep="")
      R[i]=R2
      p_v=paste("p: ",round(dat.ano$signif,3), sep="")
      # P[i]=p_v
    }
    
    gg =map$Group
    if (Micromet == "adonis") {
      ado= adonis(unif~gg,permutations=999)
      a=round(as.data.frame(ado$aov.tab[5])[1,1],3)
      R2 <- paste("adonis:R ",a, sep="")
      R[i]=R2
      b=as.data.frame(ado$aov.tab[6])[1,1]
      p_v=paste("p: ",b, sep="")
    }
    ID[i]=paste(Desep_group[1],Desep_group[2],sep="_VS_")
    P[i]=p_v
    R[i]=R2
  }
  # P
  # R
  result=data.frame(ID=ID,stat=R,p=P)
  result
  
  return(result)
}


BetaDiv=function(otu=otutab, map=metadata, tree=tree, ps=NULL,
                 group="Group", dist="bray", method="NMDS",
                 Micromet="adonis", pvalue.cutoff=0.05){
  
  # 需要的R包
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  
  # # 读取默认参数
  # otu=otutab
  # map=metadata
  # tree=tree
  # ps=NULL
  # group="Group"
  # dist="bray"
  # method ="NMDS"
  # Micromet="adonis"
  # pvalue.cutoff=0.05
  
  # 数据导入PhyloSeq#----
  if (is.null(otu)&is.null(map)) {
    ps=ps
  }else{
    otu=as.matrix(otu)
    map$Group=as.factor(map[, group])
    ps=phyloseq(otu_table(otu, taxa_are_rows=TRUE), sample_data(map))
  }
  
  # 读取树#----
  if (dist %in% c("unifrac", "wunifrac", "dpcoa")) {
    phy_tree(ps)=tree
  }
  
  #转换为相对丰度#----
  ps1_rela=transform_sample_counts(ps, function(x) x / sum(x) )
  
  #提取标准化的特征表#----
  vegan_otu= function(physeq){
    OTU=otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU=t(OTU)
    }
    return(as(OTU,"matrix"))
  }
  otu_table=as.data.frame(t(vegan_otu(ps1_rela )))
  
  #排序方法选择#----
  
  #---------DCA排序#----
  if (method == "DCA") {
    ordi=phyloseq::ordinate(ps1_rela, method=method, distance=dist)
    #提取样本坐标
    points=ordi$rproj[,1:2]
    #命名行名
    colnames(points)=c("x", "y")
    #提取特征值
    eig=ordi$evals^2
  }
  
  #---------CCA排序#----
  if (method == "CCA") {
    ordi=ordinate(ps1_rela, method=method, distance=dist)
    #样本坐标,这里可选u或者v矩阵
    points=ordi$CA$v[,1:2]
    colnames(points)=c("x", "y")
    eig=ordi$CA$eig^2
  }
  
  #---------RDA排序#----
  if (method == "RDA") {
    ordi=ordinate(ps1_rela, method=method, distance=dist)
    #样本坐标,这里可选u或者v矩阵
    points=ordi$CA$v[,1:2]
    colnames(points)=c("x", "y") #命名行名
    #提取特征值
    eig=ordi$CA$eig
  }
  
  #---------MDS排序#----
  if (method == "MDS") {
    ordi=ordinate(ps1_rela, method=method, distance=dist)
    points=ordi$vectors[,1:2]
    colnames(points)=c("x", "y")
    eig=ordi$values[,1]
  }
  
  #---------PCoA排序#----
  if (method == "PCoA") {
    unif=phyloseq::distance(ps1_rela , method=dist, type="samples")
    pcoa=cmdscale(unif, k=2, eig=T)
    # 获得坐标点get coordinate string, format to dataframme
    points=as.data.frame(pcoa$points)
    colnames(points)=c("x", "y")
    eig=pcoa$eig
  }
  
  #---------PCA排序#----
  if (method == "PCA") {
    otu.pca=prcomp(t(otu_table), scale.=TRUE)
    points=otu.pca$x[,1:2]
    colnames(points)=c("x", "y")
    # #提取荷载坐标
    # otu.pca$rotation
    # 提取解释度,这里提供的并不是特征值而是标准差，需要求其平方才是特征值
    eig=otu.pca$sdev
    eig=eig*eig
  }
  
  #---------------LDA排序#----
  if (method == "LDA") {
    #拟合模型
    library(MASS)
    data=t(otu_table)
    data=as.data.frame(data)
    data=scale(data, center=TRUE, scale=TRUE)
    dim(data)
    data1=data[,1:10]
    map=as.data.frame(sample_data(ps1_rela))
    model=lda(data, map$Group)
    # 提取坐标
    ord_in=model
    axes=c(1:2)
    points=data.frame(predict(ord_in)$x[, axes])
    colnames(points)=c("x", "y") #命名行名
    # 提取解释度
    eig= ord_in$svd^2
  }
  
  #---------------NMDS排序#----
  if (method == "NMDS") {
    ordi=ordinate(ps1_rela, method=method, distance=dist)
    #样本坐标,
    points=ordi$points[,1:2]
    colnames(points)=c("x", "y")
    #提取stress
    stress=ordi$stress
    stress=paste("stress", ": ", round(stress,2), sep="")
  }
  
  
  #---------------t-sne排序#----
  # method="t-sne"
  if (method == "t-sne") {
    data=as.data.frame(t(otu_table))
    data=scale(data, center=TRUE, scale=TRUE)
    map=as.data.frame(sample_data(ps1_rela))
    if (!requireNamespace("Rtsne", quietly=TRUE))
      install.packages("Rtsne")
    library(Rtsne)
    tsne=Rtsne(data,perplexity=3)
    
    # 提取坐标
    points=as.data.frame(tsne$Y)
    row.names(points)= row.names(map)
    colnames(points)=c("x", "y") #命名行名
    stress= NULL
  }
  
  
  #----差异分析#----
  
  #----整体差异分析#----
  title1=MicroTest(ps=ps1_rela, Micromet=Micromet, dist=dist)
  
  #----两两比较#----
  pairResult=pairMicroTest(ps=ps1_rela, Micromet=Micromet, dist=dist)
  
  
  #----绘图#----
  map=as.data.frame(sample_data(ps1_rela))
  map$Group=as.factor(map$Group)
  colbar=length(levels(map$Group))
  
  points=cbind(points, map[match(rownames(points), rownames(map)), ])
  # head(points)
  points$ID=row.names(points)
  
  #---定义配色#----
  # 改为使用Rcolorbrewer的颜色组合
  mi=colorRampPalette(c( "#CBD588", "#599861", "orange","#DA5724", "#508578", "#CD9BCD","#AD6F3B", "#673770","#D14285", "#652926", "#C84248","#8569D5", "#5E738F","#D1A33D", "#8A7C64","black"))(colbar)
  
  #----定义图形通用样式#----
  main_theme=theme(panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   plot.title=element_text(vjust=-8.5,hjust=0.1),
                   text=element_text(family="sans", size=7)
                   # axis.title.y =element_text(size=7,face="bold",colour="black"),
                   # axis.title.x =element_text(size=7,face="bold",colour="black"),
                   # axis.text=element_text(size=7,face="bold"),
                   # axis.text.x=element_text(colour="black",size=7),
                   # axis.text.y=element_text(colour="black",size=7),
                   # legend.text=element_text(size=7,face="bold")
  )
  
  if (method %in% c("DCA", "CCA", "RDA",  "MDS", "PCoA","PCA","LDA")) {
    p2 =ggplot(points, aes(x=x, y=y, fill=Group)) +
      geom_point(alpha=.7, size=2, pch=21) +
      labs(x=paste0(method," 1 (",format(100*eig[1]/sum(eig),digits=4),"%)"),
           y=paste0(method," 2 (",format(100*eig[2]/sum(eig),digits=4),"%)"),
           title=title1) +
      stat_ellipse(linetype=2,level=0.68,aes(group=Group, colour=Group))+
      scale_colour_manual(values=mi,guide=guide_legend(title=NULL))+
      scale_fill_manual(values=mi,guide=guide_legend(title=NULL))+
      guides(color=guide_legend(title=NULL),shape=guide_legend(title=NULL))
    p2=p2+theme_bw()+
      geom_hline(aes(yintercept=0), colour="black", linetype=2) +
      geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
      main_theme
    library(ggrepel)
    p3=p2+geom_text_repel(aes(label=points$ID),size=2)
  }
  
  
  if (method %in% c("NMDS","t-sne")) {
    p2 =ggplot(points, aes(x=x, y=y, fill=Group)) +
      geom_point(alpha=.7, size=2, pch=21) +
      labs(x=paste(method,"1", sep=""),
           y=paste(method,"2",sep=""),
           title=stress)+
      stat_ellipse( linetype=2,level=0.65,aes(group  =Group, colour= Group))+
      scale_colour_manual(values=mi,guide=guide_legend(title=NULL))+
      scale_fill_manual(values=mi,guide=guide_legend(title=NULL))+
      guides(color=guide_legend(title=NULL),shape=guide_legend(title=NULL))
    p2=p2+theme_bw()+
      geom_hline(aes(yintercept=0), colour="black", linetype=2) +
      geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
      main_theme
    library(ggrepel)
    p3=p2+geom_text_repel(aes(label=points$ID),size=2)
    if (method %in% c("t-sne")) {
      supp_lab=labs(x=paste(method,"1", sep=""),y=paste(method,"2",sep=""),title=title)
      p2=p2 + supp_lab
      p3=p3 + supp_lab
    }
    p2
  }
  
  # 返回结果：标准图，数据，标签图，成对比较结果，整体结果
  return(list(p2,points,p3,pairResult,title1))
  
}
library(vegan)

otutab <- read.table("data/otutab.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors=F)
metadata <- read.table("data/metadata.txt", header=T, row.names=1, sep="\t", 
                       comment.char="", stringsAsFactors=F)

### OTU数据抽平
otutab = as.data.frame(t(rrarefy(t(otutab), min(colSums(otutab)))))
# 输入抽平标准化的特征表、元数据、分组列名、距离类型、降维和统计方法
result=BetaDiv(otu=otutab, map=metadata, group="Group",
               dist="bray", method="NMDS", Micromet="adonis")
# 返回结果列表：标准图，数据，标签图，成对比较结果，整体结果

#提取排序散点图(结果列表中的1)
p <- result[[1]]
