beta_pcoa <- function(dis_mat, metadata, groupID="Group", ellipse=T, label=F, PCo=12) {
  # 依赖关系检测与安装
  p_list=c("ggplot2", "vegan", "ggrepel")
  for(p in p_list){
    if (!requireNamespace(p)){
      install.packages(p)}
    suppressWarnings(suppressMessages(library(p,character.only=T)))}
  
  # 测试默认参数
  # dis_mat=beta_unifrac
  # metadata=metadata
  # groupID="Group"
  # ellipse=T
  # label=F
  # PCo=12
  
  # 交叉筛选
  idx=rownames(metadata) %in% rownames(dis_mat)
  metadata=metadata[idx,,drop=F]
  dis_mat=dis_mat[rownames(metadata), rownames(metadata)]
  
  # 提取样品组信息,默认为group可指定
  sampFile=as.data.frame(metadata[, groupID],row.names=row.names(metadata))
  # colnames(sampFile)[1]="group"
  
  # PCoA
  pcoa=cmdscale(dis_mat, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
  points=as.data.frame(pcoa$points) # get coordinate string, format to dataframme
  eig=pcoa$eig
  points=cbind(points, sampFile[rownames(points),])
  colnames(points)=c("x", "y", "z","group")
  
  # 按1、2轴绘图
  if (PCo == 12){
    p=ggplot(points, aes(x=x, y=y, color=group))  +
      labs(x=paste("PCo 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
           y=paste("PCo 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""), color=groupID)
  }
  # 按1、3轴绘图
  if (PCo == 13){
    p=ggplot(points, aes(x=x, y=z, color=group))  +
      labs(x=paste("PCo 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
           y=paste("PCo 3 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""), color=groupID)
  }
  # 按2、3轴绘图
  if (PCo == 23){
    p=ggplot(points, aes(x=y, y=z, color=group))  +
      labs(x=paste("PCo 2 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
           y=paste("PCo 3 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""), color=groupID)
  }
  p=p + geom_point(alpha=.7, size=2) + theme_classic()
  # 是否添加置信椭圆
  if (ellipse == T){
    p=p + stat_ellipse(level=0.68)
  }
  # 是否显示样本标签
  if (label == T){
    p=p + geom_text_repel(label=paste(rownames(points)), colour="black", size=3)
  }
  p
}
