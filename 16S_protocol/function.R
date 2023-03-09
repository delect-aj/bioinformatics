library(tidyverse)
library(reshape2)
library(ggh4x)
library(patchwork) ### 图像拼接
library(RColorBrewer)
for (n in c("Ec", "KO", "pathway")){
    table <- read.table(sprintf("D:/document/老婆大人/16S分析/picrust2/function/%s.tsv", n), 
                        header = TRUE, sep = "\t", check.names = FALSE, ) %>% 
        column_to_rownames("function") %>% select(-description)
      
    metadata <- read.table("D:/document/老婆大人/16S分析/metadata.txt", header = TRUE)
    
    sum_ <- apply(table, 1, sum)
    table <- table[names(sort(sum_, decreasing = TRUE)), ] %>% apply(2, function(x) x/sum(x)) %>%
        as.data.frame() %>% rownames_to_column("description")
    names <- table$description[1 : 20]
    ## 颜色
    getPalette = colorRampPalette(brewer.pal(21, "Set3"))
    p1 <- table %>% mutate(description = if_else(description %in% names, description, "others")) %>%
        mutate(description = factor(description, levels = c(names, "others"))) %>% 
        melt(id.vars = "description", variable.name = "sample.id") %>% merge(metadata) %>% 
        ggplot() +
        geom_bar(aes(x = sample.id, y = value, fill = description), stat="identity") +   #"identity"
        scale_fill_manual(values = getPalette(21))+
        labs(x = "sample", y = "Relative Abundance", fill = "",size=20)+
        theme_bw()+
        theme( axis.text.y = element_text(color = "black",size = 20),
               axis.title.x=element_text(colour='black', size=20,face = "bold"),
               axis.text.x = element_text(color = "black",size = 20, angle = 90, vjust = 0.6, hjust = 1),
               text = element_text(size = 20),
               legend.title=element_text(size=20))+
        theme(axis.text.x = element_blank()) + 
        theme(legend.position = "bottom")
    pal <- c("#DDA0DD")
    p1 + facet_nested(.~Group, drop = T, scale = "free", space = "free", switch = "y",
                      strip = strip_nested(background_x = elem_list_rect(fill = pal), by_layer_x = T))+
        theme(panel.spacing.x = unit(0, "cm"))
    ggsave(sprintf("D:/document/老婆大人/16S分析/picrust2/function/%s/bar.png", n), width = 60, height = 30, units = "cm")
    
    ### PCA
    metadata <- read.table("D:/document/老婆大人/16S分析/metadata.txt", header = TRUE, row.names = 1)
    table <- read.table(sprintf("D:/document/老婆大人/16S分析/picrust2/function/%s.tsv", n), 
                        header = TRUE, sep = "\t", check.names = FALSE) %>% 
        column_to_rownames("function") %>% group_by(description) %>% 
        summarise_all(list(sum)) %>% column_to_rownames("description") %>%
        t()
    #特征分解
    dat_eigen <- scale(table ,scale=T)%>%cor()%>%eigen()
    #特征值提取
    eig <- dat_eigen$values 
    
    
    #主成分载荷表示各个主成分与原始变量的相关系数。
    #将中心化的变量矩阵得到每个观测值的得分
    result <- scale(table, scale=T) %*% dat_eigen$vectors
    plot_data <- data.frame(PC1 = result[,1], PC2 = result[,2], group = metadata[rownames(result), ])
    
    ggplot(plot_data, aes(x = PC1, y = PC2, color = group)) +
        geom_point(size = 3) +
        labs(x=paste("PCoA 1 (",format(100* eig[1] / sum(eig), digits=4), "%)", sep=""),
             y=paste("PCoA 2 (", format(100*eig[2] / sum(eig), digits=4), "%)", sep=""))+
        stat_ellipse(type = "norm", linetype = 2) +
        theme_bw() +
        scale_color_brewer(palette="Set3")
    
    ggsave(sprintf("D:/document/老婆大人/16S分析/picrust2/function/%s/PCA.png", n), width=20, height=18,unit="cm", dpi=400)
    
    ### 差异分析
    for (m in  c("CON", "NCUF201.1", "NCUF206.4", "NCUF207.7", "NCUF208.7", 
                 "NCUF211.1", "NCUF213.1", "NCUF214.1")){
        group <- c("HUA", sprintf("%s", m))
        metadata <- read.table("D:/document/老婆大人/16S分析/metadata.txt", header = TRUE) %>%
            filter(Group %in% group)
        table <- read.table(sprintf("D:/document/老婆大人/16S分析/picrust2/function/%s.tsv", n), 
                            header = TRUE, sep = "\t", check.names = FALSE, row.names = 1) %>% 
            select(-description) %>% apply(2, function(x) x / sum(x)) %>% as.data.frame()
        idy = colnames(table) %in% metadata$sample.id
        table <- table*100
        table <- table %>% filter(apply(table[, idy], 1, mean) > 0.1)
        table <- t(table) %>% as.data.frame() %>% rownames_to_column("sample.id") %>% merge(metadata)
        table$Group <- as.factor(table$Group)
        
        ## t-test
        diff <- table %>% 
            select_if(is.numeric) %>%
            map_df(~ broom::tidy(t.test(. ~ Group,data = table)), .id = 'var') ## t.test()
        
        diff$p.value <- p.adjust(diff$p.value,"bonferroni")
        diff <- diff %>% filter(p.value < 0.01) %>% mutate(diff = abs(estimate)) %>% arrange(desc(diff))
        write.csv(diff, sprintf("D:/document/老婆大人/16S分析/picrust2/function/%s/hua_vs_%s.csv", n, m))
        diff <- diff[1:20, ]
        
        ## 绘图数据构建
        ## 左侧条形图
        abun.bar <- table[,c(diff$var,"Group")] %>% 
            gather(variable,value,-Group) %>% 
            group_by(variable,Group) %>% 
            summarise(Mean = mean(value))
        
        ## 右侧散点图
        diff.mean <- diff[,c("var","estimate", "conf.low", "conf.high", "p.value")]
        diff.mean$Group <- c(ifelse(diff.mean$estimate >0, levels(table$Group)[1],
                                    levels(table$Group)[2]))
        diff.mean <- diff.mean[order(diff.mean$estimate, decreasing = TRUE),]
        
        ###绘图
        cbbPalette <- c("#E69F00", "#56B4E9")
        abun.bar$variable <- factor(abun.bar$variable,levels = rev(diff.mean$var))
        p1 <- ggplot(abun.bar, aes(variable,Mean, fill = Group)) +
            scale_x_discrete(limits = levels(diff.mean$var)) +
            coord_flip() +
            xlab("") +
            ylab("Mean proportion (%)") +
            theme(panel.background = element_rect(fill = 'transparent'),
                  panel.grid = element_blank(),
                  axis.ticks.length = unit(0.4,"lines"), 
                  axis.ticks = element_line(color='black'),
                  axis.line = element_line(colour = "black"),
                  axis.title.x=element_text(colour='black', size=12,face = "bold"),
                  axis.text=element_text(colour='black',size=10,face = "bold"),
                  legend.title=element_blank(),
                  legend.text=element_text(size=12,face = "bold",colour = "black",
                                           margin = margin(r = 20)),
                  legend.position = "bottom",
                  legend.direction = "horizontal",
                  legend.key.width = unit(0.8,"cm"),
                  legend.key.height = unit(0.5,"cm"))
        
        
        for (i in 1:(nrow(diff.mean) - 1)) 
            p1 <- p1 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf, 
                                fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
        
        p1 <- p1 + 
            geom_bar(stat = "identity",position = "dodge",width = 0.7, colour = "black") +
            scale_fill_manual(values=cbbPalette)
        
        
        ## 右侧散点图
        diff.mean$var <- factor(diff.mean$var,levels = levels(abun.bar$variable))
        diff.mean$p.value <- signif(diff.mean$p.value,3)
        diff.mean$p.value <- as.character(diff.mean$p.value)
        p2 <- ggplot(diff.mean,aes(var,estimate,fill = Group)) +
            theme(panel.background = element_rect(fill = 'transparent'),
                  panel.grid = element_blank(),
                  axis.ticks.length = unit(0.4,"lines"), 
                  axis.ticks = element_line(color='black'),
                  axis.line = element_line(colour = "black"),
                  axis.title.x=element_text(colour='black', size=12,face = "bold"),
                  axis.text=element_text(colour='black',size=10,face = "bold"),
                  axis.text.y = element_blank(),
                  legend.position = "none",
                  axis.line.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  plot.title = element_text(size = 15,face = "bold",colour = "black",hjust = 0.5)) +
            scale_x_discrete(limits = levels(diff.mean$var)) +
            coord_flip() +
            xlab("") +
            ylab("Difference in mean proportions (%)") +
            labs(title="95% confidence intervals") 
        
        for (i in 1:(nrow(diff.mean) - 1)) 
            p2 <- p2 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf, 
                                fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
        
        p2 <- p2 +
            geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                          position = position_dodge(0.8), width = 0.5, size = 0.5) +
            geom_point(shape = 21,size = 3) +
            scale_fill_manual(values=cbbPalette) +
            geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black')
        
        
        p3 <- ggplot(diff.mean,aes(var,estimate,fill = Group)) +
            geom_text(aes(y = 0,x = var),label = diff.mean$p.value,
                      hjust = 0,fontface = "bold",inherit.aes = FALSE,size = 3) +
            geom_text(aes(x = nrow(diff.mean)/2 +0.5,y = 0.85),label = "P-value (corrected)",
                      srt = 90,fontface = "bold",size = 5) +
            coord_flip() +
            ylim(c(0,1)) +
            theme(panel.background = element_blank(),
                  panel.grid = element_blank(),
                  axis.line = element_blank(),
                  axis.ticks = element_blank(),
                  axis.text = element_blank(),
                  axis.title = element_blank())
        
        ## 图像拼接
        p <- p1 + p2 + p3 + plot_layout(widths = c(4,6,2))
        
        ## 保存图像
        ggsave(sprintf("D:/document/老婆大人/16S分析/picrust2/function/%s/hua_vs_%s.png", n, m), width = 20,height = 10)
    }
}
