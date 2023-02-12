library(tidyverse)
library(ALDEx2)

diff_aldex2 <- function(
    otu=NULL,
    map=NULL,
    comp_group=NULL, # 指定分组名称，vector
    test=NULL, ## t, kw, glm
    paired.test=FALSE,
    file_path=NULL,
    Group="Group"){
    
    # 选择需比较的组别
    map = map[Group]
    colnames(map) = "Group"
    map = map %>% filter(Group %in% comp_group)
    
    # Extract only those ID in common between the two tables
    idx = colnames(otu) %in% rownames(map)
    otu = otu[, idx]
    map = map[colnames(otu), ]
    
    ### 差异分析
    if(test == "t"){
        x.all <- aldex(otu, map, mc.samples=126, test=test, paired.test=paired.test,
                       effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
        png(sprintf("%s/aldex2_res.png", file_path))
        aldex.plot(x.all, type="MA", test="welch", xlab="Log-ratio abundance",ylab="Difference")
        dev.off()
    } else {
      x.all <- aldex(otu, map, mc.samples=126, test=test,
                     effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
      png(sprintf("%s/aldex2_res.png", file_path))
      aldex.plot(x.all, type="MA", test="glm", xlab="Log-ratio abundance",ylab="Difference")
      dev.off()
    }
    write.csv(x.all, sprintf("%s/aldex2_res.csv", file_path))
}

otutab <- read.table("data/otutab.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors=F)
metadata <- read.table("data/metadata.txt", header=T, row.names=1, sep="\t", 
                       comment.char="", stringsAsFactors=F)
taxonomy <- read.table("data/taxonomy.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors=F)

diff_aldex2(otu = otutab, map = metadata, comp_group = c("OE", "KO"), 
            test="t", paired.test=FALSE, file_path = "difference_analysis")
