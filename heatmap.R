library(ComplexHeatmap)
library(ClusterGVis)
exps <- read.delim(
  "E:/硕士/bioinformatics/Foxg1 rnaseq/final/HD vs WT diffgene merge ALL rpm.xls",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)
exps <- read.delim(
  "E:\\硕士\\bioinformatics\\Foxg1 rnaseq\\final\\HD_diff_merge_rpm.txt",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

exps <- data.matrix(heatmap_df)
exps_log <- log2(exps + 1)

# Elbow plot
cm <- clusterData(
  obj = exps,
  scaleData = TRUE,
  cluster.method = "mfuzz",
  cluster.num = 5
)

visCluster(
  object = cm,
  plot.type = "both",
  column_names_rot = 45
)

#以下两种聚类算法都可以用
# using mfuzz for clustering
# mfuzz
cm <- clusterData(obj = exps,
                  cluster.method = "mfuzz",
                  cluster.num = 5)

# using complexheatmap row_km for clustering
# kmeans
ck <- clusterData(obj = exps,
                  cluster.method = "wgcna",
                  cluster.num = 5
                  )

#导出聚类后的数据进行通路富集
write.csv(ck$wide.res, file = "k-means clustering result.csv",quote = F)

#visCluster 绘图
#plot line only
#kmeans 结果的折线图,因为没有 membership 的信息,所以没有颜色映射
#折线图本质上是 ggplot2 对象,你可以添加其它相关的参数来进行修改细节
visCluster(object = cm,
           plot.type = "line")
#ms.col = c("royalblue","orange","gold")

# plot heatmap only
visCluster(object = cm,
           plot.type = "heatmap")

#添加注释
#   id                               term
# 1 C1              developmental process
# 2 C1   anatomical structure development
# 3 C1 multicellular organism development
# 4 C2                 system development
enrich <- enrichCluster(object = cm,
                        OrgDb = org.Mm.eg.db,
                        type = "BP",
                        idTrans = F,
                        fromType = "SYMBOL",
                        readable = F,
                        pvalueCutoff = 0.5,
                        topn = 5)
termanno <- read.csv(file = "annotation.txt" ,header = T,sep = "\t")
head(termanno)

# anno with GO terms
pdf('./HTtermCmlsrt.pdf',height = 6,width = 12)
visCluster(object = cm,
           plot.type = "both",
           add.box = T,
           add.line = F,
           mline.size = 0.05,
           go.size = 8 ,
           column_names_rot = 45,
           boxcol = ggsci::pal_npg()(9),
           annoTerm.data = enrich,
           markGenes = markGenes,
           line.side = "right",
           cluster_columns = F ,
           show_row_dend = T , #去除聚类树
           sample.group = rep(c("group1","group2","group3"),each = 3),
           #col = colorRamp2(c(-3,0,3), c("royalblue", "grey", "gold"))
           ht.col.list = list(col_range = c(-2, 0, 2), col_color = c('#36a5f4', "white", '#ec273f'))
)# + scale_color_gradient(low = "royalblue", high = "gold")

dev.off()

markGenes = rownames(exps)[sample(1:nrow(exps),30,replace = F)]
markGenes <- c(
  "Drd1", "Drd2", "Adora2a", "Ppp1r1b", "Pde10a",
  "Rgs9", "Gnal", "Tac1", "Penk", "Gpr6",
  "Bdnf", "Ntrk2", "Grin1", "Grin2b", "Dlg4",
  "Gfap", "Aif1", "C1qa", "Tyrobp", "Ppargc1a"
)
visCluster(object = cm,
           plot.type = "both",
           column_names_rot = 45,
           markGenes = markGenes)

'#26854c', "white", '#5e5b8c' "#B3A86A","white","#D4613E" '#6dead6', "white", '#c878af' '#36c5f4', "white", '#ec273f'