gsealist_HDWT <- read.table(file = "gsealist_HDWT.txt" , header = T )
gsealist_OEWT <- read.table(file = "gsealist_OEWT.txt" , header = T )
rank.HDWTlist = gsealist_HDWT[,2]
names(rank.HDWTlist) = as.character(gsealist_HDWT[,1])
rank.HDWTlist = sort(rank.HDWTlist , decreasing = T)

rank.OEWTlist = gsealist_OEWT[,2]
names(rank.OEWTlist) = as.character(gsealist_OEWT[,1])
rank.OEWTlist = sort(rank.OEWTlist , decreasing = T)

gsealist <- list(`HD vs WT`=rank.HDWTlist,`OE vs WT`=rank.OEWTlist)
compGSEAGO <- compareCluster(gsealist,
                             fun = "gseGO",
                             OrgDb = "org.Mm.eg.db",
                             keyType = "ENSEMBL",
                             ont = "BP",
                             verbose = TRUE)
write.csv(compGSEAGO,"compGSEA GO.csv")
dotplot(compGSEAGO)

library(msigdbr)
KEGG = msigdbr(species = "Mus musculus",
               category = "C2", 
               subcategory = "CP:KEGG") %>% dplyr::select(gs_name,ensembl_gene)

compGSEAKEGG <- compareCluster(gsealist,
                               TERM2GENE = KEGG,
                               fun = "GSEA")
write.csv(compGSEAKEGG,"compGSEA KEGG.csv")
dotplot(compGSEAKEGG)

egosimp <- simplify(egoMF,cutoff=0.7,by="p.adjust",
                    select_fun=min,measure="Wang")

library("genekitr")
GSEAGO_HD <- gseGO(rank.HDWTlist,
                   OrgDb = "org.Mm.eg.db",
                   keyType = "ENSEMBL",
                   ont = "BP",
                   verbose = TRUE) %>% importCP(type = "gsea")#%>% simplify(cutoff=0.7,by="p.adjust",select_fun=min) %>% importCP(type = "gsea")
GSEAGO_OE <- gseGO(rank.OEWTlist,
                   OrgDb = "org.Mm.eg.db",
                   keyType = "ENSEMBL",
                   ont = "BP",
                   verbose = TRUE) %>% importCP(type = "gsea")

library(patchwork)
pathways <- c("GO:0050803", "GO:0099054", "GO:0002377",
              "GO:0046651","GO:0106027","GO:0007269",
              "GO:0016577","GO:0002274")

p1 <- plotGSEA(GSEAGO_HD, plot_type = "volcano", show_pathway = pathways ,
               label_by = c("description"),
               main_text_size = 15 ,
               colour = c('#ec273f', "#36a5f4", '#ec273f')
               ) 

p2 <- plotGSEA(GSEAGO_OE, plot_type = "volcano", show_pathway = 3 ,
               label_by = c("description"),main_text_size = 15,
               colour = c('#ec273f', "#36a5f4", '#ec273f'))

p1 + p2# + plot_annotation(tag_levels = "A")
ggsave("GSEA volcano.pdf")
gsea_plot_HD <- data.frame(  
  NES = GSEAGO_HD$gsea_df$NES,  
  padjust = -log10(GSEAGO_HD$gsea_df$p.adjust) ,
  Description = GSEAGO_HD$gsea_df$Description
)  

gsea_plot_HD$label <- ifelse(gsea_plot_HD$padjust > 8 & gsea_plot_HD$NES >= 2, gsea_plot_HD$Description, 
                             ifelse(gsea_plot_HD$padjust > 2 & gsea_plot_HD$NES <= 2, gsea_plot_HD$Description, ""))


ggplot(data = gsea_plot_HD, aes(x = NES, y = padjust)) +  
  geom_point(aes(size = padjust,color = ifelse(NES > 0, "up", "down"))) +  
  scale_size(range = c(0.1, 5)) + # 根据需要调整点的大小范围  
  xlab("NES") +  
  ylab("-log10(padjust)") +  scale_color_manual(values=c("#36a5f4", "#ec273f", "black"))+
  theme_bw() + 
  ggrepel::geom_label_repel(aes(x = NES, y = padjust ,label = label),
                            size = 3,box.padding = unit(0.5, "lines"),point.padding = unit(0.8, "lines"),
                            segment.color = "black",show.legend = FALSE)

'#26854c', "white", '#5e5b8c' "#B3A86A","white","#D4613E" '#6dead6', "white", '#c878af' '#36c5f4', "white", '#ec273f'