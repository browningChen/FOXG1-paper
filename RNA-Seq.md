# review

## reads to counts

```bash
---shell---
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir ./out \
    -profile mamba \
    --extra_star_align_args "--alignIntronMax 1000000 --alignIntronMin 20 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --outFilterMismatchNmax 999 --outFilterMultimapNmax 20 --outFilterType BySJout --outFilterMismatchNoverLmax 0.1" \
    --star_index ./ref/GRCm39.STAR.Index \
    --fasta ./ref/Mus_musculus.GRCm39.dna.toplevel.fa \
    --gtf ./ref/Mus_musculus.GRCm39.109.gtf \
    --salmon_index ./ref/GRCm39.salmon_sa_index \
    --skip_markduplicates

find . -name quant.sf
zip quant.sf.zip `find ./* -name quant.sf`
echo "Samp" >salmon.output
cut -f 1 samplelist | xargs -i echo -e "{}\tsalmon_align/{}_quant/quant.sf" >>salmon.output
head salmon.output
sed 's/"/\t/g' ./ref/Mus_musculus.GRCm39.109.gtf | awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) print "TXname\tGene"; if($3=="transcript") print $14, $10}' > ./ref/GRCm39.tx2gene
---shell---
```

```r
---r---
rm(list=ls())
options(stringsAsFactors = F)
library(tximport) 
library(tidyverse) 
library(data.table) 
setwd("./salmon")

t2s <- fread("salmon_tx2gene.tsv", data.table = F, header = T); head(t2s)


files <- list.files(pattern="*quant.sf",recursive=T, full.names = T); files  
txi <- tximport(files, type = "salmon", tx2gene = t2s)

cn <- sapply(strsplit(files,'\\/'), function(x) x[length(x)-1]); cn
colnames(txi$counts) <- gsub('_quant','',cn); colnames(txi$counts)

counts <- as.data.frame(apply(txi$counts,2,as.integer)) 
rownames(counts) <- rownames(txi$counts) 
tpm <- as.data.frame(txi$abundance) 
colnames(tpm) <- colnames(txi$counts)
edgeR::cpm(counts)[1,1]
edgeR::cpm(counts) %>% colSums()

write.csv(counts,"raw_count.csv")
write.csv(tpm,"TPM.csv")
---r---
```

## downstream

```r
---r---
##DESeq2
library(DESeq2)
library(ggplot2)
library(ggrepel)

setwd("./nextflow_out")

count <- read.csv("raw_count.csv",header = T ,row.names = 1)
coldata <-  read.table("sampleinfo.txt",header = T)

# set `ctl` as reference level
coldata$condition <- as.factor(coldata$condition)
coldata$condition <- relevel(coldata$condition, ref = "WTControl")

dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = coldata,
                              design= ~condition)

dds <- DESeqDataSetFromTximport(txi,
                                coldata,
                                ~condition)
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
vstd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vstd, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf(file = "PCAplot.pdf", width = 12, height = 9);
ggplot(pcaData, aes(PC1, PC2,
                    color = coldata$condition,
                    label = coldata$Sample)) +
  geom_point(size=3) + coord_fixed() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  geom_text_repel(
    nudge_x = .15,
    box.padding = 0.5,
    nudge_y = 1,
    segment.curvature = -0.1,
    segment.ncp = 3,
    segment.angle = 20,
    size = 3)+
  theme_bw()
dev.off()
saveRDS(pcaData,"pcadata.rds")
write.csv(pcaData,"pcadata.csv")
#correlative plot
sampleDists <- dist(t(assay(vstd)))
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vstd$condition, vstd$Sample, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

dds <- DESeq(dds, betaPrior = FALSE)
resultsNames(dds)
#[1] "Intercept"                        "condition_FOXG1KO_vs_WTControl"   "condition_HD_vs_WTControl"       
#[4] "condition_HDFOXG1OE_vs_WTControl"
HD_vs_WT_res <- results(dds, lfcThreshold = 0,alpha = 0.1,pAdjustMethod = "BH",name="condition_HD_vs_WTControl")
#HD_vs_WT_res <- lfcShrink(dds, coef = 'condition_HD_vs_WTControl', type = 'apeglm', res = HD_vs_WT_res)
summary(HD_vs_WT_res)

OE_vs_WT_res <- results(dds,alpha = 0.1, name="condition_HDFOXG1OE_vs_WTControl")
#OE_vs_WT_res <- lfcShrink(dds, coef = 'condition_HDFOXG1OE_vs_WTControl', type = 'apeglm', res = OE_vs_WT_res)
summary(OE_vs_WT_res)

KO_vs_WT_res <- results(dds,alpha = 0.1, name="condition_FOXG1KO_vs_WTControl")
#KO_vs_WT_res <- lfcShrink(dds, coef = 'condition_FOXG1KO_vs_WTControl', type = 'apeglm', res = KO_vs_WT_res)
summary(KO_vs_WT_res)

#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(HD_vs_WT_res, plot(log2FoldChange, -log10(pvalue), pch=20, main="HD vs WT Volcano plot", ylim=c(0,30)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(HD_vs_WT_res, padj<0.1 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(HD_vs_WT_res, padj<0.01 & abs(log2FoldChange)>0), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

with(OE_vs_WT_res, plot(log2FoldChange, -log10(pvalue), pch=20, main="OE vs WT Volcano plot", ylim=c(0,30)))
with(subset(OE_vs_WT_res, padj<0.1 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(OE_vs_WT_res, padj<0.01 & abs(log2FoldChange)>0), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

with(KO_vs_WT_res, plot(log2FoldChange, -log10(pvalue), pch=20, main="KO vs WT Volcano plot", ylim=c(0,30)))
with(subset(KO_vs_WT_res, padj<0.1 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(KO_vs_WT_res, padj<0.01 & abs(log2FoldChange)>0), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
#OUTPUT
HD_vs_WT_result_data <- merge(as.data.frame(HD_vs_WT_res),
                     as.data.frame(counts(dds, normalized=TRUE)),
                     by="row.names",
                     sort=FALSE)
write.table(HD_vs_WT_result_data,file = "HD_vs_WT_result_data.xls",col.names = NA,sep = "\t",quote = F)
table(HD_vs_WT_res$padj<0.1)
diff_list_HD_vs_WT_FDR <- HD_vs_WT_res[order(HD_vs_WT_res$padj),]
diff_gene_HD_vs_WT_P <-subset(diff_list_HD_vs_WT_FDR,pvalue < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)) 
write.table(diff_gene_HD_vs_WT_P,file= "diff_gene_HD_vs_WT.xls",sep = '\t',quote = F)


OE_vs_WT_result_data <- merge(as.data.frame(OE_vs_WT_res),
                              as.data.frame(counts(dds, normalized=TRUE)),
                              by="row.names",
                              sort=FALSE)
write.table(OE_vs_WT_result_data,file = "OE_vs_WT_result_data.xls",col.names = NA,sep = "\t",quote = F) 
table(OE_vs_WT_res$padj<0.1)
diff_list_OE_vs_WT_FDR <- OE_vs_WT_res[order(OE_vs_WT_res$padj),]
diff_gene_OE_vs_WT_P <-subset(diff_list_OE_vs_WT_FDR,pvalue < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)) 
write.table(diff_gene_OE_vs_WT_P,file= "diff_gene_OE_vs_WT.xls",sep = '\t',quote = F)

KO_vs_WT_result_data <- merge(as.data.frame(KO_vs_WT_res),
                              as.data.frame(counts(dds, normalized=TRUE)),
                              by="row.names",
                              sort=FALSE)
write.table(KO_vs_WT_result_data,file = "KO_vs_WT_result_data.xls",col.names = NA,sep = "\t",quote = F) 
table(KO_vs_WT_res$padj<0.1)
diff_list_KO_vs_WT_FDR <- KO_vs_WT_res[order(KO_vs_WT_res$padj),]
diff_gene_KO_vs_WT_P <-subset(diff_list_KO_vs_WT_FDR,pvalue < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)) 
write.table(diff_gene_KO_vs_WT_P,file= "diff_gene_KO_vs_WT.xls",sep = '\t',quote = F)
#OE vs HD
coldata <-  read.table("sampleinfo.txt",header = T)
coldata$condition <- as.factor(coldata$condition)
coldata$condition <- relevel(coldata$condition, ref = "HD")
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = coldata,
                              design= ~ condition)
dds <- DESeqDataSetFromTximport(txi,
                                coldata,
                                ~condition)
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
dds <- DESeq(dds, betaPrior = FALSE)
resultsNames(dds)
OE_vs_HD_res <- results(dds,pAdjustMethod = "BH", name="condition_HDFOXG1OE_vs_HD")
KO_vs_HD_res <- results(dds,pAdjustMethod = "BH", name="condition_FOXG1KO_vs_HD")
summary(OE_vs_HD_res)
summary(KO_vs_HD_res)
with(OE_vs_HD_res, plot(log2FoldChange, -log10(pvalue), pch=20, main="OE vs HD Volcano plot", ylim=c(0,30)))
with(subset(OE_vs_HD_res, pvalue<0.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(OE_vs_HD_res, padj<0.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
OE_vs_HD_result_data <- merge(as.data.frame(OE_vs_HD_res),
                              as.data.frame(counts(dds, normalized=TRUE)),
                              by="row.names",
                              sort=FALSE)
KO_vs_HD_result_data <- merge(as.data.frame(KO_vs_HD_res),
                              as.data.frame(counts(dds, normalized=TRUE)),
                              by="row.names",
                              sort=FALSE)
write.table(OE_vs_HD_result_data,file = "OE_vs_HD_result_data.xls",col.names = NA,sep = "\t",quote = F)
write.table(KO_vs_HD_result_data,file = "KO_vs_HD_result_data.xls",col.names = NA,sep = "\t",quote = F)
table(OE_vs_HD_res$pvalue<0.05 & abs(OE_vs_HD_res$log2FoldChange) >1)
table(KO_vs_HD_res$pvalue<0.05 & abs(KO_vs_HD_res$log2FoldChange) >1)
diff_list_OE_vs_HD_P <- OE_vs_HD_res[order(OE_vs_HD_res$padj),]
diff_gene_OE_vs_HD_P <-subset(diff_list_OE_vs_HD_P,pvalue < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)) 
write.table(diff_gene_OE_vs_HD_P,file= "diff_gene_OE_vs_HD.xls",sep = '\t',quote = F)

diff_list_KO_vs_HD_P <- KO_vs_HD_res[order(KO_vs_HD_res$padj),]
diff_gene_KO_vs_HD_P <-subset(diff_list_KO_vs_HD_P,pvalue < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)) 
write.table(diff_gene_KO_vs_HD_P,file= "diff_gene_KO_vs_HD.xls",sep = '\t',quote = F)

saveRDS(HD_vs_WT_res,"HD_vs_WT_res.rds")
saveRDS(OE_vs_WT_res,"OE_vs_WT_res.rds")
saveRDS(OE_vs_HD_res,"OE_vs_HD_res.rds")
saveRDS(KO_vs_HD_res,"KO_vs_HD_res.rds")

library(clusterProfiler)
library(org.Mm.eg.db)

table(HD_vs_WT_res$pvalue<0.05 & abs(HD_vs_WT_res$log2FoldChange) >1)
table(OE_vs_WT_res$pvalue<0.05 & abs(OE_vs_WT_res$log2FoldChange) >1)
table(KO_vs_WT_res$pvalue<0.05 & abs(KO_vs_WT_res$log2FoldChange) >1)
diff_gene_HD_vs_WT_P <-subset(diff_list_HD_vs_WT_FDR,pvalue < 0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5))
diff_gene_OE_vs_WT_P <-subset(diff_list_OE_vs_WT_FDR,pvalue < 0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5))
diff_gene_KO_vs_WT_P <-subset(diff_list_KO_vs_WT_FDR,pvalue < 0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5))
diff_gene_OE_vs_HD_P <-subset(diff_list_OE_vs_HD_P,pvalue < 0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5))

list <- list(`HD vs WT`=row.names(diff_gene_HD_vs_WT_P),`OE vs WT`=row.names(diff_gene_OE_vs_WT_P),
             `KO vs WT`=row.names(diff_gene_KO_vs_WT_P),`OE vs HD`=row.names(diff_gene_OE_vs_HD_P))
compGOp005FC05 <- compareCluster(list,
                         fun = "enrichGO",
                         OrgDb = "org.Mm.eg.db",
                         ont = "ALL",
                         readable = T,
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 1,
                         keyType = "ENSEMBL")
write.csv(compGOp005FC05,"compGOp005FC05.csv")
categorys <- c("neuron migration","GABAergic neuron differentiation",
               "Wnt signaling pathway","oligodendrocyte differentiation",
               "renal system development","transforming growth factor beta receptor signaling pathway",
               "response to dopamine","dendrite membrane","axonogenesis","ERK1 and ERK2 cascade","cellular respiration")
               
dotplot(compGOp005FC05,by = "geneRatio",showCategory = categorys)+scale_color_gradient(low = "#c878af" , high = "#3388de")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("godotplot.pdf")
###HD vs OE
#DEG<- read.table("diff_gene_OE_vs_HD.xls",sep = "\t",header = T)
ego<- enrichGO(row.names(diff_gene_OE_vs_HD_P),
               OrgDb = "org.Mm.eg.db",
               ont = "ALL",
               qvalueCutoff = 1,
               readable = T,
               keyType = "ENSEMBL")
write.csv(ego,"OE vs HD GO.csv")
dotplot(ego)
---r---
```
