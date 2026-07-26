res_oe_hd <- read.csv(
  "E:/硕士/bioinformatics/Foxg1 rnaseq/nextflow/pairwise_deseq2_rerun/OE_vs_HD_pairwise_DESeq2_full.csv",
  stringsAsFactors = FALSE
)

# 阈值可改为 abs(log2FoldChange) > 1
selected_ids <- res_oe_hd$geneid[
  !is.na(res_oe_hd$pvalue) &
    res_oe_hd$pvalue < 0.05 &
    abs(res_oe_hd$log2FoldChange) > 0.5
]

# Ensembl ID -> symbol
tx2gene <- read.delim(
  "E:/硕士/bioinformatics/Foxg1 rnaseq/nextflow/salmon_tx2gene.tsv",
  header = TRUE,
  stringsAsFactors = FALSE
)
tx2gene <- tx2gene[!duplicated(tx2gene$geneid), ]

# normalized counts, then log2 transformation
expr <- log2(counts(dds, normalized = TRUE) + 1)

# Column order: HD -> OE -> WT
sample_order <- c(
  grep("^HD_", colnames(expr), value = TRUE),
  grep("^OE_", colnames(expr), value = TRUE),
  grep("^CONTROL_", colnames(expr), value = TRUE)
)

# Selected expression matrix
heatmap_expr <- expr[
  intersect(selected_ids, rownames(expr)),
  sample_order,
  drop = FALSE
]

# Set row names to gene symbols
symbols <- tx2gene$genename[
  match(rownames(heatmap_expr), tx2gene$geneid)
]
symbols[is.na(symbols) | symbols == ""] <-
  rownames(heatmap_expr)[is.na(symbols) | symbols == ""]

rownames(heatmap_expr) <- make.unique(symbols)

# Convert to data.frame if needed
heatmap_df <- as.data.frame(heatmap_expr)

head(heatmap_df)
write.csv(heatmap_df,"OE vs HD heatmap_df.csv")
