library(patchwork)
library(DESeq2)
library(ggplot2)
library(ggrepel)
setwd("E:/硕士/bioinformatics/Foxg1 rnaseq/final")

load("deseq.RData")
tx2gene <- read.delim(
  "../nextflow/salmon_tx2gene.tsv",
  header = TRUE,
  stringsAsFactors = FALSE
)
tx2gene <- tx2gene[!duplicated(tx2gene$geneid), ]
# DESeq2 normalized counts
norm_counts <- counts(dds, normalized = TRUE)

gene_symbol <- tx2gene$genename[
  match(rownames(norm_counts), tx2gene$geneid)
]
gene_symbol[is.na(gene_symbol) | gene_symbol == ""] <-
  rownames(norm_counts)[is.na(gene_symbol) | gene_symbol == ""]

# Sample groups
groups <- list(
  WT = grep("^CONTROL_", colnames(norm_counts), value = TRUE),
  HD = grep("^HD_", colnames(norm_counts), value = TRUE),
  KD = grep("^KO_", colnames(norm_counts), value = TRUE),
  OE = grep("^OE_", colnames(norm_counts), value = TRUE)
)

# Modify gene labels here
label_sets <- list(
  "HD vs WT" = c("Trpc5",
                 "Stxbp5l",
                 "Cpeb3",
                 "Gabrb1",
                 "Fis1"),
  "OE vs WT" = c("Trpc5",
                 "Stxbp5l",
                 "Cpeb3",
                 "Gabrb1",
                 "Fis1"),
  "KD vs WT" = c("Spp1", "Cd14", "Lst1",
                 "Alox5ap", "Cebpb", "Tnfaip3"),
  "OE vs HD" = c("Trpc5",
                 "Stxbp5l",
                 "Cpeb3",
                 "Gabrb1",
                 "Fis1")
)

make_scatter <- function(res, x_samples, y_samples,
                         comparison_name, x_label, y_label,
                         label_genes = character(0)) {
  
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  
  plot_df <- data.frame(
    gene_id = rownames(norm_counts),
    symbol = gene_symbol,
    x_mean = log2(rowMeans(norm_counts[, x_samples, drop = FALSE]) + 1),
    y_mean = log2(rowMeans(norm_counts[, y_samples, drop = FALSE]) + 1),
    stringsAsFactors = FALSE
  )
  
  plot_df$log2FoldChange <- res_df$log2FoldChange[
    match(plot_df$gene_id, res_df$gene_id)
  ]
  plot_df$pvalue <- res_df$pvalue[
    match(plot_df$gene_id, res_df$gene_id)
  ]
  
  # Threshold: |log2FC| > 1 and p < 0.05
  plot_df$diff <- "NS"
  plot_df$diff[
    !is.na(plot_df$pvalue) &
      plot_df$pvalue < 0.05 &
      plot_df$log2FoldChange > 0.5
  ] <- "Up"
  
  plot_df$diff[
    !is.na(plot_df$pvalue) &
      plot_df$pvalue < 0.05 &
      plot_df$log2FoldChange < -0.5
  ] <- "Down"
  
  # Draw NS first so significant points remain visible
  plot_df$diff <- factor(plot_df$diff, levels = c("NS", "Down", "Up"))
  plot_df <- plot_df[order(plot_df$diff), ]
  
  label_df <- plot_df[plot_df$symbol %in% label_genes, ]
  label_df <- label_df[!duplicated(label_df$symbol), ]
  
  n_up <- sum(plot_df$diff == "Up", na.rm = TRUE)
  n_down <- sum(plot_df$diff == "Down", na.rm = TRUE)
  
  p <- ggplot(plot_df, aes(x = x_mean, y = y_mean)) +
    geom_abline(
      intercept = 0, slope = 1,
      linewidth = 0.4, color = "black"
    ) +
    geom_abline(
      intercept = c(-1, 1), slope = 1,
      linewidth = 0.4, linetype = "dashed", color = "black"
    ) +
    geom_point(aes(color = diff), size = 1.15, alpha = 0.8) +
    scale_color_manual(
      values = c(
        "Up" = "#D64B4B",
        "Down" = "#2A6FBB",
        "NS" = "#B8B8B8"
      ),
      name = NULL
    ) +
    # Upper-left: upregulated DEG number
    annotate(
      "label",
      x = -Inf, y = Inf,
      label = paste0("Up: ", n_up),
      hjust = -0.08, vjust = 1.15,
      color = "#D64B4B",
      fill = "white",
      label.size = 0.35,
      size = 3.5,
      family = "Arial"
    ) +
    # Lower-right: downregulated DEG number
    annotate(
      "label",
      x = Inf, y = -Inf,
      label = paste0("Down: ", n_down),
      hjust = 1.08, vjust = -0.35,
      color = "#2A6FBB",
      fill = "white",
      label.size = 0.35,
      size = 3.5,
      family = "Arial"
    ) +
    coord_equal(clip = "off") +
    labs(
      title = comparison_name,
      x = x_label,
      y = y_label
    ) +
    theme_classic(base_family = "Arial") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "top",
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      plot.margin = margin(8, 12, 8, 12)
    )
  
  if (nrow(label_df) > 0) {
    p <- p +
      geom_text_repel(
        data = label_df,
        aes(label = symbol),
        size = 3,
        family = "Arial",
        max.overlaps = Inf,
        box.padding = 0.35,
        point.padding = 0.2,
        min.segment.length = 0,
        seed = 123
      )
  }
  
  p
}

p_hd_wt <- make_scatter(
  res = HD_vs_WT_res,
  x_samples = groups$WT,
  y_samples = groups$HD,
  comparison_name = "HD vs WT",
  x_label = "WT mean normalized expression",
  y_label = "HD mean normalized expression",
  label_genes = label_sets[["HD vs WT"]]
)
ggsave("p_hd_wt.png",dpi=300,width=8,height=4)
p_oe_wt <- make_scatter(
  res = OE_vs_WT_res,
  x_samples = groups$WT,
  y_samples = groups$OE,
  comparison_name = "OE vs WT",
  x_label = "WT mean normalized expression",
  y_label = "OE mean normalized expression",
  label_genes = label_sets[["OE vs WT"]]
)
ggsave("p_oe_wt.png",dpi=300,width=8,height=4)
p_kd_wt <- make_scatter(
  res = KO_vs_WT_res,
  x_samples = groups$WT,
  y_samples = groups$KD,
  comparison_name = "KD vs WT",
  x_label = "WT mean normalized expression",
  y_label = "KD mean normalized expression",
  label_genes = label_sets[["KD vs WT"]]
)
ggsave("p_kd_wt.png",dpi=300,width=8,height=4)
p_oe_hd <- make_scatter(
  res = OE_vs_HD_res,
  x_samples = groups$HD,
  y_samples = groups$OE,
  comparison_name = "OE vs HD",
  x_label = "HD mean normalized expression",
  y_label = "OE mean normalized expression",
  label_genes = label_sets[["OE vs HD"]]
)
ggsave("p_oe_hd.png",dpi=300,width=8,height=4)
final_plot <- (p_hd_wt + p_oe_wt) / (p_kd_wt + p_oe_hd) #+
  #plot_annotation(
  #  title = "FOXG1 perturbation and rescue-associated transcriptional changes"
  #)

final_plot

# Editable vector PDF
ggsave(
  "FOXG1_four_comparison_scatter.pdf",
  plot = final_plot,
  width = 14,
  height = 12,
  device = "pdf",
  bg = "white"
)


