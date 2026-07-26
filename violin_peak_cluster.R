# FOXG1 CUT&Tag peak-cluster versus RNA-seq expression
#
# This script fixes the repeated-gene weighting, raw-P filtering, and Fisher
# background problems in violin.R. It uses the existing peak clusters, so its
# biological unit is a CUT&Tag peak assigned to its nearest gene; it is NOT an
# exact reimplementation of the paper's TSS-signal k-means analysis.
#
# Rule used to obtain one observation per gene:
#   1. Keep the highest-score peak for each gene within each peak cluster.
#   2. If a gene occurs in several clusters, keep its highest-score peak across
#      all clusters. Ties are resolved deterministically by cluster name.

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(GenomicFeatures)
  library(dplyr)
  library(ggplot2)
  library(ggGenshin)
  library(GeneOverlap)
})

project_dir <- "E:/硕士/bioinformatics/Foxg1 rnaseq/final"
cuttag_cluster_file <- "E:/硕士/bioinformatics/Foxg1 CUT&TAG/重测序/cluster/FOXG1_cluster_regions.xls"
preannotated_file <- "E:/硕士/bioinformatics/Foxg1 CUT&TAG/重测序/cluster/total_peak_anno.txt"
gtf_file <- "E:/硕士/bioinformatics/GRCm39 ref/Mus_musculus.GRCm39.109.gtf"
output_dir <- file.path(project_dir, "peak_cluster_violin_corrected")

fc_cutoff <- 0.5
fdr_cutoff <- 0.05
promoter_window_bp <- 5000
plot_abs_lfc_max <- 7

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(project_dir)

stopifnot(file.exists(cuttag_cluster_file), file.exists(gtf_file))

# Reuse the annotation table created in the original CUT&Tag workflow when it
# is available. It makes repeated plotting fast. Set preannotated_file <- NULL
# to force fresh annotation from the raw deepTools BED-like table.
if (!is.null(preannotated_file) && file.exists(preannotated_file)) {
  cached <- read.delim(
    preannotated_file,
    fileEncoding = "UTF-16LE",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (ncol(cached) < 24L) stop("Unexpected format in ", preannotated_file)
  peak_anno <- data.frame(
    peak_score = suppressWarnings(as.numeric(cached[[7L]])),
    cluster = as.character(cached[[15L]]),
    geneId = as.character(cached[[22L]]),
    distanceToTSS = suppressWarnings(as.numeric(cached[[24L]])),
    stringsAsFactors = FALSE
  )
} else {
  # comment.char="" is essential: otherwise the #chrom header is skipped and
  # the first peak is misread as column names.
  peak_tbl <- read.delim(
    cuttag_cluster_file,
    header = TRUE,
    sep = "\t",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  required_peak_cols <- c("#chrom", "start", "end", "score", "deepTools_group")
  missing_peak_cols <- setdiff(required_peak_cols, names(peak_tbl))
  if (length(missing_peak_cols) > 0) {
    stop("Missing peak columns: ", paste(missing_peak_cols, collapse = ", "))
  }
  peak_gr <- GenomicRanges::GRanges(
    seqnames = peak_tbl[["#chrom"]],
    ranges = IRanges::IRanges(start = peak_tbl$start + 1L, end = peak_tbl$end),
    strand = "*",
    peak_score = as.numeric(peak_tbl$score),
    cluster = as.character(peak_tbl$deepTools_group)
  )
  txdb <- makeTxDbFromGFF(gtf_file)
  peak_anno <- as.data.frame(
    annotatePeak(
      peak_gr,
      TxDb = txdb,
      tssRegion = c(-promoter_window_bp, promoter_window_bp),
      verbose = FALSE
    )
  )
}

# Restrict to peaks within +/-5 kb of an annotated TSS. This is still a
# peak-cluster analysis, but it prevents distal peaks from being presented as
# TSS-associated FOXG1 binding.
gene_cluster_candidates <- peak_anno %>%
  transmute(
    geneId = as.character(geneId),
    cluster = as.character(cluster),
    peak_score = as.numeric(peak_score),
    distanceToTSS = as.numeric(distanceToTSS)
  ) %>%
  filter(!is.na(geneId), !is.na(cluster), abs(distanceToTSS) <= promoter_window_bp) %>%
  arrange(geneId, cluster, desc(peak_score)) %>%
  distinct(geneId, cluster, .keep_all = TRUE)

# A gene must occur once in a violin/Fisher universe. Choosing the strongest
# associated peak makes the tie-breaking rule explicit and reproducible.
cluster_gene <- gene_cluster_candidates %>%
  arrange(geneId, desc(peak_score), cluster) %>%
  distinct(geneId, .keep_all = TRUE) %>%
  mutate(cluster = factor(cluster, levels = paste0("cluster_", 1:5))) %>%
  arrange(cluster, geneId)

write.csv(peak_anno, file.path(output_dir, "all_peak_annotations.csv"), row.names = FALSE)
write.csv(gene_cluster_candidates, file.path(output_dir, "gene_cluster_candidates.csv"), row.names = FALSE)
write.csv(cluster_gene, file.path(output_dir, "unique_gene_cluster_assignment.csv"), row.names = FALSE)

read_deseq <- function(filename) {
  x <- read.delim(filename, check.names = FALSE, stringsAsFactors = FALSE)
  # DESeq2 exports from this project have an unnamed first index column.
  # dplyr rejects data frames with blank names, although that column is not
  # used below.
  blank_names <- is.na(names(x)) | names(x) == ""
  names(x)[blank_names] <- paste0("input_index_", seq_len(sum(blank_names)))
  names(x) <- make.unique(names(x))
  if (!"Row.names" %in% names(x)) {
    stop(filename, " must contain an Ensembl gene-ID column named 'Row.names'.")
  }
  needed <- c("log2FoldChange", "padj")
  missing <- setdiff(needed, names(x))
  if (length(missing) > 0) {
    stop(filename, " is missing: ", paste(missing, collapse = ", "))
  }

  x %>%
    transmute(
      geneId = as.character(Row.names),
      log2FoldChange = as.numeric(log2FoldChange),
      padj = as.numeric(padj)
    ) %>%
    filter(!is.na(geneId), is.finite(log2FoldChange)) %>%
    distinct(geneId, .keep_all = TRUE)
}

run_fisher <- function(plot_df) {
  universe <- plot_df %>%
    mutate(
      status = case_when(
        log2FoldChange >= fc_cutoff & !is.na(padj) & padj <= fdr_cutoff ~ "up",
        log2FoldChange <= -fc_cutoff & !is.na(padj) & padj <= fdr_cutoff ~ "down",
        TRUE ~ "static"
      )
    )

  statuses <- c("up", "down")
  clusters <- levels(droplevels(universe$cluster))
  result <- vector("list", length(statuses) * length(clusters))
  i <- 1L

  for (cl in clusters) {
    for (st in statuses) {
      in_cluster <- universe$cluster == cl
      in_status <- universe$status == st
      contingency <- matrix(
        c(
          sum(in_cluster & in_status),
          sum(in_cluster & !in_status),
          sum(!in_cluster & in_status),
          sum(!in_cluster & !in_status)
        ),
        nrow = 2,
        byrow = TRUE
      )
      ft <- fisher.test(contingency)
      result[[i]] <- data.frame(
        cluster = cl,
        DEG_set = st,
        cluster_genes = sum(in_cluster),
        DEG_genes = sum(in_status),
        overlap = contingency[1, 1],
        odds_ratio = unname(ft$estimate),
        p_value = ft$p.value,
        stringsAsFactors = FALSE
      )
      i <- i + 1L
    }
  }

  bind_rows(result) %>%
    mutate(padj_BH = p.adjust(p_value, method = "BH"))
}

save_geneoverlap_heatmap <- function(analysis_df, contrast_name, fisher_results) {
  prefix <- gsub("[^A-Za-z0-9]+", "_", contrast_name)
  # drawHeatmap(cutoff = 0.05) fails when every adjusted P value is above the
  # cutoff because no finite heatmap cells remain. The CSV remains the complete
  # GeneOverlap/Fisher result in that situation.
  if (!any(fisher_results$padj_BH <= 0.05, na.rm = TRUE)) {
    writeLines(
      "No cluster-by-up/down-DEG enrichment passed BH-adjusted P <= 0.05.",
      con = file.path(output_dir, paste0(prefix, "_GeneOverlap_no_significant_result.txt"))
    )
    return(invisible(NULL))
  }

  classified <- analysis_df %>%
    mutate(
      status = case_when(
        log2FoldChange >= fc_cutoff & !is.na(padj) & padj <= fdr_cutoff ~ "increased_DEG",
        log2FoldChange <= -fc_cutoff & !is.na(padj) & padj <= fdr_cutoff ~ "decreased_DEG",
        TRUE ~ "static"
      )
    )

  # All lists are unique and drawn from exactly the same DESeq2/TSS-associated
  # universe. This is the corrected equivalent of the original newGOM block.
  cluster_list <- split(classified$geneId, classified$cluster) %>%
    lapply(unique)
  deg_classes <- c("increased_DEG", "decreased_DEG", "static")
  deg_list <- setNames(
    lapply(deg_classes, function(x) unique(classified$geneId[classified$status == x])),
    deg_classes
  )

  gom <- newGOM(
    cluster_list,
    deg_list,
    genome.size = length(unique(classified$geneId))
  )

  save_one <- function(filename, device) {
    device(filename, width = 8, height = 6)
    par(mar = c(5, 5, 2, 2))
    drawHeatmap(
      gom,
      what = "odds.ratio",
      adj.p = TRUE,
      cutoff = 0.05,
      ncolused = 5,
      grid.col = "Blues",
      note.col = "Black"
    )
    dev.off()
  }
  save_one(
    file.path(output_dir, paste0(prefix, "_GeneOverlap_odds_ratio.pdf")),
    grDevices::pdf
  )
  save_one(
    file.path(output_dir, paste0(prefix, "_GeneOverlap_odds_ratio.png")),
    function(filename, width, height) grDevices::png(
      filename, width = width * 300, height = height * 300, res = 300
    )
  )
  invisible(gom)
}

plot_one_contrast <- function(de_file, contrast_name, palette) {
  de <- read_deseq(de_file)
  analysis_df <- inner_join(cluster_gene, de, by = "geneId") %>%
    filter(!is.na(cluster)) %>%
    droplevels()

  if (nrow(analysis_df) == 0L) stop("No gene IDs overlap for ", contrast_name)

  # Preserve all tested genes in Fisher's exact test, but prevent a few extreme
  # non-shrunken LFCs from making the violin unreadable.
  plot_df <- analysis_df %>%
    filter(abs(log2FoldChange) < plot_abs_lfc_max) %>%
    droplevels()
  if (nrow(plot_df) == 0L) stop("No plottable genes remain for ", contrast_name)

  significant <- plot_df %>%
    filter(abs(log2FoldChange) >= fc_cutoff, !is.na(padj), padj <= fdr_cutoff)

  # Visual styling intentionally follows the user's original violin.R:
  # ggGenshin palette, theme_light, black median, and coloured points.
  p <- ggplot(
    plot_df,
    aes(x = cluster, y = log2FoldChange, fill = cluster, colour = cluster)
  ) +
    geom_violin() +
    geom_point(
      data = significant,
      inherit.aes = TRUE,
      shape = 16,
      size = 2.5,
      alpha = 0.5,
      position = position_jitter(width = 0.2, height = 0)
    ) +
    stat_summary(fun = median, geom = "point", size = 2, colour = "black") +
    scale_fill_manual(values = palette, drop = FALSE) +
    scale_colour_manual(values = palette, drop = FALSE) +
    labs(
      x = "cluster",
      y = "Log2FC"
    ) +
    theme_light() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12, angle = 60),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 12)
    )

  fisher_results <- run_fisher(analysis_df)
  prefix <- gsub("[^A-Za-z0-9]+", "_", contrast_name)
  write.csv(plot_df, file.path(output_dir, paste0(prefix, "_violin_source_data.csv")), row.names = FALSE)
  write.csv(analysis_df, file.path(output_dir, paste0(prefix, "_fisher_universe_source_data.csv")), row.names = FALSE)
  write.csv(fisher_results, file.path(output_dir, paste0(prefix, "_fisher_enrichment.csv")), row.names = FALSE)
  ggsave(file.path(output_dir, paste0(prefix, "_violin.pdf")), p, width = 4, height = 3)
  ggsave(file.path(output_dir, paste0(prefix, "_violin.png")), p, width = 4, height = 3, dpi = 600)
  save_geneoverlap_heatmap(analysis_df, contrast_name, fisher_results)

  message(
    contrast_name, ": ", nrow(plot_df), " displayed genes; ",
    nrow(significant), " significant genes."
  )
  invisible(list(plot = p, data = plot_df, fisher_universe = analysis_df, fisher = fisher_results))
}

orange_palette <- ggGenshin::pal_xiao(alpha = 1)(5)
names(orange_palette) <- paste0("cluster_", 1:5)
green_palette <- ggGenshin::pal_nahida(alpha = 1)(5)
names(green_palette) <- paste0("cluster_", 1:5)

oe_vs_hd <- plot_one_contrast(
  de_file = "OE_vs_HD_result_data.xls",
  contrast_name = "OE vs HD",
  palette = orange_palette
)

ko_vs_wt <- plot_one_contrast(
  de_file = "KO_vs_WT_result_data.xls",
  contrast_name = "KO vs WT",
  palette = green_palette
)
