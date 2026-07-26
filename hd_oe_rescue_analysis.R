#!/usr/bin/env Rscript

# HD -> WT rescue analysis
#
# Contrast directions used here:
#   HD_vs_WT: log2(HD / WT)
#   OE_vs_HD: log2(HD FOXG1 OE / HD GFP)
# Therefore, the inferred OE-vs-WT effect is:
#   log2FC_OE_vs_WT = log2FC_HD_vs_WT + log2FC_OE_vs_HD
#
# A "true approach to WT" is defined among HD-vs-WT DEGs as a gene for which
# OE changes expression in the opposite direction to the HD-vs-WT alteration
# AND reduces its absolute distance from WT:
#   sign(HD_vs_WT) * sign(OE_vs_HD) < 0
#   |HD_vs_WT + OE_vs_HD| < |HD_vs_WT|
# A change that crosses WT is retained only if its final absolute distance from
# WT is still smaller than before OE; changes that overshoot farther away are
# excluded.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ---- Input/output paths -----------------------------------------------------
project_dir <- "E:/硕士/bioinformatics/Foxg1 rnaseq/final"
hd_file <- file.path(project_dir, "HD_vs_WT_result_data.xls")
oe_file <- file.path(project_dir, "OE_vs_HD_result_data.xls")
output_dir <- file.path(project_dir, "HD_OE_rescue_analysis")

# ---- Prespecified thresholds ------------------------------------------------
lfc_cutoff <- 0.5
fdr_cutoff <- 0.05
n_permutations <- 10000L
n_expression_bins <- 10L
seed <- 20260718L

# ---- Read DESeq2 result files robustly -------------------------------------
# The exported tables may contain a leading blank column name.  This reader
# repairs blank/NA column names before selecting the required DESeq2 columns.
read_deseq_results <- function(path, label) {
  if (!file.exists(path)) stop("Input file was not found: ", path)

  dat <- read.delim(
    path,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  bad_names <- is.na(names(dat)) | names(dat) == ""
  names(dat)[bad_names] <- paste0("input_index_", which(bad_names))
  names(dat) <- make.unique(names(dat))

  required <- c("Row.names", "baseMean", "log2FoldChange", "padj")
  absent <- setdiff(required, names(dat))
  if (length(absent) > 0) {
    stop(
      "The following required columns are absent from ", basename(path), ": ",
      paste(absent, collapse = ", ")
    )
  }

  out <- dat %>%
    transmute(
      geneId = as.character(.data$Row.names),
      baseMean = as.numeric(.data$baseMean),
      log2FC = as.numeric(.data$log2FoldChange),
      padj = as.numeric(.data$pvalue)
    ) %>%
    filter(!is.na(geneId), geneId != "", is.finite(log2FC)) %>%
    distinct(geneId, .keep_all = TRUE)

  names(out)[2:4] <- paste0(c("baseMean_", "log2FC_", "padj_"), label)
  out
}

hd_wt <- read_deseq_results(hd_file, "HD_vs_WT")
oe_hd <- read_deseq_results(oe_file, "OE_vs_HD")

tested <- inner_join(hd_wt, oe_hd, by = "geneId") %>%
  mutate(
    HD_DEG = !is.na(.data$padj_HD_vs_WT) &
      .data$padj_HD_vs_WT <= fdr_cutoff &
      abs(.data$log2FC_HD_vs_WT) >= lfc_cutoff,
    OE_responsive = !is.na(.data$padj_OE_vs_HD) &
      .data$padj_OE_vs_HD <= fdr_cutoff &
      abs(.data$log2FC_OE_vs_HD) >= lfc_cutoff,
    log2FC_OE_vs_WT = .data$log2FC_HD_vs_WT + .data$log2FC_OE_vs_HD,
    distance_to_WT_before = abs(.data$log2FC_HD_vs_WT),
    distance_to_WT_after = abs(.data$log2FC_OE_vs_WT),
    opposite_direction = sign(.data$log2FC_HD_vs_WT) *
      sign(.data$log2FC_OE_vs_HD) < 0,
    closer_to_WT = .data$distance_to_WT_after < .data$distance_to_WT_before,
    # Primary, effect-size-based rescue definition.
    true_rescue = .data$HD_DEG & .data$opposite_direction & .data$closer_to_WT,
    # More conservative version: the OE-vs-HD change also passes the DEG filter.
    supported_true_rescue = .data$true_rescue & .data$OE_responsive,
    rescue_status = case_when(
      !.data$HD_DEG ~ "not_HD_vs_WT_DEG",
      .data$supported_true_rescue ~ "supported_true_rescue",
      .data$true_rescue ~ "true_rescue_not_OE_significant",
      .data$opposite_direction & !.data$closer_to_WT ~ "opposite_but_overshot_WT",
      .data$opposite_direction ~ "opposite_without_closer_to_WT",
      TRUE ~ "not_directionally_rescued"
    )
  )

hd_deg <- tested %>% filter(.data$HD_DEG)
if (nrow(hd_deg) < 3L) {
  stop("Fewer than three HD-vs-WT DEGs met the selected thresholds; correlation is not meaningful.")
}

# ---- Spearman association ---------------------------------------------------
# Negative rho indicates that OE changes tend to oppose HD-vs-WT changes.  It
# does not alone establish rescue because it does not assess distance to WT.
spearman <- cor.test(
  hd_deg$log2FC_HD_vs_WT,
  hd_deg$log2FC_OE_vs_HD,
  method = "spearman",
  exact = FALSE
)

# ---- Recovery proportions and directional enrichment ----------------------
n_hd_deg <- nrow(hd_deg)
n_true_rescue <- sum(hd_deg$true_rescue)
n_supported_rescue <- sum(hd_deg$supported_true_rescue)
n_directional <- sum(hd_deg$opposite_direction)

true_rescue_ci <- binom.test(n_true_rescue, n_hd_deg)$conf.int
supported_rescue_ci <- binom.test(n_supported_rescue, n_hd_deg)$conf.int

# Among HD DEGs with a nonzero OE effect, ask whether opposite-direction
# changes are enriched beyond the two-sided null expectation of 50%.
directional_test_set <- hd_deg %>% filter(.data$log2FC_OE_vs_HD != 0)
directional_enrichment <- binom.test(
  sum(directional_test_set$opposite_direction),
  nrow(directional_test_set),
  p = 0.5,
  alternative = "greater"
)

# ---- Expression-stratified permutation test --------------------------------
# OE-vs-HD log2FC values are shuffled only among genes in the same HD-vs-WT
# baseMean decile.  This controls the null distribution for broad dependence
# of observed effect sizes on expression level, while preserving HD-vs-WT LFCs.
if (nrow(tested) < n_expression_bins) {
  stop("Too few genes for the requested number of expression bins.")
}

set.seed(seed)
expression_rank <- rank(
  log10(pmax(tested$baseMean_HD_vs_WT, 0) + 1),
  ties.method = "random"
)
tested$expression_bin <- pmin(
  n_expression_bins,
  ceiling(expression_rank / (length(expression_rank) / n_expression_bins))
)

hd_index <- which(tested$HD_DEG)
observed_rescue_count <- sum(tested$true_rescue[hd_index])
null_rescue_count <- numeric(n_permutations)

for (i in seq_len(n_permutations)) {
  permuted_oe_lfc <- tested$log2FC_OE_vs_HD
  for (bin in seq_len(n_expression_bins)) {
    index_in_bin <- which(tested$expression_bin == bin)
    permuted_oe_lfc[index_in_bin] <- sample(
      permuted_oe_lfc[index_in_bin],
      size = length(index_in_bin),
      replace = FALSE
    )
  }

  hd_lfc <- tested$log2FC_HD_vs_WT[hd_index]
  oe_lfc <- permuted_oe_lfc[hd_index]
  is_true_rescue <- sign(hd_lfc) * sign(oe_lfc) < 0 &
    abs(hd_lfc + oe_lfc) < abs(hd_lfc)
  null_rescue_count[i] <- sum(is_true_rescue)
}

permutation_p <- (1 + sum(null_rescue_count >= observed_rescue_count)) /
  (n_permutations + 1)
permutation_z <- if (sd(null_rescue_count) > 0) {
  (observed_rescue_count - mean(null_rescue_count)) / sd(null_rescue_count)
} else {
  NA_real_
}

# ---- Export tables and report ----------------------------------------------
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Spearman correlation plot ---------------------------------------------
# Each point is an HD-vs-WT DEG. The dashed diagonal y = -x denotes the OE
# effect that would exactly restore the gene's expression to the WT level.
scatter_df <- hd_deg %>%
  mutate(
    rescue_plot_group = case_when(
      .data$supported_true_rescue ~ "Supported true rescue",
      .data$true_rescue ~ "True rescue",
      TRUE ~ "Not directionally rescued"
    ),
    rescue_plot_group = factor(
      .data$rescue_plot_group,
      levels = c(
        "Not directionally rescued",
        "True rescue",
        "Supported true rescue"
      )
    )
  )

# A common x/y limit makes the y = -x line visually interpretable.
plot_limit <- max(
  abs(c(scatter_df$log2FC_HD_vs_WT, scatter_df$log2FC_OE_vs_HD)),
  na.rm = TRUE
) * 1.06
if (!is.finite(plot_limit) || plot_limit == 0) plot_limit <- 1

rho_label <- sprintf("Spearman rho = %.2f\nP = %s\nn = %d", 
                     unname(spearman$estimate),
                     format.pval(spearman$p.value, digits = 2, eps = 1e-300),
                     nrow(scatter_df))

rescue_scatter <- ggplot(
  scatter_df,
  aes(
    x = .data$log2FC_HD_vs_WT,
    y = .data$log2FC_OE_vs_HD,
    colour = .data$rescue_plot_group
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "#8C8C8C") +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#8C8C8C") +
  geom_abline(
    slope = -1,
    intercept = 0,
    linewidth = 0.55,
    linetype = "dashed",
    colour = "#404040"
  ) +
  geom_point(size = 1.7, alpha = 0.72) +
  annotate(
    "label",
    x = -9,
    y = 9,
    label = rho_label,
    hjust = 0,
    vjust = 1,
    size = 3.4,
    label.size = 0.2,
    fill = "white",
    colour = "black"
  ) +
  #annotate(
  #  "text",
  #  x = 0.67 * plot_limit,
  #  y = -0.72 * plot_limit,
  #  label = "complete return to WT",
  #  angle = -45,
  #  hjust = 0,
  #  size = 3.1,
  #  colour = "#404040"
  #) +
  scale_colour_manual(
    values = c(
      "Not directionally rescued" = "#8C8C8C",
      "True rescue" = "#D64B4B",
      "Supported true rescue" = "#2A6FBB"
    )
  ) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10)) +
  labs(
    x = "HD vs WT log2FC",
    y = "OE vs HD log2FC",
    colour = NULL
  ) +
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks = element_line(linewidth = 0.45, colour = "black"),
    legend.position = c(0.03, 0.03),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_blank(),
    legend.text = element_text(size = 8),
    plot.margin = margin(6, 8, 6, 6)
  )

write.csv(
  scatter_df,
  file.path(output_dir, "Spearman_correlation_plot_source_data.csv"),
  row.names = FALSE
)
ggsave(
  file.path(output_dir, "HD_vs_WT_vs_OE_vs_HD_Spearman_scatter.png"),
  rescue_scatter, width = 89, height = 82, units = "mm", dpi = 600
)
ggsave(
  file.path(output_dir, "HD_vs_WT_vs_OE_vs_HD_Spearman_scatter.pdf"),
  rescue_scatter, width = 89, height = 82, units = "mm", device = grDevices::cairo_pdf
)
ggsave(
  file.path(output_dir, "HD_vs_WT_vs_OE_vs_HD_Spearman_scatter.svg"),
  rescue_scatter, width = 89, height = 82, units = "mm", device = svglite::svglite
)

write.csv(
  tested %>% arrange(desc(.data$HD_DEG), .data$padj_HD_vs_WT),
  file.path(output_dir, "all_tested_genes_with_rescue_status.csv"),
  row.names = FALSE
)
write.csv(
  hd_deg %>% arrange(desc(.data$true_rescue), .data$padj_HD_vs_WT),
  file.path(output_dir, "HD_vs_WT_DEGs_with_rescue_status.csv"),
  row.names = FALSE
)

summary_table <- data.frame(
  metric = c(
    "tested_genes",
    "HD_vs_WT_DEGs",
    "true_rescue_genes",
    "true_rescue_proportion",
    "true_rescue_95CI_lower",
    "true_rescue_95CI_upper",
    "supported_true_rescue_genes",
    "supported_true_rescue_proportion",
    "supported_true_rescue_95CI_lower",
    "supported_true_rescue_95CI_upper",
    "Spearman_rho_HD_vs_WT_vs_OE_vs_HD",
    "Spearman_p_value",
    "directional_opposite_genes",
    "directional_test_genes",
    "directional_binomial_p_value",
    "permutation_observed_true_rescue_genes",
    "permutation_null_mean",
    "permutation_null_sd",
    "permutation_z_score",
    "permutation_one_sided_p_value"
  ),
  value = c(
    nrow(tested), n_hd_deg, n_true_rescue,
    n_true_rescue / n_hd_deg, true_rescue_ci[1], true_rescue_ci[2],
    n_supported_rescue, n_supported_rescue / n_hd_deg,
    supported_rescue_ci[1], supported_rescue_ci[2],
    unname(spearman$estimate), spearman$p.value,
    sum(directional_test_set$opposite_direction), nrow(directional_test_set),
    directional_enrichment$p.value,
    observed_rescue_count, mean(null_rescue_count), sd(null_rescue_count),
    permutation_z, permutation_p
  )
)
write.csv(summary_table, file.path(output_dir, "rescue_summary.csv"), row.names = FALSE)
write.csv(
  data.frame(permutation = seq_len(n_permutations), true_rescue_count = null_rescue_count),
  file.path(output_dir, "permutation_null_distribution.csv"),
  row.names = FALSE
)

report_lines <- c(
  "HD -> WT rescue analysis",
  paste0("Input HD-vs-WT: ", hd_file),
  paste0("Input OE-vs-HD: ", oe_file),
  paste0("DEG definition: padj <= ", fdr_cutoff, " and |log2FC| >= ", lfc_cutoff),
  "True rescue: HD DEG; OE-vs-HD effect is opposite in direction; and |HD-vs-WT + OE-vs-HD| < |HD-vs-WT|.",
  "",
  paste0("Tested genes: ", nrow(tested)),
  paste0("HD-vs-WT DEGs: ", n_hd_deg),
  paste0(
    "Spearman correlation (HD-vs-WT LFC vs OE-vs-HD LFC): rho = ",
    formatC(unname(spearman$estimate), format = "f", digits = 3),
    ", P = ", format.pval(spearman$p.value, digits = 3, eps = 1e-300)
  ),
  paste0(
    "True rescue genes: ", n_true_rescue, "/", n_hd_deg,
    " (", formatC(100 * n_true_rescue / n_hd_deg, format = "f", digits = 1),
    "%; exact 95% CI ", formatC(100 * true_rescue_ci[1], format = "f", digits = 1),
    "–", formatC(100 * true_rescue_ci[2], format = "f", digits = 1), "%)"
  ),
  paste0(
    "Supported true rescue genes (also OE-vs-HD DEGs): ", n_supported_rescue,
    "/", n_hd_deg, " (",
    formatC(100 * n_supported_rescue / n_hd_deg, format = "f", digits = 1),
    "%; exact 95% CI ",
    formatC(100 * supported_rescue_ci[1], format = "f", digits = 1), "–",
    formatC(100 * supported_rescue_ci[2], format = "f", digits = 1), "%)"
  ),
  paste0(
    "Opposite-direction OE changes: ", sum(directional_test_set$opposite_direction),
    "/", nrow(directional_test_set), "; one-sided binomial P = ",
    format.pval(directional_enrichment$p.value, digits = 3, eps = 1e-300)
  ),
  paste0(
    "Expression-stratified permutation test (", n_permutations,
    " permutations): observed = ", observed_rescue_count,
    ", null mean = ", formatC(mean(null_rescue_count), format = "f", digits = 2),
    ", null SD = ", formatC(sd(null_rescue_count), format = "f", digits = 2),
    ", Z = ", formatC(permutation_z, format = "f", digits = 2),
    ", one-sided P = ", format.pval(permutation_p, digits = 3, eps = 1e-300)
  )
)
writeLines(report_lines, file.path(output_dir, "rescue_report.txt"))

message("Analysis complete. Results written to: ", output_dir)
