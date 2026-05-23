# =============================================================================
# Phase 6: Spearman Anticorrelation Analysis — tRFs vs KIRP-Relevant Genes
# KIRP tRF Project
#
# Strategy:
#   - Use matched tumor samples (patient present in both tRF and mRNA data)
#   - Use normalized log-CPM values for correlation (not raw counts)
#   - Test all significant tRF × significant gene pairs
#   - Filter: Spearman r < 0, BH-adjusted p < 0.05
#   - Prioritize: tRF up + gene down (or tRF down + gene up) with known RCC biology
#
# Input:
#   data_processed/trf_count_matrix.tsv
#   data_processed/trf_sample_metadata.tsv
#   data_processed/mrna_count_matrix.tsv
#   data_processed/mrna_sample_metadata.tsv
#   results/tables/trf_DE_significant.tsv
#   results/tables/gene_DE_significant.tsv
#
# Output:
#   results/tables/anticorrelation_full.tsv
#   results/tables/anticorrelation_significant.tsv
#   results/tables/anticorrelation_prioritized.tsv
#   results/figures/anticorrelation_heatmap.pdf
#   results/figures/anticorrelation_scatterplots_top.pdf
# =============================================================================

library(edgeR)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)

# ── Paths ─────────────────────────────────────────────────────────────────────
# Edit scripts/config.R to set ROOT for your machine, or change the fallback below:
if (file.exists("scripts/config.R")) source("scripts/config.R") else if (file.exists("config.R")) source("config.R") else ROOT <- "/Users/nam.tnguyen2022/KIRP_tRF_project"
TRF_COUNT   <- file.path(ROOT, "data_processed", "trf_count_matrix.tsv")
TRF_META    <- file.path(ROOT, "data_processed", "trf_sample_metadata.tsv")
MRNA_COUNT  <- file.path(ROOT, "data_processed", "mrna_count_matrix.tsv")
MRNA_META   <- file.path(ROOT, "data_processed", "mrna_sample_metadata.tsv")
SIG_TRF     <- file.path(ROOT, "results", "tables", "trf_DE_significant.tsv")
SIG_GENE    <- file.path(ROOT, "results", "tables", "gene_DE_significant.tsv")
TABLE_DIR   <- file.path(ROOT, "results", "tables")
FIG_DIR     <- file.path(ROOT, "results", "figures")

# ── Helper: normalize a count matrix to log-CPM ───────────────────────────────
norm_lcpm <- function(counts_df) {
  dge  <- DGEList(counts = counts_df)
  dge  <- calcNormFactors(dge, method = "TMM")
  cpm(dge, log = TRUE, prior.count = 1)
}

# ── 1. Load metadata and find matched tumor patients ─────────────────────────
cat("Loading metadata ...\n")
trf_meta  <- read.table(TRF_META,  sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE,
                        colClasses = c(sample_type_code = "character"))
mrna_meta <- read.table(MRNA_META, sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE,
                        colClasses = c(sample_type_code = "character"))

# 12-char patient ID
trf_meta$patient_id  <- substr(trf_meta$sample_id,  1, 12)
mrna_meta$patient_id <- substr(mrna_meta$sample_id, 1, 12)

# Matched tumor patients
trf_tumor  <- trf_meta[trf_meta$sample_type_code == "01", ]
mrna_tumor <- mrna_meta[mrna_meta$sample_type_code == "01", ]
matched_patients <- intersect(trf_tumor$patient_id, mrna_tumor$patient_id)

cat(sprintf("Matched tumor patients: %d\n", length(matched_patients)))

# One sample per patient (take first if duplicates)
trf_matched  <- trf_tumor[match(matched_patients, trf_tumor$patient_id), ]
mrna_matched <- mrna_tumor[match(matched_patients, mrna_tumor$patient_id), ]

# ── 2. Load count matrices and extract matched samples ────────────────────────
cat("Loading count matrices ...\n")
trf_counts  <- read.table(TRF_COUNT,  sep = "\t", header = TRUE,
                          row.names = 1, check.names = FALSE)
mrna_counts <- read.table(MRNA_COUNT, sep = "\t", header = TRUE,
                          row.names = 1, check.names = FALSE)

# Subset to matched samples
trf_mat  <- trf_counts[,  trf_matched$sample_id[trf_matched$sample_id  %in% colnames(trf_counts)]]
mrna_mat <- mrna_counts[, mrna_matched$sample_id[mrna_matched$sample_id %in% colnames(mrna_counts)]]

cat(sprintf("tRF matrix:  %d tRFs x %d samples\n",  nrow(trf_mat),  ncol(trf_mat)))
cat(sprintf("mRNA matrix: %d genes x %d samples\n", nrow(mrna_mat), ncol(mrna_mat)))

# ── 3. Load significant DE results ───────────────────────────────────────────
sig_trf  <- read.table(SIG_TRF,  sep = "\t", header = TRUE)
sig_gene <- read.table(SIG_GENE, sep = "\t", header = TRUE)

cat(sprintf("Significant tRFs: %d  |  Significant genes: %d\n",
            nrow(sig_trf), nrow(sig_gene)))
cat(sprintf("Total pairs to test: %d\n", nrow(sig_trf) * nrow(sig_gene)))

# ── 4. Normalize to log-CPM (matched samples only) ───────────────────────────
cat("Normalizing to log-CPM ...\n")

# tRF: subset to significant tRFs present in matrix
sig_trfs_present <- sig_trf$trf_id[sig_trf$trf_id %in% rownames(trf_mat)]
trf_sub  <- trf_mat[sig_trfs_present, ]

# Remove all-zero rows before normalization
trf_sub  <- trf_sub[rowSums(trf_sub) > 0, ]
trf_lcpm <- norm_lcpm(trf_sub)

# mRNA: subset to significant genes, align to same samples
# Use Ensembl gene_id → need to check if rownames are Ensembl or gene names
# (mrna_count_matrix uses Ensembl IDs as rownames)
# gene_DE_significant uses gene_name; we need the Ensembl ID for lookup
# Load gene map from the DE results (we have logFC etc but need Ensembl IDs)
# Workaround: load mrna_count_matrix column names and search by gene_name stored in meta
# We actually stored gene_names as rownames of the sig gene table, so we need to
# look up Ensembl IDs from the full mRNA matrix (which uses Ensembl IDs)
# The mrna_sample_metadata does not have gene info, so re-use the STAR counts file.

GDC_DIR <- file.path(ROOT, "data_raw", "gdc_rnaseq",
                     "gdc_download_20260506_091015.911939")
uuid_dirs <- list.dirs(GDC_DIR, full.names = TRUE, recursive = FALSE)
ref_tsv   <- NULL
for (d in uuid_dirs) {
  tsvs <- list.files(d, pattern = "\\.tsv$", full.names = TRUE)
  if (length(tsvs) > 0) { ref_tsv <- tsvs[1]; break }
}
gene_map <- read.table(ref_tsv, sep = "\t", header = TRUE, comment.char = "#",
                       stringsAsFactors = FALSE)
gene_map <- gene_map[grepl("^ENSG", gene_map$gene_id),
                     c("gene_id", "gene_name")]
gene_map <- gene_map[!duplicated(gene_map$gene_name), ]

# Map sig gene names -> Ensembl IDs
target_ensembl <- gene_map$gene_id[gene_map$gene_name %in% sig_gene$gene_name]
names(target_ensembl) <- gene_map$gene_name[gene_map$gene_name %in% sig_gene$gene_name]

mrna_sub <- mrna_mat[target_ensembl[target_ensembl %in% rownames(mrna_mat)], ]

# Replace Ensembl IDs with gene names
gene_name_lookup <- setNames(names(target_ensembl), target_ensembl)
rownames(mrna_sub) <- gene_name_lookup[rownames(mrna_sub)]

mrna_lcpm <- norm_lcpm(mrna_sub)

cat(sprintf("tRF normalized matrix:  %d x %d\n", nrow(trf_lcpm),  ncol(trf_lcpm)))
cat(sprintf("mRNA normalized matrix: %d x %d\n", nrow(mrna_lcpm), ncol(mrna_lcpm)))

# ── 5. Align samples between tRF and mRNA matrices ───────────────────────────
# Samples are named by TCGA barcode (full aliquot for tRF, 16-char for mRNA)
# Match by 12-char patient ID
trf_patients  <- substr(colnames(trf_lcpm),  1, 12)
mrna_patients <- substr(colnames(mrna_lcpm), 1, 12)
common_patients <- intersect(trf_patients, mrna_patients)

trf_aligned  <- trf_lcpm[,  match(common_patients, trf_patients)]
mrna_aligned <- mrna_lcpm[, match(common_patients, mrna_patients)]

# Ensure column count matches
stopifnot(ncol(trf_aligned) == ncol(mrna_aligned))
cat(sprintf("Aligned samples for correlation: %d\n", ncol(trf_aligned)))

# ── 6. Compute Spearman correlations for all tRF × gene pairs ─────────────────
cat("Computing Spearman correlations for all pairs ...\n")
n_trf  <- nrow(trf_aligned)
n_gene <- nrow(mrna_aligned)
n_pairs <- n_trf * n_gene

results <- data.frame(
  trf      = character(n_pairs),
  gene     = character(n_pairs),
  spearman_r = numeric(n_pairs),
  p_value  = numeric(n_pairs),
  stringsAsFactors = FALSE
)

idx <- 1
for (g in rownames(mrna_aligned)) {
  gene_expr <- as.numeric(mrna_aligned[g, ])
  for (t in rownames(trf_aligned)) {
    trf_expr <- as.numeric(trf_aligned[t, ])
    ct <- suppressWarnings(cor.test(trf_expr, gene_expr,
                                    method = "spearman", exact = FALSE))
    results[idx, ] <- list(t, g, ct$estimate, ct$p.value)
    idx <- idx + 1
  }
  cat(sprintf("  Completed gene: %s (%d / %d pairs done)\n",
              g, idx - 1, n_pairs))
}

results$spearman_r <- as.numeric(results$spearman_r)
results$p_value    <- as.numeric(results$p_value)

# BH FDR correction
results$FDR <- p.adjust(results$p_value, method = "BH")

# ── 7. Merge DE statistics ────────────────────────────────────────────────────
trf_lfc  <- setNames(sig_trf$logFC,  sig_trf$trf_id)
gene_lfc <- setNames(sig_gene$logFC, sig_gene$gene_name)
trf_fdr  <- setNames(sig_trf$FDR,   sig_trf$trf_id)
gene_fdr <- setNames(sig_gene$FDR,  sig_gene$gene_name)

results$trf_logFC  <- trf_lfc[results$trf]
results$gene_logFC <- gene_lfc[results$gene]
results$trf_FDR    <- trf_fdr[results$trf]
results$gene_FDR   <- gene_fdr[results$gene]

# Direction classification
results$trf_direction  <- ifelse(results$trf_logFC  > 0, "Up", "Down")
results$gene_direction <- ifelse(results$gene_logFC > 0, "Up", "Down")
results$pair_type <- paste0("tRF_", results$trf_direction,
                            "__gene_", results$gene_direction)

# Sort by Spearman r (most negative first)
results <- results[order(results$spearman_r), ]

# Save full results
write.table(results,
            file.path(TABLE_DIR, "anticorrelation_full.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nSaved full anticorrelation table -> %s\n",
            file.path(TABLE_DIR, "anticorrelation_full.tsv")))

# ── 8. Filter for significant anticorrelations ────────────────────────────────
anti_sig <- results[results$spearman_r < 0 & results$FDR < 0.05, ]
anti_sig <- anti_sig[order(anti_sig$spearman_r), ]

cat(sprintf("Significant anticorrelations (r < 0, FDR < 0.05): %d pairs\n",
            nrow(anti_sig)))
cat(sprintf("  tRF up + gene down: %d\n",
            sum(anti_sig$pair_type == "tRF_Up__gene_Down")))
cat(sprintf("  tRF down + gene up: %d\n",
            sum(anti_sig$pair_type == "tRF_Down__gene_Up")))

write.table(anti_sig,
            file.path(TABLE_DIR, "anticorrelation_significant.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Saved significant anticorrelations -> %s\n",
            file.path(TABLE_DIR, "anticorrelation_significant.tsv")))

# ── 9. Prioritized pairs ──────────────────────────────────────────────────────
# Priority 1: tRF up + gene down (tRF may suppress gene)
# Priority 2: tRF down + gene up (loss of tRF may allow gene overexpression)
# Priority 3: strongest r regardless of direction
anti_sig$priority <- ifelse(anti_sig$pair_type == "tRF_Up__gene_Down",   1,
                     ifelse(anti_sig$pair_type == "tRF_Down__gene_Up",   2, 3))
prioritized <- anti_sig[order(anti_sig$priority, anti_sig$spearman_r), ]

write.table(prioritized,
            file.path(TABLE_DIR, "anticorrelation_prioritized.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Saved prioritized pairs -> %s\n",
            file.path(TABLE_DIR, "anticorrelation_prioritized.tsv")))

# Print top 20
cat("\nTop anticorrelated pairs (prioritized):\n")
print(head(prioritized[, c("trf","gene","spearman_r","FDR",
                            "trf_logFC","gene_logFC","pair_type")], 20))

# ── 10. Anticorrelation heatmap (top pairs per gene) ─────────────────────────
cat("\nGenerating anticorrelation heatmap ...\n")
if (nrow(anti_sig) >= 2) {
  # Build r matrix: top 50 anticorrelated tRFs x all sig genes
  top_trfs <- head(unique(prioritized$trf), 50)
  r_mat <- matrix(NA, nrow = length(top_trfs), ncol = nrow(mrna_aligned),
                  dimnames = list(top_trfs, rownames(mrna_aligned)))
  for (g in rownames(mrna_aligned)) {
    for (t in top_trfs) {
      val <- results$spearman_r[results$trf == t & results$gene == g]
      if (length(val) > 0) r_mat[t, g] <- val
    }
  }
  r_mat[is.na(r_mat)] <- 0

  ann_row <- data.frame(
    Direction = ifelse(trf_lfc[top_trfs] > 0, "Up in Tumor", "Down in Tumor")
  )
  rownames(ann_row) <- top_trfs

  ann_colors <- list(
    Direction = c("Up in Tumor"   = "#D73027",
                  "Down in Tumor" = "#4575B4")
  )

  pdf(file.path(FIG_DIR, "anticorrelation_heatmap.pdf"), width = 7, height = 14)
  pheatmap(r_mat,
           annotation_row   = ann_row,
           annotation_colors = ann_colors,
           color            = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
           breaks           = seq(-1, 1, length.out = 101),
           main             = "Spearman r: Top Anticorrelated tRF-Gene Pairs",
           fontsize_row     = 6,
           fontsize_col     = 10,
           cluster_cols     = FALSE,
           cluster_rows     = TRUE,
           border_color     = NA,
           display_numbers  = TRUE,
           number_format    = "%.2f",
           fontsize_number  = 6)
  dev.off()
  cat("  Saved: anticorrelation_heatmap.pdf\n")
}

# ── 11. Scatterplots for top anticorrelated pairs ─────────────────────────────
cat("Generating scatterplots for top anticorrelated pairs ...\n")

n_scatter <- min(12, nrow(prioritized))
top_pairs <- head(prioritized, n_scatter)

scatter_list <- lapply(seq_len(nrow(top_pairs)), function(i) {
  trf_id   <- top_pairs$trf[i]
  gene_id  <- top_pairs$gene[i]
  r_val    <- round(top_pairs$spearman_r[i], 3)
  fdr_val  <- signif(top_pairs$FDR[i], 3)
  trf_fc   <- round(top_pairs$trf_logFC[i], 2)
  gene_fc  <- round(top_pairs$gene_logFC[i], 2)
  ptype    <- top_pairs$pair_type[i]

  df_plot <- data.frame(
    trf_expr  = as.numeric(trf_aligned[trf_id, ]),
    gene_expr = as.numeric(mrna_aligned[gene_id, ])
  )

  ggplot(df_plot, aes(x = trf_expr, y = gene_expr)) +
    geom_point(alpha = 0.5, size = 1.2, color = "#555555") +
    geom_smooth(method = "lm", se = TRUE, color = "#D73027",
                linewidth = 0.8, fill = "#FFCCCC") +
    labs(
      title    = sprintf("%s  x  %s", trf_id, gene_id),
      subtitle = sprintf("r = %s | FDR = %s | tRF logFC = %s | gene logFC = %s\n%s",
                         r_val, fdr_val, trf_fc, gene_fc, ptype),
      x = sprintf("%s (log2 CPM)", trf_id),
      y = sprintf("%s (log2 CPM)", gene_id)
    ) +
    theme_bw(base_size = 8) +
    theme(plot.title    = element_text(size = 7, face = "bold"),
          plot.subtitle = element_text(size = 6))
})

n_cols_s <- min(3, n_scatter)
n_rows_s <- ceiling(n_scatter / n_cols_s)
pdf(file.path(FIG_DIR, "anticorrelation_scatterplots_top.pdf"),
    width = n_cols_s * 4, height = n_rows_s * 3.5)
grid.arrange(grobs = scatter_list, ncol = n_cols_s)
dev.off()
cat(sprintf("  Saved: anticorrelation_scatterplots_top.pdf  (%d pairs)\n", n_scatter))

cat("\n✅ Phase 6 (anticorrelation analysis) complete.\n")
