# =============================================================================
# Phase 4: Differential Expression of tRFs using edgeR
# KIRP tRF Project
#
# Input:
#   data_processed/trf_count_matrix.tsv
#   data_processed/trf_sample_metadata.tsv
#
# Output:
#   results/tables/trf_DE_full.tsv         (all tRFs)
#   results/tables/trf_DE_significant.tsv  (FDR < 0.05, |logFC| >= 1)
#   results/figures/trf_volcano.pdf
#   results/figures/trf_boxplots_top10.pdf
#   results/figures/trf_heatmap_top50.pdf
# =============================================================================

library(edgeR)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)

# ── Paths ─────────────────────────────────────────────────────────────────────
# Edit scripts/config.R to set ROOT for your machine, or change the fallback below:
if (file.exists("scripts/config.R")) source("scripts/config.R") else if (file.exists("config.R")) source("config.R") else ROOT <- "/Users/nam.tnguyen2022/KIRP_tRF_project"

COUNT_FILE  <- file.path(ROOT, "data_processed", "trf_count_matrix.tsv")
META_FILE   <- file.path(ROOT, "data_processed", "trf_sample_metadata.tsv")
TABLE_DIR   <- file.path(ROOT, "results", "tables")
FIG_DIR     <- file.path(ROOT, "results", "figures")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,   recursive = TRUE, showWarnings = FALSE)

# ── 1. Load data ──────────────────────────────────────────────────────────────
cat("Loading count matrix ...\n")
counts <- read.table(COUNT_FILE, sep = "\t", header = TRUE,
                     row.names = 1, check.names = FALSE)

meta <- read.table(META_FILE, sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE,
                   colClasses = c(sample_type_code = "character"))

# Align samples
common_samples <- intersect(colnames(counts), meta$sample_id)
counts <- counts[, common_samples]
meta   <- meta[meta$sample_id %in% common_samples, ]
meta   <- meta[match(colnames(counts), meta$sample_id), ]

cat(sprintf("Samples: %d tumor, %d normal\n",
            sum(meta$sample_type_code == "01"),
            sum(meta$sample_type_code == "11")))
cat(sprintf("tRFs before filtering: %d\n", nrow(counts)))

# ── 2. Create DGEList and filter low-count tRFs ───────────────────────────────
group <- factor(ifelse(meta$sample_type_code == "01", "Tumor", "Normal"),
                levels = c("Normal", "Tumor"))  # Normal = reference

dge <- DGEList(counts = counts, group = group)

# Filter: keep tRFs with CPM > 1 in at least (min group size) samples
min_n     <- min(table(group))
keep      <- rowSums(cpm(dge) > 1) >= min_n
dge       <- dge[keep, , keep.lib.sizes = FALSE]
cat(sprintf("tRFs after CPM > 1 filter: %d\n", nrow(dge)))

# ── 3. Normalize (TMM) ────────────────────────────────────────────────────────
dge <- calcNormFactors(dge, method = "TMM")

# ── 4. Design matrix & dispersion ─────────────────────────────────────────────
design <- model.matrix(~ group)
dge    <- estimateDisp(dge, design, robust = TRUE)

cat(sprintf("Common dispersion: %.4f\n", dge$common.dispersion))

# ── 5. Fit GLM and test (Tumor vs Normal) ─────────────────────────────────────
fit  <- glmQLFit(dge, design, robust = TRUE)
qlf  <- glmQLFTest(fit, coef = 2)   # coef 2 = groupTumor

# ── 6. Extract results with BH FDR correction ─────────────────────────────────
res_full <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "PValue")$table
res_full$trf_id <- rownames(res_full)
res_full <- res_full[, c("trf_id", "logFC", "logCPM", "F", "PValue", "FDR")]

cat(sprintf("\nResults summary:\n"))
cat(sprintf("  Total tested:             %d\n",    nrow(res_full)))
cat(sprintf("  FDR < 0.05:               %d\n",    sum(res_full$FDR < 0.05)))
cat(sprintf("  FDR < 0.05 & |logFC|>=1:  %d\n",
            sum(res_full$FDR < 0.05 & abs(res_full$logFC) >= 1)))
cat(sprintf("  Upregulated in tumor:     %d\n",
            sum(res_full$FDR < 0.05 & res_full$logFC >= 1)))
cat(sprintf("  Downregulated in tumor:   %d\n",
            sum(res_full$FDR < 0.05 & res_full$logFC <= -1)))

# ── 7. Save tables ─────────────────────────────────────────────────────────────
write.table(res_full,
            file = file.path(TABLE_DIR, "trf_DE_full.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

res_sig <- res_full[res_full$FDR < 0.05 & abs(res_full$logFC) >= 1, ]
res_sig <- res_sig[order(res_sig$FDR), ]
write.table(res_sig,
            file = file.path(TABLE_DIR, "trf_DE_significant.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\nSaved full results    -> %s\n", file.path(TABLE_DIR, "trf_DE_full.tsv")))
cat(sprintf("Saved sig. results    -> %s\n",  file.path(TABLE_DIR, "trf_DE_significant.tsv")))

# ── 8. Volcano plot ───────────────────────────────────────────────────────────
cat("Generating volcano plot ...\n")
vol_df <- res_full
vol_df$sig <- "Not significant"
vol_df$sig[vol_df$FDR < 0.05 & vol_df$logFC >= 1]  <- "Up in Tumor"
vol_df$sig[vol_df$FDR < 0.05 & vol_df$logFC <= -1] <- "Down in Tumor"
vol_df$sig <- factor(vol_df$sig,
                     levels = c("Up in Tumor", "Down in Tumor", "Not significant"))

# Label top 10 by FDR
top_labels <- head(res_sig[order(res_sig$FDR), "trf_id"], 10)
vol_df$label <- ifelse(vol_df$trf_id %in% top_labels, vol_df$trf_id, NA)

p_volcano <- ggplot(vol_df, aes(x = logFC, y = -log10(FDR),
                                 color = sig, label = label)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_text_repel(size = 2.8, max.overlaps = 20, na.rm = TRUE,
                  segment.color = "grey50") +
  scale_color_manual(values = c("Up in Tumor"    = "#D73027",
                                "Down in Tumor"  = "#4575B4",
                                "Not significant"= "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(title    = "Differential tRF Expression: KIRP Tumor vs Normal",
       subtitle = sprintf("FDR < 0.05, |logFC| >= 1  |  Up: %d  Down: %d",
                          sum(vol_df$sig == "Up in Tumor"),
                          sum(vol_df$sig == "Down in Tumor")),
       x = "log2 Fold Change (Tumor / Normal)",
       y = "-log10(FDR)",
       color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

pdf(file.path(FIG_DIR, "trf_volcano.pdf"), width = 8, height = 6)
print(p_volcano)
dev.off()
cat("  Saved: trf_volcano.pdf\n")

# ── 9. Boxplots for top 10 significant tRFs ───────────────────────────────────
cat("Generating boxplots for top 10 tRFs ...\n")

# Get normalized log-CPM values
lcpm <- cpm(dge, log = TRUE, prior.count = 1)

top10 <- head(res_sig[order(res_sig$FDR), "trf_id"], 10)

plot_list <- lapply(top10, function(trf) {
  expr_vals <- lcpm[trf, ]
  df_plot <- data.frame(
    sample   = names(expr_vals),
    logCPM   = as.numeric(expr_vals),
    group    = group
  )
  direction <- ifelse(res_sig[trf, "logFC"] > 0, "Up", "Down")
  fdr_val   <- signif(res_sig[trf, "FDR"], 3)
  lfc_val   <- round(res_sig[trf, "logFC"], 2)

  ggplot(df_plot, aes(x = group, y = logCPM, fill = group)) +
    geom_boxplot(outlier.size = 0.8, width = 0.5) +
    geom_jitter(width = 0.15, size = 0.6, alpha = 0.5) +
    scale_fill_manual(values = c("Normal" = "#4575B4", "Tumor" = "#D73027")) +
    labs(title    = trf,
         subtitle = sprintf("%s  |  logFC = %s  |  FDR = %s", direction, lfc_val, fdr_val),
         x = NULL, y = "log2 CPM") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none")
})

pdf(file.path(FIG_DIR, "trf_boxplots_top10.pdf"), width = 12, height = 8)
gridExtra_available <- requireNamespace("gridExtra", quietly = TRUE)
if (gridExtra_available) {
  gridExtra::grid.arrange(grobs = plot_list, ncol = 5)
} else {
  for (p in plot_list) print(p)
}
dev.off()
cat("  Saved: trf_boxplots_top10.pdf\n")

# ── 10. Heatmap for top 50 significant tRFs ───────────────────────────────────
cat("Generating heatmap for top 50 tRFs ...\n")

top50 <- head(res_sig[order(res_sig$FDR), "trf_id"], 50)
mat   <- lcpm[top50, ]

# Z-score rows
mat_z <- t(scale(t(mat)))
mat_z[mat_z >  3] <-  3   # cap extremes
mat_z[mat_z < -3] <- -3

ann_col <- data.frame(Group = group)
rownames(ann_col) <- colnames(mat_z)

ann_colors <- list(Group = c(Normal = "#4575B4", Tumor = "#D73027"))

pdf(file.path(FIG_DIR, "trf_heatmap_top50.pdf"), width = 14, height = 12)
pheatmap(mat_z,
         annotation_col  = ann_col,
         annotation_colors = ann_colors,
         show_colnames   = FALSE,
         fontsize_row    = 7,
         color           = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         main            = "Top 50 Differentially Expressed tRFs (Z-score)",
         cluster_cols    = TRUE,
         cluster_rows    = TRUE,
         border_color    = NA)
dev.off()
cat("  Saved: trf_heatmap_top50.pdf\n")

cat("\n✅ Phase 4 (tRF differential expression) complete.\n")
