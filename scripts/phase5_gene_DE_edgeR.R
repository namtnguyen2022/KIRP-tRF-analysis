# =============================================================================
# Phase 5: Differential Expression of KIRP-Relevant Genes using edgeR
# KIRP tRF Project
#
# Gene list: MET, FH, CDKN2A, SETD2, NF2, BAP1, PBRM1, VHL, MTOR, TFE3,
#            TFEB, NFE2L2, KEAP1, TERT, PTEN, TP53
#
# Input:
#   data_processed/mrna_count_matrix.tsv
#   data_processed/mrna_sample_metadata.tsv
#   One STAR counts TSV (to extract gene_id <-> gene_name mapping)
#
# Output:
#   results/tables/gene_DE_full.tsv
#   results/tables/gene_DE_significant.tsv
#   results/figures/gene_boxplots_significant.pdf
#   results/figures/gene_heatmap.pdf
#   results/figures/gene_volcano.pdf
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
COUNT_FILE <- file.path(ROOT, "data_processed", "mrna_count_matrix.tsv")
META_FILE  <- file.path(ROOT, "data_processed", "mrna_sample_metadata.tsv")
GDC_DIR    <- file.path(ROOT, "data_raw", "gdc_rnaseq",
                         "gdc_download_20260506_091015.911939")
TABLE_DIR  <- file.path(ROOT, "results", "tables")
FIG_DIR    <- file.path(ROOT, "results", "figures")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,   recursive = TRUE, showWarnings = FALSE)

# ── KIRP candidate gene list ──────────────────────────────────────────────────
KIRP_GENES <- c(
  "MET",     # Type 1 papillary RCC driver (proto-oncogene)
  "FH",      # Type 2 papillary RCC driver (fumarate hydratase)
  "CDKN2A",  # Tumor suppressor; cell cycle
  "SETD2",   # Chromatin remodeling; frequently mutated in KIRP
  "NF2",     # Tumor suppressor
  "BAP1",    # Chromatin; mutated in aggressive RCC
  "PBRM1",   # SWI/SNF chromatin remodeling
  "VHL",     # Clear cell RCC suppressor; occasionally in KIRP
  "MTOR",    # mTOR pathway; therapeutic target in RCC
  "TFE3",    # Transcription factor; translocation-associated KIRP subtype
  "TFEB",    # MiT/TFE family; amplified in RCC
  "NFE2L2",  # Oxidative stress response (NRF2)
  "KEAP1",   # NFE2L2 regulator
  "TERT",    # Telomerase; promoter mutations in many cancers
  "PTEN",    # PI3K/AKT pathway suppressor
  "TP53"     # Master tumor suppressor
)

# ── 1. Build gene_id <-> gene_name map from one STAR counts file ──────────────
cat("Building gene ID <-> name map from STAR counts file ...\n")
uuid_dirs  <- list.dirs(GDC_DIR, full.names = TRUE, recursive = FALSE)
ref_tsv    <- NULL
for (d in uuid_dirs) {
  tsvs <- list.files(d, pattern = "\\.tsv$", full.names = TRUE)
  if (length(tsvs) > 0) { ref_tsv <- tsvs[1]; break }
}

gene_map <- read.table(ref_tsv, sep = "\t", header = TRUE, comment.char = "#",
                       stringsAsFactors = FALSE)
# Remove STAR summary rows
gene_map <- gene_map[grepl("^ENSG", gene_map$gene_id), c("gene_id", "gene_name")]
gene_map <- gene_map[!duplicated(gene_map$gene_name), ]

# Find Ensembl IDs for our candidate genes
target_map <- gene_map[gene_map$gene_name %in% KIRP_GENES, ]
found_genes <- target_map$gene_name
missing     <- setdiff(KIRP_GENES, found_genes)

cat(sprintf("  Candidate genes found in matrix: %d / %d\n",
            length(found_genes), length(KIRP_GENES)))
if (length(missing) > 0)
  cat(sprintf("  Not found: %s\n", paste(missing, collapse = ", ")))

# ── 2. Load mRNA count matrix ─────────────────────────────────────────────────
cat("Loading mRNA count matrix ...\n")
counts_all <- read.table(COUNT_FILE, sep = "\t", header = TRUE,
                         row.names = 1, check.names = FALSE)

# Subset to candidate gene Ensembl IDs
target_ids   <- target_map$gene_id
target_ids   <- target_ids[target_ids %in% rownames(counts_all)]
counts_genes <- counts_all[target_ids, ]

# Replace Ensembl IDs with gene names as rownames
rownames(counts_genes) <- target_map$gene_name[
  match(rownames(counts_genes), target_map$gene_id)]

cat(sprintf("  Genes in count matrix: %d\n", nrow(counts_genes)))

# ── 3. Load and align metadata ────────────────────────────────────────────────
meta <- read.table(META_FILE, sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE,
                   colClasses = c(sample_type_code = "character"))

common_samples <- intersect(colnames(counts_genes), meta$sample_id)
counts_genes   <- counts_genes[, common_samples]
meta           <- meta[match(common_samples, meta$sample_id), ]

cat(sprintf("  Samples: %d tumor, %d normal\n",
            sum(meta$sample_type_code == "01"),
            sum(meta$sample_type_code == "11")))

# ── 4. edgeR DE analysis ──────────────────────────────────────────────────────
group  <- factor(ifelse(meta$sample_type_code == "01", "Tumor", "Normal"),
                 levels = c("Normal", "Tumor"))

dge    <- DGEList(counts = counts_genes, group = group)

# Filter: keep genes with CPM > 1 in at least min(group size) samples
min_n  <- min(table(group))
keep   <- rowSums(cpm(dge) > 1) >= min_n
dge    <- dge[keep, , keep.lib.sizes = FALSE]
cat(sprintf("  Genes after CPM filter: %d / %d\n", nrow(dge), nrow(counts_genes)))

dge    <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(~ group)
dge    <- estimateDisp(dge, design, robust = TRUE)
fit    <- glmQLFit(dge, design, robust = TRUE)
qlf    <- glmQLFTest(fit, coef = 2)

res <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "PValue")$table
res$gene_name <- rownames(res)
res <- res[, c("gene_name", "logFC", "logCPM", "F", "PValue", "FDR")]

cat("\nResults:\n")
print(res[, c("gene_name", "logFC", "FDR")])

cat(sprintf("\n  FDR < 0.05: %d genes\n",   sum(res$FDR < 0.05)))
cat(sprintf("  FDR < 0.05 & |logFC|>=1: %d genes\n",
            sum(res$FDR < 0.05 & abs(res$logFC) >= 1)))

# ── 5. Save tables ─────────────────────────────────────────────────────────────
write.table(res,
            file.path(TABLE_DIR, "gene_DE_full.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

res_sig <- res[res$FDR < 0.05 & abs(res$logFC) >= 1, ]
res_sig <- res_sig[order(res_sig$FDR), ]
write.table(res_sig,
            file.path(TABLE_DIR, "gene_DE_significant.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\nSaved full results  -> %s\n", file.path(TABLE_DIR, "gene_DE_full.tsv")))
cat(sprintf("Saved sig. results  -> %s\n",  file.path(TABLE_DIR, "gene_DE_significant.tsv")))

# ── 6. Volcano plot ───────────────────────────────────────────────────────────
cat("Generating volcano plot ...\n")
vol_df        <- res
vol_df$status <- "Not significant"
vol_df$status[vol_df$FDR < 0.05 & vol_df$logFC >= 1]  <- "Up in Tumor"
vol_df$status[vol_df$FDR < 0.05 & vol_df$logFC <= -1] <- "Down in Tumor"
vol_df$status <- factor(vol_df$status,
                        levels = c("Up in Tumor", "Down in Tumor", "Not significant"))

p_vol <- ggplot(vol_df, aes(x = logFC, y = -log10(FDR),
                             color = status, label = gene_name)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(size = 3.2, max.overlaps = 20,
                  segment.color = "grey50", fontface = "italic") +
  scale_color_manual(values = c("Up in Tumor"     = "#D73027",
                                "Down in Tumor"   = "#4575B4",
                                "Not significant" = "grey60")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(title    = "KIRP-Relevant Gene Expression: Tumor vs Normal",
       subtitle = sprintf("FDR < 0.05 & |logFC| >= 1  |  n = %d genes", nrow(res)),
       x = "log2 Fold Change (Tumor / Normal)",
       y = "-log10(FDR)",
       color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

pdf(file.path(FIG_DIR, "gene_volcano.pdf"), width = 8, height = 6)
print(p_vol)
dev.off()
cat("  Saved: gene_volcano.pdf\n")

# ── 7. Boxplots for significant genes ─────────────────────────────────────────
cat("Generating boxplots ...\n")
lcpm <- cpm(dge, log = TRUE, prior.count = 1)

# If no significant genes, plot all; otherwise only significant
plot_genes <- if (nrow(res_sig) > 0) res_sig$gene_name else res$gene_name

plot_list <- lapply(plot_genes, function(gene) {
  expr_vals <- lcpm[gene, ]
  df_plot <- data.frame(
    logCPM = as.numeric(expr_vals),
    Group  = group
  )
  direction <- ifelse(res[gene, "logFC"] > 0, "Up in Tumor", "Down in Tumor")
  fdr_val   <- signif(res[gene, "FDR"], 3)
  lfc_val   <- round(res[gene, "logFC"], 2)

  ggplot(df_plot, aes(x = Group, y = logCPM, fill = Group)) +
    geom_boxplot(outlier.size = 0.8, width = 0.5, alpha = 0.85) +
    geom_jitter(width = 0.15, size = 0.7, alpha = 0.4) +
    scale_fill_manual(values = c("Normal" = "#4575B4", "Tumor" = "#D73027")) +
    labs(title    = gene,
         subtitle = sprintf("%s | logFC = %s | FDR = %s",
                            direction, lfc_val, fdr_val),
         x = NULL, y = "log2 CPM") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none",
          plot.title    = element_text(face = "italic", size = 11),
          plot.subtitle = element_text(size = 7.5))
})

# Arrange in a grid: up to 4 per row
n_genes  <- length(plot_list)
n_cols   <- min(4, n_genes)
n_rows   <- ceiling(n_genes / n_cols)
pdf_w    <- n_cols * 3.2
pdf_h    <- n_rows * 3.2

pdf(file.path(FIG_DIR, "gene_boxplots_significant.pdf"),
    width = pdf_w, height = pdf_h)
if (n_genes > 1) {
  grid.arrange(grobs = plot_list, ncol = n_cols)
} else {
  print(plot_list[[1]])
}
dev.off()
cat(sprintf("  Saved: gene_boxplots_significant.pdf  (%d genes)\n", n_genes))

# ── 8. Heatmap of all candidate genes ─────────────────────────────────────────
cat("Generating heatmap ...\n")
mat_z <- t(scale(t(lcpm)))
mat_z[mat_z >  2.5] <-  2.5
mat_z[mat_z < -2.5] <- -2.5

# Annotate columns by group, ordered Normal first
col_order <- order(group)
mat_z_ord <- mat_z[, col_order]

ann_col <- data.frame(Group = group[col_order])
rownames(ann_col) <- colnames(mat_z_ord)
ann_colors <- list(Group = c(Normal = "#4575B4", Tumor = "#D73027"))

# Annotate rows with direction
direction_col <- ifelse(res[rownames(mat_z_ord), "FDR"] < 0.05 &
                          res[rownames(mat_z_ord), "logFC"] > 0, "Up",
                 ifelse(res[rownames(mat_z_ord), "FDR"] < 0.05 &
                          res[rownames(mat_z_ord), "logFC"] < 0, "Down", "NS"))
ann_row <- data.frame(Direction = factor(direction_col,
                                          levels = c("Up","Down","NS")))
rownames(ann_row) <- rownames(mat_z_ord)
ann_colors$Direction <- c(Up = "#D73027", Down = "#4575B4", NS = "grey80")

pdf(file.path(FIG_DIR, "gene_heatmap.pdf"), width = 14, height = 6)
pheatmap(mat_z_ord,
         annotation_col   = ann_col,
         annotation_row   = ann_row,
         annotation_colors = ann_colors,
         show_colnames    = FALSE,
         fontsize_row     = 9,
         color            = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         main             = "KIRP-Relevant Gene Expression (Z-score, Tumor vs Normal)",
         cluster_cols     = FALSE,
         cluster_rows     = TRUE,
         border_color     = NA,
         gaps_col         = sum(group[col_order] == "Normal"))
dev.off()
cat("  Saved: gene_heatmap.pdf\n")

cat("\n✅ Phase 5 (gene differential expression) complete.\n")
