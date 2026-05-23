## build_final_figures.R  (v4)
## 5 publication-ready figures (Fig 1 cohort removed per user request)
## Fig 1 = tRF DE  |  Fig 2 = Gene DE  |  Fig 3 = tRF-MET axis
## Fig 4 = ROC/AUC  |  Fig 5 = Stage/Survival

suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(pheatmap)
  library(patchwork); library(cowplot); library(scales)
  library(dplyr); library(tidyr); library(grid)
  library(RColorBrewer); library(pROC)
  library(survival); library(survminer)
})

ROOT <- "/Users/nam.tnguyen2022/KIRP_tRF_project"
DATA <- file.path(ROOT, "data_processed")
TABS <- file.path(ROOT, "results/tables")
RAW  <- file.path(ROOT, "data_raw")
OUT  <- file.path(ROOT, "results/figures")

## ── Publication theme (Nature-style) ────────────────────────────────────────
pub_theme <- theme_classic(base_size = 9, base_family = "Helvetica") +
  theme(
    ## axes
    axis.line         = element_line(color = "black", linewidth = 0.4),
    axis.ticks        = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(0.15, "cm"),
    axis.title        = element_text(size = 8, color = "black"),
    axis.text         = element_text(size = 7.5, color = "black"),
    ## panels — completely white, no border
    panel.background  = element_rect(fill = "white", color = NA),
    panel.border      = element_blank(),
    panel.grid        = element_blank(),
    ## plot area
    plot.background   = element_rect(fill = "white", color = NA),
    plot.title        = element_text(size = 9, face = "bold", color = "black",
                                     hjust = 0, margin = margin(b = 4)),
    plot.margin       = margin(4, 6, 4, 6),
    ## legend — no box by default, small text
    legend.background = element_blank(),
    legend.key        = element_blank(),
    legend.text       = element_text(size = 7.5),
    legend.title      = element_text(size = 8, face = "bold"),
    legend.key.size   = unit(0.32, "cm"),
    ## strips for facets
    strip.background  = element_blank(),
    strip.text        = element_text(size = 7.5, face = "bold", color = "black"),
    ## panel tags (a, b, c)
    plot.tag          = element_text(size = 10, face = "bold", color = "black")
  )
theme_set(pub_theme)

GROUP_COLS <- c(Tumor = "#D62728", Normal = "#1F77B4")
STAGE_COLS <- c("Stage I"   = "#2166AC", "Stage II"  = "#74ADD1",
                "Stage III" = "#F46D43", "Stage IV"  = "#A50026")

wrap_caption <- function(x, width = 190) {
  paste(strwrap(x, width = width), collapse = "\n")
}

## ── Load all data ─────────────────────────────────────────────────────────────
cat("Loading data...\n")

trf_mat   <- read.table(file.path(DATA, "trf_count_matrix.tsv"),
                        sep="\t", header=TRUE, row.names=1, check.names=FALSE)
trf_meta  <- read.table(file.path(DATA, "trf_sample_metadata.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE,
                        colClasses=c(sample_type_code="character"))
trf_de    <- read.table(file.path(TABS, "trf_DE_significant.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE)
trf_full  <- read.table(file.path(TABS, "trf_DE_full.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE)
gene_full <- read.table(file.path(TABS, "gene_DE_full.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE)
anticorr  <- read.table(file.path(TABS, "anticorrelation_prioritized.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE)
roc_tab   <- read.table(file.path(TABS, "roc_auc_top_trfs.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE)
cox_tab   <- read.table(file.path(TABS, "cox_univariate_results.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE)
kw_tab    <- read.table(file.path(TABS, "stage_kruskal_results.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE)
rna22     <- read.table(file.path(TABS, "rna22_binding_summary.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE)
mrna_mat  <- read.table(file.path(DATA, "mrna_count_matrix.tsv"),
                        sep="\t", header=TRUE, row.names=1, check.names=FALSE)
mrna_meta <- read.table(file.path(DATA, "mrna_sample_metadata.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE,
                        colClasses=c(sample_type_code="character"))

trf_lcpm  <- log2(sweep(trf_mat,  2, colSums(trf_mat),  "/") * 1e6 + 1)
mrna_lcpm <- log2(sweep(mrna_mat, 2, colSums(mrna_mat), "/") * 1e6 + 1)
trf_meta$group  <- ifelse(trf_meta$sample_type_code  == "01", "Tumor", "Normal")
mrna_meta$group <- ifelse(mrna_meta$sample_type_code == "01", "Tumor", "Normal")
rownames(trf_meta)  <- trf_meta$sample_id
rownames(mrna_meta) <- mrna_meta$sample_id

star_dir <- list.files(
  file.path(RAW, "gdc_rnaseq/gdc_download_20260506_091015.911939"),
  pattern = "augmented_star_gene_counts.tsv$", recursive = TRUE, full.names = TRUE)
gmap <- read.table(star_dir[1], sep="\t", header=TRUE,
                   comment.char="#")[, c("gene_id","gene_name")]

clin_raw <- read.table(file.path(RAW, "gdc_clinical/clinical.tsv"),
                       sep="\t", header=TRUE, stringsAsFactors=FALSE,
                       fill=TRUE, quote="")
clin <- clin_raw[, c("cases.submitter_id","diagnoses.ajcc_pathologic_stage",
                     "demographic.vital_status","demographic.days_to_death",
                     "diagnoses.days_to_last_follow_up")]
colnames(clin) <- c("patient_id","stage","vital_status","days_death","days_fu")
clin <- clin[!duplicated(clin$patient_id), ]
primary_stage <- clin_raw[
  as.character(clin_raw$diagnoses.diagnosis_is_primary_disease) %in% c("TRUE", "True", "true"),
  c("cases.submitter_id", "diagnoses.submitter_id",
    "diagnoses.ajcc_pathologic_stage")
]
colnames(primary_stage) <- c("patient_id", "diagnosis_id", "stage")
primary_stage <- primary_stage[order(primary_stage$patient_id,
                                     primary_stage$diagnosis_id), ]
primary_stage <- primary_stage[!duplicated(primary_stage$patient_id), ]
primary_stage$stage[primary_stage$stage == "'--"] <- NA
primary_stage$stage_simple <- sub(" [ABC]$", "", primary_stage$stage)
clin$stage <- NULL
clin <- merge(clin, primary_stage[, c("patient_id", "stage", "stage_simple")],
              by = "patient_id", all.x = TRUE, sort = FALSE)
cn <- function(x) { x[x %in% c("'--","--","")] <- NA; as.numeric(x) }
clin$days_death <- cn(clin$days_death)
clin$days_fu    <- cn(clin$days_fu)
clin$event      <- as.integer(clin$vital_status == "Dead")
clin$os_days    <- ifelse(!is.na(clin$days_death), clin$days_death, clin$days_fu)

tumor_trf  <- trf_meta[trf_meta$group == "Tumor", ]
tumor_mrna <- mrna_meta[mrna_meta$group == "Tumor", ]
tumor_trf$patient_id  <- substr(tumor_trf$sample_id,  1, 12)
tumor_mrna$patient_id <- substr(tumor_mrna$sample_id, 1, 12)
cat("Data loaded.\n\n")

## ════════════════════════════════════════════════════════════════════════════
## FIGURE 1 — tRF Differential Expression
## ════════════════════════════════════════════════════════════════════════════
cat("Figure 1: tRF DE...\n")

trf_full$sig <- case_when(
  trf_full$FDR < 0.05 & trf_full$logFC >  1 ~ "Up",
  trf_full$FDR < 0.05 & trf_full$logFC < -1 ~ "Down",
  TRUE ~ "NS")
trf_full$sig <- factor(trf_full$sig, levels = c("Up","Down","NS"))
n_up  <- sum(trf_full$sig == "Up")
n_dn  <- sum(trf_full$sig == "Down")

## Label top 5 up + top 5 down by FDR
top_label <- bind_rows(
  trf_full %>% filter(sig == "Up")   %>% arrange(FDR) %>% slice_head(n=5),
  trf_full %>% filter(sig == "Down") %>% arrange(FDR) %>% slice_head(n=5))
## Shorten label: strip "tRF-NN-" prefix, keep the ID part only
top_label$short_id <- sub("^tRF-[0-9]+-", "", top_label$trf_id)

## Actual data x/y range — don't pad beyond the data
x_lo  <- floor(min(trf_full$logFC))   - 0.3
x_hi  <- ceiling(max(trf_full$logFC)) + 0.3
y_max <- max(-log10(trf_full$FDR[trf_full$FDR > 0]), na.rm = TRUE)

## ── Panel A: Volcano ─────────────────────────────────────────────────────────
p1A <- ggplot(trf_full, aes(x = logFC, y = -log10(FDR), color = sig)) +
  ## NS cloud first (behind)
  geom_point(data = ~filter(.x, sig == "NS"),
             size = 0.22, alpha = 0.14, color = "grey75") +
  ## Significant points
  geom_point(data = ~filter(.x, sig != "NS"),
             size = 0.75, alpha = 0.75) +
  ## Labels for top hits
  geom_text_repel(
    data          = top_label,
    aes(label     = short_id),
    size          = 2.3, segment.size = 0.28, segment.color = "grey50",
    segment.alpha = 0.8, min.segment.length = 0.15, max.overlaps = 20,
    show.legend   = FALSE, fontface = "italic",
    box.padding   = 0.45, point.padding = 0.2, force = 3,
    direction     = "y") +
  ## Threshold lines
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             color = "grey45", linewidth = 0.38) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey45", linewidth = 0.38) +
  scale_color_manual(
    values = c(Up = "#C0392B", Down = "#2980B9", NS = "grey75"),
    labels = c(paste0("Up (n=", n_up, ")"),
               paste0("Down (n=", n_dn, ")"), "NS"),
    name   = NULL) +
  scale_x_continuous(
    limits = c(x_lo, x_hi),
    breaks = seq(ceiling(x_lo), floor(x_hi), by = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.07))) +
  labs(x = expression(log[2]~"fold change (Tumor / Normal)"),
       y = expression(-log[10](FDR))) +
  guides(color = guide_legend(
    override.aes = list(size = 2.5, alpha = 1))) +
  theme(legend.position        = "inside",
        legend.position.inside = c(0.5, 0.93),
        legend.direction       = "horizontal",
        legend.background      = element_blank())

## ── Panel B: Heatmap (top 30 by |logFC|) ─────────────────────────────────────
top30  <- trf_de %>% arrange(desc(abs(logFC))) %>%
  slice_head(n = 30) %>% pull(trf_id)
top30  <- top30[top30 %in% rownames(trf_lcpm)]
## Order columns: all Normal first, then all Tumor
col_ord <- c(which(trf_meta$group == "Normal"),
             which(trf_meta$group == "Tumor"))
mat30   <- trf_lcpm[top30, trf_meta$sample_id[col_ord]]
## Shorten row names for readability
rownames(mat30) <- sub("^tRF-[0-9]+-", "", rownames(mat30))
ann_df  <- data.frame(Group = trf_meta$group[col_ord],
                      row.names = colnames(mat30))
hm_png  <- tempfile(fileext = ".png")
pheatmap(mat30,
  annotation_col    = ann_df,
  annotation_colors = list(Group = GROUP_COLS),
  show_colnames     = FALSE,
  scale             = "row",
  cluster_cols      = FALSE,
  cluster_rows      = TRUE,
  color             = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
  fontsize_row      = 7.5,
  fontsize          = 8.5,
  border_color      = NA,
  annotation_legend = TRUE,
  annotation_names_col = FALSE,
  main              = "Top 30 tRFs by |log2FC|\n(row z-score)",
  filename          = hm_png,
  width = 4.4, height = 6.0, res = 280)
p1B <- ggdraw() + draw_image(hm_png, scale = 1)

## ── Panel C: Boxplots — top 5 Down + top 4 Up by FDR (balanced) ─────────────
box_ids <- c(
  trf_full %>% filter(sig == "Down") %>% arrange(FDR) %>%
    slice_head(n=5) %>% pull(trf_id),
  trf_full %>% filter(sig == "Up")   %>% arrange(FDR) %>%
    slice_head(n=4) %>% pull(trf_id))
box_ids <- box_ids[box_ids %in% rownames(trf_lcpm)]

## Build long data frame; add direction for strip color
box_df <- bind_rows(lapply(box_ids, function(t) {
  dir <- as.character(trf_full$sig[trf_full$trf_id == t])
  data.frame(
    tRF   = sub("^tRF-[0-9]+-", "", t),   # shortened label
    full_id = t,
    dir   = dir,
    expr  = as.numeric(trf_lcpm[t, trf_meta$sample_id]),
    Group = factor(trf_meta$group, levels = c("Normal","Tumor")))
}))
box_df$tRF <- factor(box_df$tRF,
  levels = sub("^tRF-[0-9]+-", "", box_ids))

## FDR + FC annotation per facet
box_annot <- bind_rows(lapply(box_ids, function(t) {
  row   <- trf_full[trf_full$trf_id == t, ]
  fdr_s <- sub("e-0","e-", formatC(row$FDR, format="e", digits=1))
  data.frame(
    tRF   = sub("^tRF-[0-9]+-", "", t),
    dir   = as.character(row$sig),
    lbl   = paste0("FC=", sprintf("%+.2f", row$logFC),
                   "\nFDR=", fdr_s))
}))
box_annot$tRF <- factor(box_annot$tRF, levels = levels(box_df$tRF))
## strip direction colors: Down=blue, Up=red (applied via labeller)
dir_map <- setNames(
  trf_full$sig[match(box_ids, trf_full$trf_id)],
  sub("^tRF-[0-9]+-", "", box_ids))
strip_col_vec <- ifelse(dir_map[levels(box_df$tRF)] == "Down",
                        "#2980B9", "#C0392B")

p1C <- ggplot(box_df, aes(x = Group, y = expr, fill = Group, color = Group)) +
  ## filled box with whiskers, no outlier points (shown by jitter)
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.4,
               alpha = 0.35, color = "black") +
  ## individual points overlaid — height=0 prevents vertical jitter below zero
  geom_jitter(width = 0.15, height = 0, size = 0.55, alpha = 0.55, shape = 16) +
  ## Annotation: top-right inside each facet
  geom_text(data = box_annot,
            aes(x = Inf, y = Inf, label = lbl),
            inherit.aes = FALSE,
            hjust = 1.05, vjust = 1.3,
            size = 2.0, color = "grey30", fontface = "italic") +
  facet_wrap(~tRF, nrow = 3, scales = "free_y") +
  scale_fill_manual(values = GROUP_COLS, guide = "none") +
  scale_color_manual(values = GROUP_COLS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.32))) +
  labs(x = NULL, y = expression(log[2]~"(CPM + 1)")) +
  theme(axis.text.x  = element_text(size = 7.5),
        axis.text.y  = element_text(size = 7),
        strip.text   = element_text(size = 6.5, face = "bold"))

## ── Assemble ─────────────────────────────────────────────────────────────────
fig1_top <- cowplot::plot_grid(
  p1A, p1B,
  ncol = 2, rel_widths = c(1, 1),
  labels = c("A", "B"),
  label_size = 10, label_fontface = "bold",
  align = "hv", axis = "tblr")
fig1_body <- cowplot::plot_grid(
  fig1_top, p1C,
  nrow = 2, rel_heights = c(1, 1),
  labels = c("", "C"),
  label_size = 10, label_fontface = "bold")
fig1_caption <- wrap_caption(
  paste0(
    "a: Volcano plot (dashed lines: |log2FC|=1, FDR=0.05).  ",
    "b: Heatmap of top 30 tRFs by |log2FC| (row z-score).  ",
    "c: Expression of top 5 down- and top 4 up-regulated tRFs."),
  width = 188)
fig1_header <- cowplot::ggdraw() +
  cowplot::draw_label(
    "Figure 1. Differential expression of tRFs in KIRP tumor vs. normal tissue.",
    fontface = "bold", size = 9,
    x = 0.003, y = 0.80, hjust = 0, vjust = 0.5) +
  cowplot::draw_label(
    fig1_caption,
    size = 7.2, color = "grey40",
    x = 0.003, y = 0.38, hjust = 0, vjust = 0.5,
    lineheight = 1.05)
fig1 <- cowplot::ggdraw() +
  cowplot::draw_grob(rectGrob(gp = gpar(fill = "white", col = NA))) +
  cowplot::draw_plot(
    cowplot::plot_grid(fig1_header, fig1_body,
                       nrow = 2, rel_heights = c(0.072, 1)),
    x = 0, y = 0, width = 1, height = 1)

ggsave(file.path(OUT, "fig2_trf_DE.pdf"),
       fig1, width = 13, height = 12,
       device = "pdf", useDingbats = FALSE)
ggsave(file.path(OUT, "fig2_trf_DE.png"),
       fig1, width = 13, height = 12,
       device = "png", dpi = 600, bg = "white")
cat("  Done.\n")

## ════════════════════════════════════════════════════════════════════════════
## FIGURE 2 — Candidate Gene Expression
## ════════════════════════════════════════════════════════════════════════════
cat("Figure 2: Gene DE...\n")

## Significance tiers
gene_full$sig_tier <- case_when(
  gene_full$FDR < 0.05 & gene_full$logFC >  1  ~ "Up (|FC|>2)",
  gene_full$FDR < 0.05 & gene_full$logFC < -1  ~ "Down (|FC|>2)",
  gene_full$FDR < 0.05                          ~ "Sig (small FC)",
  TRUE                                           ~ "NS")
gene_full$sig_tier <- factor(gene_full$sig_tier,
  levels = c("Up (|FC|>2)", "Down (|FC|>2)", "Sig (small FC)", "NS"))
TIER_COLS <- c("Up (|FC|>2)"    = "#C0392B",
               "Down (|FC|>2)"  = "#2980B9",
               "Sig (small FC)" = "#8E44AD",
               "NS"             = "grey65")

## ── Panel A: ranked lollipop ─────────────────────────────────────────────────
lollipop_df <- gene_full %>%
  arrange(logFC) %>%
  mutate(
    gene_f   = factor(gene_name, levels = gene_name),
    log10fdr = pmin(-log10(FDR), 30),
    sig_lbl  = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE        ~ ""),
    lbl_x = ifelse(logFC >= 0,
                   pmax(logFC + 0.22, 0.5),
                   pmin(logFC - 0.22, -0.5)))

x_lo_g <- floor(min(lollipop_df$logFC))   - 0.5
x_hi_g <- ceiling(max(lollipop_df$logFC)) + 1.2

p2A <- ggplot(lollipop_df,
              aes(y = gene_f, x = logFC, color = sig_tier)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.35) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             color = "grey80", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = logFC, yend = gene_f),
               linewidth = 0.7, alpha = 0.6) +
  geom_point(aes(size = log10fdr), alpha = 0.9) +
  geom_text(aes(x = lbl_x, label = sig_lbl),
            hjust = ifelse(lollipop_df$logFC >= 0, 0, 1),
            vjust = 0.4, size = 3.2, color = "grey20", fontface = "bold") +
  scale_color_manual(values = TIER_COLS, name = NULL) +
  scale_size_continuous(name  = expression(-log[10](FDR)),
                        range = c(2, 8), breaks = c(5, 15, 25)) +
  scale_x_continuous(limits = c(x_lo_g, x_hi_g), expand = c(0, 0),
                     breaks = seq(floor(x_lo_g)+1, floor(x_hi_g), by=1)) +
  labs(x = expression(log[2]~"fold change (Tumor / Normal)"), y = NULL) +
  guides(color = guide_legend(order=1, override.aes=list(size=3)),
         size  = guide_legend(order=2)) +
	  theme(axis.text.y         = element_text(size=8.5, face="italic",
	                                            color="black"),
	        axis.line.y         = element_blank(),
	        axis.ticks.y        = element_blank(),
	        panel.grid.major.x  = element_line(color="grey93", linewidth=0.3),
	        legend.position     = "bottom",
	        legend.direction    = "horizontal",
	        legend.box          = "horizontal",
	        legend.justification = "center",
	        legend.title        = element_text(size = 6.6, face = "bold"),
	        legend.text         = element_text(size = 6.5),
	        legend.key.size     = unit(0.25, "cm"),
	        legend.margin       = margin(t = -2, b = 0),
	        plot.margin         = margin(4, 6, 0, 6))

## ── Panel B: KIRP-relevant gene heatmap ─────────────────────────────────────
gene_order <- lollipop_df$gene_name
gene_ids <- gmap$gene_id[match(gene_order, gmap$gene_name)]
names(gene_ids) <- gene_order
gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids %in% rownames(mrna_lcpm)]
col_order_gene <- c(which(mrna_meta$group == "Normal"),
                    which(mrna_meta$group == "Tumor"))
gene_hm_mat <- mrna_lcpm[gene_ids, mrna_meta$sample_id[col_order_gene]]
rownames(gene_hm_mat) <- names(gene_ids)
gene_hm_z <- t(scale(t(gene_hm_mat)))
gene_hm_z[!is.finite(gene_hm_z)] <- 0
gene_hm_df <- as.data.frame(as.table(gene_hm_z), stringsAsFactors = FALSE)
colnames(gene_hm_df) <- c("gene", "sample_id", "z")
gene_hm_df$gene <- factor(gene_hm_df$gene, levels = rev(rownames(gene_hm_z)))
gene_hm_df$sample_index <- match(gene_hm_df$sample_id, colnames(gene_hm_z))

gene_hm_ann <- data.frame(
  Group = ifelse(mrna_meta$group[col_order_gene] == "Normal",
                 "Solid normal", "Tumor"),
  row.names = colnames(gene_hm_mat))
gene_hm_ann$Group <- factor(gene_hm_ann$Group,
                            levels = c("Solid normal", "Tumor"))
gene_hm_cols <- list(Group = c("Solid normal" = unname(GROUP_COLS["Normal"]),
                               "Tumor" = unname(GROUP_COLS["Tumor"])))
gene_hm_breaks <- seq(-5, 5, length.out = 101)
gene_hm <- pheatmap(
  gene_hm_mat,
  annotation_col = gene_hm_ann,
  annotation_colors = gene_hm_cols,
  show_colnames = FALSE,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
  scale = "row",
  breaks = gene_hm_breaks,
  fontsize_row = 7.4,
  fontsize = 7.2,
  border_color = NA,
  annotation_legend = TRUE,
  annotation_names_col = FALSE,
  legend = TRUE,
  legend_breaks = c(-5, 0, 5),
  legend_labels = c("-5", "0", "5"),
  main = "KIRP candidate gene expression heatmap (row z-score)",
  margins = c(5, 5, 3, 2),
  silent = TRUE)
p2B <- cowplot::ggdraw() +
  cowplot::draw_grob(gene_hm$gtable) +
  theme(plot.margin = margin(0, 2, 2, 2))

## ── Panels C-E: strip plots for MET, CDKN2A, TERT ───────────────────────────
key_genes <- c("MET", "CDKN2A", "TERT")

make_gene_strip <- function(gname, show_ylab = FALSE) {
  gid   <- gmap$gene_id[gmap$gene_name == gname][1]
  if (is.na(gid) || !gid %in% rownames(mrna_lcpm)) return(NULL)
  row_g   <- gene_full[gene_full$gene_name == gname, ]
  tier    <- as.character(row_g$sig_tier[1])
  dot_col <- TIER_COLS[tier]

  df <- data.frame(
    expr  = as.numeric(mrna_lcpm[gid, mrna_meta$sample_id]),
    Group = factor(mrna_meta$group, levels = c("Normal","Tumor")))
  n_N <- sum(df$Group == "Normal")
  n_T <- sum(df$Group == "Tumor")

  ## Wilcoxon p — just show stars above the strip
  pv    <- wilcox.test(expr ~ Group, data = df)$p.value
  ## Format p-value as "P = 0.001" or "P < 0.0001" like the Nature reference
  pv_lbl <- if (pv < 0.0001) "P < 0.0001" else paste0("P = ", signif(pv, 2))

  ## Bracket: horizontal bar above the highest point in either group,
  ## with short downward feet — matching the Nature bracket style
  y_N    <- max(df$expr[df$Group == "Normal"])
  y_T    <- max(df$expr[df$Group == "Tumor"])
  y_rng  <- diff(range(df$expr))
  y_bar  <- max(y_N, y_T) + 0.10 * y_rng   # horizontal crossbar
  y_foot <- y_bar - 0.04 * y_rng            # foot tick drops slightly below bar
  y_txt  <- y_bar + 0.06 * y_rng            # P-value label above bar

  ggplot(df, aes(x = Group, y = expr, fill = Group, color = Group)) +
    ## filled box, no outlier points
    geom_boxplot(outlier.shape = NA, width = 0.52, linewidth = 0.4,
                 alpha = 0.30, color = "black") +
    ## individual points overlaid — height=0 prevents vertical jitter below zero
    geom_jitter(width = 0.13, height = 0, size = 0.9, alpha = 0.55, shape = 16) +
    ## Nature-style bracket: horizontal bar + short downward feet
    annotate("segment", x=1,   xend=2,   y=y_bar,  yend=y_bar,
             linewidth=0.45, color="black") +
    annotate("segment", x=1,   xend=1,   y=y_foot, yend=y_bar,
             linewidth=0.45, color="black") +
    annotate("segment", x=2,   xend=2,   y=y_foot, yend=y_bar,
             linewidth=0.45, color="black") +
    annotate("text", x=1.5, y=y_txt, label=pv_lbl,
             size=2.8, hjust=0.5, vjust=0, color="black", fontface="italic") +
    scale_x_discrete(labels = c(
      "Normal" = paste0("Normal\n(n=", n_N, ")"),
      "Tumor"  = paste0("Tumor\n(n=", n_T, ")"))) +
    scale_fill_manual(values  = GROUP_COLS, guide = "none") +
    scale_color_manual(values = GROUP_COLS, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.07, 0.30))) +
    labs(x = NULL,
         y = if (show_ylab) expression(log[2]~"(CPM + 1)") else NULL,
         title = gname) +
    theme(plot.title  = element_text(hjust=0.5, face="bold.italic",
                                     size=9, color=dot_col),
          axis.text.x = element_text(size=7.5),
          axis.text.y = element_text(size=7.5),
          plot.margin = margin(4, 8, 4, 8))
}

strip_plots <- Filter(Negate(is.null),
  lapply(seq_along(key_genes), function(i)
    make_gene_strip(key_genes[i], show_ylab = (i == 1))))

## ── Assemble revised Figure 2: A/B on top, C-E on bottom ───────────────────
p2_top <- cowplot::plot_grid(
  p2A, p2B,
  ncol = 2, rel_widths = c(1.05, 1.15),
  labels = c("A", "B"),
  label_size = 10, label_fontface = "bold",
  align = "hv", axis = "tblr")
p2_bottom <- cowplot::plot_grid(
  plotlist = strip_plots,
  ncol = 3,
  labels = c("C", "D", "E"),
  label_size = 10, label_fontface = "bold",
  align = "hv", axis = "tblr")
fig2_caption <- wrap_caption(
  paste0(
    "(A) Sixteen KIRP-relevant genes ranked by log2 fold change. Dot size indicates -log10(FDR), ",
    "and dashed lines indicate |log2FC| = 1. (B) Heatmap of the 16 KIRP-relevant genes across tumor ",
    "and solid normal samples, scaled by row (z-score). (C-E) Individual sample expression of MET, CDKN2A, ",
    "and TERT in tumor versus solid normal tissue."),
  width = 188)
fig2_body <- cowplot::plot_grid(
  p2_top, p2_bottom,
  nrow = 2, rel_heights = c(1.18, 0.82))
fig2_header <- cowplot::ggdraw() +
  cowplot::draw_label(
    "Figure 2. Differential expression of KIRP candidate target genes.",
    fontface = "bold", size = 9,
    x = 0.003, y = 0.80, hjust = 0, vjust = 0.5) +
  cowplot::draw_label(
    fig2_caption,
    size = 7.2, color = "grey40",
    x = 0.003, y = 0.38, hjust = 0, vjust = 0.5,
    lineheight = 1.05)
fig2_full <- cowplot::ggdraw() +
  cowplot::draw_grob(rectGrob(gp = gpar(fill = "white", col = NA))) +
  cowplot::draw_plot(
    cowplot::plot_grid(fig2_header, fig2_body,
                       nrow = 2, rel_heights = c(0.072, 1)),
    x = 0, y = 0, width = 1, height = 1)

ggsave(file.path(OUT, "fig3_gene_DE.pdf"),
       fig2_full, width = 13, height = 11,
       device = "pdf", useDingbats = FALSE)
ggsave(file.path(OUT, "fig3_gene_DE.png"),
       fig2_full, width = 13, height = 11,
       device = "png", dpi = 600, bg = "white")
ggsave(file.path(OUT, "fig3_gene_DE_revised.pdf"),
       fig2_full, width = 13, height = 11,
       device = "pdf", useDingbats = FALSE)
ggsave(file.path(OUT, "fig3_gene_DE_revised.png"),
       fig2_full, width = 13, height = 11,
       device = "png", dpi = 600, bg = "white")
cat("  Done.\n")

## ════════════════════════════════════════════════════════════════════════════
## FIGURE 3 — tRF-MET Regulatory Axis
## ════════════════════════════════════════════════════════════════════════════
cat("Figure 3: tRF-MET axis...\n")

met_id  <- gmap$gene_id[gmap$gene_name == "MET"][1]
shared  <- intersect(tumor_trf$patient_id, tumor_mrna$patient_id)
tp_trf  <- tumor_trf[match(shared, tumor_trf$patient_id), ]
tp_mrna <- tumor_mrna[match(shared, tumor_mrna$patient_id), ]
cat("  Shared tumor patients:", length(shared), "\n")

top5 <- anticorr %>% filter(gene == "MET") %>%
  arrange(FDR) %>% slice_head(n = 5)

## Scatter: one panel per tRF
make_scatter <- function(i) {
  trf <- top5$trf[i]
  if (!trf %in% rownames(trf_lcpm) || !met_id %in% rownames(mrna_lcpm)) return(NULL)
  df <- data.frame(
    x = as.numeric(trf_lcpm[trf,     tp_trf$sample_id]),
    y = as.numeric(mrna_lcpm[met_id, tp_mrna$sample_id]))
  r_val <- round(top5$spearman_r[i], 3)
  fdr_s <- formatC(top5$FDR[i], format = "e", digits = 1)
  fdr_s <- sub("e-0", "e-", fdr_s)
  nm    <- sub("tRF-", "", trf)
  ann   <- paste0("r = ", r_val, ", FDR = ", fdr_s)
  ggplot(df, aes(x = x, y = y)) +
    geom_point(size = 0.75, alpha = 0.65, shape = 16, color = "black") +
    geom_smooth(method = "lm", se = TRUE, color = "#C0392B",
                fill = "#CCCCCC", alpha = 0.20, linewidth = 0.7) +
    annotate("text",
             x = -Inf, y = Inf,
             hjust = -0.05, vjust = 1.4,
             label = ann,
             size = 2.3, color = "grey25", fontface = "italic") +
    labs(x     = expression(log[2]~"(CPM+1)  tRF"),
         y     = expression(log[2]~"(CPM+1)  MET"),
         title = nm) +
    theme(plot.title  = element_text(size = 8, face = "bold", hjust = 0.5),
          axis.title  = element_text(size = 7),
          axis.text   = element_text(size = 6.5),
          plot.margin = margin(5, 7, 5, 7))
}
scatters <- Filter(Negate(is.null), lapply(seq_len(nrow(top5)), make_scatter))

## RNA22 lollipop — MET predictions only
## Dots at the tip (left end, negative deltaG); text to the RIGHT of 0
rna22_met <- rna22 %>%
  filter(gene == "MET") %>%
  arrange(p_value) %>%
  mutate(
    tRF_short = sub("tRF-", "", tRF_id),
    tRF_short = factor(tRF_short, levels = rev(tRF_short)),
    log10p    = -log10(p_value),
    sig       = p_value < 0.05,
    pos_lbl   = paste0("3'UTR pos. ", format(target_position, big.mark=",")))

# x range: pad 15% on left for labels, small gap on right of 0
x_min <- min(rna22_met$folding_energy_kcal) * 1.18
p3F <- ggplot(rna22_met,
              aes(y = tRF_short, x = folding_energy_kcal)) +
  ## stem (lollipop stick) from 0 to the deltaG value
  geom_segment(aes(x = 0, xend = folding_energy_kcal,
                   yend = tRF_short, color = sig),
               linewidth = 0.85, alpha = 0.65) +
  ## dot at deltaG tip
  geom_point(aes(color = sig, size = log10p), alpha = 0.92) +
  ## 3'UTR position label — placed at x=0.8 (just right of the 0 baseline)
  geom_text(aes(x = 1.2, label = pos_lbl),
            hjust = 0, size = 2.5, color = "grey30") +
  geom_vline(xintercept = 0, linewidth = 0.45, color = "grey30") +
  scale_color_manual(
    values = c("TRUE" = "#C0392B", "FALSE" = "#95A5A6"),
    labels = c("TRUE" = "p < 0.05", "FALSE" = "p >= 0.05"),
    name   = NULL) +
  scale_size_continuous(
    name   = "-log10(p)",
    range  = c(3, 9),
    breaks = c(1, 2, 3)) +
  scale_x_continuous(
    limits = c(x_min, 14),
    expand = c(0, 0),
    breaks = seq(-30, 0, by = 5)) +
  labs(x = "Predicted folding energy, dG (kcal/mol)",
       y = NULL,
       title = "RNA22 v2: tRF-MET 3'UTR binding predictions (red = p<0.05)") +
  theme(axis.text.y         = element_text(size = 8, face = "bold"),
        axis.title.x        = element_text(size = 8),
        legend.position     = "none",
        panel.grid.major.y  = element_line(color = "grey92", linewidth = 0.35))

## Layout: row1 = top 3 scatters, row2 = remaining 2 scatters + RNA22
n_sc <- length(scatters)
row1_panels <- scatters[seq_len(min(3, n_sc))]
row2_panels <- if (n_sc > 3) scatters[4:n_sc] else list()
row2_panels <- c(row2_panels, list(p3F))

row1 <- wrap_plots(row1_panels, ncol = 3)
row2 <- wrap_plots(row2_panels, ncol = 3)

fig3 <- (row1 / row2) +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(
    tag_levels = "A",
    title      = "Figure 3. The tRF-MET regulatory axis in KIRP.",
    subtitle   = paste0(
      "a-e: Spearman anticorrelation (BH-FDR) across ", length(shared),
      " matched KIRP tumor samples (top 5 anti-correlated tRFs).  ",
      "f: RNA22 v2 predicted tRF-MET 3'UTR binding; dot size = -log10(p), red = p<0.05."),
    theme      = theme(
      plot.tag      = element_text(size = 10, face = "bold"),
      plot.title    = element_text(size = 9, face = "bold", hjust = 0,
                                   margin = margin(b = 3)),
      plot.subtitle = element_text(size = 7.5, color = "grey40", hjust = 0,
                                   margin = margin(b = 6)))
  ) & theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(OUT, "fig4_trf_met_axis.pdf"),
       fig3, width = 14, height = 9,
       device = "pdf", useDingbats = FALSE)
ggsave(file.path(OUT, "fig4_trf_met_axis.png"),
       fig3, width = 14, height = 9,
       device = "png", dpi = 600, bg = "white")
cat("  Done.\n")

## ════════════════════════════════════════════════════════════════════════════
## FIGURE 4 — ROC / AUC
## ════════════════════════════════════════════════════════════════════════════
cat("Figure 4: ROC...\n")

response  <- as.numeric(trf_meta$group == "Tumor")
top5_roc  <- roc_tab$tRF[1:min(5, nrow(roc_tab))]
pal5      <- c("#E74C3C","#2980B9","#27AE60","#8E44AD","#E67E22")

roc_df <- bind_rows(Filter(Negate(is.null), lapply(seq_along(top5_roc), function(i) {
  t <- top5_roc[i]
  if (!t %in% rownames(trf_lcpm)) return(NULL)
  r   <- roc(response, as.numeric(trf_lcpm[t, ]), quiet = TRUE)
  a   <- round(as.numeric(auc(r)), 3)
  lbl <- paste0(sub("tRF-","",t),"  (AUC = ",a,")")
  data.frame(FPR=1-r$specificities, TPR=r$sensitivities,
             label=lbl, col=pal5[i], stringsAsFactors=FALSE)
})))
roc_df$label <- factor(roc_df$label, levels = unique(roc_df$label))
pal_named    <- setNames(unique(roc_df$col), levels(roc_df$label))

p4A <- ggplot(roc_df, aes(x = FPR, y = TPR, color = label)) +
  geom_abline(slope=1, intercept=0, linetype="dashed",
              color="grey70", linewidth=0.4) +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  scale_color_manual(values = pal_named, name = NULL) +
  scale_x_continuous(labels = percent_format(), limits = c(0,1),
                     breaks = c(0,.25,.5,.75,1)) +
  scale_y_continuous(labels = percent_format(), limits = c(0,1),
                     breaks = c(0,.25,.5,.75,1)) +
  coord_fixed() +
  labs(x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme(legend.position        = "inside",
        legend.position.inside = c(0.62, 0.18),
        legend.text            = element_text(size = 7),
        legend.key.height      = unit(0.45, "cm"))

auc_df <- roc_tab %>%
  arrange(AUC) %>%
  mutate(tRF_short = sub("tRF-","",tRF),
         tRF_short = factor(tRF_short, levels = tRF_short),
         hi        = AUC >= 0.9)

p4B <- ggplot(auc_df, aes(x = AUC, y = tRF_short, color = hi)) +
  geom_segment(aes(x = 0.5, xend = AUC, y = tRF_short, yend = tRF_short),
               linewidth = 0.9, alpha = 0.60) +
  geom_point(size = 2.5, shape = 16) +
  ## white-background label so it masks the lollipop line behind it
  geom_label(aes(label = sprintf("%.3f", AUC)),
             hjust = -0.18, size = 2.4, color = "black",
             fill = "white", label.size = 0,
             label.padding = unit(0.08, "lines")) +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             color = "grey50", linewidth = 0.35) +
  geom_vline(xintercept = 0.9, linetype = "dotted",
             color = "#C0392B", linewidth = 0.4) +
  scale_color_manual(values = c("TRUE"="#C0392B","FALSE"="#2980B9"),
                     guide = "none") +
  scale_x_continuous(limits = c(0.25, 1.16), expand = c(0,0),
                     breaks = c(0.25, 0.5, 0.7, 0.9, 1.0)) +
  labs(x = "AUC", y = NULL) +
  theme(axis.text.y       = element_text(size = 7),
        axis.line.y       = element_blank(),
        axis.ticks.y      = element_blank(),
        panel.grid.major.x = element_line(color = "grey93", linewidth = 0.35))

fig4 <- (p4A | p4B) +
  plot_annotation(
    tag_levels = "A",
    title      = "Figure 4. Internal tumor-normal discrimination by candidate tRFs.",
    subtitle   = paste0(
      "a: ROC curves for the top 5 tRF classifiers (dashed diagonal = random classifier).  ",
      "b: Ranked AUC values for the top candidate tRFs (dotted red line = AUC 0.90)."),
    theme      = theme(
      plot.tag      = element_text(size = 10, face = "bold"),
      plot.title    = element_text(size = 9, face = "bold", hjust = 0,
                                   margin = margin(b = 3)),
      plot.subtitle = element_text(size = 7.5, color = "grey40", hjust = 0,
                                   margin = margin(b = 6)))
  ) & theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(OUT, "fig5_roc.pdf"),
       fig4, width = 12, height = 5.5,
       device = "pdf", useDingbats = FALSE)
ggsave(file.path(OUT, "fig5_roc.png"),
       fig4, width = 12, height = 5.5,
       device = "png", dpi = 600, bg = "white")
cat("  Done.\n")

## ════════════════════════════════════════════════════════════════════════════
## FIGURE 5 — Stage Association + Survival
## ════════════════════════════════════════════════════════════════════════════
cat("Figure 5: Stage + Survival...\n")

stage_trfs <- c(
  "tRF-22-7SIRMM121",
  "tRF-22-7SERML921",
  "tRF-23-7SIR3DR2DV",
  "tRF-23-7SIRMM12V",
  "tRF-25-RNLNKSEK51")
stage_trfs <- stage_trfs[stage_trfs %in% kw_tab$tRF]

stage_clin <- merge(tumor_trf, clin[, c("patient_id","stage_simple")],
                    by = "patient_id")
stage_clin <- stage_clin[
  stage_clin$stage_simple %in%
    c("Stage I","Stage II","Stage III","Stage IV"), ]
stage_clin$stage_simple <- factor(stage_clin$stage_simple,
  levels = c("Stage I","Stage II","Stage III","Stage IV"))
n_stg  <- as.data.frame(table(stage_clin$stage_simple))
x_labs <- setNames(
  paste0(sub("Stage ", "", levels(stage_clin$stage_simple)),
         "\n(n=", n_stg$Freq, ")"),
  levels(stage_clin$stage_simple))
fmt_p <- function(x) {
  if (is.na(x)) return("NA")
  if (x < 0.001) {
    exponent <- floor(log10(x))
    mantissa <- x / 10^exponent
    return(paste0(formatC(mantissa, format = "f", digits = 1),
                  " × 10^", exponent))
  }
  sprintf("%.3f", x)
}

stage_long <- bind_rows(lapply(stage_trfs, function(trf) {
  fdr_v <- kw_tab$kruskal_fdr[kw_tab$tRF == trf]
  data.frame(
    tRF_label = paste0(sub("tRF-","",trf),
                       "\nFDR = ", fmt_p(fdr_v)),
    stage     = stage_clin$stage_simple,
    expr      = as.numeric(trf_lcpm[trf, stage_clin$sample_id]))
}))
stage_long$tRF_label <- factor(
  stage_long$tRF_label,
  levels = unique(stage_long$tRF_label))

p5A <- ggplot(stage_long, aes(x = stage, y = expr, fill = stage, color = stage)) +
  geom_boxplot(outlier.shape = NA, width = 0.58, linewidth = 0.42,
               alpha = 0.32, color = "black") +
  geom_jitter(width = 0.13, height = 0, size = 0.35, alpha = 0.34, shape = 16) +
  scale_fill_manual(values  = STAGE_COLS, guide = "none") +
  scale_color_manual(values = STAGE_COLS, guide = "none") +
  scale_x_discrete(labels = x_labs) +
  facet_wrap(~tRF_label, ncol = 5, scales = "free_y") +
  labs(x = NULL, y = expression(log[2]~"(CPM + 1)")) +
  theme(panel.grid.major.y = element_line(color = "grey93", linewidth = 0.35),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.line.x        = element_line(color = "black", linewidth = 0.4),
        axis.line.y        = element_line(color = "black", linewidth = 0.4),
        axis.ticks.x       = element_line(color = "black", linewidth = 0.4),
        axis.text.x        = element_text(size = 6.2, color = "black",
                                          lineheight = 0.9),
        axis.text.y        = element_text(size = 6.2, color = "black"),
        strip.text         = element_text(size = 6.4, face = "bold",
                                          color = "black", lineheight = 0.95),
        strip.background   = element_rect(fill = "grey96", color = NA),
        panel.spacing.x    = unit(0.48, "lines"),
        panel.border       = element_blank())

cox_plot_df <- cox_tab %>%
  arrange(HR) %>%
  mutate(
    tRF_short = sub("tRF-", "", tRF),
    tRF_short = factor(tRF_short, levels = tRF_short),
    significant = FDR < 0.05)

p5B <- ggplot(cox_plot_df, aes(y = tRF_short)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey45", linewidth = 0.35) +
  geom_segment(aes(x = HR_lo, xend = HR_hi, yend = tRF_short,
                   color = significant),
               linewidth = 0.52, alpha = 0.9) +
  geom_point(aes(x = HR, fill = significant),
             shape = 21, size = 1.85, color = "black", stroke = 0.25) +
  scale_color_manual(values = c("FALSE" = "grey55", "TRUE" = "#C0392B"),
                     guide = "none") +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#C0392B"),
                    labels = c("FALSE" = "FDR >= 0.05", "TRUE" = "FDR < 0.05"),
                    name = NULL) +
  scale_x_continuous(limits = c(0.45, 1.36),
                     breaks = c(0.6, 0.8, 1.0, 1.2),
                     expand = c(0, 0)) +
  labs(x = "Hazard ratio (95% CI)", y = NULL,
       title = "Univariable Cox regression for top candidate tRFs") +
  theme(panel.grid.major.x = element_line(color = "grey92", linewidth = 0.35),
        panel.grid.major.y = element_blank(),
        axis.text.y        = element_text(size = 5.9, color = "black"),
        axis.text.x        = element_text(size = 6.7, color = "black"),
        axis.title.x       = element_text(size = 7.3),
        plot.title         = element_text(size = 8.4, face = "bold",
                                          hjust = 0.5, margin = margin(b = 4)),
        legend.position    = "top",
        legend.text        = element_text(size = 6.6),
        legend.key.size    = unit(0.25, "cm"),
        plot.margin        = margin(4, 7, 2, 5))

surv_trf <- cox_tab %>% arrange(FDR) %>% slice_head(n = 1) %>% pull(tRF)
cat("  Survival tRF:", surv_trf, "\n")
surv_d <- merge(tumor_trf, clin[, c("patient_id","os_days","event")],
                by = "patient_id")
surv_d <- surv_d[!is.na(surv_d$os_days) & surv_d$os_days > 0 &
                   surv_d$sample_id %in% colnames(trf_lcpm), ]
cat("  n =", nrow(surv_d), "\n")

p5C <- ggplot() + theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  annotate("text", x=.5, y=.5,
           label="Survival data unavailable", size=9, color="grey55")

if (surv_trf %in% rownames(trf_lcpm) && nrow(surv_d) > 20) {
  expr_s <- as.numeric(trf_lcpm[surv_trf, surv_d$sample_id])
  med    <- median(expr_s, na.rm = TRUE)
  if (med == 0) {
    surv_d$grp <- factor(
      ifelse(expr_s > 0, "Detected", "Not detected"),
      levels = c("Detected","Not detected"))
  } else {
    surv_d$grp <- factor(
      ifelse(expr_s >= med, "High", "Low"),
      levels = c("High","Low"))
  }
  cat("  Groups:", paste(names(table(surv_d$grp)), "n=",
                         as.integer(table(surv_d$grp)), collapse = "; "), "\n")
  cox_f <- coxph(Surv(os_days, event) ~ expr_s, data = surv_d)
  hr_v  <- round(exp(coef(cox_f)), 2)
  pv_v  <- summary(cox_f)$coefficients[, "Pr(>|z|)"]
  fdr_v <- cox_tab$FDR[cox_tab$tRF == surv_trf][1]
  lr    <- survdiff(Surv(os_days, event) ~ grp, data = surv_d)
  lr_p  <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
  stats_label <- paste0(
    "Univariable Cox HR = ", sprintf("%.2f", hr_v),
    "; Cox p = ", fmt_p(pv_v),
    "; Cox FDR = ", fmt_p(fdr_v),
    "; log-rank p = ", fmt_p(lr_p))

  ## Build legend labels that clearly describe the two groups
  grp_levels <- levels(surv_d$grp)
  grp_ns     <- as.integer(table(surv_d$grp))
  legend_lbs <- paste0(grp_levels, " (n=", grp_ns, ")")

  sf <- survfit(Surv(os_days, event) ~ grp, data = surv_d)
  km <- ggsurvplot(sf, data = surv_d,
    palette             = c("#C0392B","#2980B9"),
    pval                = FALSE,
    censor.shape        = 124,
    censor.size         = 2.4,
    conf.int            = TRUE,  conf.int.alpha = 0.12,
    risk.table          = "nrisk_cumcensor",
    risk.table.col      = "black",
    risk.table.height   = 0.24,
    risk.table.y.text   = FALSE,
    risk.table.y.text.col = TRUE,
    risk.table.fontsize = 2.25,
    risk.table.title    = "No. at risk",
    legend.labs         = legend_lbs,
    xlab                = "Follow-up time (days)",
    ylab                = "Overall survival probability",
    break.time.by       = 1000,
    legend.title        = paste0(sub("tRF-","",surv_trf), " expression"),
    title               = "tRF-30-81PV6RRNLNK8 expression and overall survival",
    font.title          = c(9, "bold"),
    font.x = 8, font.y = 8, font.tickslab = 7.4,
    ggtheme             = pub_theme)

  km$plot <- km$plot +
    labs(subtitle = stats_label) +
    theme(
      legend.position      = "top",
      legend.direction     = "horizontal",
      legend.background    = element_blank(),
      legend.key           = element_blank(),
      legend.text          = element_text(size = 7),
      legend.title         = element_text(size = 7.5, face = "bold"),
      plot.title           = element_text(size = 9, face = "bold",
                                          hjust = 0.5, margin = margin(b = 3)),
      plot.subtitle        = element_text(size = 6.8, color = "grey30",
                                          hjust = 0.5, margin = margin(b = 5)),
      plot.background      = element_rect(fill = "white", color = NA),
      panel.background     = element_rect(fill = "white", color = NA),
      plot.margin          = margin(4, 6, 2, 6))

  km$table <- km$table +
    labs(x = NULL, y = NULL) +
    theme(
      legend.position  = "none",
      axis.title.y     = element_blank(),
      axis.title.x     = element_blank(),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title       = element_text(size = 6.5, face = "bold", hjust = 0),
      axis.text.x      = element_text(size = 6.1),
      plot.margin      = margin(0, 6, 0, 6))

  p5C <- cowplot::plot_grid(
    km$plot, km$table,
    nrow = 2, rel_heights = c(1, 0.30),
    align = "v", axis = "lr")
  p5C <- cowplot::plot_grid(NULL, p5C, ncol = 2, rel_widths = c(0.035, 1))
}

p5_bottom <- cowplot::plot_grid(
  p5B, p5C, ncol = 2,
  labels = c("B","C"),
  label_size = 10, label_fontface = "bold",
  rel_widths = c(0.43, 0.57))

fig5_body <- cowplot::plot_grid(
  p5A, p5_bottom, nrow = 2,
  rel_heights = c(0.80, 1.05),
  labels = c("A",""),
  label_size = 10, label_fontface = "bold")
fig5_caption_txt <- wrap_caption(
  paste0(
    "(A) Expression of candidate tRFs across AJCC pathologic stages. FDR values are from Kruskal-Wallis tests with ",
    "Benjamini-Hochberg correction. (B) Univariable Cox proportional-hazards regression for the top 20 candidate ",
    "tRFs; points show hazard ratios and horizontal lines show 95% confidence intervals. (C) Kaplan-Meier overall ",
    "survival curves stratified by detected versus not detected tRF-30-81PV6RRNLNK8 expression. Survival analyses are exploratory."),
  width = 238)
fig5_header <- cowplot::ggdraw() +
  cowplot::draw_label(
    "Figure 5. Clinical associations of candidate tRF expression in TCGA-KIRP.",
    fontface = "bold", size = 9,
    x = 0.003, y = 0.76, hjust = 0, vjust = 0.5) +
  cowplot::draw_label(
    fig5_caption_txt,
    size = 6.25, color = "grey40",
    x = 0.003, y = 0.43, hjust = 0, vjust = 0.5,
    lineheight = 1.02)
fig5_layout <- cowplot::plot_grid(fig5_header, fig5_body,
  nrow = 2, rel_heights = c(0.092, 1))
fig5_full <- cowplot::ggdraw() +
  cowplot::draw_grob(rectGrob(gp = gpar(fill = "white", col = NA))) +
  cowplot::draw_plot(fig5_layout, x = 0, y = 0, width = 1, height = 1)

ggsave(file.path(OUT, "fig6_stage_survival.pdf"),
       fig5_full, width = 14, height = 10,
       device = "pdf", useDingbats = FALSE, bg = "white")
ggsave(file.path(OUT, "fig6_stage_survival.png"),
       fig5_full, width = 14, height = 10,
       device = "png", dpi = 600, bg = "white")
ggsave(file.path(OUT, "fig6_stage_survival_revised.pdf"),
       fig5_full, width = 14, height = 10,
       device = "pdf", useDingbats = FALSE, bg = "white")
ggsave(file.path(OUT, "fig6_stage_survival_revised.png"),
       fig5_full, width = 14, height = 10,
       device = "png", dpi = 600, bg = "white")
cat("  Done.\n")

## ════════════════════════════════════════════════════════════════════════════
## Summary
## ════════════════════════════════════════════════════════════════════════════
cat("\n===========================================\n")
figs   <- c("fig2_trf_DE.png","fig3_gene_DE.png","fig3_gene_DE_revised.png",
            "fig4_trf_met_axis.png","fig5_roc.png","fig6_stage_survival.png",
            "fig6_stage_survival_revised.png")
labels <- c("Fig 1 (tRF DE)","Fig 2 (Gene DE)","Fig 2 revised",
            "Fig 3 (tRF-MET)","Fig 4 (ROC/AUC)","Fig 5 (Stage/Cox/Surv)",
            "Fig 5 revised")
for (i in seq_along(figs)) {
  fp <- file.path(OUT, figs[i])
  ok <- file.exists(fp)
  cat(sprintf("  %s  %-22s  %-34s  %s KB\n",
    if (ok) "[OK]" else "[!!]", labels[i], figs[i],
    if (ok) round(file.info(fp)$size / 1024) else "MISSING"))
}
cat("===========================================\n")
