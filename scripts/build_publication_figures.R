## Publication-Ready Main Figures
## Uses patchwork for layout, cowplot for panel labels, consistent theme

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(patchwork)
  library(cowplot)
  library(ggplotify)
  library(scales)
  library(dplyr)
  library(tidyr)
  library(grid)
  library(gridExtra)
  library(RColorBrewer)
  library(pROC)
  library(survival)
  library(survminer)
})

# Edit scripts/config.R to set ROOT for your machine, or change the fallback below:
if (file.exists("scripts/config.R")) source("scripts/config.R") else if (file.exists("config.R")) source("config.R") else ROOT <- "/Users/nam.tnguyen2022/KIRP_tRF_project"
DATA <- file.path(ROOT, "data_processed")
TABS <- "/Users/nam.tnguyen2022/KIRP_tRF_project/results/tables"
RAW  <- "/Users/nam.tnguyen2022/KIRP_tRF_project/data_raw"
OUT  <- "/Users/nam.tnguyen2022/KIRP_tRF_project/results/figures"

# ── Global publication theme ─────────────────────────────────────────────────
pub_theme <- theme_classic(base_size = 11, base_family = "Helvetica") +
  theme(
    plot.title        = element_text(face = "bold", size = 12, hjust = 0),
    axis.title        = element_text(size = 10),
    axis.text         = element_text(size = 9, color = "black"),
    axis.line         = element_line(color = "black", linewidth = 0.4),
    axis.ticks        = element_line(color = "black", linewidth = 0.4),
    legend.text       = element_text(size = 9),
    legend.title      = element_text(size = 9, face = "bold"),
    legend.key.size   = unit(0.4, "cm"),
    strip.background  = element_rect(fill = "#F0F0F0", color = NA),
    strip.text        = element_text(face = "bold", size = 9),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    plot.margin       = margin(8, 8, 8, 8)
  )

theme_set(pub_theme)

GROUP_COLS <- c(Tumor = "#C0392B", Normal = "#2980B9")
STAGE_COLS <- c("Stage I"="#2166AC","Stage II"="#74ADD1",
                "Stage III"="#F46D43","Stage IV"="#A50026")

# ── Load data ────────────────────────────────────────────────────────────────
cat("Loading data...\n")

trf_mat   <- read.table(file.path(DATA,"trf_count_matrix.tsv"),
                        sep="\t", header=TRUE, row.names=1, check.names=FALSE)
trf_meta  <- read.table(file.path(DATA,"trf_sample_metadata.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE,
                        colClasses=c(sample_type_code="character"))
trf_de    <- read.table(file.path(TABS,"trf_DE_significant.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE)
trf_full  <- read.table(file.path(TABS,"trf_DE_full.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE)
gene_full <- read.table(file.path(TABS,"gene_DE_full.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE)
gene_sig  <- read.table(file.path(TABS,"gene_DE_significant.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE)
anticorr  <- read.table(file.path(TABS,"anticorrelation_prioritized.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE)
roc_tab   <- read.table(file.path(TABS,"roc_auc_top_trfs.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE)
cox_tab   <- read.table(file.path(TABS,"cox_univariate_results.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE)
kw_tab    <- read.table(file.path(TABS,"stage_kruskal_results.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE)
rna22_sum <- read.table(file.path(TABS,"rna22_binding_summary.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE)
mrna_mat  <- read.table(file.path(DATA,"mrna_count_matrix.tsv"),
                        sep="\t", header=TRUE, row.names=1, check.names=FALSE)
mrna_meta <- read.table(file.path(DATA,"mrna_sample_metadata.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE,
                        colClasses=c(sample_type_code="character"))

# CPM
trf_lcpm  <- log2(sweep(trf_mat, 2, colSums(trf_mat), "/") * 1e6 + 1)
mrna_lcpm <- log2(sweep(mrna_mat, 2, colSums(mrna_mat), "/") * 1e6 + 1)
trf_meta$group  <- ifelse(trf_meta$sample_type_code == "01","Tumor","Normal")
mrna_meta$group <- ifelse(mrna_meta$sample_type_code == "01","Tumor","Normal")
rownames(trf_meta)  <- trf_meta$sample_id
rownames(mrna_meta) <- mrna_meta$sample_id

# Gene id map
star_files <- list.files(
  file.path(RAW,"gdc_rnaseq/gdc_download_20260506_091015.911939"),
  pattern="*.rna_seq.augmented_star_gene_counts.tsv",
  recursive=TRUE, full.names=TRUE)
gmap <- read.table(star_files[1], sep="\t", header=TRUE,
                   comment.char="#")[, c("gene_id","gene_name")]

# Clinical
clin_raw <- read.table(file.path(RAW,"gdc_clinical/clinical.tsv"),
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
primary_stage$stage_simple <- sub(" [ABC]$","",primary_stage$stage)
clin$stage <- NULL
clin <- merge(clin, primary_stage[, c("patient_id","stage","stage_simple")],
              by="patient_id", all.x=TRUE, sort=FALSE)
cn <- function(x) { x[x %in% c("'--","--","")] <- NA; as.numeric(x) }
clin$days_death <- cn(clin$days_death); clin$days_fu <- cn(clin$days_fu)
clin$event   <- as.integer(clin$vital_status == "Dead")
clin$os_days <- ifelse(!is.na(clin$days_death), clin$days_death, clin$days_fu)

tumor_trf  <- trf_meta[trf_meta$group == "Tumor",]
tumor_mrna <- mrna_meta[mrna_meta$group == "Tumor",]
tumor_trf$patient_id  <- substr(tumor_trf$sample_id, 1, 12)
tumor_mrna$patient_id <- substr(tumor_mrna$sample_id, 1, 12)

cat("Data loaded.\n\n")

# ════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Cohort Overview + Pipeline
# ════════════════════════════════════════════════════════════════════════════
cat("Figure 1...\n")

cohort_df <- data.frame(
  Platform = factor(rep(c("tRF\n(MINTbase v2.0)","mRNA\n(GDC STAR)"), each=2),
                    levels=c("tRF\n(MINTbase v2.0)","mRNA\n(GDC STAR)")),
  Group    = factor(rep(c("Tumor","Normal"), 2), levels=c("Tumor","Normal")),
  N        = c(291, 34, 290, 32)
)

p1A <- ggplot(cohort_df, aes(x=Platform, y=N, fill=Group)) +
  geom_col(position=position_dodge(0.65), width=0.6,
           color="white", linewidth=0.3) +
  geom_text(aes(label=N), position=position_dodge(0.65),
            vjust=-0.5, size=3.5, fontface="bold") +
  scale_fill_manual(values=GROUP_COLS, name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0,0.15)),
                     breaks=c(0,50,100,150,200,250,300)) +
  labs(x=NULL, y="Number of Samples") +
  theme(legend.position=c(0.88,0.88),
        legend.background=element_rect(fill="white", color="grey80", linewidth=0.3))

pipeline_df <- data.frame(
  y    = 7:1,
  step = paste0("Step ", 7:1),
  phase = c("Manuscript",
            "ROC / Stage / Survival",
            "RNA22 Binding Prediction",
            "Spearman Anticorrelation",
            "Differential Expression",
            "Count Matrix Construction",
            "Data Acquisition"),
  result = c("",
             "AUC=0.949 | HR=1.16 (FDR=0.033)",
             "18 tRF-gene binding pairs",
             "104 anticorrelated pairs",
             "1,585 tRFs | 3 genes",
             "30,076 tRFs | 60,660 genes",
             "MINTbase v2.0 + GDC STAR counts")
)

p1B <- ggplot(pipeline_df) +
  geom_rect(aes(xmin=0, xmax=10, ymin=y-0.42, ymax=y+0.42),
            fill="#EAF4FB", color="#AED6F1", linewidth=0.4) +
  geom_text(aes(x=0.3, y=y, label=phase), hjust=0, size=3.3, fontface="bold") +
  geom_text(aes(x=9.8, y=y, label=result), hjust=1, size=2.9,
            color="grey30", fontface="italic") +
  scale_x_continuous(limits=c(0,10), expand=c(0,0)) +
  scale_y_continuous(limits=c(0.4, 7.6)) +
  labs(x=NULL, y=NULL) +
  theme_void() +
  theme(plot.margin=margin(4,4,4,4))

fig1 <- (p1A | p1B) +
  plot_annotation(
    title   = "Figure 1. Study cohort and analytical pipeline.",
    caption = "TCGA-KIRP | May 2026",
    theme   = theme(
      plot.title   = element_text(size=11, face="bold", hjust=0),
      plot.caption = element_text(size=8, color="grey50", hjust=1)
    )
  ) +
  plot_layout(widths=c(1,1.4)) &
  theme(plot.background=element_rect(fill="white", color=NA))

ggsave(file.path(OUT,"fig1_cohort_workflow.pdf"),
       fig1, width=10, height=4.5, device=cairo_pdf)
cat("  Done.\n")

# ════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — tRF Differential Expression
# ════════════════════════════════════════════════════════════════════════════
cat("Figure 2...\n")

# Panel A: Volcano
trf_full$sig <- case_when(
  trf_full$FDR < 0.05 & trf_full$logFC >  1 ~ "Up",
  trf_full$FDR < 0.05 & trf_full$logFC < -1 ~ "Down",
  TRUE ~ "NS")
trf_full$sig <- factor(trf_full$sig, levels=c("Up","Down","NS"))

n_up   <- sum(trf_full$sig == "Up")
n_down <- sum(trf_full$sig == "Down")
top_lab <- trf_full %>% filter(sig != "NS") %>%
  arrange(FDR) %>% slice_head(n=10)

p2A <- ggplot(trf_full, aes(x=logFC, y=-log10(FDR), color=sig)) +
  geom_point(data=~filter(.x, sig=="NS"),
             size=0.4, alpha=0.3, color="grey70") +
  geom_point(data=~filter(.x, sig!="NS"),
             size=0.7, alpha=0.7) +
  geom_text_repel(data=top_lab, aes(label=trf_id),
                  size=2.4, segment.size=0.3, segment.color="grey50",
                  min.segment.length=0.2, max.overlaps=12,
                  show.legend=FALSE, fontface="italic") +
  geom_vline(xintercept=c(-1,1), linetype="dashed",
             color="grey40", linewidth=0.4) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed",
             color="grey40", linewidth=0.4) +
  scale_color_manual(values=c(Up="#C0392B", Down="#2980B9", NS="grey70"),
                     labels=c(paste0("Up (",n_up,")"),
                               paste0("Down (",n_down,")"), "NS"),
                     name=NULL) +
  scale_x_continuous(limits=c(-6, 8)) +
  annotate("text", x=-5.5, y=max(-log10(trf_full$FDR[is.finite(trf_full$FDR)]),na.rm=TRUE)*0.95,
           label=paste0(n_down," down"), color="#2980B9", size=3.2, hjust=0, fontface="bold") +
  annotate("text", x=5.5,  y=max(-log10(trf_full$FDR[is.finite(trf_full$FDR)]),na.rm=TRUE)*0.95,
           label=paste0(n_up," up"), color="#C0392B", size=3.2, hjust=1, fontface="bold") +
  labs(x=expression(log[2]~"Fold Change (Tumor/Normal)"),
       y=expression(-log[10]~"(FDR)")) +
  theme(legend.position="none")

# Panel B: Heatmap (top 30 by |logFC|, made with pheatmap -> ggplotify)
top30_ids <- trf_de %>% arrange(desc(abs(logFC))) %>% slice_head(n=30) %>% pull(trf_id)
top30_ids <- top30_ids[top30_ids %in% rownames(trf_lcpm)]
col_order  <- c(which(trf_meta$group=="Normal"), which(trf_meta$group=="Tumor"))
mat30      <- trf_lcpm[top30_ids, trf_meta$sample_id[col_order]]
ann_col    <- data.frame(Group=trf_meta$group[col_order],
                         row.names=colnames(mat30))
ann_colors <- list(Group=GROUP_COLS)

p2B <- as.ggplot(pheatmap(mat30,
  annotation_col  = ann_col,
  annotation_colors = ann_colors,
  show_colnames   = FALSE,
  color           = colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
  scale           = "row",
  cluster_cols    = FALSE,
  fontsize_row    = 7,
  fontsize        = 8,
  border_color    = NA,
  annotation_legend = TRUE,
  main            = "",
  silent          = TRUE))

# Panel C: Top 9 tRFs boxplots
top9 <- trf_de %>% arrange(FDR) %>% slice_head(n=9) %>% pull(trf_id)
top9 <- top9[top9 %in% rownames(trf_lcpm)]
box_df <- bind_rows(lapply(top9, function(t)
  data.frame(tRF=t, expr=as.numeric(trf_lcpm[t, trf_meta$sample_id]),
             Group=trf_meta$group)))
box_df$tRF <- factor(box_df$tRF, levels=top9)
# Short labels
box_df$tRF_short <- factor(sub("tRF-","",box_df$tRF), levels=sub("tRF-","",top9))

p2C <- ggplot(box_df, aes(x=Group, y=expr, fill=Group)) +
  geom_boxplot(outlier.size=0.3, outlier.alpha=0.4,
               linewidth=0.4, width=0.65) +
  stat_summary(fun=median, geom="point", shape=18, size=1.5,
               color="white", show.legend=FALSE) +
  facet_wrap(~tRF_short, nrow=3, scales="free_y") +
  scale_fill_manual(values=GROUP_COLS, name=NULL) +
  labs(x=NULL, y=expression(log[2]~"(CPM+1)")) +
  theme(axis.text.x=element_text(angle=30, hjust=1, size=8),
        legend.position="none",
        strip.text=element_text(size=7.5))

# Assemble Figure 2
fig2 <- (p2A | p2B) / p2C +
  plot_annotation(
    tag_levels  = "A",
    title       = "Figure 2. Differential expression of tRFs in KIRP tumor vs. normal tissue.",
    theme       = theme(
      plot.tag   = element_text(size=13, face="bold"),
      plot.title = element_text(size=11, face="bold", hjust=0))
  ) +
  plot_layout(heights=c(1.2, 1)) &
  theme(plot.background=element_rect(fill="white", color=NA))

ggsave(file.path(OUT,"fig2_trf_DE.pdf"),
       fig2, width=12, height=11, device=cairo_pdf)
cat("  Done.\n")

# ════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Candidate Gene Expression
# ════════════════════════════════════════════════════════════════════════════
cat("Figure 3...\n")

# Panel A: gene volcano
gene_full$sig <- case_when(
  gene_full$FDR < 0.05 & gene_full$logFC >  1 ~ "Up",
  gene_full$FDR < 0.05 & gene_full$logFC < -1 ~ "Down",
  gene_full$FDR < 0.05 ~ "Sig (|FC|<2)",
  TRUE ~ "NS")
gene_full$sig <- factor(gene_full$sig,
                        levels=c("Up","Down","Sig (|FC|<2)","NS"))

p3A <- ggplot(gene_full, aes(x=logFC, y=-log10(FDR), color=sig)) +
  geom_point(size=3, alpha=0.85) +
  geom_text_repel(aes(label=gene_name), size=2.8,
                  segment.size=0.3, show.legend=FALSE,
                  max.overlaps=20, fontface="italic",
                  box.padding=0.4, min.segment.length=0.1) +
  geom_vline(xintercept=c(-1,1), linetype="dashed",
             color="grey40", linewidth=0.4) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed",
             color="grey40", linewidth=0.4) +
  scale_color_manual(
    values=c(Up="#C0392B", Down="#2980B9",
             "Sig (|FC|<2)"="#8E44AD", NS="grey70"),
    name=NULL) +
  labs(x=expression(log[2]~"Fold Change"),
       y=expression(-log[10]~"(FDR)")) +
  theme(legend.position=c(0.78,0.88),
        legend.background=element_rect(fill="white",color="grey80",linewidth=0.3))

# Panels B–D: MET, CDKN2A, TERT
make_gene_box <- function(gname) {
  gid  <- gmap$gene_id[gmap$gene_name == gname][1]
  if (is.na(gid) || !gid %in% rownames(mrna_lcpm)) return(NULL)
  lfc  <- round(gene_full$logFC[gene_full$gene_name == gname], 2)
  fdr  <- formatC(gene_full$FDR[gene_full$gene_name == gname],
                  format="e", digits=1)
  df   <- data.frame(
    expr  = as.numeric(mrna_lcpm[gid, mrna_meta$sample_id]),
    Group = mrna_meta$group)
  ggplot(df, aes(x=Group, y=expr, fill=Group)) +
    geom_boxplot(outlier.size=0.4, outlier.alpha=0.5,
                 linewidth=0.4, width=0.6) +
    stat_summary(fun=median, geom="point", shape=18,
                 size=2, color="white", show.legend=FALSE) +
    scale_fill_manual(values=GROUP_COLS) +
    annotate("text", x=1.5, y=max(df$expr)*0.98,
             label=paste0("FC=",lfc,"\nFDR=",fdr),
             size=2.8, hjust=0.5, vjust=1, color="grey20") +
    labs(x=NULL, y=expression(log[2]~"(CPM+1)"),
         title=bquote(italic(.(gname)))) +
    theme(legend.position="none",
          plot.title=element_text(face="bold.italic", hjust=0.5))
}

p3B <- make_gene_box("MET")
p3C <- make_gene_box("CDKN2A")
p3D <- make_gene_box("TERT")

fig3 <- p3A / (p3B | p3C | p3D) +
  plot_annotation(
    tag_levels  = "A",
    title       = "Figure 3. Differential expression of KIRP candidate genes.",
    theme       = theme(
      plot.tag   = element_text(size=13, face="bold"),
      plot.title = element_text(size=11, face="bold", hjust=0))
  ) +
  plot_layout(heights=c(1.3, 1)) &
  theme(plot.background=element_rect(fill="white", color=NA))

ggsave(file.path(OUT,"fig3_gene_DE.pdf"),
       fig3, width=10, height=9, device=cairo_pdf)
cat("  Done.\n")

# ════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — tRF–MET Regulatory Axis
# ════════════════════════════════════════════════════════════════════════════
cat("Figure 4...\n")

met_id <- gmap$gene_id[gmap$gene_name == "MET"][1]
shared  <- intersect(tumor_trf$patient_id, tumor_mrna$patient_id)
tp_trf  <- tumor_trf[match(shared, tumor_trf$patient_id),]
tp_mrna <- tumor_mrna[match(shared, tumor_mrna$patient_id),]

# Top 5 anticorrelated tRF-MET pairs
top5 <- anticorr %>% filter(gene=="MET") %>%
  arrange(FDR) %>% slice_head(n=5)

scatter_list <- lapply(seq_len(nrow(top5)), function(i) {
  trf <- top5$trf[i]
  if (!trf %in% rownames(trf_lcpm) || !met_id %in% rownames(mrna_lcpm))
    return(NULL)
  df <- data.frame(
    x = as.numeric(trf_lcpm[trf,  tp_trf$sample_id]),
    y = as.numeric(mrna_lcpm[met_id, tp_mrna$sample_id]))
  r   <- round(top5$spearman_r[i], 3)
  fdr <- formatC(top5$FDR[i], format="e", digits=1)
  short_name <- sub("tRF-","",trf)
  ggplot(df, aes(x=x, y=y)) +
    geom_point(size=0.9, alpha=0.45, color="#5D6D7E") +
    geom_smooth(method="lm", se=TRUE, color="#C0392B",
                fill="#F1948A", alpha=0.2, linewidth=0.9) +
    annotate("text", x=Inf, y=Inf, hjust=1.05, vjust=1.6,
             size=2.8, color="grey20",
             label=paste0("r = ",r,"\nFDR = ",fdr)) +
    labs(x=expression(log[2]~"(CPM+1)"),
         y=expression(italic(MET)~log[2]~"(CPM+1)"),
         title=paste0("tRF-",short_name)) +
    theme(plot.title=element_text(size=8.5, face="bold"))
})
scatter_list <- Filter(Negate(is.null), scatter_list)

# RNA22 binding table panel
rna22_met <- rna22_sum %>%
  filter(gene == "MET") %>%
  arrange(p_value) %>%
  slice_head(n=8) %>%
  transmute(tRF     = sub("tRF-","",tRF_id),
            Position = target_position,
            `ΔG`    = paste0(folding_energy_kcal," kcal/mol"),
            `p-value`= signif(p_value, 3))

tbl_p <- ggplot() +
  theme_void() +
  annotation_custom(
    tableGrob(rna22_met, rows=NULL,
              theme=ttheme_minimal(
                base_size=8,
                core=list(
                  fg_params=list(hjust=0, x=0.04, fontsize=7.5),
                  bg_params=list(fill=c("white","#F2F3F4"), col=NA)),
                colhead=list(
                  fg_params=list(fontsize=8.5, fontface="bold", hjust=0, x=0.04),
                  bg_params=list(fill="#D5E8D4", col=NA))))) +
  labs(title="RNA22 Predicted MET Binding Sites") +
  theme(plot.title=element_text(size=9.5, face="bold", margin=margin(b=4)))

# Assemble with scatter 2x2 top + RNA22 table bottom right
fig4_top  <- wrap_plots(scatter_list[1:4], ncol=2)
fig4_bot  <- wrap_plots(scatter_list[5], tbl_p, ncol=2, widths=c(1,1.6))

fig4 <- fig4_top / fig4_bot +
  plot_annotation(
    tag_levels  = "A",
    title       = "Figure 4. The tRF-MET regulatory axis in KIRP.",
    subtitle    = "Top anticorrelated tRF-MET pairs across 290 matched tumor samples (Spearman, BH-adjusted)",
    theme       = theme(
      plot.tag      = element_text(size=13, face="bold"),
      plot.title    = element_text(size=11, face="bold", hjust=0),
      plot.subtitle = element_text(size=9, color="grey40", hjust=0))
  ) +
  plot_layout(heights=c(1, 0.8)) &
  theme(plot.background=element_rect(fill="white", color=NA))

ggsave(file.path(OUT,"fig4_trf_met_axis.pdf"),
       fig4, width=11, height=10, device=cairo_pdf)
cat("  Done.\n")

# ════════════════════════════════════════════════════════════════════════════
# FIGURE 5 — ROC / AUC
# ════════════════════════════════════════════════════════════════════════════
cat("Figure 5...\n")

response <- as.numeric(trf_meta$group == "Tumor")
top5_roc <- roc_tab$tRF[1:5]
pal5     <- c("#E74C3C","#2980B9","#27AE60","#8E44AD","#E67E22")

# Build ROC curves as data frames
roc_data <- bind_rows(lapply(seq_along(top5_roc), function(i) {
  trf  <- top5_roc[i]
  if (!trf %in% rownames(trf_lcpm)) return(NULL)
  expr <- as.numeric(trf_lcpm[trf,])
  r    <- roc(response, expr, quiet=TRUE)
  auc_val <- round(as.numeric(auc(r)), 3)
  data.frame(
    FPR   = 1 - r$specificities,
    TPR   = r$sensitivities,
    tRF   = paste0(sub("tRF-","",trf),"\n(AUC=",auc_val,")"),
    color = pal5[i])
}))

p5A <- ggplot(roc_data, aes(x=FPR, y=TPR, color=tRF)) +
  geom_abline(slope=1, intercept=0, linetype="dashed",
              color="grey50", linewidth=0.5) +
  geom_line(linewidth=0.9, alpha=0.9) +
  scale_color_manual(values=setNames(pal5,
    unique(roc_data$tRF[!duplicated(roc_data$tRF)])),
    name=NULL) +
  scale_x_continuous(labels=percent_format(), limits=c(0,1)) +
  scale_y_continuous(labels=percent_format(), limits=c(0,1)) +
  labs(x="False Positive Rate",
       y="True Positive Rate",
       title="ROC Curves — Top 5 tRFs") +
  theme(legend.position=c(0.62,0.25),
        legend.text=element_text(size=7.5),
        legend.key.height=unit(0.55,"cm"),
        legend.background=element_rect(fill="white",color="grey80",linewidth=0.3))

# Panel B: AUC barplot (all 20)
roc_plot <- roc_tab %>% arrange(AUC) %>%
  mutate(tRF_short = sub("tRF-","",tRF),
         tRF_short = factor(tRF_short, levels=tRF_short))

p5B <- ggplot(roc_plot, aes(x=AUC, y=tRF_short, fill=AUC)) +
  geom_col(width=0.75, color=NA) +
  geom_text(aes(label=round(AUC,3)), hjust=-0.1, size=2.8) +
  geom_vline(xintercept=0.5, linetype="dashed", color="grey40", linewidth=0.4) +
  geom_vline(xintercept=0.9, linetype="dotted", color="#C0392B", linewidth=0.4) +
  scale_fill_gradientn(
    colors=c("#AED6F1","#2980B9","#1A5276"),
    values=rescale(c(0.5, 0.8, 1.0)), guide="none") +
  scale_x_continuous(limits=c(0, 1.05), expand=c(0,0),
                     breaks=c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(x="AUC", y=NULL, title="AUC — Top 20 tRFs") +
  theme(axis.text.y=element_text(size=7.5))

fig5 <- p5A | p5B +
  plot_annotation(
    tag_levels  = "A",
    title       = "Figure 5. Diagnostic performance of tRFs distinguishing tumor from normal tissue.",
    theme       = theme(
      plot.tag   = element_text(size=13, face="bold"),
      plot.title = element_text(size=11, face="bold", hjust=0))
  ) &
  theme(plot.background=element_rect(fill="white", color=NA))

ggsave(file.path(OUT,"fig5_roc.pdf"),
       fig5, width=11, height=5.5, device=cairo_pdf)
cat("  Done.\n")

# ════════════════════════════════════════════════════════════════════════════
# FIGURE 6 — Stage + Survival
# ════════════════════════════════════════════════════════════════════════════
cat("Figure 6...\n")

# Panel A: Stage boxplots for the 2 significant tRFs
stage_trfs <- kw_tab %>% filter(kruskal_fdr < 0.05) %>%
  arrange(kruskal_fdr) %>% slice_head(n=2) %>% pull(tRF)

stage_clin <- merge(tumor_trf, clin[,c("patient_id","stage_simple")],
                    by="patient_id")
stage_clin <- stage_clin[stage_clin$stage_simple %in%
                           c("Stage I","Stage II","Stage III","Stage IV"),]
stage_clin$stage_simple <- factor(stage_clin$stage_simple,
  levels=c("Stage I","Stage II","Stage III","Stage IV"))

stage_long <- bind_rows(lapply(stage_trfs, function(trf) {
  fdr_val <- kw_tab$kruskal_fdr[kw_tab$tRF == trf]
  data.frame(
    tRF_label = paste0(sub("tRF-","",trf),"\n(KW FDR=",round(fdr_val,3),")"),
    stage     = stage_clin$stage_simple,
    expr      = as.numeric(trf_lcpm[trf, stage_clin$sample_id]))
}))
n_stage <- stage_clin %>% count(stage_simple)
x_lab <- paste0(levels(stage_clin$stage_simple), "\n(n=",n_stage$n,")")

p6A <- ggplot(stage_long, aes(x=stage, y=expr, fill=stage)) +
  geom_boxplot(outlier.size=0.4, outlier.alpha=0.5,
               linewidth=0.4, width=0.65) +
  stat_summary(fun=median, geom="point", shape=18,
               size=2, color="white", show.legend=FALSE) +
  scale_fill_manual(values=STAGE_COLS, name=NULL) +
  scale_x_discrete(labels=x_lab) +
  facet_wrap(~tRF_label, nrow=2, scales="free_y") +
  labs(x=NULL, y=expression(log[2]~"(CPM+1)"),
       title="Stage Association") +
  theme(axis.text.x=element_text(angle=30, hjust=1, size=8),
        legend.position="none")

# Panel B: KM curve
surv_trf <- "tRF-30-81PV6RRNLNK8"
surv_data <- merge(tumor_trf, clin[,c("patient_id","os_days","event")],
                   by="patient_id")
surv_data <- surv_data[!is.na(surv_data$os_days) & surv_data$os_days > 0,]

km_p <- NULL
if (surv_trf %in% rownames(trf_lcpm) && nrow(surv_data) > 10) {
  expr_s <- as.numeric(trf_lcpm[surv_trf, surv_data$sample_id])
  surv_data$km_group <- factor(
    ifelse(expr_s >= median(expr_s, na.rm=TRUE), "High","Low"),
    levels=c("High","Low"))
  sf  <- survfit(Surv(os_days, event) ~ km_group, data=surv_data)
  cox <- coxph(Surv(os_days, event) ~ expr_s, data=surv_data)
  hr  <- round(exp(coef(cox)), 2)
  pv  <- round(summary(cox)$coefficients[,"Pr(>|z|)"], 3)

  km_p <- ggsurvplot(sf, data=surv_data,
    palette    = c("#C0392B","#2980B9"),
    pval       = TRUE,
    pval.size  = 3.5,
    conf.int   = TRUE,
    conf.int.alpha = 0.1,
    risk.table = TRUE,
    risk.table.height = 0.25,
    risk.table.fontsize = 3,
    xlab       = "Time (days)",
    ylab       = "Overall Survival",
    legend.title = "Expression",
    title      = paste0("tRF-30-81PV6RRNLNK8\n(HR=",hr,", p=",pv,")"),
    font.title = c(9, "bold"),
    font.x = 9, font.y = 9, font.tickslab = 8,
    ggtheme    = pub_theme)
}

# Save Figure 6
pdf(file.path(OUT,"fig6_stage_survival.pdf"), width=11, height=7)
if (!is.null(km_p)) {
  # Left: stage panel, Right: KM
  pushViewport(viewport(layout=grid.layout(1, 2, widths=unit(c(0.42,0.58),"npc"))))
  pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
  print(p6A + plot_annotation(tag_levels=list(c("A")),
    theme=theme(plot.tag=element_text(size=13,face="bold"))), newpage=FALSE)
  popViewport()
  pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
  print(km_p, newpage=FALSE)
  popViewport()
  popViewport()
  grid.text("B", x=0.44, y=0.97, gp=gpar(fontsize=13, fontface="bold"))
  grid.text("Figure 6. Clinical associations of tRF expression in KIRP.",
            x=0.5, y=0.005, hjust=0.5, vjust=0,
            gp=gpar(fontsize=9, fontface="bold"))
} else {
  print(p6A)
}
dev.off()
cat("  Done.\n")

cat("\n=====================================\n")
cat("All 6 publication-ready figures saved:\n")
for (f in c("fig1_cohort_workflow.pdf","fig2_trf_DE.pdf","fig3_gene_DE.pdf",
            "fig4_trf_met_axis.pdf","fig5_roc.pdf","fig6_stage_survival.pdf")) {
  fp <- file.path(OUT,f)
  if (file.exists(fp))
    cat(sprintf("  [OK] %-36s %s KB\n", f,
                round(file.info(fp)$size/1024)))
}
cat("=====================================\n")
