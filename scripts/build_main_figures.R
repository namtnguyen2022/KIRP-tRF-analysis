## Build 6 composite main figures for the manuscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(gridExtra)
  library(grid)
  library(pROC)
  library(survival)
  library(survminer)
  library(dplyr)
})

DATA  <- "/Users/nam.tnguyen2022/KIRP_tRF_project/data_processed"
TABS  <- "/Users/nam.tnguyen2022/KIRP_tRF_project/results/tables"
RAW   <- "/Users/nam.tnguyen2022/KIRP_tRF_project/data_raw"
OUT   <- "/Users/nam.tnguyen2022/KIRP_tRF_project/results/figures"

# ── Shared helpers ────────────────────────────────────────────────────────────
label_panel <- function(label, x=0.01, y=0.97) {
  annotation_custom(
    grob = textGrob(label, x=x, y=y, just=c("left","top"),
                    gp=gpar(fontsize=14, fontface="bold")))
}

# ── Load count matrices and metadata ─────────────────────────────────────────
cat("Loading data...\n")
trf_mat  <- read.table(file.path(DATA,"trf_count_matrix.tsv"),
                       sep="\t", header=TRUE, row.names=1, check.names=FALSE)
trf_meta <- read.table(file.path(DATA,"trf_sample_metadata.tsv"),
                       sep="\t", header=TRUE, stringsAsFactors=FALSE,
                       colClasses=c(sample_type_code="character"))
trf_de   <- read.table(file.path(TABS,"trf_DE_significant.tsv"),
                       sep="\t", header=TRUE, stringsAsFactors=FALSE)
trf_de_full <- read.table(file.path(TABS,"trf_DE_full.tsv"),
                          sep="\t", header=TRUE, stringsAsFactors=FALSE)
gene_de  <- read.table(file.path(TABS,"gene_DE_significant.tsv"),
                       sep="\t", header=TRUE, stringsAsFactors=FALSE)
gene_de_full <- read.table(file.path(TABS,"gene_DE_full.tsv"),
                           sep="\t", header=TRUE, stringsAsFactors=FALSE)
anticorr <- read.table(file.path(TABS,"anticorrelation_prioritized.tsv"),
                       sep="\t", header=TRUE, stringsAsFactors=FALSE)
roc_df   <- read.table(file.path(TABS,"roc_auc_top_trfs.tsv"),
                       sep="\t", header=TRUE, stringsAsFactors=FALSE)
cox_df   <- read.table(file.path(TABS,"cox_univariate_results.tsv"),
                       sep="\t", header=TRUE, stringsAsFactors=FALSE)
kw_df    <- read.table(file.path(TABS,"stage_kruskal_results.tsv"),
                       sep="\t", header=TRUE, stringsAsFactors=FALSE)

mrna_mat  <- read.table(file.path(DATA,"mrna_count_matrix.tsv"),
                        sep="\t", header=TRUE, row.names=1, check.names=FALSE)
mrna_meta <- read.table(file.path(DATA,"mrna_sample_metadata.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE,
                        colClasses=c(sample_type_code="character"))

# CPM normalise
trf_cpm  <- sweep(trf_mat, 2, colSums(trf_mat), "/") * 1e6
trf_lcpm <- log2(trf_cpm + 1)
mrna_cpm  <- sweep(mrna_mat, 2, colSums(mrna_mat), "/") * 1e6
mrna_lcpm <- log2(mrna_cpm + 1)

trf_meta$group  <- ifelse(trf_meta$sample_type_code == "01", "Tumor", "Normal")
mrna_meta$group <- ifelse(mrna_meta$sample_type_code == "01", "Tumor", "Normal")
rownames(trf_meta)  <- trf_meta$sample_id
rownames(mrna_meta) <- mrna_meta$sample_id

GROUP_COLS <- c(Tumor="#D6604D", Normal="#4393C3")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Study Workflow + Cohort Summary
# ═══════════════════════════════════════════════════════════════════════════════
cat("Building Figure 1...\n")

# Panel A: cohort bar chart (sample counts by type and platform)
cohort_df <- data.frame(
  Platform  = rep(c("tRF (MINTbase)", "mRNA (GDC)"), each=2),
  Group     = rep(c("Tumor","Normal"), 2),
  N         = c(291, 34, 290, 32)
)
cohort_df$Platform <- factor(cohort_df$Platform,
                             levels=c("tRF (MINTbase)","mRNA (GDC)"))
cohort_df$Group    <- factor(cohort_df$Group, levels=c("Tumor","Normal"))

pA <- ggplot(cohort_df, aes(x=Platform, y=N, fill=Group)) +
  geom_col(position="dodge", width=0.6) +
  geom_text(aes(label=N), position=position_dodge(0.6),
            vjust=-0.4, size=4, fontface="bold") +
  scale_fill_manual(values=GROUP_COLS) +
  scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
  labs(title="A  Cohort Sample Counts",
       x=NULL, y="Number of Samples", fill=NULL) +
  theme_bw(base_size=12) +
  theme(legend.position="top",
        plot.title=element_text(face="bold"))

# Panel B: pipeline text table
steps <- data.frame(
  Step  = paste0("Step ", 1:7),
  Phase = c("Data acquisition",
            "Count matrix construction",
            "Differential expression (edgeR)",
            "Spearman anticorrelation",
            "RNA22 binding prediction",
            "ROC / stage / survival",
            "Manuscript"),
  Output = c("MINTbase tRF counts\nGDC STAR mRNA counts",
             "30,076 tRFs × 325\n60,660 genes × 322",
             "1,585 sig. tRFs\n3 sig. genes",
             "104 anticorrelated pairs",
             "18 tRF–gene binding pairs",
             "AUC up to 0.949\nSurvival HR=1.16",
             "This manuscript")
)
pB <- ggplot(steps, aes(y=reorder(Step, -as.integer(Step)))) +
  geom_tile(aes(x=1), fill="#E8F4FD", color="grey70", width=0.9, height=0.85) +
  geom_text(aes(x=0.58, label=Phase), hjust=0, size=3.5, fontface="bold") +
  geom_text(aes(x=1.42, label=Output), hjust=1, size=3, color="grey30") +
  scale_x_continuous(limits=c(0.1, 1.9)) +
  labs(title="B  Analytical Pipeline", x=NULL, y=NULL) +
  theme_void(base_size=11) +
  theme(plot.title=element_text(face="bold", size=12),
        axis.text.y=element_blank())

pdf(file.path(OUT,"fig1_cohort_workflow.pdf"), width=10, height=5)
grid.arrange(pA, pB, ncol=2, widths=c(1,1.2))
dev.off()
cat("  Figure 1 saved.\n")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — tRF Differential Expression (volcano + heatmap + boxplots)
# ═══════════════════════════════════════════════════════════════════════════════
cat("Building Figure 2...\n")

# Panel A: volcano
trf_de_full$sig <- ifelse(trf_de_full$FDR < 0.05 & abs(trf_de_full$logFC) >= 1,
                          ifelse(trf_de_full$logFC > 0, "Up", "Down"), "NS")
trf_de_full$sig <- factor(trf_de_full$sig, levels=c("Up","Down","NS"))
top_lab <- trf_de_full[order(trf_de_full$FDR)[1:8], ]

pA_vol <- ggplot(trf_de_full, aes(x=logFC, y=-log10(FDR), color=sig)) +
  geom_point(size=0.6, alpha=0.5) +
  geom_point(data=top_lab, size=1.5) +
  geom_text_repel(data=top_lab, aes(label=trf_id), size=2.5,
                  max.overlaps=10, segment.size=0.3, show.legend=FALSE) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="grey50") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey50") +
  scale_color_manual(values=c(Up="#D6604D", Down="#4393C3", NS="grey70"),
                     name=NULL) +
  labs(title="A  tRF Differential Expression",
       x="log\u2082 Fold Change (Tumor/Normal)",
       y="-log\u2081\u2080(FDR)") +
  theme_bw(base_size=11) +
  theme(legend.position=c(0.85,0.85),
        plot.title=element_text(face="bold"))

# Panel B: heatmap — top 30 by |logFC|, saved separately then read back
top30 <- trf_de[order(abs(trf_de$logFC), decreasing=TRUE), "trf_id"][1:30]
top30 <- top30[top30 %in% rownames(trf_lcpm)]
mat30 <- trf_lcpm[top30, trf_meta$sample_id]
ann_col <- data.frame(Group=trf_meta$group, row.names=trf_meta$sample_id)
ann_colors <- list(Group=GROUP_COLS)
col_order  <- order(trf_meta$group)

hmap_file <- file.path(OUT, "fig2B_heatmap_tmp.pdf")
pdf(hmap_file, width=7, height=5.5)
pheatmap(mat30[, col_order],
         annotation_col=ann_col,
         annotation_colors=ann_colors,
         show_colnames=FALSE,
         color=colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
         scale="row",
         cluster_cols=FALSE,
         fontsize_row=7,
         main="B  Top 30 tRFs by |log\u2082FC|",
         border_color=NA)
dev.off()

# Panel C: boxplots top 6
top6 <- trf_de[order(trf_de$FDR)[1:6], "trf_id"]
top6 <- top6[top6 %in% rownames(trf_lcpm)]
box_df <- do.call(rbind, lapply(top6, function(t) {
  data.frame(tRF=t,
             expr=as.numeric(trf_lcpm[t, trf_meta$sample_id]),
             Group=trf_meta$group)
}))
box_df$tRF <- factor(box_df$tRF, levels=top6)

pC_box <- ggplot(box_df, aes(x=Group, y=expr, fill=Group)) +
  geom_boxplot(outlier.size=0.4, alpha=0.85) +
  facet_wrap(~tRF, scales="free_y", ncol=3) +
  scale_fill_manual(values=GROUP_COLS) +
  labs(title="C  Top 6 tRFs by FDR",
       x=NULL, y="log\u2082(CPM+1)") +
  theme_bw(base_size=9) +
  theme(legend.position="none",
        strip.text=element_text(size=7),
        plot.title=element_text(face="bold"))

# Save volcano and boxplot panels individually (avoids grid.arrange + ggrepel conflict)
pdf(file.path(OUT,"fig2A_volcano.pdf"), width=7, height=6)
print(pA_vol)
dev.off()

pdf(file.path(OUT,"fig2C_boxplots.pdf"), width=10, height=5)
print(pC_box)
dev.off()
cat("  Figure 2 panels saved.\n")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Candidate Gene Expression (MET, CDKN2A, TERT + volcano)
# ═══════════════════════════════════════════════════════════════════════════════
cat("Building Figure 3...\n")

# Panel A: gene volcano
gene_de_full$sig <- ifelse(gene_de_full$FDR < 0.05 & abs(gene_de_full$logFC) >= 1,
                           ifelse(gene_de_full$logFC > 0, "Up", "Down"), "NS")
pA_gv <- ggplot(gene_de_full, aes(x=logFC, y=-log10(FDR), color=sig)) +
  geom_point(size=2.5) +
  geom_text_repel(aes(label=gene_name), size=3, max.overlaps=20,
                  segment.size=0.3, show.legend=FALSE) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="grey50") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey50") +
  scale_color_manual(values=c(Up="#D6604D", Down="#4393C3", NS="grey60"),
                     name=NULL) +
  labs(title="A  KIRP Candidate Gene DE",
       x="log\u2082 Fold Change", y="-log\u2081\u2080(FDR)") +
  theme_bw(base_size=11) +
  theme(legend.position=c(0.85,0.15),
        plot.title=element_text(face="bold"))

# Panels B-D: boxplots for MET, CDKN2A, TERT
# Need gene_name → row mapping in mrna_lcpm (rownames are gene_id)
# Read one star file to get gene_id->gene_name
star_files <- list.files(
  "/Users/nam.tnguyen2022/KIRP_tRF_project/data_raw/gdc_rnaseq/gdc_download_20260506_091015.911939",
  pattern="*.rna_seq.augmented_star_gene_counts.tsv",
  recursive=TRUE, full.names=TRUE)
gmap <- read.table(star_files[1], sep="\t", header=TRUE,
                   comment.char="#", stringsAsFactors=FALSE)[, c("gene_id","gene_name")]

make_gene_box <- function(gname, panel_label) {
  gid <- gmap$gene_id[gmap$gene_name == gname][1]
  if (is.na(gid) || !gid %in% rownames(mrna_lcpm)) return(NULL)
  df <- data.frame(
    expr  = as.numeric(mrna_lcpm[gid, mrna_meta$sample_id]),
    Group = mrna_meta$group
  )
  lfc  <- round(gene_de_full$logFC[gene_de_full$gene_name == gname], 2)
  fdr  <- signif(gene_de_full$FDR[gene_de_full$gene_name == gname], 2)
  ggplot(df, aes(x=Group, y=expr, fill=Group)) +
    geom_boxplot(outlier.size=0.5, alpha=0.85) +
    scale_fill_manual(values=GROUP_COLS) +
    labs(title=bquote(bold(.(panel_label)) ~ ~ italic(.(gname))),
         subtitle=paste0("log\u2082FC=", lfc, ", FDR=", fdr),
         x=NULL, y="log\u2082(CPM+1)") +
    theme_bw(base_size=11) +
    theme(legend.position="none",
          plot.title=element_text(face="bold"),
          plot.subtitle=element_text(size=9))
}

pB <- make_gene_box("MET",    "B")
pC <- make_gene_box("CDKN2A", "C")
pD <- make_gene_box("TERT",   "D")

pdf(file.path(OUT,"fig3A_gene_volcano.pdf"), width=6, height=5)
print(pA_gv)
dev.off()

plots_bcd <- Filter(Negate(is.null), list(pB, pC, pD))
pdf(file.path(OUT,"fig3BCD_gene_boxes.pdf"), width=9, height=4)
grid.arrange(grobs=plots_bcd, ncol=3)
dev.off()
cat("  Figure 3 saved.\n")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — tRF–MET axis: top scatterplot + RNA22 binding table
# ═══════════════════════════════════════════════════════════════════════════════
cat("Building Figure 4...\n")

# Get MET gene_id
met_id <- gmap$gene_id[gmap$gene_name == "MET"][1]

# Tumor samples matched
tumor_trf  <- trf_meta[trf_meta$group == "Tumor", ]
tumor_mrna <- mrna_meta[mrna_meta$group == "Tumor", ]
tumor_trf$patient  <- substr(tumor_trf$sample_id,  1, 12)
tumor_mrna$patient <- substr(tumor_mrna$sample_id, 1, 12)
shared_pts <- intersect(tumor_trf$patient, tumor_mrna$patient)

trf_pt  <- tumor_trf[tumor_trf$patient   %in% shared_pts, ]
mrna_pt <- tumor_mrna[tumor_mrna$patient %in% shared_pts, ]
trf_pt  <- trf_pt[order(trf_pt$patient), ]
mrna_pt <- mrna_pt[order(mrna_pt$patient), ]

top5_trfs <- anticorr$trf[1:5]
scatter_plots <- lapply(seq_along(top5_trfs), function(i) {
  trf <- top5_trfs[i]
  if (!trf %in% rownames(trf_lcpm) || !met_id %in% rownames(mrna_lcpm)) return(NULL)
  df <- data.frame(
    tRF_expr = as.numeric(trf_lcpm[trf,  trf_pt$sample_id]),
    MET_expr = as.numeric(mrna_lcpm[met_id, mrna_pt$sample_id])
  )
  r   <- anticorr$spearman_r[anticorr$trf == trf][1]
  fdr <- anticorr$FDR[anticorr$trf == trf][1]
  panel_lbl <- LETTERS[i]
  ggplot(df, aes(x=tRF_expr, y=MET_expr)) +
    geom_point(size=0.8, alpha=0.5, color="#4393C3") +
    geom_smooth(method="lm", se=TRUE, color="#D6604D", linewidth=0.8) +
    annotate("text", x=Inf, y=Inf, hjust=1.1, vjust=1.5, size=3,
             label=paste0("r=", round(r,3), "\nFDR=", signif(fdr,2))) +
    labs(title=paste0(panel_lbl, "  ", trf),
         x="tRF log\u2082(CPM+1)", y="MET log\u2082(CPM+1)") +
    theme_bw(base_size=9) +
    theme(plot.title=element_text(face="bold", size=8))
})
scatter_plots <- Filter(Negate(is.null), scatter_plots)

# Panel F: RNA22 top binding table as ggplot grob
rna22_sum <- read.table(file.path(TABS,"rna22_binding_summary.tsv"),
                        sep="\t", header=TRUE, stringsAsFactors=FALSE)
rna22_met <- rna22_sum[rna22_sum$gene == "MET",
                        c("tRF_id","target_position","folding_energy_kcal","p_value")]
colnames(rna22_met) <- c("tRF","Position","ΔG (kcal/mol)","p-value")
rna22_met$`p-value` <- signif(rna22_met$`p-value`, 3)
rna22_met$`ΔG (kcal/mol)` <- round(rna22_met$`ΔG (kcal/mol)`, 1)
rna22_met <- head(rna22_met, 8)

tbl_grob <- tableGrob(rna22_met, rows=NULL,
                      theme=ttheme_minimal(base_size=8,
                        core=list(fg_params=list(hjust=0, x=0.05)),
                        colhead=list(fg_params=list(fontface="bold"))))

pF_title <- textGrob("F  RNA22 Predicted MET Binding Sites (p < 0.05)",
                     gp=gpar(fontface="bold", fontsize=10), x=0.05, hjust=0)
pF <- arrangeGrob(pF_title, tbl_grob, ncol=1, heights=c(0.12, 0.88))

pdf(file.path(OUT,"fig4_trf_met_axis.pdf"), width=14, height=9)
grid.arrange(
  scatter_plots[[1]], scatter_plots[[2]], scatter_plots[[3]],
  scatter_plots[[4]], scatter_plots[[5]], pF,
  ncol=3
)
dev.off()
cat("  Figure 4 saved.\n")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 5 — ROC curves (top 5 tRFs)
# ═══════════════════════════════════════════════════════════════════════════════
cat("Building Figure 5...\n")

response <- as.numeric(trf_meta$group == "Tumor")
top5_roc <- roc_df$tRF[1:5]
pal5 <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00")

# Panel A: overlaid ROC curves
pdf(file.path(OUT,"fig5_roc.pdf"), width=11, height=5.5)
par(mfrow=c(1,2), mar=c(4,4,3,1))

# Left: overlaid curves
plot(0:1, 0:1, type="l", lty=2, col="grey60",
     xlab="False Positive Rate", ylab="True Positive Rate",
     main="A  ROC Curves — Top 5 tRFs", cex.main=1.1, font.main=2)
for (i in seq_along(top5_roc)) {
  trf  <- top5_roc[i]
  expr <- as.numeric(trf_lcpm[trf, ])
  r    <- roc(response, expr, quiet=TRUE)
  lines(1-r$specificities, r$sensitivities, col=pal5[i], lwd=2)
}
legend("bottomright",
       legend=paste0(top5_roc, "  AUC=",
                     round(roc_df$AUC[match(top5_roc, roc_df$tRF)], 3)),
       col=pal5, lwd=2, cex=0.65, bty="n")

# Right: AUC barplot
roc_ord <- roc_df[order(roc_df$AUC), ]
barplot(roc_ord$AUC, names.arg=roc_ord$tRF,
        horiz=TRUE, las=1, cex.names=0.6,
        col=colorRampPalette(c("#fdd49e","#d7301f"))(nrow(roc_ord)),
        xlab="AUC", main="B  AUC — Top 20 tRFs",
        cex.main=1.1, font.main=2,
        xlim=c(0.4, 1.0))
abline(v=0.5, lty=2, col="grey40")

dev.off()
cat("  Figure 5 saved.\n")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 6 — Stage + Survival
# ═══════════════════════════════════════════════════════════════════════════════
cat("Building Figure 6...\n")

# Clinical data
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
primary_stage$stage_simple <- sub(" [ABC]$","", primary_stage$stage)
clin$stage <- NULL
clin <- merge(clin, primary_stage[, c("patient_id","stage","stage_simple")],
              by="patient_id", all.x=TRUE, sort=FALSE)
clean_val <- function(x) { x[x %in% c("'--","--","")] <- NA; as.numeric(x) }
clin$days_death <- clean_val(clin$days_death)
clin$days_fu    <- clean_val(clin$days_fu)
clin$event      <- as.integer(clin$vital_status == "Dead")
clin$os_days    <- ifelse(!is.na(clin$days_death), clin$days_death, clin$days_fu)

tumor_trf2 <- trf_meta[trf_meta$group == "Tumor", ]
tumor_trf2$patient_id <- substr(tumor_trf2$sample_id, 1, 12)

# Panel A: stage boxplot for top 2 sig tRFs
stage_trfs <- kw_df$tRF[kw_df$kruskal_fdr < 0.05][1:2]
stage_clin <- merge(tumor_trf2, clin[, c("patient_id","stage_simple")],
                    by="patient_id")
stage_clin <- stage_clin[stage_clin$stage_simple %in%
                           c("Stage I","Stage II","Stage III","Stage IV"), ]
stage_clin$stage_simple <- factor(stage_clin$stage_simple,
                                   levels=c("Stage I","Stage II","Stage III","Stage IV"))
stage_long <- do.call(rbind, lapply(stage_trfs, function(trf) {
  data.frame(tRF=trf, stage=stage_clin$stage_simple,
             expr=as.numeric(trf_lcpm[trf, stage_clin$sample_id]))
}))
fdr_labs <- setNames(paste0("KW FDR=",
                             signif(kw_df$kruskal_fdr[kw_df$tRF %in% stage_trfs], 2)),
                     stage_trfs)

pA_stage <- ggplot(stage_long, aes(x=stage, y=expr, fill=stage)) +
  geom_boxplot(outlier.size=0.5, alpha=0.85) +
  facet_wrap(~tRF, scales="free_y", ncol=1,
             labeller=labeller(tRF=function(x) paste0(x, "\n", fdr_labs[x]))) +
  scale_fill_manual(values=c("Stage I"="#2166AC","Stage II"="#92C5DE",
                              "Stage III"="#F4A582","Stage IV"="#D6604D")) +
  labs(title="A  Stage Association",
       x="AJCC Pathologic Stage", y="log\u2082(CPM+1)") +
  theme_bw(base_size=10) +
  theme(axis.text.x=element_text(angle=30, hjust=1),
        legend.position="none",
        plot.title=element_text(face="bold"))

# Panel B: KM for tRF-30-81PV6RRNLNK8
surv_trf <- "tRF-30-81PV6RRNLNK8"
surv_data <- merge(tumor_trf2, clin[, c("patient_id","os_days","event")],
                   by="patient_id")
surv_data <- surv_data[!is.na(surv_data$os_days) & surv_data$os_days > 0, ]

if (surv_trf %in% rownames(trf_lcpm)) {
  expr_s <- as.numeric(trf_lcpm[surv_trf, surv_data$sample_id])
  surv_data$km_group <- ifelse(expr_s >= median(expr_s, na.rm=TRUE), "High","Low")
  sf <- survfit(Surv(os_days, event) ~ km_group, data=surv_data)
  km_plot <- ggsurvplot(sf, data=surv_data,
                        title="B  Overall Survival\ntRF-30-81PV6RRNLNK8",
                        palette=c("#D6604D","#4393C3"),
                        pval=TRUE, conf.int=FALSE,
                        risk.table=TRUE, risk.table.height=0.28,
                        xlab="Days", ylab="Survival Probability",
                        legend.title="Expression",
                        font.title=c(10,"bold"),
                        font.x=9, font.y=9, font.tickslab=8)
}

pdf(file.path(OUT,"fig6_stage_survival.pdf"), width=11, height=7)
if (exists("km_plot")) {
  print(km_plot, newpage=FALSE)
  vp <- viewport(x=0.02, y=0, width=0.35, height=1, just=c("left","bottom"))
  pushViewport(vp)
  print(pA_stage, newpage=FALSE)
  popViewport()
} else {
  grid.arrange(pA_stage, ncol=1)
}
dev.off()
cat("  Figure 6 saved.\n")

cat("\n✓ All 6 main figures saved to:", OUT, "\n")
