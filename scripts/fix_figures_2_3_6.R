## Fix figures 2, 3, and 6 — robust approach using temp PNG panels

suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(pheatmap); library(patchwork)
  library(cowplot); library(scales); library(dplyr); library(tidyr)
  library(grid); library(gridExtra); library(RColorBrewer)
  library(pROC); library(survival); library(survminer)
})

DATA <- "/Users/nam.tnguyen2022/KIRP_tRF_project/data_processed"
TABS <- "/Users/nam.tnguyen2022/KIRP_tRF_project/results/tables"
RAW  <- "/Users/nam.tnguyen2022/KIRP_tRF_project/data_raw"
OUT  <- "/Users/nam.tnguyen2022/KIRP_tRF_project/results/figures"
TMP  <- tempdir()

pub_theme <- theme_classic(base_size=11, base_family="Helvetica") +
  theme(plot.title=element_text(face="bold",size=12,hjust=0),
        axis.title=element_text(size=10), axis.text=element_text(size=9,color="black"),
        axis.line=element_line(color="black",linewidth=0.4),
        axis.ticks=element_line(color="black",linewidth=0.4),
        legend.text=element_text(size=9), legend.title=element_text(size=9,face="bold"),
        legend.key.size=unit(0.4,"cm"), strip.background=element_rect(fill="#F0F0F0",color=NA),
        strip.text=element_text(face="bold",size=9),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        plot.margin=margin(8,8,8,8))
theme_set(pub_theme)
GROUP_COLS <- c(Tumor="#C0392B", Normal="#2980B9")
STAGE_COLS <- c("Stage I"="#2166AC","Stage II"="#74ADD1","Stage III"="#F46D43","Stage IV"="#A50026")

cat("Loading data...\n")
trf_mat   <- read.table(file.path(DATA,"trf_count_matrix.tsv"), sep="\t",header=TRUE,row.names=1,check.names=FALSE)
trf_meta  <- read.table(file.path(DATA,"trf_sample_metadata.tsv"),sep="\t",header=TRUE,stringsAsFactors=FALSE,colClasses=c(sample_type_code="character"))
trf_de    <- read.table(file.path(TABS,"trf_DE_significant.tsv"),sep="\t",header=TRUE,stringsAsFactors=FALSE)
trf_full  <- read.table(file.path(TABS,"trf_DE_full.tsv"),sep="\t",header=TRUE,stringsAsFactors=FALSE)
gene_full <- read.table(file.path(TABS,"gene_DE_full.tsv"),sep="\t",header=TRUE,stringsAsFactors=FALSE)
kw_tab    <- read.table(file.path(TABS,"stage_kruskal_results.tsv"),sep="\t",header=TRUE,stringsAsFactors=FALSE)
cox_tab   <- read.table(file.path(TABS,"cox_univariate_results.tsv"),sep="\t",header=TRUE,stringsAsFactors=FALSE)
mrna_mat  <- read.table(file.path(DATA,"mrna_count_matrix.tsv"),sep="\t",header=TRUE,row.names=1,check.names=FALSE)
mrna_meta <- read.table(file.path(DATA,"mrna_sample_metadata.tsv"),sep="\t",header=TRUE,stringsAsFactors=FALSE,colClasses=c(sample_type_code="character"))

trf_lcpm  <- log2(sweep(trf_mat,2,colSums(trf_mat),"/")*1e6+1)
mrna_lcpm <- log2(sweep(mrna_mat,2,colSums(mrna_mat),"/")*1e6+1)
trf_meta$group  <- ifelse(trf_meta$sample_type_code=="01","Tumor","Normal")
mrna_meta$group <- ifelse(mrna_meta$sample_type_code=="01","Tumor","Normal")
rownames(trf_meta) <- trf_meta$sample_id
rownames(mrna_meta) <- mrna_meta$sample_id

star_files <- list.files(
  file.path(RAW,"gdc_rnaseq/gdc_download_20260506_091015.911939"),
  pattern="*.rna_seq.augmented_star_gene_counts.tsv",recursive=TRUE,full.names=TRUE)
gmap <- read.table(star_files[1],sep="\t",header=TRUE,comment.char="#")[,c("gene_id","gene_name")]

clin_raw <- read.table(file.path(RAW,"gdc_clinical/clinical.tsv"),sep="\t",header=TRUE,stringsAsFactors=FALSE,fill=TRUE,quote="")
clin <- clin_raw[,c("cases.submitter_id","diagnoses.ajcc_pathologic_stage","demographic.vital_status","demographic.days_to_death","diagnoses.days_to_last_follow_up")]
colnames(clin) <- c("patient_id","stage","vital_status","days_death","days_fu")
clin <- clin[!duplicated(clin$patient_id),]
primary_stage <- clin_raw[
  as.character(clin_raw$diagnoses.diagnosis_is_primary_disease) %in% c("TRUE", "True", "true"),
  c("cases.submitter_id", "diagnoses.submitter_id",
    "diagnoses.ajcc_pathologic_stage")
]
colnames(primary_stage) <- c("patient_id", "diagnosis_id", "stage")
primary_stage <- primary_stage[order(primary_stage$patient_id,
                                     primary_stage$diagnosis_id),]
primary_stage <- primary_stage[!duplicated(primary_stage$patient_id),]
primary_stage$stage[primary_stage$stage=="'--"] <- NA
primary_stage$stage_simple <- sub(" [ABC]$","",primary_stage$stage)
clin$stage <- NULL
clin <- merge(clin, primary_stage[, c("patient_id","stage","stage_simple")],
              by="patient_id", all.x=TRUE, sort=FALSE)
cn <- function(x){x[x %in% c("'--","--","")] <- NA; as.numeric(x)}
clin$days_death <- cn(clin$days_death); clin$days_fu <- cn(clin$days_fu)
clin$event   <- as.integer(clin$vital_status=="Dead")
clin$os_days <- ifelse(!is.na(clin$days_death),clin$days_death,clin$days_fu)

tumor_trf  <- trf_meta[trf_meta$group=="Tumor",]
tumor_mrna <- mrna_meta[mrna_meta$group=="Tumor",]
tumor_trf$patient_id  <- substr(tumor_trf$sample_id,1,12)
tumor_mrna$patient_id <- substr(tumor_mrna$sample_id,1,12)
cat("Done.\n\n")

# ─── Helper: save panel as PNG then return as cowplot image ─────────────────
panel_png <- function(p, file, width, height, dpi=300) {
  ggsave(file, p, width=width, height=height, dpi=dpi, bg="white")
  ggdraw() + draw_image(file, scale=1)
}

# ════════════════════════════════════════════════════════════════════════════
# FIG 2 — tRF Differential Expression (3-panel composite)
# ════════════════════════════════════════════════════════════════════════════
cat("Figure 2...\n")

trf_full$sig <- case_when(
  trf_full$FDR<0.05 & trf_full$logFC>1  ~ "Up",
  trf_full$FDR<0.05 & trf_full$logFC<(-1) ~ "Down",
  TRUE ~ "NS")
trf_full$sig <- factor(trf_full$sig, levels=c("Up","Down","NS"))
n_up <- sum(trf_full$sig=="Up"); n_down <- sum(trf_full$sig=="Down")
top_lab <- trf_full %>% filter(sig!="NS") %>% arrange(FDR) %>% slice_head(n=10)
ylim_max <- max(-log10(trf_full$FDR[is.finite(trf_full$FDR) & trf_full$FDR>0]),na.rm=TRUE)

p2A <- ggplot(trf_full, aes(x=logFC, y=-log10(FDR), color=sig)) +
  geom_point(data=~filter(.x,sig=="NS"),size=0.4,alpha=0.25,color="grey75") +
  geom_point(data=~filter(.x,sig!="NS"),size=0.7,alpha=0.7) +
  geom_text_repel(data=top_lab,aes(label=trf_id),size=2.3,segment.size=0.3,
    segment.color="grey50",min.segment.length=0.2,max.overlaps=12,
    show.legend=FALSE,fontface="italic") +
  geom_vline(xintercept=c(-1,1),linetype="dashed",color="grey40",linewidth=0.4) +
  geom_hline(yintercept=-log10(0.05),linetype="dashed",color="grey40",linewidth=0.4) +
  scale_color_manual(values=c(Up="#C0392B",Down="#2980B9",NS="grey75"),
    labels=c(paste0("Up (",n_up,")"),paste0("Down (",n_down,")"),"NS"),name=NULL) +
  scale_x_continuous(limits=c(-6,8)) +
  annotate("text",x=-5.8,y=ylim_max*0.96,label=paste0(n_down," down"),
    color="#2980B9",size=3.2,hjust=0,fontface="bold") +
  annotate("text",x=6,y=ylim_max*0.96,label=paste0(n_up," up"),
    color="#C0392B",size=3.2,hjust=1,fontface="bold") +
  labs(x=expression(log[2]~"Fold Change (Tumor/Normal)"),
       y=expression(-log[10]~"(FDR)")) +
  theme(legend.position="none")

# Heatmap panel: save to PNG, load back
top30_ids <- trf_de %>% arrange(desc(abs(logFC))) %>% slice_head(n=30) %>% pull(trf_id)
top30_ids <- top30_ids[top30_ids %in% rownames(trf_lcpm)]
col_order  <- c(which(trf_meta$group=="Normal"), which(trf_meta$group=="Tumor"))
mat30 <- trf_lcpm[top30_ids, trf_meta$sample_id[col_order]]
ann_col <- data.frame(Group=trf_meta$group[col_order], row.names=colnames(mat30))
hm_png <- file.path(TMP,"fig2B_heatmap.png")
pheatmap(mat30, annotation_col=ann_col,
  annotation_colors=list(Group=GROUP_COLS),
  show_colnames=FALSE,
  color=colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
  scale="row", cluster_cols=FALSE, fontsize_row=7, fontsize=8,
  border_color=NA, annotation_legend=TRUE, main="Top 30 DE tRFs (by |FC|)",
  filename=hm_png, width=4, height=6)
p2B_img <- ggdraw() + draw_image(hm_png, scale=0.97)

# Boxplots
top9 <- trf_de %>% arrange(FDR) %>% slice_head(n=9) %>% pull(trf_id)
top9 <- top9[top9 %in% rownames(trf_lcpm)]
box_df <- bind_rows(lapply(top9, function(t)
  data.frame(tRF=t, expr=as.numeric(trf_lcpm[t, trf_meta$sample_id]),
             Group=trf_meta$group)))
box_df$tRF_short <- factor(sub("tRF-","",box_df$tRF), levels=sub("tRF-","",top9))

p2C <- ggplot(box_df, aes(x=Group, y=expr, fill=Group)) +
  geom_boxplot(outlier.size=0.3, outlier.alpha=0.4, linewidth=0.4, width=0.65) +
  stat_summary(fun=median,geom="point",shape=18,size=1.5,color="white",show.legend=FALSE) +
  facet_wrap(~tRF_short, nrow=3, scales="free_y") +
  scale_fill_manual(values=GROUP_COLS, name=NULL) +
  labs(x=NULL, y=expression(log[2]~"(CPM+1)"),
       title="Top 9 DE tRFs (by FDR)") +
  theme(axis.text.x=element_text(angle=30,hjust=1,size=8), legend.position="none",
        strip.text=element_text(size=7))

# Save Fig 2 panels as PNG then assemble
p2A_f <- file.path(TMP,"p2A.png")
p2C_f <- file.path(TMP,"p2C.png")
ggsave(p2A_f, p2A, width=5.5, height=5, dpi=250, bg="white")
ggsave(p2C_f, p2C, width=5.5, height=5, dpi=250, bg="white")
p2A_i <- ggdraw() + draw_image(p2A_f)
p2C_i <- ggdraw() + draw_image(p2C_f)

top_row   <- plot_grid(p2A_i, p2B_img, ncol=2, labels=c("A","B"),
                       label_size=14, label_fontface="bold")
bot_row   <- plot_grid(p2C_i, ncol=1, labels="C",
                       label_size=14, label_fontface="bold")
fig2_body <- plot_grid(top_row, bot_row, nrow=2, rel_heights=c(1.1,1))
title2    <- ggdraw() + draw_label(
  "Figure 2. Differential expression of tRFs in KIRP tumor vs. normal tissue.",
  fontface="bold", size=11, x=0.01, hjust=0)
fig2 <- plot_grid(title2, fig2_body, nrow=2, rel_heights=c(0.06, 1))
ggsave(file.path(OUT,"fig2_trf_DE.pdf"), fig2, width=12, height=11)
cat("  Done.\n")

# ════════════════════════════════════════════════════════════════════════════
# FIG 3 — Candidate Gene Expression
# ════════════════════════════════════════════════════════════════════════════
cat("Figure 3...\n")

cat("  gene_full dims:", nrow(gene_full), "cols:", paste(colnames(gene_full),collapse=","), "\n")
cat("  gene_full FDR range:", range(gene_full$FDR, na.rm=TRUE), "\n")
cat("  sig genes:", sum(gene_full$FDR<0.05, na.rm=TRUE), "\n")

gene_full$sig <- case_when(
  gene_full$FDR<0.05 & gene_full$logFC>1  ~ "Up",
  gene_full$FDR<0.05 & gene_full$logFC<(-1) ~ "Down",
  gene_full$FDR<0.05  ~ "Sig (|FC|<2)",
  TRUE ~ "NS")
gene_full$sig <- factor(gene_full$sig, levels=c("Up","Down","Sig (|FC|<2)","NS"))

# Show all genes (small set — don't hide labels)
p3A <- ggplot(gene_full, aes(x=logFC, y=-log10(FDR), color=sig)) +
  geom_point(data=~filter(.x,sig=="NS"),size=1.2,alpha=0.3,color="grey75") +
  geom_point(data=~filter(.x,sig!="NS"),size=3,alpha=0.85) +
  geom_text_repel(data=~filter(.x,sig!="NS"),
    aes(label=gene_name),size=3,segment.size=0.3,show.legend=FALSE,
    max.overlaps=30,fontface="italic",box.padding=0.5,min.segment.length=0.1) +
  geom_vline(xintercept=c(-1,1),linetype="dashed",color="grey40",linewidth=0.4) +
  geom_hline(yintercept=-log10(0.05),linetype="dashed",color="grey40",linewidth=0.4) +
  scale_color_manual(values=c(Up="#C0392B",Down="#2980B9","Sig (|FC|<2)"="#8E44AD",NS="grey75"),name=NULL) +
  labs(x=expression(log[2]~"Fold Change"),y=expression(-log[10]~"(FDR)"),
       title="Gene DE Volcano") +
  theme(legend.position=c(0.78,0.88),
        legend.background=element_rect(fill="white",color="grey80",linewidth=0.3))

make_gene_box <- function(gname) {
  rows <- which(gmap$gene_name==gname)
  cat("    ",gname,"- matching rows:", length(rows), "\n")
  if (length(rows)==0) return(NULL)
  gid <- gmap$gene_id[rows[1]]
  cat("    ",gname,"- gene_id:", gid, "in mrna_lcpm:", gid %in% rownames(mrna_lcpm), "\n")
  if (!gid %in% rownames(mrna_lcpm)) return(NULL)
  df_g <- gene_full[gene_full$gene_name==gname,]
  lfc  <- if(nrow(df_g)>0) round(df_g$logFC[1],2) else "?"
  fdr  <- if(nrow(df_g)>0) formatC(df_g$FDR[1],format="e",digits=1) else "?"
  df   <- data.frame(
    expr  = as.numeric(mrna_lcpm[gid, mrna_meta$sample_id]),
    Group = mrna_meta$group)
  ggplot(df, aes(x=Group, y=expr, fill=Group)) +
    geom_boxplot(outlier.size=0.4, outlier.alpha=0.5, linewidth=0.4, width=0.6) +
    stat_summary(fun=median,geom="point",shape=18,size=2,color="white",show.legend=FALSE) +
    scale_fill_manual(values=GROUP_COLS) +
    annotate("text",x=1.5,y=quantile(df$expr,0.98),
      label=paste0("FC=",lfc,"\nFDR=",fdr),
      size=2.8,hjust=0.5,vjust=1,color="grey20") +
    labs(x=NULL,y=expression(log[2]~"(CPM+1)"),
         title=bquote(italic(.(gname)))) +
    theme(legend.position="none",
          plot.title=element_text(face="bold.italic",hjust=0.5,size=11))
}

p3B <- make_gene_box("MET")
p3C <- make_gene_box("CDKN2A")
p3D <- make_gene_box("TERT")

# Fallback: if gene box returns NULL, show placeholder
placeholder <- ggplot() + theme_void() +
  annotate("text",x=0.5,y=0.5,label="No data",size=5,color="grey50")
if (is.null(p3B)) p3B <- placeholder
if (is.null(p3C)) p3C <- placeholder
if (is.null(p3D)) p3D <- placeholder

p3A_f <- file.path(TMP,"p3A.png")
p3B_f <- file.path(TMP,"p3B.png")
p3C_f <- file.path(TMP,"p3C.png")
p3D_f <- file.path(TMP,"p3D.png")
ggsave(p3A_f, p3A, width=6, height=5, dpi=250, bg="white")
ggsave(p3B_f, p3B, width=3.2, height=3.5, dpi=250, bg="white")
ggsave(p3C_f, p3C, width=3.2, height=3.5, dpi=250, bg="white")
ggsave(p3D_f, p3D, width=3.2, height=3.5, dpi=250, bg="white")
p3A_i <- ggdraw()+draw_image(p3A_f)
p3BCD  <- plot_grid(ggdraw()+draw_image(p3B_f),
                    ggdraw()+draw_image(p3C_f),
                    ggdraw()+draw_image(p3D_f),
                    ncol=3, labels=c("B","C","D"),
                    label_size=14,label_fontface="bold")
top3   <- plot_grid(p3A_i, ncol=1, labels="A", label_size=14, label_fontface="bold")
fig3_body <- plot_grid(top3, p3BCD, nrow=2, rel_heights=c(1.3,1))
title3    <- ggdraw() + draw_label(
  "Figure 3. Differential expression of KIRP candidate genes.",
  fontface="bold",size=11,x=0.01,hjust=0)
fig3 <- plot_grid(title3, fig3_body, nrow=2, rel_heights=c(0.06,1))
ggsave(file.path(OUT,"fig3_gene_DE.pdf"), fig3, width=10, height=9)
cat("  Done.\n")

# ════════════════════════════════════════════════════════════════════════════
# FIG 5 — Stage + Survival
# ════════════════════════════════════════════════════════════════════════════
cat("Figure 5...\n")

stage_trfs <- c(
  "tRF-22-7SIRMM121",
  "tRF-22-7SERML921",
  "tRF-23-7SIR3DR2DV",
  "tRF-23-7SIRMM12V",
  "tRF-25-RNLNKSEK51")
stage_trfs <- stage_trfs[stage_trfs %in% kw_tab$tRF]

stage_clin <- merge(tumor_trf, clin[,c("patient_id","stage_simple")], by="patient_id")
stage_clin <- stage_clin[stage_clin$stage_simple %in% c("Stage I","Stage II","Stage III","Stage IV"),]
stage_clin$stage_simple <- factor(stage_clin$stage_simple,
  levels=c("Stage I","Stage II","Stage III","Stage IV"))

n_stage <- stage_clin %>% count(stage_simple)
x_lab <- setNames(paste0(levels(stage_clin$stage_simple),"\n(n=",n_stage$n,")"),
                  levels(stage_clin$stage_simple))
fmt_p <- function(x) {
  if (is.na(x)) return("NA")
  if (x < 0.001) {
    exponent <- floor(log10(x))
    mantissa <- x / 10^exponent
    return(paste0(formatC(mantissa, format="f", digits=1), " × 10^", exponent))
  }
  sprintf("%.3f", x)
}

stage_long <- bind_rows(lapply(stage_trfs, function(trf) {
  fdr_val <- kw_tab$kruskal_fdr[kw_tab$tRF==trf]
  data.frame(tRF_label=paste0(sub("tRF-","",trf),"\nFDR = ",fmt_p(fdr_val)),
             stage=stage_clin$stage_simple,
             expr=as.numeric(trf_lcpm[trf, stage_clin$sample_id]))
}))
stage_long$tRF_label <- factor(stage_long$tRF_label, levels=unique(stage_long$tRF_label))

p6A <- ggplot(stage_long, aes(x=stage, y=expr, fill=stage, color=stage)) +
  geom_boxplot(outlier.shape=NA, linewidth=0.4, width=0.58, alpha=0.32, color="black") +
  geom_jitter(width=0.13, height=0, size=0.35, alpha=0.34, shape=16) +
  scale_fill_manual(values=STAGE_COLS, name=NULL) +
  scale_color_manual(values=STAGE_COLS, guide="none") +
  scale_x_discrete(labels=x_lab) +
  facet_wrap(~tRF_label, ncol=5, scales="free_y") +
  labs(x=NULL, y=expression(log[2]~"(CPM + 1)")) +
  theme(axis.text.x=element_text(size=6.2, lineheight=0.9),
        axis.text.y=element_text(size=6.2),
        strip.text=element_text(size=6.4, lineheight=0.95),
        panel.spacing.x=unit(0.48, "lines"),
        legend.position="none")

# Survival: KM curves for top Cox tRF
surv_trf <- cox_tab %>% arrange(FDR) %>% slice_head(n=1) %>% pull(tRF)
cat("  Survival tRF:", surv_trf, "\n")

surv_data <- merge(tumor_trf, clin[,c("patient_id","os_days","event")], by="patient_id")
surv_data <- surv_data[!is.na(surv_data$os_days) & surv_data$os_days>0,]
cat("  Survival data rows:", nrow(surv_data), "\n")

km_png_path <- file.path(TMP,"km_plot.png")
if (surv_trf %in% rownames(trf_lcpm) && nrow(surv_data)>10) {
  # Use only samples whose IDs are in trf_lcpm
  matched <- surv_data$sample_id[surv_data$sample_id %in% colnames(trf_lcpm)]
  surv_s  <- surv_data[surv_data$sample_id %in% matched, ]
  expr_s  <- as.numeric(trf_lcpm[surv_trf, surv_s$sample_id])
  med     <- median(expr_s, na.rm=TRUE)
  # If median == 0, split at detected (>0) vs not detected (=0).
  if (med == 0) {
    surv_s$km_group <- factor(ifelse(expr_s > 0, "Detected", "Not detected"),
                              levels=c("Detected","Not detected"))
    split_label <- paste0(sub("tRF-","",surv_trf), " expression")
  } else {
    surv_s$km_group <- factor(ifelse(expr_s >= med, "High", "Low"),
                              levels=c("High","Low"))
    split_label <- paste0(sub("tRF-","",surv_trf), " expression")
  }
  cat("  KM split -", levels(surv_s$km_group)[1], ":",
      sum(surv_s$km_group==levels(surv_s$km_group)[1]),
      levels(surv_s$km_group)[2], ":",
      sum(surv_s$km_group==levels(surv_s$km_group)[2]), "\n")
  sf  <- survfit(Surv(os_days,event) ~ km_group, data=surv_s)
  cox_fit <- coxph(Surv(os_days,event) ~ expr_s, data=surv_s)
  hr  <- round(exp(coef(cox_fit)),2)
  pv  <- round(summary(cox_fit)$coefficients[,"Pr(>|z|)"],3)
  km_palette <- c("#C0392B","#2980B9")

  km_out <- ggsurvplot(sf, data=surv_s,
    palette=km_palette, pval=TRUE, pval.size=3.5,
    conf.int=TRUE, conf.int.alpha=0.1,
    risk.table=TRUE, risk.table.height=0.24, risk.table.fontsize=2.25,
    xlab="Time (days)", ylab="Survival Probability",
    legend.title=split_label,
	    title="tRF-30-81PV6RRNLNK8 expression and overall survival",
	    font.title=c(9,"bold"), font.x=9, font.y=9, font.tickslab=8,
	    ggtheme=pub_theme)

  km_out$table <- km_out$table +
    labs(x=NULL, y=NULL) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.margin=margin(0, 6, 0, 6))

  # Save KM plot using survminer helper
  km_combined <- cowplot::plot_grid(km_out$plot, km_out$table,
    nrow=2, rel_heights=c(1, 0.30), align="v", axis="lr")
  ggsave(km_png_path, km_combined, width=4.7, height=4.4, dpi=250, bg="white")
  p6B_i <- ggdraw() + draw_image(km_png_path, scale=0.97)
} else {
  p6B_i <- ggdraw() + draw_label("Survival data unavailable", size=12, color="grey50")
}

p6A_f <- file.path(TMP,"p6A.png")
ggsave(p6A_f, p6A, width=11, height=3.6, dpi=250, bg="white")
p6A_i <- ggdraw() + draw_image(p6A_f)

p6B_small <- plot_grid(NULL, p6B_i, NULL, ncol=3,
  rel_widths=c(0.10, 0.80, 0.10))
fig6_body <- plot_grid(p6A_i, p6B_small, nrow=2,
  labels=c("a","b"), label_size=10, label_fontface="bold",
  rel_heights=c(0.94, 0.86))
title6 <- ggdraw() + draw_label(
  "Figure 5. Clinical associations of candidate tRF expression in TCGA-KIRP.",
  fontface="bold", size=11, x=0.01, hjust=0)
fig6 <- plot_grid(title6, fig6_body, nrow=2, rel_heights=c(0.05,1))
ggsave(file.path(OUT,"fig6_stage_survival.pdf"), fig6, width=12, height=8.2)
ggsave(file.path(OUT,"fig6_stage_survival.png"), fig6, width=12, height=8.2, dpi=600, bg="white")
cat("  Done.\n")

cat("\n=====================================\n")
cat("Fixed figures:\n")
for (f in c("fig2_trf_DE.pdf","fig3_gene_DE.pdf","fig6_stage_survival.pdf")) {
  fp <- file.path(OUT,f)
  cat(sprintf("  %s  %s  %s KB\n", if(file.exists(fp)) "[OK]" else "[MISSING]",
              f, if(file.exists(fp)) round(file.info(fp)$size/1024) else "N/A"))
}
cat("=====================================\n")
