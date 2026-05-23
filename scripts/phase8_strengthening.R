## Phase 8: Strengthening Analysis
## 1. ROC/AUC — top tRFs distinguishing tumor vs normal
## 2. Tumor stage association (Kruskal-Wallis + boxplots)
## 3. Kaplan-Meier survival analysis

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(pROC)
  library(survival)
  library(survminer)
})

# Edit scripts/config.R to set ROOT for your machine, or change the fallback below:
if (file.exists("scripts/config.R")) source("scripts/config.R") else if (file.exists("config.R")) source("config.R") else ROOT <- "/Users/nam.tnguyen2022/KIRP_tRF_project"
DATA  <- file.path(ROOT, "data_processed")
RAW   <- file.path(ROOT, "data_raw")
FIGS  <- file.path(ROOT, "results", "figures")
TABS  <- file.path(ROOT, "results", "tables")

# ── Load data ──────────────────────────────────────────────────────────────────
cat("Loading data...\n")

trf_mat  <- read.table(file.path(DATA, "trf_count_matrix.tsv"),
                       sep="\t", header=TRUE, row.names=1,
                       check.names=FALSE)
trf_meta <- read.table(file.path(DATA, "trf_sample_metadata.tsv"),
                       sep="\t", header=TRUE, stringsAsFactors=FALSE,
                       colClasses=c(sample_type_code="character"))
trf_de   <- read.table(file.path(TABS, "trf_DE_significant.tsv"),
                       sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Top 20 tRFs by absolute logFC (from significant set)
top_trfs <- trf_de[order(abs(trf_de$logFC), decreasing=TRUE), "trf_id"][1:20]

# CPM-normalise
library_sizes <- colSums(trf_mat)
cpm_mat <- sweep(trf_mat, 2, library_sizes, "/") * 1e6
log_cpm  <- log2(cpm_mat + 1)

# Align metadata to matrix
trf_meta <- trf_meta[trf_meta$sample_id %in% colnames(log_cpm), ]
rownames(trf_meta) <- trf_meta$sample_id
trf_meta <- trf_meta[colnames(log_cpm), ]
trf_meta$group <- ifelse(trf_meta$sample_type_code == "01", "Tumor", "Normal")

# ── 1. ROC / AUC ───────────────────────────────────────────────────────────────
cat("Running ROC analysis...\n")

response <- as.numeric(trf_meta$group == "Tumor")  # 1=Tumor, 0=Normal

roc_res <- lapply(top_trfs, function(trf) {
  if (!trf %in% rownames(log_cpm)) return(NULL)
  expr <- as.numeric(log_cpm[trf, ])
  r <- roc(response, expr, quiet=TRUE)
  data.frame(tRF=trf, AUC=as.numeric(auc(r)),
             stringsAsFactors=FALSE)
})
roc_df <- do.call(rbind, Filter(Negate(is.null), roc_res))
roc_df  <- roc_df[order(roc_df$AUC, decreasing=TRUE), ]

# Save AUC table
write.table(roc_df, file.path(TABS, "roc_auc_top_trfs.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("AUC table: %d tRFs → %s\n", nrow(roc_df),
            file.path(TABS, "roc_auc_top_trfs.tsv")))

# Plot ROC curves for top 5
pdf(file.path(FIGS, "roc_curves_top5.pdf"), width=6, height=6)
top5 <- roc_df$tRF[1:min(5, nrow(roc_df))]
colors <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00")
plot(0:1, 0:1, type="l", lty=2, col="grey60",
     xlab="False Positive Rate", ylab="True Positive Rate",
     main="ROC Curves — Top 5 tRFs (Tumor vs Normal)")
for (i in seq_along(top5)) {
  trf <- top5[i]
  expr <- as.numeric(log_cpm[trf, ])
  r <- roc(response, expr, quiet=TRUE)
  lines(1 - r$specificities, r$sensitivities,
        col=colors[i], lwd=2)
}
legend("bottomright",
       legend=paste0(top5, " (AUC=",
                     round(roc_df$AUC[1:length(top5)], 3), ")"),
       col=colors[1:length(top5)], lwd=2, cex=0.7, bty="n")
dev.off()
cat("ROC curves saved.\n")

# Bar plot of AUC values
p_auc <- ggplot(roc_df, aes(x=reorder(tRF, AUC), y=AUC, fill=AUC)) +
  geom_col() +
  geom_hline(yintercept=0.5, linetype="dashed", color="grey40") +
  scale_fill_gradient(low="#fdd49e", high="#d7301f") +
  coord_flip() +
  labs(title="AUC — Top 20 tRFs (Tumor vs Normal)",
       x=NULL, y="AUC") +
  theme_bw(base_size=11) +
  theme(legend.position="none")
ggsave(file.path(FIGS, "roc_auc_barplot.pdf"), p_auc,
       width=7, height=6)
cat("AUC barplot saved.\n")

# ── 2. Tumor Stage Association ──────────────────────────────────────────────────
cat("\nRunning stage association...\n")

clin_raw <- read.table(file.path(RAW, "gdc_clinical/clinical.tsv"),
                       sep="\t", header=TRUE, stringsAsFactors=FALSE,
                       fill=TRUE, quote="")

# Keep unique rows with stage info
# Use the primary KIRP diagnosis row for AJCC staging. The GDC clinical export
# can include additional non-primary diagnoses for the same patient.
clin <- clin_raw[
  as.character(clin_raw$diagnoses.diagnosis_is_primary_disease) %in% c("TRUE", "True", "true"),
  c("cases.submitter_id", "diagnoses.submitter_id",
    "diagnoses.ajcc_pathologic_stage")
]
colnames(clin) <- c("patient_id", "diagnosis_id", "stage")
clin <- clin[order(clin$patient_id, clin$diagnosis_id), ]
clin <- clin[!duplicated(clin$patient_id), ]
clin$stage[clin$stage == "'--"] <- NA
clin <- clin[!is.na(clin$stage) & clin$stage != "", ]

# Simplify stage: Stage I/II/III/IV
clin$stage_simple <- sub(" [ABC]$", "", clin$stage)
clin <- clin[clin$stage_simple %in% c("Stage I","Stage II","Stage III","Stage IV"), ]

cat(sprintf("Stage distribution:\n"))
print(table(clin$stage_simple))

# Tumor samples only
tumor_meta <- trf_meta[trf_meta$group == "Tumor", ]
# Extract 12-char patient ID from barcode
tumor_meta$patient_id <- substr(tumor_meta$sample_id, 1, 12)

# Merge with clinical
tumor_clin <- merge(tumor_meta, clin, by="patient_id")
cat(sprintf("Samples with stage info: %d\n", nrow(tumor_clin)))

if (nrow(tumor_clin) > 20) {
  # Top tRF by AUC for stage plot
  stage_trfs <- roc_df$tRF[1:min(6, nrow(roc_df))]

  stage_long <- lapply(stage_trfs, function(trf) {
    if (!trf %in% rownames(log_cpm)) return(NULL)
    data.frame(
      sample_id = tumor_clin$sample_id,
      stage     = tumor_clin$stage_simple,
      expr      = as.numeric(log_cpm[trf, tumor_clin$sample_id]),
      tRF       = trf
    )
  })
  stage_df <- do.call(rbind, Filter(Negate(is.null), stage_long))
  stage_df$stage <- factor(stage_df$stage,
                           levels=c("Stage I","Stage II","Stage III","Stage IV"))

  p_stage <- ggplot(stage_df, aes(x=stage, y=expr, fill=stage)) +
    geom_boxplot(outlier.size=0.8, alpha=0.8) +
    facet_wrap(~tRF, scales="free_y", ncol=3) +
    scale_fill_manual(values=c("#2166AC","#92C5DE","#F4A582","#D6604D")) +
    labs(title="tRF Expression by Tumor Stage",
         x="AJCC Pathologic Stage", y="log2(CPM+1)") +
    theme_bw(base_size=10) +
    theme(axis.text.x=element_text(angle=35, hjust=1),
          legend.position="none")
  ggsave(file.path(FIGS, "stage_association_boxplots.pdf"), p_stage,
         width=9, height=6)
  cat("Stage boxplots saved.\n")

  # Kruskal-Wallis p-values per tRF
  kw_res <- stage_df %>%
    group_by(tRF) %>%
    summarise(kruskal_p = kruskal.test(expr ~ stage)$p.value,
              .groups="drop") %>%
    mutate(kruskal_fdr = p.adjust(kruskal_p, method="BH"))
  write.table(kw_res, file.path(TABS, "stage_kruskal_results.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("Kruskal-Wallis results:\n")
  print(kw_res)
}

# ── 3. Kaplan–Meier Survival ────────────────────────────────────────────────────
cat("\nRunning survival analysis...\n")

clin_surv <- clin_raw[, c("cases.submitter_id",
                           "demographic.vital_status",
                           "demographic.days_to_death",
                           "diagnoses.days_to_last_follow_up")]
colnames(clin_surv) <- c("patient_id","vital_status","days_to_death","days_to_last_fu")
clin_surv <- clin_surv[!duplicated(clin_surv$patient_id), ]

# Clean up
clean_val <- function(x) {
  x[x %in% c("'--", "--", "")] <- NA
  as.numeric(x)
}
clin_surv$days_to_death  <- clean_val(clin_surv$days_to_death)
clin_surv$days_to_last_fu <- clean_val(clin_surv$days_to_last_fu)
clin_surv$event <- as.integer(clin_surv$vital_status == "Dead")
clin_surv$os_days <- ifelse(!is.na(clin_surv$days_to_death),
                             clin_surv$days_to_death,
                             clin_surv$days_to_last_fu)
clin_surv <- clin_surv[!is.na(clin_surv$os_days) & clin_surv$os_days > 0, ]

# Merge tumor samples with survival
surv_data <- merge(tumor_meta, clin_surv, by="patient_id")
cat(sprintf("Samples with survival data: %d (events: %d)\n",
            nrow(surv_data), sum(surv_data$event)))

if (nrow(surv_data) > 20 && sum(surv_data$event) >= 5) {
  # KM for top 3 tRFs by AUC — split by median expression
  km_trfs <- roc_df$tRF[1:min(3, nrow(roc_df))]

  pdf(file.path(FIGS, "km_survival_top3_trfs.pdf"), width=7, height=6)
  for (trf in km_trfs) {
    if (!trf %in% rownames(log_cpm)) next
    expr <- as.numeric(log_cpm[trf, surv_data$sample_id])
    med  <- median(expr, na.rm=TRUE)
    surv_data$group_km <- ifelse(expr >= med, "High", "Low")

    sf <- survfit(Surv(os_days, event) ~ group_km, data=surv_data)
    p_km <- ggsurvplot(sf, data=surv_data,
                       title=paste0("KM — ", trf),
                       palette=c("#D6604D","#4393C3"),
                       pval=TRUE, conf.int=FALSE,
                       risk.table=TRUE, risk.table.height=0.25,
                       xlab="Days", ylab="Overall Survival",
                       legend.labs=c("High","Low"),
                       legend.title="Expression")
    print(p_km)
  }
  dev.off()
  cat("KM curves saved.\n")

  # Cox univariate for all top 20 tRFs
  cox_res <- lapply(top_trfs, function(trf) {
    if (!trf %in% rownames(log_cpm)) return(NULL)
    expr <- as.numeric(log_cpm[trf, surv_data$sample_id])
    fit  <- coxph(Surv(os_days, event) ~ expr, data=surv_data)
    s    <- summary(fit)
    data.frame(tRF=trf,
               HR=s$coefficients[,"exp(coef)"],
               HR_lo=s$conf.int[,"lower .95"],
               HR_hi=s$conf.int[,"upper .95"],
               p_value=s$coefficients[,"Pr(>|z|)"],
               stringsAsFactors=FALSE)
  })
  cox_df <- do.call(rbind, Filter(Negate(is.null), cox_res))
  cox_df$FDR <- p.adjust(cox_df$p_value, method="BH")
  cox_df <- cox_df[order(cox_df$p_value), ]
  write.table(cox_df, file.path(TABS, "cox_univariate_results.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("Cox results:\n")
  print(cox_df[cox_df$p_value < 0.1, ])

  # Forest plot
  cox_plot <- cox_df[order(cox_df$HR), ]
  cox_plot$tRF <- factor(cox_plot$tRF, levels=cox_plot$tRF)
  p_forest <- ggplot(cox_plot, aes(x=HR, xmin=HR_lo, xmax=HR_hi,
                                    y=tRF, color=FDR < 0.05)) +
    geom_pointrange(size=0.5) +
    geom_vline(xintercept=1, linetype="dashed", color="grey40") +
    scale_color_manual(values=c("grey60","#D6604D"),
                       labels=c("FDR≥0.05","FDR<0.05"),
                       name=NULL) +
    labs(title="Cox Univariate HR — Top 20 tRFs",
         x="Hazard Ratio (95% CI)", y=NULL) +
    theme_bw(base_size=10)
  ggsave(file.path(FIGS, "cox_forest_plot.pdf"), p_forest,
         width=7, height=6)
  cat("Forest plot saved.\n")
} else {
  cat("Not enough survival events for KM/Cox — skipping.\n")
}

cat("\n✓ Phase 8 complete.\n")
