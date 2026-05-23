#!/usr/bin/env Rscript
# generate_table1.R
# Generates Table 1 (cohort characteristics) and Table 2 (prioritised tRF candidates)
# from TCGA-KIRP metadata and analysis results.
# Outputs: results/tables/table1_cohort.tsv, results/tables/table2_candidates.tsv

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

# Edit scripts/config.R to set ROOT for your machine, or change the fallback below:
if (file.exists("scripts/config.R")) source("scripts/config.R") else if (file.exists("config.R")) source("config.R") else ROOT <- "/Users/nam.tnguyen2022/KIRP_tRF_project"
TABS  <- file.path(ROOT, "results/tables")

# ── Load data ─────────────────────────────────────────────────────────────────
clin   <- read_tsv(file.path(ROOT, "data_raw/gdc_clinical/clinical.tsv"),
                   show_col_types = FALSE, na = c("'--", "--", "NA", ""))
fup    <- read_tsv(file.path(ROOT, "data_raw/gdc_clinical/follow_up.tsv"),
                   show_col_types = FALSE, na = c("'--", "--", "NA", ""))
smeta  <- read_tsv(file.path(ROOT, "data_processed/sample_metadata.tsv"),
                   show_col_types = FALSE)
trf_meta  <- read_tsv(file.path(ROOT, "data_processed/trf_sample_metadata.tsv"),
                      show_col_types = FALSE)
mrna_meta <- read_tsv(file.path(ROOT, "data_processed/mrna_sample_metadata.tsv"),
                      show_col_types = FALSE)

# One row per patient (deduplicate clinical)
clin1 <- clin %>%
  distinct(`cases.submitter_id`, .keep_all = TRUE)

# Primary KIRP diagnosis row per patient. Some patients have additional
# non-primary diagnoses in the GDC clinical export; those should not contribute
# to KIRP AJCC stage summaries.
clin_primary <- clin %>%
  filter(`diagnoses.diagnosis_is_primary_disease` == TRUE) %>%
  arrange(`cases.submitter_id`, `diagnoses.submitter_id`) %>%
  distinct(`cases.submitter_id`, .keep_all = TRUE)

# ── Sample counts — read directly from the analysis-ready metadata files ──────
trf_tumor   <- trf_meta  %>% filter(sample_type_code == "01") %>% nrow()
trf_normal  <- trf_meta  %>% filter(sample_type_code == "11") %>% nrow()
mrna_tumor  <- mrna_meta %>% filter(sample_type_code == "01") %>% nrow()
mrna_normal <- mrna_meta %>% filter(sample_type_code == "11") %>% nrow()
both_tumor  <- length(intersect(
  trf_meta  %>% filter(sample_type_code == "01") %>% pull(patient_id),
  mrna_meta %>% filter(sample_type_code == "01") %>% pull(patient_id)
))
both_normal <- length(intersect(
  trf_meta  %>% filter(sample_type_code == "11") %>% pull(patient_id),
  mrna_meta %>% filter(sample_type_code == "11") %>% pull(patient_id)
))

cat(sprintf("tRF tumor=%d  normal=%d\n", trf_tumor, trf_normal))
cat(sprintf("mRNA tumor=%d  normal=%d\n", mrna_tumor, mrna_normal))
cat(sprintf("Matched both: tumor=%d  normal=%d\n", both_tumor, both_normal))

# ── Age ───────────────────────────────────────────────────────────────────────
ages <- clin1 %>%
  filter(!is.na(`demographic.age_at_index`)) %>%
  pull(`demographic.age_at_index`) %>%
  as.numeric()
age_med <- median(ages, na.rm = TRUE)
age_iqr <- quantile(ages, c(0.25, 0.75), na.rm = TRUE)
age_str <- sprintf("%.0f (IQR: %.0f–%.0f)", age_med, age_iqr[1], age_iqr[2])
cat("Age:", age_str, "\n")

# ── Sex ───────────────────────────────────────────────────────────────────────
sex_tab <- clin1 %>%
  count(`demographic.gender`) %>%
  filter(!is.na(`demographic.gender`))
n_total_sex <- sum(sex_tab$n)
sex_rows <- sex_tab %>%
  mutate(label = paste0(`demographic.gender`, ", n (%)"),
         val   = sprintf("%d (%.1f%%)", n, 100*n/n_total_sex))
cat("Sex:\n"); print(sex_rows)

# ── Race ──────────────────────────────────────────────────────────────────────
race_tab <- clin1 %>%
  mutate(race = case_when(
    grepl("white", `demographic.race`, ignore.case=TRUE) ~ "White",
    grepl("black", `demographic.race`, ignore.case=TRUE) ~ "Black or African American",
    grepl("asian", `demographic.race`, ignore.case=TRUE) ~ "Asian",
    grepl("not reported|unknown", `demographic.race`, ignore.case=TRUE) | is.na(`demographic.race`) ~ "Not reported/Unknown",
    TRUE ~ "Other"
  )) %>%
  count(race)
n_total_race <- sum(race_tab$n)
cat("Race:\n"); print(race_tab)

# ── AJCC stage ────────────────────────────────────────────────────────────────
stage_tab <- clin_primary %>%
  mutate(stage = `diagnoses.ajcc_pathologic_stage`) %>%
  filter(!is.na(stage), stage != "not reported") %>%
  mutate(stage = sub("Stage ", "Stage ", stage)) %>%
  count(stage) %>%
  arrange(stage)
n_staged <- sum(stage_tab$n)
cat("Stage:\n"); print(stage_tab)

# ── Survival ──────────────────────────────────────────────────────────────────
vital <- clin1 %>%
  mutate(
    vital_status = `demographic.vital_status`,
    days_to_death = as.numeric(`demographic.days_to_death`)
  )

n_surv_avail  <- vital %>% filter(!is.na(vital_status)) %>% nrow()
n_dead        <- vital %>% filter(grepl("dead", vital_status, ignore.case=TRUE)) %>% nrow()
pct_dead      <- 100 * n_dead / n_surv_avail

# Median follow-up: use days_to_death for deceased, days to last follow-up for alive
followup <- fup %>%
  group_by(`cases.submitter_id`) %>%
  summarise(max_followup = max(as.numeric(`follow_ups.days_to_follow_up`), na.rm=TRUE),
            .groups="drop") %>%
  filter(is.finite(max_followup))

dead_days <- vital %>%
  filter(grepl("dead", vital_status, ignore.case=TRUE), !is.na(days_to_death)) %>%
  select(patient_id = `cases.submitter_id`, days = days_to_death)

alive_days <- followup %>%
  rename(patient_id = `cases.submitter_id`, days = max_followup)

all_os <- bind_rows(dead_days, alive_days) %>%
  group_by(patient_id) %>%
  summarise(os_days = max(days, na.rm=TRUE), .groups="drop") %>%
  filter(is.finite(os_days), os_days > 0)

med_fu_days   <- median(all_os$os_days, na.rm=TRUE)
med_fu_months <- round(med_fu_days / 30.44, 1)
cat(sprintf("Median follow-up: %.1f months (%d days)\n", med_fu_months, round(med_fu_days)))
cat(sprintf("Deaths: %d / %d (%.1f%%)\n", n_dead, n_surv_avail, pct_dead))

# ── Build Table 1 ─────────────────────────────────────────────────────────────
# Helper for race
race_white  <- race_tab %>% filter(race=="White") %>% pull(n) %>% sum()
race_black  <- race_tab %>% filter(race=="Black or African American") %>% pull(n) %>% sum()
race_asian  <- race_tab %>% filter(race=="Asian") %>% pull(n) %>% sum()
race_other  <- race_tab %>% filter(!race %in% c("White","Black or African American","Asian","Not reported/Unknown")) %>% pull(n) %>% sum()
race_nr     <- race_tab %>% filter(race=="Not reported/Unknown") %>% pull(n) %>% sum()

fmt_pct <- function(n, tot) sprintf("%d (%.1f%%)", n, 100*n/tot)

stage_i   <- stage_tab %>% filter(stage %in% c("Stage I","Stage IA","Stage IB")) %>% pull(n) %>% sum()
stage_ii  <- stage_tab %>% filter(stage %in% c("Stage II","Stage IIA","Stage IIB")) %>% pull(n) %>% sum()
stage_iii <- stage_tab %>% filter(stage %in% c("Stage III","Stage IIIA","Stage IIIB","Stage IIIC")) %>% pull(n) %>% sum()
stage_iv  <- stage_tab %>% filter(stage %in% c("Stage IV","Stage IVA","Stage IVB","Stage IVC")) %>% pull(n) %>% sum()

male_n   <- sex_rows %>% filter(grepl("male", `demographic.gender`, ignore.case=TRUE) &
                                  !grepl("female", `demographic.gender`, ignore.case=TRUE)) %>% pull(n) %>% sum()
female_n <- sex_rows %>% filter(grepl("female", `demographic.gender`, ignore.case=TRUE)) %>% pull(n) %>% sum()
sex_tot  <- male_n + female_n

table1 <- tibble(
  Characteristic = c(
    "Sample availability",
    "  tRF data: tumor / solid normal",
    "  mRNA data: tumor / solid normal",
    "  Matched tRF+mRNA: tumor / solid normal",
    "Patients with clinical data, n",
    "Age at diagnosis, median (IQR), years",
    "Sex, n (%)",
    "  Male",
    "  Female",
    "Race/ethnicity, n (%)",
    "  White",
    "  Black or African American",
    "  Asian",
    "  Other",
    "  Not reported / Unknown",
    "AJCC pathologic stage, n (%)",
    "  Stage I",
    "  Stage II",
    "  Stage III",
    "  Stage IV",
    "  Not reported",
    "Overall survival data available, n",
    "Death events, n (%)",
    "Median follow-up, months"
  ),
  Value = c(
    "",
    sprintf("%d / %d", trf_tumor, trf_normal),
    sprintf("%d / %d", mrna_tumor, mrna_normal),
    sprintf("%d / %d", both_tumor, both_normal),
    as.character(nrow(clin1)),
    age_str,
    "",
    fmt_pct(male_n, sex_tot),
    fmt_pct(female_n, sex_tot),
    "",
    fmt_pct(race_white, n_total_race),
    fmt_pct(race_black, n_total_race),
    fmt_pct(race_asian, n_total_race),
    fmt_pct(race_other, n_total_race),
    fmt_pct(race_nr,    n_total_race),
    "",
    fmt_pct(stage_i,   n_staged),
    fmt_pct(stage_ii,  n_staged),
    fmt_pct(stage_iii, n_staged),
    fmt_pct(stage_iv,  n_staged),
    fmt_pct(nrow(clin1) - n_staged, nrow(clin1)),
    as.character(n_surv_avail),
    sprintf("%d (%.1f%%)", n_dead, pct_dead),
    as.character(med_fu_months)
  )
)

write_tsv(table1, file.path(TABS, "table1_cohort.tsv"))
cat("\n✓ Table 1 saved:", file.path(TABS, "table1_cohort.tsv"), "\n")
print(table1, n=Inf)

# ── Build Table 2 (5-column evidence-based format) ────────────────────────────
table2 <- tibble(
  `Candidate tRF`         = c(
    "tRF-23-79MP9P9MDD",
    "tRF-23-7SIR3DR2DV",
    "tRF-22-7SIRMM121",
    "tRF-22-7SERML921",
    "tRF-30-81PV6RRNLNK8"
  ),
  `Evidence category`     = c(
    "Regulatory candidate",
    "Internal discrimination / stage",
    "Internal discrimination / stage",
    "Internal discrimination / stage",
    "Exploratory survival"
  ),
  `Tumor direction`       = c(
    "Down",
    "Down",
    "Down",
    "Down",
    "Down"
  ),
  `Key supporting evidence` = c(
    "Anticorrelated with MET: r = \u22120.380, FDR = 9.5\u00d710\u207b\u2078; RNA22 MET 3\u02b9UTR pos. 522, \u0394G = \u221219.0 kcal/mol, p = 0.036",
    "AUC = 0.949; stage FDR = 0.018; requires external validation",
    "AUC = 0.929; stage FDR = 7.7\u00d710\u207b\u2074",
    "AUC = 0.901; stage FDR = 0.011",
    "Univariable Cox HR = 1.16, 95% CI 1.06\u20131.28; FDR = 0.033\u1d43"
  ),
  `Interpretation`        = c(
    "Candidate tRF\u2013MET regulatory interaction",
    "Strong internal tumor-normal discrimination; stage-associated in primary-diagnosis analysis",
    "Candidate internal discrimination / stage-associated tRF",
    "Candidate internal discrimination / stage-associated tRF",
    "Exploratory survival-associated tRF"
  )
)

write_tsv(table2, file.path(TABS, "table2_candidates.tsv"))
cat("\n\u2713 Table 2 saved:", file.path(TABS, "table2_candidates.tsv"), "\n")
print(table2, n=Inf, width=Inf)
cat("\nFootnote: \u1d43 Exploratory univariable association only; no adjustment for age, sex, or stage. Not equivalent to independent prognostic evidence.\n")
