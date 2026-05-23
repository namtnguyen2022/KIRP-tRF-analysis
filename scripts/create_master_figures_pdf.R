## Create one master PDF containing main and supplemental/supporting figures.

suppressPackageStartupMessages({
  library(grid)
})

ROOT <- "/Users/nam.tnguyen2022/KIRP_tRF_project"
FIGS <- file.path(ROOT, "results", "figures")

OUT_PDF <- file.path(FIGS, "master_figures_with_supplemental.pdf")
OUT_MANIFEST <- file.path(FIGS, "master_figures_with_supplemental_manifest.tsv")

main_figures <- c(
  "fig2_trf_DE.pdf",
  "fig3_gene_DE.pdf",
  "fig4_trf_met_axis.pdf",
  "fig5_roc.pdf",
  "fig6_stage_survival.pdf"
)

supplemental_figures <- c(
  "fig2A_volcano.pdf",
  "fig2B_heatmap_tmp.pdf",
  "fig2C_boxplots.pdf",
  "fig3A_gene_volcano.pdf",
  "fig3BCD_gene_boxes.pdf",
  "trf_volcano.pdf",
  "trf_heatmap_top50.pdf",
  "trf_boxplots_top10.pdf",
  "gene_volcano.pdf",
  "gene_heatmap.pdf",
  "gene_boxplots_significant.pdf",
  "anticorrelation_heatmap.pdf",
  "anticorrelation_scatterplots_top.pdf",
  "roc_curves_top5.pdf",
  "roc_auc_barplot.pdf",
  "stage_association_boxplots.pdf",
  "km_survival_top3_trfs.pdf",
  "cox_forest_plot.pdf"
)

make_title_page <- function(path, title, subtitle = NULL, lines = character()) {
  pdf(path, width = 8.5, height = 11, useDingbats = FALSE)
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  grid.text(title, x = 0.08, y = 0.73, just = "left",
            gp = gpar(fontsize = 22, fontface = "bold"))
  if (!is.null(subtitle)) {
    grid.text(subtitle, x = 0.08, y = 0.68, just = "left",
              gp = gpar(fontsize = 12, col = "grey35"))
  }
  if (length(lines) > 0) {
    y <- 0.58
    for (line in lines) {
      grid.text(line, x = 0.08, y = y, just = "left",
                gp = gpar(fontsize = 10.5, col = "grey25"))
      y <- y - 0.032
    }
  }
  dev.off()
}

figure_path <- function(x) file.path(FIGS, x)

main_paths <- figure_path(main_figures)
supp_paths <- figure_path(supplemental_figures)

missing_main <- main_figures[!file.exists(main_paths)]
missing_supp <- supplemental_figures[!file.exists(supp_paths)]

if (length(missing_main) > 0) {
  stop("Missing required main figure(s): ", paste(missing_main, collapse = ", "))
}

if (length(missing_supp) > 0) {
  warning("Skipping missing supplemental/supporting figure(s): ",
          paste(missing_supp, collapse = ", "))
  supplemental_figures <- supplemental_figures[file.exists(supp_paths)]
  supp_paths <- figure_path(supplemental_figures)
}

tmpdir <- tempfile("master_figures_pages_")
dir.create(tmpdir)
on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

title_pdf <- file.path(tmpdir, "00_title.pdf")
main_divider <- file.path(tmpdir, "01_main_figures.pdf")
supp_divider <- file.path(tmpdir, "02_supplemental_figures.pdf")

make_title_page(
  title_pdf,
  "KIRP tRF Project",
  "Master figures PDF",
  c(
    paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
    paste0("Main figure PDFs: ", length(main_figures)),
    paste0("Supplemental/supporting figure PDFs: ", length(supplemental_figures)),
    "Source directory: results/figures"
  )
)

make_title_page(
  main_divider,
  "Main Figures",
  "Current composite manuscript figures",
  main_figures
)

make_title_page(
  supp_divider,
  "Supplemental and Supporting Figures",
  "Standalone figure outputs and supporting panels",
  supplemental_figures
)

ordered_files <- c(title_pdf, main_divider, main_paths,
                   supp_divider, supp_paths)

manifest <- data.frame(
  order = seq_along(ordered_files),
  section = c(
    "cover",
    "main_divider",
    rep("main", length(main_paths)),
    "supplemental_divider",
    rep("supplemental", length(supp_paths))
  ),
  file = c(
    basename(title_pdf),
    basename(main_divider),
    main_figures,
    basename(supp_divider),
    supplemental_figures
  ),
  path = c(
    "generated cover page",
    "generated section divider",
    main_paths,
    "generated section divider",
    supp_paths
  ),
  stringsAsFactors = FALSE
)

write.table(manifest, OUT_MANIFEST, sep = "\t",
            row.names = FALSE, quote = FALSE)

pdfunite <- Sys.which("pdfunite")
if (pdfunite == "") {
  stop("pdfunite was not found on PATH.")
}

if (file.exists(OUT_PDF)) {
  file.remove(OUT_PDF)
}

status <- system2(pdfunite, c(ordered_files, OUT_PDF))
if (!identical(status, 0L)) {
  stop("pdfunite failed with status ", status)
}

cat("Created master figures PDF:\n  ", OUT_PDF, "\n", sep = "")
cat("Created manifest:\n  ", OUT_MANIFEST, "\n", sep = "")
