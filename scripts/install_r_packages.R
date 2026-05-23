pkgs <- c('edgeR','ggplot2','ggrepel','pheatmap','RColorBrewer','gridExtra')
missing <- pkgs[!pkgs %in% installed.packages()[,'Package']]
if (length(missing) > 0) {
  cat('Installing:', paste(missing, collapse=', '), '\n')
  if (!requireNamespace('BiocManager', quietly=TRUE))
    install.packages('BiocManager', repos='https://cloud.r-project.org')
  bioc_pkgs <- intersect(missing, c('edgeR','limma'))
  cran_pkgs <- setdiff(missing, bioc_pkgs)
  if (length(bioc_pkgs) > 0) BiocManager::install(bioc_pkgs, ask=FALSE)
  if (length(cran_pkgs) > 0) install.packages(cran_pkgs, repos='https://cloud.r-project.org')
} else {
  cat('All packages already installed.\n')
}
