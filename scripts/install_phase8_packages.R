pkgs <- c('pROC', 'survival', 'survminer')
for (p in pkgs) {
  if (!requireNamespace(p, quietly=TRUE)) {
    install.packages(p, repos='https://cloud.r-project.org')
  } else {
    cat(p, "already installed\n")
  }
}
