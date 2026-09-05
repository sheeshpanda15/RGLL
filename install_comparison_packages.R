repos <- "https://cloud.r-project.org"
needed <- c("mclust")
for (pkg in needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = repos)
}
cat("Optional comparison packages installed/available.\n")
