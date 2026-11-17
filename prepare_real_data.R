#!/usr/bin/env Rscript
# Recover the original input tables without restoring legacy experiment code.
prepare_real_data_inputs <- function(directory = "real_data/data", revision = "0f45646",
                                     datasets = c("student", "wage")) {
  files <- c(
    "student-mat.csv" = "HPC test/true dataset/student/student-mat.csv",
    "student-por.csv" = "HPC test/true dataset/student/student-por.csv",
    "2004.csv" = "HPC test/true dataset/wage/2004.csv"
  )
  stopifnot(length(datasets) > 0L, all(datasets %in% c("student", "wage")))
  files <- files[c("student", "student", "wage") %in% datasets]
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  for (name in names(files)) {
    target <- file.path(directory, name)
    if (file.exists(target) && file.info(target)$size > 0) next
    status <- system2("git", c("show", shQuote(paste0(revision, ":", files[[name]]))),
                      stdout = target)
    if (status != 0L || !file.exists(target) || file.info(target)$size == 0) {
      unlink(target)
      stop("Cannot recover ", name, ". Copy the original real_data/data directory to this project.")
    }
  }
  list(revision = revision, original_paths = files,
       files = normalizePath(file.path(directory, names(files))),
       md5 = tools::md5sum(file.path(directory, names(files))))
}

if (sys.nframe() == 0L) {
  inputs <- prepare_real_data_inputs()
  print(inputs$md5)
}
