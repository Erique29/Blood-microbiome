#!/usr/bin/env Rscript
# Shared configuration; run scripts 01-12 in order.
# Set SKCCC_PRIVATE_DIR to a private directory outside this folder.
# Inputs: <SKCCC_PRIVATE_DIR>/input; outputs: <SKCCC_PRIVATE_DIR>/analysis_output.

if (!exists("SKCCC_SCRIPTS_DIR", inherits = FALSE)) {
  .cli <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .srcs <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
  .src <- if (length(.srcs)) .srcs[[length(.srcs)]] else NULL
  SKCCC_SCRIPTS_DIR <- if (!is.null(.src)) dirname(normalizePath(.src)) else if (length(.cli)) {
    dirname(normalizePath(sub("^--file=", "", .cli[1])))
  } else normalizePath(getwd())
}
.private <- Sys.getenv("SKCCC_PRIVATE_DIR", "")
if (!nzchar(.private)) stop("Set SKCCC_PRIVATE_DIR to an approved private working directory before running.", call. = FALSE)
PRIVATE_DIR <- normalizePath(path.expand(.private), mustWork = TRUE)
.public <- normalizePath(SKCCC_SCRIPTS_DIR, mustWork = TRUE)
if (identical(PRIVATE_DIR, .public) || startsWith(PRIVATE_DIR, paste0(.public, .Platform$file.sep))) {
  stop("The private workspace must be outside the public script directory.", call. = FALSE)
}
INPUT_DIR <- file.path(PRIVATE_DIR, "input")
PROJECT_DIR <- file.path(PRIVATE_DIR, "analysis_output")
project <- PROJECT_DIR
private_input <- function(...) {
  f <- file.path(INPUT_DIR, ...)
  if (!file.exists(f)) stop("Required private input is missing: ", basename(f), call. = FALSE)
  f
}
save_session <- function(directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(utils::sessionInfo()), file.path(directory, "sessionInfo.txt"))
}
