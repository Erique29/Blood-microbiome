#!/usr/bin/env Rscript
# Calculate alpha-diversity metrics from corrected profiles.
# Upstream of Fig. 3 and Fig. S2; input: script 01.

local({
  .files <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
  .file <- if (length(.files)) .files[[length(.files)]] else NULL
  .cli <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .dir <- if (!is.null(.file)) dirname(normalizePath(.file)) else if (length(.cli)) {
    dirname(normalizePath(sub("^--file=", "", .cli[1])))
  } else getwd()
  assign("SKCCC_SCRIPTS_DIR", .dir, envir = .GlobalEnv)
  sys.source(file.path(.dir, "00_config_and_analysis_framework.R"), envir = .GlobalEnv)
})

suppressPackageStartupMessages({library(phyloseq);library(microbiome)})
physeq <- readRDS(file.path(PROJECT_DIR,"Output/physeq.rds"))
physeq_alpha_div <- microbiome::alpha(physeq, index="all")
physeq_meta <- data.frame(sample_data(physeq), check.names=FALSE)
physeq_alpha_div$Sample_id <- rownames(physeq_alpha_div)
physeq_alpha_div <- merge(physeq_alpha_div,physeq_meta,by="Sample_id")
physeq_alpha_div$Race <- physeq_alpha_div$race
write.csv(physeq_alpha_div,file.path(PROJECT_DIR,"Output/physeq_alpha_div.csv"))
save_session(file.path(PROJECT_DIR,"Output"))
