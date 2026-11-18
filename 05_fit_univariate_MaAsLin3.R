#!/usr/bin/env Rscript
# Fit separate MaAsLin3 models for age, race, ER and CTC.
# Upstream of Fig. 4; Data S5. Input: script 01.

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

suppressPackageStartupMessages({library(phyloseq);library(dplyr);library(maaslin3)})
physeq <- readRDS(file.path(PROJECT_DIR,"Output/physeq.rds"))
meta <- data.frame(sample_data(physeq),check.names=FALSE)
meta <- meta %>% mutate(ctc2=case_when(CTC_count_cat %in% c("0","0-5") ~ "0–5",
                                     CTC_count_cat==">5" ~ ">5", TRUE ~ NA_character_))
meta_ctc <- meta %>% filter(!is.na(ctc2))
physeq_ctc <- prune_samples(rownames(meta_ctc),physeq)
meta_ctc$ctc2 <- relevel(factor(meta_ctc$ctc2),ref="0–5")
run_single_factor <- function(physeq_obj, metadata_df, factor_name, out_dir) {

  abund <- as.data.frame(t(phyloseq::otu_table(physeq_obj)))
  meta  <- metadata_df[, factor_name, drop = FALSE]

  maaslin3(
    input_data    = abund,
    input_metadata = meta,
    output        = out_dir,
    formula       = paste0("~ ", factor_name),
    min_prevalence = 0.5,
    normalization = "TSS",
    transform = "LOG",
    augment = TRUE,
    standardize = TRUE,
    median_comparison_abundance = TRUE,
    median_comparison_prevalence = FALSE,
    max_significance = 0.25,
    cores = 4
    )

  cat("Finished:", factor_name, "\n")
}

base <- file.path(PROJECT_DIR,"Output/Maaslin3/Paper used analysis")
for (v in c("age_group","race","ER","ctc2")) {
  dest <- file.path(base,paste0(v,"_only_0.5"))
  if (file.exists(file.path(dest,"all_results.tsv"))) stop("Use a fresh private output directory for new fits.")
  if (v=="ctc2") run_single_factor(physeq_ctc,meta_ctc,v,dest) else run_single_factor(physeq,meta,v,dest)
}
save_session(base)
