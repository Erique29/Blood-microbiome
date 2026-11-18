#!/usr/bin/env Rscript
# Prepare healthy-control profiles using prevalence-based decontam.
# Upstream of Figs. 1 and 5F.
# Inputs: healthy_abundance.csv, healthy_metadata.csv and the raw object from script 01.

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

suppressPackageStartupMessages({library(dplyr);library(tibble);library(phyloseq);library(decontam)})
old_raw_controls_path <- file.path(PROJECT_DIR, "Output", "physeq_raw.RDS")
out_dir <- file.path(PROJECT_DIR, "Output", "healthy_controls")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tax_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

make_taxa_id <- function(tax_df) {
  tax_df[] <- lapply(tax_df, function(x) ifelse(is.na(x) | x == "", "Unclassified", as.character(x)))
  make.unique(apply(tax_df, 1, paste, collapse = ";"))
}

otu_mat <- function(ps) {
  mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) mat <- t(mat)
  storage.mode(mat) <- "numeric"
  mat
}

tax_mat <- function(ps) {
  tx <- as(tax_table(ps), "matrix")
  colnames(tx) <- colnames(tax_table(ps))
  tx
}

combine_phyloseq <- function(ps_list) {
  all_taxa <- Reduce(union, lapply(ps_list, taxa_names))
  all_ranks <- Reduce(union, lapply(ps_list, function(ps) colnames(tax_table(ps))))
  all_samples <- unlist(lapply(ps_list, sample_names), use.names = FALSE)

  otu_out <- matrix(0, nrow = length(all_taxa), ncol = length(all_samples),
                    dimnames = list(all_taxa, all_samples))
  tax_out <- matrix(NA_character_, nrow = length(all_taxa), ncol = length(all_ranks),
                    dimnames = list(all_taxa, all_ranks))
  sample_frames <- list()

  sample_start <- 1
  for (i in seq_along(ps_list)) {
    ps <- ps_list[[i]]
    mat <- otu_mat(ps)
    tx <- tax_mat(ps)
    smp <- data.frame(sample_data(ps), check.names = FALSE)
    smp$.sample_name <- rownames(smp)
    sample_frames[[i]] <- smp

    sample_idx <- sample_start:(sample_start + nsamples(ps) - 1)
    otu_out[rownames(mat), sample_idx] <- mat
    tx_new <- tx[, intersect(colnames(tx), all_ranks), drop = FALSE]
    missing_tax <- is.na(tax_out[rownames(tx_new), colnames(tx_new), drop = FALSE])
    existing <- tax_out[rownames(tx_new), colnames(tx_new), drop = FALSE]
    existing[missing_tax] <- tx_new[missing_tax]
    tax_out[rownames(tx_new), colnames(tx_new)] <- existing
    sample_start <- sample_start + nsamples(ps)
  }

  sample_df <- bind_rows(sample_frames)
  rownames(sample_df) <- sample_df$.sample_name
  sample_df$.sample_name <- NULL
  sample_df <- sample_df[all_samples, , drop = FALSE]

  phyloseq(
    otu_table(otu_out, taxa_are_rows = TRUE),
    tax_table(tax_out),
    sample_data(sample_df)
  )
}

read_health_abundance <- function() {
  abund <- read.csv(private_input("healthy_abundance.csv"), check.names = FALSE)
  stopifnot(identical(names(abund)[1:7], tax_cols))
  list(abundance = abund)
}

build_health_ps <- function(abund_info) {
  health_meta <- read.csv(private_input("healthy_metadata.csv"), check.names = FALSE)
  stopifnot(!anyDuplicated(health_meta$Sample_id))
  health_meta$sampleID <- health_meta$Sample_id
  health_meta$Sample_or_Control <- "True Sample"
  health_meta$case_control <- "control"
  rownames(health_meta) <- health_meta$Sample_id

  abund <- abund_info$abundance
  keep_samples <- intersect(health_meta$Sample_id, colnames(abund))
  stopifnot(setequal(keep_samples, health_meta$Sample_id), length(keep_samples) > 0)
  health_meta <- health_meta[keep_samples, , drop = FALSE]

  tax_df <- abund[, tax_cols, drop = FALSE]
  taxa_id <- make_taxa_id(tax_df)
  rownames(tax_df) <- taxa_id
  otu_df <- abund[, keep_samples, drop = FALSE]
  rownames(otu_df) <- taxa_id

  phyloseq(
    otu_table(as.matrix(otu_df), taxa_are_rows = TRUE),
    tax_table(as.matrix(tax_df)),
    sample_data(health_meta)
  )
}

build_old_controls_ps <- function() {
  old_raw <- readRDS(old_raw_controls_path)
  old_meta <- data.frame(sample_data(old_raw), check.names = FALSE)
  control_names <- sample_names(old_raw)[old_meta$Sample_or_Control == "Control Sample"]
  controls <- prune_samples(control_names, old_raw)
  controls <- prune_taxa(taxa_sums(controls) > 0, controls)

  meta <- data.frame(sample_data(controls), check.names = FALSE)
  meta$Sample_id <- rownames(meta)
  meta$sampleID <- rownames(meta)
  meta$case_control <- "negative_control"
  meta$age <- NA_real_
  meta$age_MBC_dx <- NA_real_
  meta$interviewyr <- NA_real_
  meta$collectionyr <- NA_real_
  sample_data(controls) <- sample_data(meta)
  controls
}

run_decontam <- function(ps_health, ps_controls) {
  ps_raw <- combine_phyloseq(list(ps_health, ps_controls))
  sample_data(ps_raw)$is.neg <- sample_data(ps_raw)$Sample_or_Control == "Control Sample"
  contam <- isContaminant(ps_raw, method = "prevalence", neg = "is.neg", threshold = 0.1)
  contam$FeatureID <- rownames(contam)
  contam <- contam %>% relocate(FeatureID)

  contam_taxa <- contam$FeatureID[contam$contaminant %in% TRUE]
  ps_clean_with_controls <- prune_taxa(!(taxa_names(ps_raw) %in% contam_taxa), ps_raw)
  ps_health_clean <- subset_samples(ps_clean_with_controls, Sample_or_Control == "True Sample")
  ps_health_clean <- prune_taxa(taxa_sums(ps_health_clean) > 0, ps_health_clean)

  list(raw = ps_raw, clean_with_controls = ps_clean_with_controls,
       health_clean = ps_health_clean, contaminant_table = contam)
}

health_info <- read_health_abundance()
ps_health <- build_health_ps(health_info)
ps_controls <- build_old_controls_ps()
decontam_result <- run_decontam(ps_health, ps_controls)
saveRDS(ps_health, file.path(out_dir, "phyloseq_health_raw.rds"))
saveRDS(decontam_result$raw, file.path(out_dir, "phyloseq_health_plus_negative_raw.rds"))
saveRDS(decontam_result$health_clean, file.path(out_dir, "phyloseq_health_decontam_clean.rds"))
write.csv(decontam_result$contaminant_table, file.path(out_dir, "decontam_prevalence_results.csv"), row.names = FALSE)
save_session(out_dir)
