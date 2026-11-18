#!/usr/bin/env Rscript
# Prepare case profiles: sample inclusion, decontam and MMUPHin batch correction.
# Upstream of Figs. 1-6 and S1-S6.
# Inputs: case_abundance.csv, case_metadata.csv, sample_inclusion.csv and final_case_metadata.csv.

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

suppressPackageStartupMessages({library(phyloseq); library(decontam); library(MMUPHin)})
tax_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
abund <- read.csv(private_input("case_abundance.csv"), check.names = FALSE)
meta <- read.csv(private_input("case_metadata.csv"), check.names = FALSE)
manifest <- read.csv(private_input("sample_inclusion.csv"), check.names = FALSE)
stopifnot(identical(names(abund)[1:7], tax_cols), !anyDuplicated(meta$Sample_id),
          !anyDuplicated(manifest$Sample_id), all(meta$Sample_id %in% names(abund)))
include <- as.logical(manifest$include_for_decontam)
stopifnot(!anyNA(include), setequal(manifest$Sample_id, meta$Sample_id))
ids <- as.character(manifest$Sample_id)
meta <- meta[match(ids, meta$Sample_id), , drop = FALSE]
rownames(meta) <- ids
tax <- as.matrix(abund[, tax_cols, drop = FALSE])
features <- apply(tax, 1, paste, collapse = ";")
stopifnot(!anyDuplicated(features))
rownames(tax) <- features
x <- as.matrix(abund[, ids, drop = FALSE]); storage.mode(x) <- "double"
rownames(x) <- features
stopifnot(all(is.finite(x)), all(x >= 0), all(colSums(x) > 0))
ps_raw <- phyloseq(otu_table(x, taxa_are_rows = TRUE), tax_table(tax), sample_data(meta))
out <- file.path(PROJECT_DIR, "Output")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
saveRDS(ps_raw, file.path(out, "physeq_raw.RDS"))

ps <- prune_samples(ids[include], ps_raw)
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"
stopifnot(any(sample_data(ps)$is.neg), any(!sample_data(ps)$is.neg))
contam <- decontam::isContaminant(ps, method = "prevalence", neg = "is.neg", threshold = 0.1)
contam_taxa <- rownames(contam)[contam$contaminant %in% TRUE]
ps_clean <- prune_taxa(!(taxa_names(ps) %in% contam_taxa), ps)
saveRDS(ps_clean, file.path(out, "phyloseq_cleaned.rds"))
flags <- data.frame(taxon = rownames(contam), contam, check.names = FALSE)
flag_out <- file.path(PROJECT_DIR, "revision_analysis", "02_contamination_control")
dir.create(flag_out, recursive = TRUE, showWarnings = FALSE)
write.csv(flags[flags$contaminant %in% TRUE, ], file.path(flag_out, "decontam_flagged_taxa_details.csv"), row.names = FALSE)

ps_case <- subset_samples(ps_clean, Sample_or_Control != "Control Sample")
ps_case <- transform_sample_counts(ps_case, function(x) x / sum(x))
ps_case <- prune_taxa(taxa_sums(ps_case) > 0, ps_case)
feature_abd <- as.data.frame(otu_table(ps_case))
md <- data.frame(sample_data(ps_case), check.names = FALSE)
md <- md[colnames(feature_abd), , drop = FALSE]
md$batch <- as.factor(md$batch)
stopifnot(identical(rownames(md), colnames(feature_abd)))
species_adj <- MMUPHin::adjust_batch(feature_abd = feature_abd, batch = "batch",
                                    covariates = NULL, data = md,
                                    control = list(verbose = TRUE))
corrected <- as.matrix(species_adj$feature_abd_adj)
stopifnot(identical(dimnames(corrected), dimnames(as.matrix(feature_abd))),
          all(is.finite(corrected)), all(corrected >= 0),
          max(abs(colSums(corrected) - 1)) < 1e-6)
otu_table(ps_case) <- otu_table(corrected, taxa_are_rows = TRUE)
saveRDS(ps_case, file.path(out, "phyloseq.rds"))

final_meta <- read.csv(private_input("final_case_metadata.csv"), check.names = FALSE)
stopifnot(!anyDuplicated(final_meta$Sample_id), !anyDuplicated(final_meta$Patient_id),
          setequal(final_meta$Sample_id, sample_names(ps_case)),
          all(is.finite(final_meta$survival_months)), all(final_meta$survival_months > 0),
          all(final_meta$vital_status %in% c(0, 1)))
final_meta$TNBC <- ifelse(final_meta$ER == "N" & final_meta$PR == "N" & final_meta$HER2 == "N", "Y",
                          ifelse(final_meta$ER == "Y" | final_meta$PR == "Y" | final_meta$HER2 == "Y", "N", NA_character_))
if ("BMI" %in% names(final_meta)) {
  final_meta$BMI_cat <- as.character(cut(final_meta$BMI, breaks = c(-Inf, 18.5, 25, 30, Inf),
                                       labels = c("Underweight", "Normal", "Overweight", "Obesity"), right = FALSE))
}
final_meta$plot_name <- paste0("P", seq_len(nrow(final_meta)))
rownames(final_meta) <- final_meta$Sample_id
mat <- as(otu_table(ps_case), "matrix")[, final_meta$Sample_id, drop = FALSE]
ps_final <- phyloseq(otu_table(mat, taxa_are_rows = TRUE), tax_table(ps_case), sample_data(final_meta))
stopifnot(identical(sample_names(ps_final), final_meta$Sample_id))
saveRDS(ps_final, file.path(out, "physeq.rds"))
save_session(out)
