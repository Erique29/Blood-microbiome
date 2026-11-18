#!/usr/bin/env Rscript
# Summarize fastp QC and contamination screening.
# Data S1; upstream of Fig. 1.
# Inputs: fastp_report_map.csv, fastp_reports/ and outputs from script 01.

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

suppressPackageStartupMessages({library(phyloseq);library(decontam);library(dplyr);library(tidyr);library(readr);library(jsonlite);library(tibble)})
root <- file.path(PROJECT_DIR, "revision_analysis")
for (d in c("01_sequencing_qc", "02_contamination_control")) dir.create(file.path(root,d), recursive=TRUE, showWarnings=FALSE)
write_csv2 <- function(x,path) readr::write_csv(x,path,na="")
safe_ratio <- function(a,b) ifelse(is.na(b) | b == 0, NA_real_, a/b)
get_otu_taxa_by_samples <- function(ps) {
  x <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) x <- t(x)
  as.data.frame(x, check.names=FALSE)
}
ps_raw0 <- readRDS(file.path(PROJECT_DIR,"Output/physeq_raw.RDS"))
ps_clean <- readRDS(file.path(PROJECT_DIR,"Output/phyloseq_cleaned.rds"))
ps_final <- readRDS(file.path(PROJECT_DIR,"Output/physeq.rds"))
manifest <- read.csv(private_input("sample_inclusion.csv"),check.names=FALSE)
stopifnot(!anyNA(as.logical(manifest$include_for_decontam)))
excluded_nonanalytic <- manifest$Sample_id[!as.logical(manifest$include_for_decontam)]
report_map <- read.csv(private_input("fastp_report_map.csv"), check.names=FALSE)
stopifnot(!anyDuplicated(report_map$report_file), !anyDuplicated(report_map$Sample_id),
          all(basename(report_map$report_file)==report_map$report_file))
read_fastp <- function(path) {
  x <- fromJSON(path, simplifyVector = FALSE)
  sample <- sub("\\.fastp\\.json$", "", basename(path))
  before <- x$summary$before_filtering
  after <- x$summary$after_filtering
  filtering <- x$filtering_result
  adapter <- x$adapter_cutting
  dup <- x$duplication
  tibble(
    fastp_sample = sample,
    reads_before = before$total_reads,
    reads_after = after$total_reads,
    read_pairs_before = before$total_reads / 2,
    read_pairs_after = after$total_reads / 2,
    bases_before = before$total_bases,
    bases_after = after$total_bases,
    q20_rate_after = after$q20_rate,
    q30_rate_after = after$q30_rate,
    gc_content_after = after$gc_content,
    mean_length_after = after$read1_mean_length,
    passed_filter_reads = filtering$passed_filter_reads,
    low_quality_reads = filtering$low_quality_reads,
    too_many_N_reads = filtering$too_many_N_reads,
    too_short_reads = filtering$too_short_reads,
    adapter_trimmed_reads = adapter$adapter_trimmed_reads,
    adapter_trimmed_bases = adapter$adapter_trimmed_bases,
    duplication_rate = dup$rate
  )
}
make_summary <- function(dat, metrics, group_var = NULL) {
  if (is.null(group_var)) {
    dat$.group <- "All included samples"
    group_var <- ".group"
  }
  dat %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      across(
        all_of(metrics),
        list(
          n = ~sum(!is.na(.x)),
          mean = ~mean(.x, na.rm = TRUE),
          median = ~median(.x, na.rm = TRUE),
          sd = ~sd(.x, na.rm = TRUE),
          min = ~min(.x, na.rm = TRUE),
          max = ~max(.x, na.rm = TRUE)
        ),
        .names = "{.col}__{.fn}"
      ),
      .groups = "drop"
    ) %>%
    rename(group = all_of(group_var)) %>%
    pivot_longer(-group, names_to = c("metric", "stat"), names_sep = "__", values_to = "value") %>%
    pivot_wider(names_from = stat, values_from = value)
}

reports <- lapply(report_map$report_file, function(n) read_fastp(private_input("fastp_reports",n)))
fastp <- bind_rows(reports); fastp$Sample_id <- report_map$Sample_id
raw_meta <- data.frame(sample_data(ps_raw0),check.names=FALSE)
idx <- match(fastp$Sample_id, rownames(raw_meta)); stopifnot(!anyNA(idx))
fastp$analysis_group <- ifelse(raw_meta$Sample_or_Control[idx]=="Control Sample", "Control",
                               ifelse(fastp$Sample_id %in% sample_names(ps_final),"Article true sample","Excluded"))
fastp_included <- fastp %>% filter(analysis_group %in% c("Article true sample","Control")) %>%
  mutate(read_retention_rate=reads_after/reads_before,
         read_filter_removed_rate=1-read_retention_rate,
         adapter_trimmed_rate=adapter_trimmed_reads/reads_before) %>% arrange(analysis_group,Sample_id)
qc_metrics <- c("reads_before","reads_after","read_pairs_after","read_retention_rate",
                "read_filter_removed_rate","q20_rate_after","q30_rate_after","gc_content_after",
                "adapter_trimmed_rate","duplication_rate")
write_csv2(fastp_included,file.path(root,"01_sequencing_qc/fastp_qc_per_sample_included.csv"))
write_csv2(bind_rows(make_summary(fastp_included,qc_metrics,"analysis_group"),
                     make_summary(fastp_included,qc_metrics)),file.path(root,"01_sequencing_qc/fastp_qc_summary_included.csv"))
ps_raw <- prune_samples(!(sample_names(ps_raw0) %in% excluded_nonanalytic), ps_raw0)
sample_data(ps_raw)$is.neg <- sample_data(ps_raw)$Sample_or_Control == "Control Sample"
contamdf <- decontam::isContaminant(ps_raw, method = "prevalence", neg = "is.neg", threshold = 0.1)
contamdf$taxon <- rownames(contamdf)
contam_taxa <- contamdf$taxon[contamdf$contaminant]

tax_raw <- as.data.frame(tax_table(ps_raw), check.names = FALSE)
tax_raw$taxon <- rownames(tax_raw)
meta_raw <- data.frame(sample_data(ps_raw), check.names = FALSE)
raw_otu <- get_otu_taxa_by_samples(ps_raw)
article_samples_raw <- intersect(sample_names(ps_final), colnames(raw_otu))
control_samples_raw <- rownames(meta_raw)[meta_raw$Sample_or_Control == "Control Sample"]

all_taxa <- tax_raw %>%
  mutate(
    contaminant = taxon %in% contam_taxa,
    retained_after_decontam = taxon %in% taxa_names(ps_clean),
    retained_final_article_object = taxon %in% taxa_names(ps_final),
    article_prevalence_n_raw = rowSums(raw_otu[taxon, article_samples_raw, drop = FALSE] > 0, na.rm = TRUE),
    control_prevalence_n_raw = rowSums(raw_otu[taxon, control_samples_raw, drop = FALSE] > 0, na.rm = TRUE),
    article_prevalence_raw = article_prevalence_n_raw / length(article_samples_raw),
    control_prevalence_raw = control_prevalence_n_raw / length(control_samples_raw),
    article_mean_abundance_raw_pct = 100 * rowMeans(raw_otu[taxon, article_samples_raw, drop = FALSE], na.rm = TRUE),
    control_mean_abundance_raw_pct = 100 * rowMeans(raw_otu[taxon, control_samples_raw, drop = FALSE], na.rm = TRUE),
    article_max_abundance_raw_pct = 100 * apply(raw_otu[taxon, article_samples_raw, drop = FALSE], 1, max, na.rm = TRUE),
    control_max_abundance_raw_pct = 100 * apply(raw_otu[taxon, control_samples_raw, drop = FALSE], 1, max, na.rm = TRUE),
    control_to_article_mean_ratio = safe_ratio(control_mean_abundance_raw_pct, article_mean_abundance_raw_pct)
  ) %>%
  left_join(contamdf %>% select(taxon, p, prev), by = "taxon") %>%
  arrange(desc(contaminant), desc(control_prevalence_raw), desc(control_mean_abundance_raw_pct))
write_csv2(all_taxa, file.path(root, "02_contamination_control", "all_taxa_article_samples_vs_controls_decontam.csv"))

per_sample_contam <- tibble(
  Sample_id = c(article_samples_raw, control_samples_raw),
  analysis_group = c(rep("Article true sample", length(article_samples_raw)), rep("Control", length(control_samples_raw))),
  contaminant_abundance_pct = 100 * colSums(raw_otu[contam_taxa, c(article_samples_raw, control_samples_raw), drop = FALSE], na.rm = TRUE),
  total_abundance_pct = 100 * colSums(raw_otu[, c(article_samples_raw, control_samples_raw), drop = FALSE], na.rm = TRUE),
  contaminant_fraction_of_profile_pct = 100 * colSums(raw_otu[contam_taxa, c(article_samples_raw, control_samples_raw), drop = FALSE], na.rm = TRUE) /
    colSums(raw_otu[, c(article_samples_raw, control_samples_raw), drop = FALSE], na.rm = TRUE)
)
write_csv2(per_sample_contam, file.path(root, "02_contamination_control", "per_sample_contaminant_abundance_burden_article_samples_and_controls.csv"))

decontam_summary <- tibble(
  metric = c(
    "Article true samples used for decontam/control comparison",
    "Controls used for decontam/control comparison",
    "Excluded non-analytic true samples",
    "Raw taxa evaluated by decontam",
    "Taxa flagged as contaminants",
    "Percent taxa removed",
    "Mean contaminant abundance burden in article samples (% normalized MAP2B profile)",
    "Median contaminant abundance burden in article samples (% normalized MAP2B profile)",
    "Range contaminant abundance burden in article samples (% normalized MAP2B profile)",
    "Mean contaminant abundance burden in controls (% normalized MAP2B profile)",
    "Total contaminant abundance fraction in article samples (%)",
    "Total contaminant abundance fraction in controls (%)",
    "Final analytic taxa in article object"
  ),
  value = c(
    length(article_samples_raw),
    length(control_samples_raw),
    length(excluded_nonanalytic),
    ntaxa(ps_raw),
    length(contam_taxa),
    round(100 * length(contam_taxa) / ntaxa(ps_raw), 2),
    round(mean(per_sample_contam$contaminant_fraction_of_profile_pct[per_sample_contam$analysis_group == "Article true sample"], na.rm = TRUE), 4),
    round(median(per_sample_contam$contaminant_fraction_of_profile_pct[per_sample_contam$analysis_group == "Article true sample"], na.rm = TRUE), 4),
    paste0(
      round(min(per_sample_contam$contaminant_fraction_of_profile_pct[per_sample_contam$analysis_group == "Article true sample"], na.rm = TRUE), 4),
      "-",
      round(max(per_sample_contam$contaminant_fraction_of_profile_pct[per_sample_contam$analysis_group == "Article true sample"], na.rm = TRUE), 4)
    ),
    round(mean(per_sample_contam$contaminant_fraction_of_profile_pct[per_sample_contam$analysis_group == "Control"], na.rm = TRUE), 4),
    round(100 * sum(raw_otu[contam_taxa, article_samples_raw, drop = FALSE], na.rm = TRUE) / sum(raw_otu[, article_samples_raw, drop = FALSE], na.rm = TRUE), 4),
    round(100 * sum(raw_otu[contam_taxa, control_samples_raw, drop = FALSE], na.rm = TRUE) / sum(raw_otu[, control_samples_raw, drop = FALSE], na.rm = TRUE), 4),
    ntaxa(ps_final)
  )
)
write_csv2(decontam_summary, file.path(root, "02_contamination_control", "decontam_removal_summary_article_samples_and_controls.csv"))

top_contam <- all_taxa %>%
  filter(contaminant) %>%
  select(taxon, Phylum, Genus, Species, p, prev, article_prevalence_raw, control_prevalence_raw,
         article_mean_abundance_raw_pct, control_mean_abundance_raw_pct, control_to_article_mean_ratio) %>%
  arrange(desc(control_prevalence_raw), desc(control_mean_abundance_raw_pct), desc(article_mean_abundance_raw_pct))
write_csv2(top_contam, file.path(root, "02_contamination_control", "decontam_flagged_taxa_details.csv"))
write_csv2(head(top_contam, 30), file.path(root, "02_contamination_control", "top30_decontam_flagged_taxa.csv"))

save_session(file.path(root,"01_sequencing_qc"))
