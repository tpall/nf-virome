#!/usr/bin/env Rscript
# Build an nf-virome samplesheet from the eluring data_paths.csv manifest.
# Resolves per-sample HPC paths for assembled contigs, processed reads, and
# refined bins, so one shared vOTU catalog is built across all selected samples.
#
# Selection is by an explicit sample-id list (--samples), a project list
# (--projects), or both (intersection). The 620-sample aging cohort is driven
# off assets/aging_samples_620.txt (written from shared/output/ps_aging.rds).
#
# Usage:
#   Rscript assets/build_samplesheet.R <output.csv> \
#     [--samples <file>] [--projects p1,p2,...] \
#     [--eluring <path>] [--hpc-root <path>] \
#     [--check-files] [--single-end <ids|file>]
#
# Examples:
#   # 620-sample aging cohort (8 projects, one shared catalog)
#   Rscript assets/build_samplesheet.R assets/samplesheet_aging620.csv \
#     --samples assets/aging_samples_620.txt --eluring ../eluring
#
#   # A single project
#   Rscript assets/build_samplesheet.R assets/samplesheet_sajased.csv \
#     --projects sajased --eluring ../eluring
#
# HPC directory:
#   The data_paths.csv `folder` column is the authoritative per-sample HPC
#   directory (e.g. sajased -> "sajased-aftekas"). Paths are built under
#   <hpc-root>/<folder>/results/... so no project->dir mapping can drift.
#
# Single-end handling:
#   Paired-end reads are <sample>_processed_{1,2}.fq.gz; single-end is
#   <sample>_processed.fq.gz. For single-end rows reads_1 holds the unsuffixed
#   file and reads_2 is empty.
#   --check-files   stat each row's reads_1; if absent, fall back to the
#                   single-end name. Only meaningful where the HPC paths are
#                   reachable (i.e. run this on the HPC before launching).
#   --single-end    comma-separated ids, OR a file with one id per line; forces
#                   those rows single-end without a file check (for local builds).
#
# Output columns: sample, project, cohort, contigs, reads_1, reads_2, bins_dir,
# single_end. nf-virome reads sample/contigs/reads_1/reads_2/bins_dir/single_end
# by header name; project/cohort ride along as provenance.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || startsWith(args[1], "--")) {
  stop("Usage: build_samplesheet.R <output.csv> [--samples <file>] [--projects p1,p2,...] [--eluring <path>] [--hpc-root <path>] [--check-files] [--single-end <ids|file>]")
}
output_path <- args[1]

opt_val <- function(flag) {
  if (flag %in% args) args[which(args == flag) + 1L] else NA_character_
}
read_list <- function(x) {
  if (is.na(x)) return(character(0))
  v <- if (file.exists(x)) readLines(x, warn = FALSE) else unlist(strsplit(x, ",", fixed = TRUE))
  trimws(v[nzchar(trimws(v))])
}

eluring_root <- if ("--eluring" %in% args) opt_val("--eluring") else "../eluring"
hpc_root     <- if ("--hpc-root" %in% args) opt_val("--hpc-root") else "/gpfs/helios/home/taavi74/Projects"
check_files  <- "--check-files" %in% args
sample_filter  <- read_list(opt_val("--samples"))
project_filter <- read_list(opt_val("--projects"))
se_override    <- read_list(opt_val("--single-end"))

if (length(sample_filter) == 0 && length(project_filter) == 0) {
  stop("Provide at least one of --samples <file> or --projects p1,p2,...")
}

# Project -> cohort label. hpc_dir comes from the data_paths.csv `folder`
# column (authoritative); only the cohort tag lives here.
cohort_of <- c(
  sajased = "aging", wu_rampelli = "aging", chulenbayeva = "aging",
  curated_metagenomic = "aging", samples_65_74_new = "aging",
  samples_12_18_new = "aging", samples_middle_fill = "aging", dengl2025 = "aging",
  extraves = "newborn-term", newborn = "newborn-term"
)

data_paths_file <- file.path(eluring_root, "data/data_paths.csv")
if (!file.exists(data_paths_file)) {
  stop(sprintf("data_paths.csv not found at %s — pass --eluring <path>", data_paths_file))
}

dp <- read_csv(data_paths_file, show_col_types = FALSE) |>
  mutate(sample = trimws(sample), project = trimws(project), folder = trimws(folder))

# Warn on requested ids/projects that aren't in the manifest.
if (length(sample_filter)) {
  miss <- setdiff(sample_filter, dp$sample)
  if (length(miss)) warning(sprintf("--samples: %d id(s) not in data_paths.csv: %s",
                                     length(miss), paste(head(miss, 10), collapse = ", ")))
}
if (length(project_filter)) {
  missp <- setdiff(project_filter, dp$project)
  if (length(missp)) stop(sprintf("--projects: unknown project(s): %s. Known: %s",
                                   paste(missp, collapse = ", "),
                                   paste(sort(unique(dp$project)), collapse = ", ")))
}

samples <- dp
if (length(sample_filter))  samples <- samples |> filter(sample %in% sample_filter)
if (length(project_filter)) samples <- samples |> filter(project %in% project_filter)

if (nrow(samples) == 0) stop("No samples matched the given --samples / --projects filters.")

samples <- samples |>
  mutate(
    cohort    = unname(ifelse(project %in% names(cohort_of), cohort_of[project], NA_character_)),
    reads_dir = file.path(hpc_root, folder, "results/processed_reads"),
    contigs   = file.path(hpc_root, folder, "results/assembly/contigs",
                          paste0("MEGAHIT_", sample, "_contigs.fa.gz")),
    reads_pe1 = file.path(reads_dir, paste0(sample, "_processed_1.fq.gz")),
    reads_pe2 = file.path(reads_dir, paste0(sample, "_processed_2.fq.gz")),
    reads_se  = file.path(reads_dir, paste0(sample, "_processed.fq.gz")),
    bins_dir  = file.path(hpc_root, folder, "results/binrefine/final_bins"),
    single_end = sample %in% se_override
  )

if (any(is.na(samples$cohort))) {
  warning(sprintf("%d sample(s) from unmapped project(s) have cohort=NA: %s",
                  sum(is.na(samples$cohort)),
                  paste(sort(unique(samples$project[is.na(samples$cohort)])), collapse = ", ")))
}

if (check_files) {
  flipped <- 0L; missing <- character(0)
  for (i in seq_len(nrow(samples))) {
    if (samples$single_end[i]) next
    if (file.exists(samples$reads_pe1[i])) next
    if (file.exists(samples$reads_se[i])) { samples$single_end[i] <- TRUE; flipped <- flipped + 1L }
    else missing <- c(missing, samples$sample[i])
  }
  if (flipped) message(sprintf("--check-files: detected %d single-end samples", flipped))
  if (length(missing)) warning(sprintf("--check-files: %d samples have neither paired nor single-end reads on disk: %s",
                                        length(missing), paste(missing, collapse = ", ")))
}

samples <- samples |>
  mutate(
    reads_1 = if_else(single_end, reads_se, reads_pe1),
    reads_2 = if_else(single_end, NA_character_, reads_pe2)
  ) |>
  arrange(project, sample) |>
  select(sample, project, cohort, contigs, reads_1, reads_2, bins_dir, single_end)

unknown_se <- setdiff(se_override, samples$sample)
if (length(unknown_se)) warning(sprintf("--single-end: %d id(s) not in samplesheet: %s",
                                         length(unknown_se), paste(unknown_se, collapse = ", ")))

n_se <- sum(samples$single_end)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write_csv(samples, output_path, na = "")
message(sprintf("Wrote %d samples (%d paired, %d single-end; %d projects) to %s",
                nrow(samples), nrow(samples) - n_se, n_se,
                dplyr::n_distinct(samples$project), output_path))
if (!check_files && n_se == 0) {
  message("Note: all rows marked paired-end. Re-run with --check-files on the HPC ",
          "(where read paths are reachable) to auto-detect single-end samples before launching.")
}
