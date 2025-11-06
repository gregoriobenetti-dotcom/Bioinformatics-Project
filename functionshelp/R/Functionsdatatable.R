############################################################
# datatable_functions.R
# All project tasks (1–12 + final) — data.table self-contained versions
############################################################

#' Task 1 — Filter and summarize bulk counts
#'
#' Loads bulk RNA counts and metadata, filters treated samples,
#' and computes mean/median per gene and per condition.
#' @return A list of two data.tables: per-gene and per-condition summaries.
#' @export
task1_filter_counts <- function() {
  suppressPackageStartupMessages(library(data.table))

  # ---- Load data ----
  counts <- fread("~/workspace/bulk_counts_long.csv")
  meta   <- fread("~/workspace/sample_metadata.csv")

  # ---- Standardize column names ----
  setnames(counts, tolower(trimws(names(counts))))
  setnames(meta,   tolower(trimws(names(meta))))
  #trim spaces and lowercase all column names

  # ---- Filter genes ----
  filtered <- counts[grepl("^GENE_00", gene)]
  #keeps genes whose name starts with GENE_00

  # ---- Merge metadata to add condition info ----
  filtered <- merge(filtered, meta, by = "sample_id", all.x = TRUE)

  # Keep treated only
  filtered <- filtered[condition == "treated"]

  # ---- Mean and median count by gene ----
  summary_gene <- filtered[, .(
    mean_count = mean(count),
    median_count = median(count)
  ), by = gene]
  #in summary gene creates the median and mean of the counts

  # ---- Per-condition mean count ----
  merged <- merge(counts, meta, by = "sample_id", all.x = TRUE)
  per_condition_summary <- merged[, .(
    mean_count = mean(count)
  ), by = .(gene, condition)]

  # ---- Return ----
  list(
    filtered = filtered, #filtered treated samples
    summary_gene = summary_gene, #mean/median per gene
    per_condition_summary = per_condition_summary # mean count per gene and condition
  )
}

#' Task 2 — Add QC-style columns
#'
#' Loads bulk counts, adds log2-transformed counts and flags high values.
#' @return A data.table with log2_count and high flag columns.
#' @export
task2_add_qc_flags <- function() {

  suppressPackageStartupMessages(library(data.table))

  # ---- Load Data ----
  counts <- fread("~/workspace/bulk_counts_long.csv")

  # ---- Trim spaces and lowercase all names (standard column names) ----
  setnames(counts, tolower(trimws(names(counts))))

  counts[, log2_count := log2(count + 1)] # new log2 count column
  counts[, high := count > median(count), by = gene] # flag 'high' for counts above the median

  return(counts) # table with extra columns

}


#' Task 3 — Join metadata and index
#'
#' Joins metadata into counts and benchmarks indexing efficiency.
#' @return A joined data.table with sample metadata.
#' @export

task3_join_and_index <- function() {

  suppressPackageStartupMessages(library(data.table))

  # ---- Load data ----
  meta   <- fread("~/workspace/sample_metadata.csv")
  counts <- fread("~/workspace/bulk_counts_long.csv")

  # ---- Trims and lower case ----
  setnames(meta, tolower(trimws(names(meta))))
  setnames(counts, tolower(trimws(names(counts))))

  # ---- Joins metadata and count data by sample_id ----
  setkey(meta, sample_id)
  joined <- meta[counts, on = "sample_id"]

  # ---- Pick one gene + sample to test lookup ----
  gene_to_test <- joined$gene[1]
  sample_to_test <- joined$sample_id[1]

  # ----Benchmark before index ----
  time_before <- system.time({
    res_before <- joined[gene == gene_to_test & sample_id == sample_to_test]
  })

  # ---- Add secondary index ----
  setindex(joined, gene, sample_id)

  # ---- Benchmark after index ----
  time_after <- system.time({
    res_after <- joined[gene == gene_to_test & sample_id == sample_to_test]
  })

  # ---- Return only small useful output ----
  list(
    preview = head(joined, 10),
    benchmark = data.frame(
      stage = c("before_index", "after_index"),
      user_time = c(time_before["user.self"], time_after["user.self"])
    )
  )
}

#' Task 4 — Annotate counts with patient info
#'
#' Merges metadata into counts, computes totals and top genes.
#' @return A list of two data.tables: patient_totals and top_genes.
#' @export
task4_annotate_counts <- function() {

  suppressPackageStartupMessages(library(data.table))

  # ---- Load Data ----
  meta   <- fread("~/workspace/sample_metadata.csv")
  counts <- fread("~/workspace/bulk_counts_long.csv")

  # ---- Trims & Lower Case ----
  setnames(meta, tolower(trimws(names(meta))))
  setnames(counts, tolower(trimws(names(counts))))

  # ---- Adds metadata columns to counts ----
  setkey(meta, sample_id)
  joined <- meta[counts, on = "sample_id"] #merges meta onto counts

  patient_totals <- joined[, .(total_count = sum(count)), by = patient_id] #sums count per patient

  top_genes <- joined[, .(avg_count = mean(count)), by = .(condition, gene)
  ][order(condition, -avg_count), head(.SD, 10), by = condition] #finds 10 top genes

  return(list(patient_totals = patient_totals, top_genes = top_genes)) # total counts per patient and top genes by average count
}

#' Task 5 — Classify labs as normal or out_of_range
#'
#' Classifies lab values based on reference ranges and computes abnormal rates.
#' @return A list with abnormal rates by patient and by lab.
#' @export
task5_classify_labs <- function() {

  suppressPackageStartupMessages(library(data.table))

  # ---- Load Data ----
  labs   <- fread("~/workspace/clinical_labs.csv")
  ranges <- fread("~/workspace/lab_reference_ranges.csv")

  # ---- Trim & Lower case ----
  setnames(labs,   tolower(trimws(names(labs))))
  setnames(ranges, tolower(trimws(names(ranges))))

  # ---- Row selection and classification ----
  setkey(ranges, lab)
  classified <- ranges[labs, on = .(lab, lower <= value, upper >= value), nomatch = 0L]
  classified[, status := ifelse(!is.na(lower), "normal", "out_of_range")]
  #joins lab results and reference and labels values

  out_all <- merge(labs, classified, all.x = TRUE)
  out_all[is.na(status), status := "out_of_range"]

  abnormal_rates_patient <- out_all[, .(abnormal_rate = mean(status == "out_of_range")), by = patient_id]
  abnormal_rates_lab     <- out_all[, .(abnormal_rate = mean(status == "out_of_range")), by = lab]
  #abnormal rate per patient and per lab test

  return(list(patient = abnormal_rates_patient, lab = abnormal_rates_lab))
}

#' Task 6 — Match labs to nearest vitals
#'
#' Finds the nearest heart rate and SBP values for each lab draw.
#' @return Correlation summary between CRP and vitals per patient.
#' @export
task6_match_labs_vitals <- function() {

  suppressPackageStartupMessages(library(data.table))

  # ---- Load Data ----
  labs   <- fread("~/workspace/clinical_labs.csv")
  vitals <- fread("~/workspace/vitals_time_series.csv")

  # ---- Name Adjustments ----
  setnames(labs, tolower(trimws(names(labs))))
  setnames(vitals, tolower(trimws(names(vitals))))

  setnames(labs, "time_iso", "time")
  setnames(vitals, "time_iso", "time")

  # ---- Makes timestamps comparable ----
  labs[, time := as.POSIXct(time, tz = "UTC")]
  vitals[, time := as.POSIXct(time, tz = "UTC")]

  vitals_sub <- vitals[vital %in% c("HR", "SBP")]
  hr  <- vitals_sub[vital == "HR",  .(patient_id, time, HR = value)]
  sbp <- vitals_sub[vital == "SBP", .(patient_id, time, SBP = value)]

  setkey(hr, patient_id, time)
  setkey(sbp, patient_id, time)
  setkey(labs, patient_id, time)

  labs_hr  <- hr[labs,  roll = "nearest"]
  labs_sbp <- sbp[labs, roll = "nearest"]

  merged <- merge(
    labs_hr[, .(patient_id, time, lab, value, HR, time_HR = time)],
    labs_sbp[, .(patient_id, time, SBP, time_SBP = time)],
    by = c("patient_id", "time"),
    all.x = TRUE
  )

  corr_summary <- merged[lab == "CRP",
                         .(cor_HR = cor(value, HR, use = "pairwise.complete.obs"),
                           cor_SBP = cor(value, SBP, use = "pairwise.complete.obs")),
                         by = patient_id]

  return(corr_summary)
}

#' Task 7 — Filter genomic peaks
#' @return Top 50 peaks on chr2 between 2–4 Mb by score.
#' @export
task7_filter_genomic_peaks <- function() {

  suppressPackageStartupMessages(library(data.table))

  # ---- Load Data ----
  peaks <- fread("~/workspace/atac_peaks.bed.csv")

  setnames(peaks, tolower(trimws(names(peaks))))

  # ---- Region Selection ----
  filtered <- peaks[chr == "chr2" & start >= 2e6 & start <= 4e6]

  # ---- Sort by Score ----
  setorder(filtered, -score)

  return(filtered[1:50])
}

#' Task 8 — Summarize gene stats
#' @return Data.table of summary statistics for genes.
#' @export
task8_summarize_gene_stats <- function() {

  suppressPackageStartupMessages(library(data.table))

  # ---- Load Data ----
  counts <- fread("~/workspace/bulk_counts_long.csv")
  meta   <- fread("~/workspace/sample_metadata.csv")

  # ---- Trim & Lowercase ----
  setnames(counts, tolower(trimws(names(counts))))
  setnames(meta,   tolower(trimws(names(meta))))

  setkey(meta, sample_id)
  dt <- meta[counts, on = "sample_id"]

  # ---- CAlculates mean median and quartiles of expression ----
  summ_stats <- dt[, .(
    mean   = mean(count),
    median = median(count),
    q1     = quantile(count, 0.25),
    q3     = quantile(count, 0.75)
  ), by = .(gene, condition)]

  # ---- Selects genes where treated mean >= 2x ----
  wide <- dcast(summ_stats, gene ~ condition, value.var = "mean")
  filtered_genes <- wide[treated >= 2 * control, gene] # Keeps genes with >= 2x mean in treated

  result <- summ_stats[gene %in% filtered_genes]
  return(result)
}

#' Task 9 — Reshape counts wide→long→wide
#' @return Data.table reshaped for plotting.
#' @export
task9_reshape_counts <- function() {

  suppressPackageStartupMessages(library(data.table))


  wide <- fread("~/workspace/bulk_counts_wide.csv")


  setnames(wide, tolower(trimws(names(wide))))

  # ---- Converts wide gene table ito long format ----
  long <- melt(wide, id.vars = "gene", variable.name = "sample_id", value.name = "count")

  # ---- Adds conditions labels ----
  long[, condition := fifelse(grepl("treated", sample_id, ignore.case = TRUE), "treated", "control")]

  # ---- Computes mean count per gene and condition ----
  summary_wide <- long[, .(mean_count = mean(count)), by = .(gene, condition)]

  # ---- Converts back to the wide format ----
  final <- dcast(summary_wide, gene ~ condition, value.var = "mean_count")

  return(final)
}

#' Task 10 — Map ATAC peaks to genes
#' @return Top 20 genes by total overlap length.
#' @export
task10_map_peaks_to_genes <- function() {

  suppressPackageStartupMessages(library(data.table))

  # ---- Load Data ----
  peaks <- fread("~/workspace/atac_peaks.bed.csv")
  genes <- fread("~/workspace/gene_annotation.bed.csv")

  # ---- Trim & Lowercase ----
  setnames(peaks, tolower(trimws(names(peaks))))
  setnames(genes, tolower(trimws(names(genes))))

  setkey(peaks, chr, start, end)
  setkey(genes, chr, start, end)

  # ---- Finds Overlaps ----
  overlaps <- foverlaps(peaks, genes, by.x = c("chr", "start", "end"),
                        by.y = c("chr", "start", "end"), type = "any", nomatch = 0L)

  # ---- Calculates total overlap length ----
  overlaps[, overlap_bp := pmax(0, pmin(end, i.end) - pmax(start, i.start))]


  summary <- overlaps[, .(total_overlap_bp = sum(overlap_bp)), by = gene]

  # ---- Selects top 20 genes ----
  setorder(summary, -total_overlap_bp)

  return(summary[1:20])
}

#' Task 11 — Map variants to genes
#' @return List of variant counts and genes with high-impact variants.
#' @export
task11_map_snps <- function() {

  suppressPackageStartupMessages(library(data.table))

  # ---- Load Data ----
  variants <- fread("~/workspace/variants.csv")
  genes    <- fread("~/workspace/gene_annotation.bed.csv")


  setnames(variants, tolower(trimws(names(variants))))
  setnames(genes, tolower(trimws(names(genes))))


  variants[, start := pos]
  variants[, end := pos]


  setkey(variants, chr, start, end)
  setkey(genes, chr, start, end)

  # ---- Finds overlaps between genes and variants ----
  overlaps <- foverlaps(variants, genes, by.x = c("chr", "start", "end"),
                        by.y = c("chr", "start", "end"), type = "any", nomatch = 0L)

  # ---- Finds high impact variants and counts them
  high_imp <- overlaps[tolower(impact) == "high"]
  counts_gene_sample <- high_imp[, .N, by = .(gene, sample_id)]

  # ---- Detects genes with high impact variants in all samples ----
  setnames(counts_gene_sample, "N", "high_variant_count")
  samples_with_high <- unique(high_imp[, .(gene, sample_id)])
  gene_sample_counts <- samples_with_high[, .N, by = gene]
  total_samples <- uniqueN(variants$sample_id)
  genes_all_samples <- gene_sample_counts[N == total_samples, gene]

  return(list(counts = counts_gene_sample, all_samples = genes_all_samples))
}

#' Task 12 — Combine cohorts
#' @return Summary of top variable genes by cohort and condition.
#' @export
task12_combine_cohorts <- function() {

  suppressPackageStartupMessages(library(data.table))

  # ---- Load Data ----
  A <- fread("~/workspace/cohortA_samples.csv")
  B <- fread("~/workspace/cohortB_samples.csv")
  counts <- fread("~/workspace/bulk_counts_long.csv")

  # ---- Name settings ----
  setnames(A, tolower(trimws(names(A))))
  setnames(B, tolower(trimws(names(B))))
  setnames(counts, tolower(trimws(names(counts))))

  # ---- Combines metadata from two cohorts ----
  combined_meta <- rbindlist(list(A, B), use.names = TRUE, fill = TRUE)
  setorder(combined_meta, cohort, condition, sample_id)

  # ---- Merges with count metadata ----
  setkey(combined_meta, sample_id)
  dt <- combined_meta[counts, on = "sample_id"]

  # ---- Calculates gene variance ----
  gene_var <- dt[, .(var_count = var(count)), by = gene]

  # ---- Selects top 100 most variable genes ----
  top100_genes <- gene_var[order(-var_count)][1:100, gene]
  dt_top <- dt[gene %in% top100_genes]

  # ---- Calculates statistics ----
  summary_stats <- dt_top[, .(mean_count = mean(count)),
                          by = .(cohort, condition, gene)]
  return(summary_stats)
}
