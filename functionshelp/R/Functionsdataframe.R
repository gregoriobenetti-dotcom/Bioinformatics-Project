suppressPackageStartupMessages(library(data.table))

#' Task 1 (data.frame): Filter, summarize, and group bulk counts
#' @return A list with per-gene and per-condition summaries.
task1_df_filter_counts <- function() {
  suppressPackageStartupMessages(library(data.table))

  counts <- as.data.frame(fread("~/workspace/bulk_counts_long.csv"))
  meta   <- as.data.frame(fread("~/workspace/sample_metadata.csv"))

  names(counts) <- tolower(trimws(names(counts)))
  names(meta)   <- tolower(trimws(names(meta)))

  # Filter genes and join metadata
  filtered <- subset(counts, grepl("^GENE_00", gene))
  filtered <- merge(filtered, meta, by = "sample_id", all.x = TRUE)
  filtered <- subset(filtered, condition == "treated")

  # Compute mean and median counts by gene
  mean_counts <- aggregate(count ~ gene, data = filtered, FUN = mean)
  median_counts <- aggregate(count ~ gene, data = filtered, FUN = median)

  # Merge them together
  summary_gene <- merge(mean_counts, median_counts, by = "gene", suffixes = c("_mean", "_median"))
  names(summary_gene) <- c("gene", "mean_count", "median_count")

  # Compute per-condition mean counts
  merged <- merge(counts, meta, by = "sample_id", all.x = TRUE)
  per_condition_summary <- aggregate(count ~ gene + condition, data = merged, FUN = mean)
  names(per_condition_summary)[3] <- "mean_count"

  list(
    summary_gene = summary_gene,
    per_condition_summary = per_condition_summary
  )
}


#' Task 2 (data.frame): Add QC-style derived columns
#' @return Data.frame with log2 and high flag columns.
task2_df_add_qc_flags <- function() {
  counts <- as.data.frame(fread("~/workspace/bulk_counts_long.csv"))
  names(counts) <- tolower(trimws(names(counts)))

  counts$log2_count <- log2(counts$count + 1)
  gene_medians <- aggregate(count ~ gene, data = counts, median)
  names(gene_medians)[2] <- "median_count"
  counts <- merge(counts, gene_medians, by = "gene", all.x = TRUE)
  counts$high <- counts$count > counts$median_count
  counts
}

#' Task 3 (data.frame): Join metadata and counts
#' @return Merged data.frame of counts and metadata.
task3_df_join_and_index <- function() {
  meta   <- as.data.frame(fread("~/workspace/sample_metadata.csv"))
  counts <- as.data.frame(fread("~/workspace/bulk_counts_long.csv"))
  names(meta) <- tolower(trimws(names(meta)))
  names(counts) <- tolower(trimws(names(counts)))
  merge(counts, meta, by = "sample_id", all.x = TRUE)
}

#' Task 4 (data.frame): Annotate counts with sample/patient info
#' @return A list of per-patient totals and top genes per condition.
task4_df_annotate_counts <- function() {
  meta   <- as.data.frame(fread("~/workspace/sample_metadata.csv"))
  counts <- as.data.frame(fread("~/workspace/bulk_counts_long.csv"))
  names(meta) <- tolower(trimws(names(meta)))
  names(counts) <- tolower(trimws(names(counts)))

  joined <- merge(counts, meta, by = "sample_id", all.x = TRUE)
  patient_totals <- aggregate(count ~ patient_id, data = joined, sum)
  names(patient_totals)[2] <- "total_count"

  avg_counts <- aggregate(count ~ condition + gene, data = joined, mean)
  names(avg_counts)[3] <- "avg_count"
  avg_counts <- avg_counts[order(avg_counts$condition, -avg_counts$avg_count), ]
  top_genes <- do.call(rbind, lapply(split(avg_counts, avg_counts$condition), function(x) head(x, 10)))

  list(patient_totals = patient_totals, top_genes = top_genes)
}

#' Task 5 (data.frame): Classify labs as normal or out_of_range
#' @return A list with abnormal rates by patient and lab.
task5_df_classify_labs <- function() {
  labs   <- as.data.frame(fread("~/workspace/clinical_labs.csv"))
  ranges <- as.data.frame(fread("~/workspace/lab_reference_ranges.csv"))
  names(labs) <- tolower(trimws(names(labs)))
  names(ranges) <- tolower(trimws(names(ranges)))

  joined <- merge(labs, ranges, by = "lab", all.x = TRUE)
  joined$status <- ifelse(joined$value >= joined$lower & joined$value <= joined$upper,
                          "normal", "out_of_range")

  abnormal_rates_patient <- aggregate(status == "out_of_range" ~ patient_id, data = joined, mean)
  names(abnormal_rates_patient)[2] <- "abnormal_rate"
  abnormal_rates_lab <- aggregate(status == "out_of_range" ~ lab, data = joined, mean)
  names(abnormal_rates_lab)[2] <- "abnormal_rate"

  list(by_patient = abnormal_rates_patient, by_lab = abnormal_rates_lab)
}

#' Task 6 (data.frame): Match vitals to nearest lab times
#' @return Merged data.frame with HR/SBP and correlation summary.
task6_df_match_labs_vitals <- function() {
  labs   <- as.data.frame(fread("~/workspace/clinical_labs.csv"))
  vitals <- as.data.frame(fread("~/workspace/vitals_time_series.csv"))
  names(labs)   <- tolower(trimws(names(labs)))
  names(vitals) <- tolower(trimws(names(vitals)))
  if ("time_iso" %in% names(labs))   names(labs)[names(labs) == "time_iso"] <- "time"
  if ("time_iso" %in% names(vitals)) names(vitals)[names(vitals) == "time_iso"] <- "time"
  labs$time   <- as.POSIXct(labs$time, tz = "UTC")
  vitals$time <- as.POSIXct(vitals$time, tz = "UTC")

  vitals_sub <- subset(vitals, vital %in% c("HR", "SBP"))
  hr  <- subset(vitals_sub, vital == "HR",  select = c("patient_id", "time", "value"))
  names(hr)[3] <- "HR"
  sbp <- subset(vitals_sub, vital == "SBP", select = c("patient_id", "time", "value"))
  names(sbp)[3] <- "SBP"

  nearest_match <- function(lab_df, vital_df, vital_name) {
    res <- lab_df
    res[[vital_name]] <- NA
    for (i in seq_len(nrow(res))) {
      pid <- res$patient_id[i]; t <- res$time[i]
      sub_vital <- subset(vital_df, patient_id == pid)
      if (nrow(sub_vital) > 0) {
        idx <- which.min(abs(difftime(sub_vital$time, t, units = "mins")))
        res[[vital_name]][i] <- sub_vital[[vital_name]][idx]
      }
    }
    res
  }

  labs_hr  <- nearest_match(labs, hr, "HR")
  labs_sbp <- nearest_match(labs, sbp, "SBP")
  merge(labs_hr, labs_sbp[, c("patient_id", "time", "SBP")],
        by = c("patient_id", "time"), all.x = TRUE)
}

#' Task 7 (data.frame): Slice genomics windows
#' @return Top 50 peaks in chr2 between 2–4 Mb.
task7_df_filter_genomic_peaks <- function() {
  peaks <- as.data.frame(fread("~/workspace/atac_peaks.bed.csv"))
  names(peaks) <- tolower(trimws(names(peaks)))
  filtered <- subset(peaks, chr == "chr2" & start >= 2e6 & start <= 4e6)
  filtered <- filtered[order(-filtered$score), ]
  head(filtered, 50)
}

#' Task 8 (data.frame): Multi-column stats per group
#' @return Data.frame of summary stats for selected genes.
task8_df_summarize_gene_stats <- function() {
  counts <- as.data.frame(fread("~/workspace/bulk_counts_long.csv"))
  meta   <- as.data.frame(fread("~/workspace/sample_metadata.csv"))
  names(counts) <- tolower(trimws(names(counts)))
  names(meta)   <- tolower(trimws(names(meta)))

  dt <- merge(counts, meta, by = "sample_id", all.x = TRUE)
  summ_stats <- aggregate(count ~ gene + condition, data = dt,
                          FUN = function(x) c(mean = mean(x), median = median(x),
                                              q1 = quantile(x, 0.25), q3 = quantile(x, 0.75)))
  summ_stats_df <- data.frame(
    gene = summ_stats$gene, condition = summ_stats$condition,
    mean = sapply(summ_stats$count, `[`, 1),
    median = sapply(summ_stats$count, `[`, 2)
  )

  summary_wide <- reshape(summ_stats_df[, c("gene", "condition", "mean")],
                          timevar = "condition", idvar = "gene", direction = "wide")
  names(summary_wide) <- gsub("mean\\.", "", names(summary_wide))
  summary_wide[is.na(summary_wide)] <- 0
  filtered_genes <- subset(summary_wide, treated >= 2 * control)$gene
  subset(summ_stats_df, gene %in% filtered_genes)
}

#' Task 9 (data.frame): Wide → long → wide reshape
#' @return Data.frame reshaped for plotting.
task9_df_reshape_counts <- function() {
  wide <- as.data.frame(fread("~/workspace/bulk_counts_wide.csv"))
  names(wide) <- tolower(trimws(names(wide)))

  long <- reshape(wide, varying = names(wide)[-1], v.names = "count",
                  timevar = "sample_id", times = names(wide)[-1],
                  idvar = "gene", direction = "long")
  row.names(long) <- NULL

  sample_totals <- aggregate(count ~ sample_id, data = long, sum)
  names(sample_totals)[2] <- "total_count"
  long <- merge(long, sample_totals, by = "sample_id", all.x = TRUE)
  long$condition <- ifelse(grepl("treated", long$sample_id, ignore.case = TRUE),
                           "treated", "control")
  summary_df <- aggregate(count ~ gene + condition, data = long, mean)
  final <- reshape(summary_df, timevar = "condition", idvar = "gene", direction = "wide")
  names(final) <- gsub("count\\.", "", names(final))
  final
}

#' Task 10 (data.frame): ATAC peaks → genes mapping
#' @return Top 20 genes by total overlap.
task10_df_map_peaks_to_genes <- function() {
  peaks <- as.data.frame(fread("~/workspace/atac_peaks.bed.csv"))
  genes <- as.data.frame(fread("~/workspace/gene_annotation.bed.csv"))
  names(peaks) <- tolower(trimws(names(peaks)))
  names(genes) <- tolower(trimws(names(genes)))

  overlaps <- do.call(rbind, lapply(split(peaks, peaks$chr), function(pk_chr) {
    gn_chr <- subset(genes, chr == unique(pk_chr$chr))
    res <- list()
    for (i in seq_len(nrow(pk_chr))) {
      p <- pk_chr[i, ]
      overlap_idx <- which(gn_chr$start <= p$end & gn_chr$end >= p$start)
      if (length(overlap_idx) > 0) {
        tmp <- gn_chr[overlap_idx, ]
        tmp$peak_start <- p$start; tmp$peak_end <- p$end; tmp$peak_score <- p$score
        res[[length(res) + 1]] <- tmp
      }
    }
    if (length(res) > 0) do.call(rbind, res) else NULL
  }))

  overlaps$overlap_bp <- pmax(0, pmin(overlaps$end, overlaps$peak_end) -
                                pmax(overlaps$start, overlaps$peak_start))
  peaks_per_gene <- aggregate(peak_start ~ gene, data = overlaps, FUN = length)
  names(peaks_per_gene)[2] <- "n_peaks"
  bp_per_gene <- aggregate(overlap_bp ~ gene, data = overlaps, sum)
  names(bp_per_gene)[2] <- "total_overlap_bp"
  summary_df <- merge(peaks_per_gene, bp_per_gene, by = "gene")
  head(summary_df[order(-summary_df$total_overlap_bp), ], 20)
}

#' Task 11 (data.frame): SNPs → genes mapping
#' @return Data.frame of high-impact variants per gene and sample.
task11_df_map_snps <- function() {
  variants <- as.data.frame(fread("~/workspace/variants.csv"))
  genes    <- as.data.frame(fread("~/workspace/gene_annotation.bed.csv"))
  names(variants) <- tolower(trimws(names(variants)))
  names(genes)    <- tolower(trimws(names(genes)))
  variants$start <- variants$pos; variants$end <- variants$pos

  overlaps <- do.call(rbind, lapply(split(variants, variants$chr), function(v_chr) {
    g_chr <- subset(genes, chr == unique(v_chr$chr))
    res <- list()
    for (i in seq_len(nrow(v_chr))) {
      v <- v_chr[i, ]
      overlap_idx <- which(g_chr$start <= v$end & g_chr$end >= v$start)
      if (length(overlap_idx) > 0) {
        tmp <- g_chr[overlap_idx, ]
        tmp$sample_id <- v$sample_id; tmp$impact <- v$impact
        res[[length(res) + 1]] <- tmp
      }
    }
    if (length(res) > 0) do.call(rbind, res) else NULL
  }))

  high_imp <- subset(overlaps, tolower(impact) == "high")
  aggregate(impact ~ gene + sample_id, data = high_imp, FUN = length)
}

#' Task 12 (data.frame): Combine cohorts and summarize
#' @return Summary of top 100 variable genes per cohort/condition.
task12_df_combine_cohorts <- function() {
  A <- as.data.frame(fread("~/workspace/cohortA_samples.csv"))
  B <- as.data.frame(fread("~/workspace/cohortB_samples.csv"))
  counts <- as.data.frame(fread("~/workspace/bulk_counts_long.csv"))
  names(A) <- tolower(trimws(names(A)))
  names(B) <- tolower(trimws(names(B)))
  names(counts) <- tolower(trimws(names(counts)))

  combined_meta <- rbind(A, B)
  dt <- merge(counts, combined_meta, by = "sample_id", all.x = TRUE)
  gene_var <- aggregate(count ~ gene, data = dt, var)
  names(gene_var)[2] <- "var_count"
  top100_genes <- head(gene_var[order(-gene_var$var_count), "gene"], 100)
  dt_top <- subset(dt, gene %in% top100_genes)
  aggregate(count ~ cohort + condition + gene, data = dt_top, mean)
}
