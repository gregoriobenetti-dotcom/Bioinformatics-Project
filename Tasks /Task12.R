############################################################
# task12.R — Combine cohorts safely
# Data: cohortA_samples.csv, cohortB_samples.csv, bulk_counts_long.csv
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
A <- fread("~/workspace/cohortA_samples.csv")
B <- fread("~/workspace/cohortB_samples.csv")
counts <- fread("~/workspace/bulk_counts_long.csv")

# ---- Standardize column names ----
setnames(A, tolower(trimws(names(A))))
setnames(B, tolower(trimws(names(B))))
setnames(counts, tolower(trimws(names(counts))))

# Expect:
# cohorts: cohort, sample_id, condition, patient_id
# counts:  gene, sample_id, count

# ---- Combine cohorts safely ----
combined_meta <- rbindlist(list(A, B), use.names = TRUE, fill = TRUE)

# Verify alignment (optional sanity check)
stopifnot(all(c("sample_id", "condition", "cohort") %in% names(combined_meta)))

# ---- Order by cohort, condition, sample_id ----
setorder(combined_meta, cohort, condition, sample_id)

# ---- Join back to counts ----
setkey(combined_meta, sample_id)
dt <- combined_meta[counts, on = "sample_id"]

# ---- Compute gene variability ----
# Variance across all samples to pick top 100 most variable genes
gene_var <- dt[, .(var_count = var(count)), by = gene]
top100_genes <- gene_var[order(-var_count)][1:100, gene]

# ---- Filter for those genes ----
dt_top <- dt[gene %in% top100_genes]

# ---- Per-cohort, per-condition mean counts ----
summary_stats <- dt_top[, .(mean_count = mean(count)),
                        by = .(cohort, condition, gene)]

# ---- Output preview ----
print(head(summary_stats))
cat("\n✅ Combined cohorts and summarized top 100 most variable genes.\n")
