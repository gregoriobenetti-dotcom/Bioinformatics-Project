############################################################
# task8.R — Multi-column operations per group
# Data: bulk_counts_long.csv, sample_metadata.csv
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
counts <- fread("~/workspace/bulk_counts_long.csv")
meta   <- fread("~/workspace/sample_metadata.csv")

# ---- Standardize column names ----
setnames(counts, tolower(trimws(names(counts))))
setnames(meta,   tolower(trimws(names(meta))))

# Expect: counts has gene, sample_id, count
#         meta has sample_id, condition

# ---- Join metadata to counts ----
setkey(meta, sample_id)
dt <- meta[counts, on = "sample_id"]

# ---- Compute summary stats by gene and condition ----
# Q1 = 25th percentile, Q3 = 75th percentile
summ_stats <- dt[, .(
  mean   = mean(count),
  median = median(count),
  q1     = quantile(count, 0.25),
  q3     = quantile(count, 0.75)
), by = .(gene, condition)]

# ---- Reshape: wide form for treated/control comparison ----
wide <- dcast(summ_stats, gene ~ condition, value.var = "mean")

# ---- Filter: treated mean ≥ 2× control mean ----
# (use na.rm to handle genes missing in one condition)
filtered_genes <- wide[treated >= 2 * control, gene]

# ---- Final results (stats only for selected genes) ----
result <- summ_stats[gene %in% filtered_genes]

# ---- Output preview ----
print(head(result))
cat("\n✅ Genes where treated mean ≥ 2× control mean:", length(filtered_genes), "\n")
