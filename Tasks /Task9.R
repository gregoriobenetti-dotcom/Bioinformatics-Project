############################################################
# task9.R — Go wide → long → wide for downstream plotting
# Data: bulk_counts_wide.csv
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
wide <- fread("~/workspace/bulk_counts_wide.csv")

# ---- Standardize column names ----
setnames(wide, tolower(trimws(names(wide))))

# Expect: first column = gene, remaining columns = samples
# Example columns: gene, sample_1, sample_2, sample_3, ...

# ---- Convert wide → long ----
long <- melt(
  wide,
  id.vars = "gene",
  variable.name = "sample_id",
  value.name = "count"
)

# ---- Add per-sample totals ----
sample_totals <- long[, .(total_count = sum(count)), by = sample_id]
long <- merge(long, sample_totals, by = "sample_id")

# ---- Derive condition from sample name (if encoded) ----
# e.g., "treated_01", "control_03" → "treated", "control"
# Modify this regex if your sample names differ
long[, condition := fifelse(grepl("treated", sample_id, ignore.case = TRUE),
                            "treated", "control")]

# ---- Compute mean counts by gene × condition ----
summary_wide <- long[, .(mean_count = mean(count)), by = .(gene, condition)]

# ---- Convert back to wide format for plotting ----
final <- dcast(summary_wide, gene ~ condition, value.var = "mean_count")

# ---- Output preview ----
print(head(final))
cat("\n✅ Reshaped successfully: wide → long → wide.\n")
