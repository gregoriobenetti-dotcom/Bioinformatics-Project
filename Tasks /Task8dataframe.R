############################################################
# task8_df.R — Multi-column operations per group
# Data: bulk_counts_long.csv, sample_metadata.csv
# Base R (data.frame) version using fread() for input
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
counts <- as.data.frame(fread("~/workspace/bulk_counts_long.csv"))
meta   <- as.data.frame(fread("~/workspace/sample_metadata.csv"))

# ---- Standardize column names ----
names(counts) <- tolower(trimws(names(counts)))
names(meta)   <- tolower(trimws(names(meta)))

# Expect: counts has gene, sample_id, count
#         meta has sample_id, condition

# ---- 1️⃣ Join metadata to counts ----
dt <- merge(counts, meta, by = "sample_id", all.x = TRUE)

# ---- 2️⃣ Compute summary stats by gene and condition ----
# Helper for mean/median/quantiles
summ_stats <- aggregate(
  count ~ gene + condition,
  data = dt,
  FUN = function(x) c(
    mean = mean(x),
    median = median(x),
    q1 = quantile(x, 0.25),
    q3 = quantile(x, 0.75)
  )
)

# Convert list columns to separate numeric columns
summ_stats_df <- data.frame(
  gene = summ_stats$gene,
  condition = summ_stats$condition,
  mean = sapply(summ_stats$count, function(x) x[1]),
  median = sapply(summ_stats$count, function(x) x[2]),
  q1 = sapply(summ_stats$count, function(x) x[3]),
  q3 = sapply(summ_stats$count, function(x) x[4])
)

# ---- 3️⃣ Reshape to wide format for treated/control comparison ----
summary_wide <- reshape(
  summ_stats_df[, c("gene", "condition", "mean")],
  timevar = "condition",
  idvar = "gene",
  direction = "wide"
)
# Rename columns for clarity
names(summary_wide) <- gsub("mean\\.", "", names(summary_wide))

# ---- 4️⃣ Filter: treated mean ≥ 2× control mean ----
summary_wide$treated[is.na(summary_wide$treated)] <- 0
summary_wide$control[is.na(summary_wide$control)] <- 0
filtered_genes <- subset(summary_wide, treated >= 2 * control)$gene

# ---- 5️⃣ Keep only selected genes ----
result <- subset(summ_stats_df, gene %in% filtered_genes)

# ---- Output preview ----
cat("\n✅ Genes where treated mean ≥ 2× control mean:\n")
print(head(result))
cat("\nNumber of selected genes:", length(filtered_genes), "\n")
