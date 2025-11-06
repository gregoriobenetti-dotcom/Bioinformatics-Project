############################################################
# task9_df.R — Go wide → long → wide for downstream plotting
# Data: bulk_counts_wide.csv
# Base R (data.frame) version using fread() for input
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
wide <- as.data.frame(fread("~/workspace/bulk_counts_wide.csv"))

# ---- Standardize column names ----
names(wide) <- tolower(trimws(names(wide)))

# Expect: first column = gene, remaining columns = samples
# Example: gene, sample_1, sample_2, ...

# ---- 1️⃣ Convert wide → long ----
# Use base R reshape()
long <- reshape(
  wide,
  varying = names(wide)[-1],
  v.names = "count",
  timevar = "sample_id",
  times = names(wide)[-1],
  idvar = "gene",
  direction = "long"
)
row.names(long) <- NULL

# ---- 2️⃣ Add per-sample totals ----
sample_totals <- aggregate(count ~ sample_id, data = long, sum)
names(sample_totals)[2] <- "total_count"

long <- merge(long, sample_totals, by = "sample_id", all.x = TRUE)

# ---- 3️⃣ Derive condition from sample name ----
# e.g., "treated_01", "control_03" → "treated", "control"
long$condition <- ifelse(grepl("treated", long$sample_id, ignore.case = TRUE),
                         "treated", "control")

# ---- 4️⃣ Compute mean counts by gene × condition ----
summary_df <- aggregate(count ~ gene + condition, data = long, mean)
names(summary_df)[3] <- "mean_count"

# ---- 5️⃣ Convert back to wide format ----
final <- reshape(
  summary_df,
  timevar = "condition",
  idvar = "gene",
  direction = "wide"
)
names(final) <- gsub("mean_count\\.", "", names(final))

# ---- Output preview ----
cat("\n✅ Reshaped successfully: wide → long → wide.\n")
print(head(final))
