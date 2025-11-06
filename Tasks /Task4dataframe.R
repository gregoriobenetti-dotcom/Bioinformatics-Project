############################################################
# task4_df.R — Annotate counts with sample/patient info
# Data: sample_metadata.csv, bulk_counts_long.csv
# Base R (data.frame) version using fread() for input
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
meta   <- as.data.frame(fread("~/workspace/sample_metadata.csv"))
counts <- as.data.frame(fread("~/workspace/bulk_counts_long.csv"))

# ---- Standardize column names ----
names(meta)   <- tolower(trimws(names(meta)))
names(counts) <- tolower(trimws(names(counts)))

# Expect:
# meta: sample_id, patient_id, condition
# counts: gene, sample_id, count

# ---- 1️⃣ Join sample metadata ----
joined <- merge(counts, meta, by = "sample_id", all.x = TRUE)

# ---- 2️⃣ Per-patient total counts ----
patient_totals <- aggregate(count ~ patient_id, data = joined, sum)
names(patient_totals)[2] <- "total_count"

# ---- 3️⃣ Top 10 genes by average count within each condition ----
# Compute mean count per gene × condition
avg_counts <- aggregate(count ~ condition + gene, data = joined, mean)
names(avg_counts)[3] <- "avg_count"

# Order by condition and descending avg_count
avg_counts <- avg_counts[order(avg_counts$condition, -avg_counts$avg_count), ]

# Take top 10 genes per condition
top_genes <- do.call(rbind, lapply(split(avg_counts, avg_counts$condition), function(x) head(x, 10)))

# ---- 4️⃣ Output preview ----
cat("\n✅ Per-patient total counts:\n")
print(head(patient_totals))

cat("\n✅ Top 10 genes by average count within each condition:\n")
print(head(top_genes))
