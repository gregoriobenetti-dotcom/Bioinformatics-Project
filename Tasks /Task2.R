############################################################
# task2.R — Add QC-style derived columns without copying
# Data: bulk_counts_long.csv
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
counts <- fread("~/workspace/bulk_counts_long.csv")

# ---- Standardize column names ----
setnames(counts, tolower(trimws(names(counts))))
# Expect columns: gene, sample_id, count

# ---- 1️⃣ Add log2(count) column ----
# Add a small offset (e.g. +1) to avoid log(0) errors
counts[, log2_count := log2(count + 1)]

# ---- 2️⃣ Add a binary flag 'high' if count > 100 ----
counts[, high := count > 100]

# ---- 3️⃣ Overwrite 'high' to use a gene-wise threshold ----
# For each gene, compute its median count, then mark high if above that median
counts[, high := count > median(count), by = gene]

# ---- Output preview ----
cat("✅ Added derived columns:\n")
print(head(counts))

