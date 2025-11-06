############################################################
# task3_df.R — Join metadata and counts (base R version)
# Data: sample_metadata.csv, bulk_counts_long.csv
# Using fread() for input, base R for operations
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
meta   <- as.data.frame(fread("~/workspace/sample_metadata.csv"))
counts <- as.data.frame(fread("~/workspace/bulk_counts_long.csv"))

# ---- Standardize column names ----
names(meta)   <- tolower(trimws(names(meta)))
names(counts) <- tolower(trimws(names(counts)))

# Expect:
# meta: sample_id, condition, patient_id
# counts: gene, sample_id, count

# ---- 1️⃣ Merge metadata with counts ----
joined <- merge(counts, meta, by = "sample_id", all.x = TRUE)

# ---- 2️⃣ Define a simple lookup (simulating index-based query) ----
query <- function(df, gene_name, sample_name) {
  subset(df, gene == gene_name & sample_id == sample_name)
}

# ---- 3️⃣ Benchmark the lookup ----
cat("\nBenchmarking query performance (base R)...\n")
system.time(query(joined, "GENE_005", "SAMPLE_01"))
system.time(query(joined, "GENE_005", "SAMPLE_01"))

# ---- 4️⃣ Output preview ----
cat("\n✅ Joined metadata (head):\n")
print(head(joined))
