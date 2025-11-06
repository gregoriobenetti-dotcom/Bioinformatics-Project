############################################################
# task2_df.R — Add QC-style derived columns (base R version)
# Data: bulk_counts_long.csv
# Uses fread() for input, base R for operations
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
counts <- as.data.frame(fread("~/workspace/bulk_counts_long.csv"))

# ---- Standardize column names ----
names(counts) <- tolower(trimws(names(counts)))

# Expect: gene, sample_id, count

# ---- 1️⃣ Add log2(count) column ----
# Add a small offset (+1) to avoid log(0) errors
counts$log2_count <- log2(counts$count + 1)

# ---- 2️⃣ Add a binary flag 'high' if count > 100 ----
counts$high <- counts$count > 100

# ---- 3️⃣ Overwrite 'high' using a gene-wise threshold ----
# For each gene, compute its median count and flag if above median
gene_medians <- aggregate(count ~ gene, data = counts, median)
names(gene_medians)[2] <- "median_count"

# Merge back to main table
counts <- merge(counts, gene_medians, by = "gene", all.x = TRUE)

# Update 'high' flag
counts$high <- counts$count > counts$median_count

# ---- Output preview ----
cat("\n✅ Added derived columns (log2_count and high):\n")
print(head(counts))
