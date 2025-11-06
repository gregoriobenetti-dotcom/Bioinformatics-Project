############################################################
# task1_df.R — Filter, summarize, and group bulk counts
# Data: bulk_counts_long.csv, sample_metadata.csv
# Base R (data.frame) version — fixed to match data.table logic
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
counts <- as.data.frame(fread("~/workspace/bulk_counts_long.csv"))
meta   <- as.data.frame(fread("~/workspace/sample_metadata.csv"))

# ---- Standardize column names ----
names(counts) <- tolower(trimws(names(counts)))
names(meta)   <- tolower(trimws(names(meta)))

# Expect:
# counts: gene, sample_id, count
# meta:   sample_id, condition

# ---- 1️⃣ Filter: keep only treated samples and genes starting with "GENE_00" ----
filtered <- subset(counts, grepl("^GENE_00", gene))

# Join metadata (same direction as data.table: meta[filtered, on="sample_id"])
filtered <- merge(meta, filtered, by = "sample_id", all.y = TRUE)

# Keep only treated samples
filtered <- subset(filtered, condition == "treated")

# ---- 2️⃣ Compute mean and median count by gene ----
summary_gene <- aggregate(filtered$count,
                          by = list(gene = filtered$gene),
                          FUN = function(x) c(mean = mean(x), median = median(x)))
summary_gene <- data.frame(
  gene = summary_gene$gene,
  mean_count = summary_gene$x[, "mean"],
  median_count = summary_gene$x[, "median"]
)

# ---- 3️⃣ Join metadata and compute per-condition mean counts by gene ----
merged <- merge(meta, counts, by = "sample_id", all.y = TRUE)
per_condition_summary <- aggregate(merged$count,
                                   by = list(gene = merged$gene, condition = merged$condition),
                                   FUN = mean)
names(per_condition_summary)[3] <- "mean_count"

# ---- Output preview ----
cat("\n✅ Summary of treated GENE_00 genes:\n")
print(head(summary_gene))

cat("\n✅ Per-condition mean counts by gene:\n")
print(head(per_condition_summary))

