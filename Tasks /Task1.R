############################################################
# task1.R — Filter, summarize, and group bulk counts
# Data: bulk_counts_long.csv, sample_metadata.CSV
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
counts <- fread("~/workspace/bulk_counts_long.csv")
meta   <- fread("~/workspace/sample_metadata.csv")

# ---- Standardize column names ----
setnames(counts, tolower(trimws(names(counts))))
setnames(meta,   tolower(trimws(names(meta))))

# Expect:
# counts: gene, sample_id, count
# meta:   sample_id, condition

# ---- 1️⃣ Filter: keep only treated samples and genes starting with "GENE_00" ----
filtered <- counts[grepl("^GENE_00", gene)]

# Join to get condition info
setkey(meta, sample_id)
filtered <- meta[filtered, on = "sample_id"]

filtered <- filtered[condition == "treated"]

# ---- 2️⃣ Compute mean and median count by gene ----
summary_gene <- filtered[, .(
  mean_count   = mean(count),
  median_count = median(count)
), by = gene]

# ---- 3️⃣ Join metadata and compute per-condition mean counts by gene (in one pipeline) ----
per_condition_summary <- meta[counts, on = "sample_id"
][, .(mean_count = mean(count)), by = .(gene, condition)]

