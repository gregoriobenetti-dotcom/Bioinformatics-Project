############################################################
# task12_df.R — Combine cohorts safely and compute variability
# Data: cohortA_samples.csv, cohortB_samples.csv, bulk_counts_long.csv
# Base R (data.frame) version using fread() for input
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
A <- as.data.frame(fread("~/workspace/cohortA_samples.csv"))
B <- as.data.frame(fread("~/workspace/cohortB_samples.csv"))
counts <- as.data.frame(fread("~/workspace/bulk_counts_long.csv"))

# ---- Standardize column names ----
names(A) <- tolower(trimws(names(A)))
names(B) <- tolower(trimws(names(B)))
names(counts) <- tolower(trimws(names(counts)))

# Expect:
# cohorts: cohort, sample_id, condition, patient_id
# counts:  gene, sample_id, count

# ---- 1️⃣ Combine cohorts safely ----
combined_meta <- rbind(A, B)

# ---- 2️⃣ Ensure key columns exist ----
required_cols <- c("sample_id", "condition", "cohort")
stopifnot(all(required_cols %in% names(combined_meta)))

# ---- 3️⃣ Order by cohort, condition, sample_id ----
combined_meta <- combined_meta[order(combined_meta$cohort, combined_meta$condition, combined_meta$sample_id), ]

# ---- 4️⃣ Join back to counts ----
dt <- merge(counts, combined_meta, by = "sample_id", all.x = TRUE)

# ---- 5️⃣ Compute variance per gene ----
gene_var <- aggregate(count ~ gene, data = dt, var)
names(gene_var)[2] <- "var_count"

# ---- 6️⃣ Identify top 100 most variable genes ----
top100_genes <- head(gene_var[order(-gene_var$var_count), "gene"], 100)

# ---- 7️⃣ Filter for top genes ----
dt_top <- subset(dt, gene %in% top100_genes)

# ---- 8️⃣ Compute mean counts per cohort and condition ----
summary_stats <- aggregate(count ~ cohort + condition + gene, data = dt_top, mean)
names(summary_stats)[4] <- "mean_count"

# ---- Output preview ----
cat("\n✅ Combined cohorts and summarized top 100 most variable genes:\n")
print(head(summary_stats))
