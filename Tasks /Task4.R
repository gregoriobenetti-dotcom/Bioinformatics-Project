# task4.R — Annotate counts with sample/patient info
suppressPackageStartupMessages(library(data.table))

# Load
meta  <- fread("~/workspace/sample_metadata.csv")
counts <- fread("~/workspace/bulk_counts_long.csv")

#--- data.table version ---
dt_meta <- copy(meta)
dt_counts <- copy(counts)

# Join sample metadata
setkey(dt_meta, sample_id)
joined <- dt_meta[dt_counts, on = "sample_id"]

# Per-patient total counts
patient_totals <- joined[, .(total_count = sum(count)), by = patient_id]

# Top 10 genes by average count within each condition
top_genes <- joined[, .(avg_count = mean(count)), by = .(condition, gene)
][order(condition, -avg_count), head(.SD, 10), by = condition]

#--- data.frame version ---
df_meta <- as.data.frame(meta)
df_counts <- as.data.frame(counts)

df_joined <- merge(df_counts, df_meta, by = "sample_id", all.x = TRUE)
df_patient_total <- aggregate(count ~ patient_id, df_joined, sum)
df_top_genes <- do.call(rbind, lapply(split(df_joined, df_joined$condition), function(x) {
  t <- aggregate(count ~ gene, x, mean)
  head(t[order(-t$count), ], 10)
}))

#--- Output small checks ---
print(head(patient_totals))
print(head(top_genes))

