############################################################
# final_task_df.R — Integration cluster–cell type association
# Data: SeuratIntegration.csv + nt_combined_clustering.output.csv
# Approach: base R (data.frame)
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
integration <- as.data.frame(fread("~/workspace/annotated_GSM3516673_normal_annotated_GSM3516672_tumor_SeuratIntegration.csv"))
nt_cluster  <- as.data.frame(fread("~/workspace/nt_combined_clustering.output.csv"))

# ---- Standardize column names ----
names(integration) <- tolower(trimws(names(integration)))
names(nt_cluster)  <- tolower(trimws(names(nt_cluster)))

# Expect:
# integration: cell, integration_cluster
# nt_cluster:  cell, cell_type, sample_type

# ---- Merge by cell ----
merged <- merge(integration, nt_cluster, by = "cell", all.x = TRUE)

# ---- Summarize counts per integration_cluster × sample_type × cell_type ----
summary_counts <- aggregate(cell ~ integration_cluster + sample_type + cell_type,
                            data = merged, FUN = length)
names(summary_counts)[4] <- "n_cells"

# ---- Optional: compute proportions per cluster ----
cluster_totals <- aggregate(n_cells ~ integration_cluster, data = summary_counts, sum)
names(cluster_totals)[2] <- "cluster_total"
summary_counts <- merge(summary_counts, cluster_totals, by = "integration_cluster", all.x = TRUE)
summary_counts$prop <- summary_counts$n_cells / summary_counts$cluster_total

# ---- Sort results ----
summary_counts <- summary_counts[order(summary_counts$integration_cluster, -summary_counts$n_cells), ]

# ---- Output preview ----
cat("\n✅ Cluster–CellType associations:\n")
print(head(summary_counts))

