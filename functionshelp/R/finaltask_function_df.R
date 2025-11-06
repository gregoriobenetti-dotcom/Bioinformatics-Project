#' final_task_df — Associate cell types to integration clusters (data.frame)
#'
#' This function merges Seurat integration results with cell type annotations
#' and summarizes the number and proportion of cells of each type
#' per integration cluster and tissue type (Normal/Tumor).
#'
#' @return data.frame with columns:
#' integration_cluster, sample_type, cell_type, n_cells, prop
#' @export
final_task_df <- function() {
  suppressPackageStartupMessages(library(data.table))

  # ---- Load data ----
  integration <- as.data.frame(fread("~/workspace/annotated_GSM3516673_normal_annotated_GSM3516672_tumor_SeuratIntegration.csv"))
  nt_cluster  <- as.data.frame(fread("~/workspace/nt_combined_clustering.output.csv"))

  # ---- Standardize column names ----
  names(integration) <- tolower(trimws(names(integration)))
  names(nt_cluster)  <- tolower(trimws(names(nt_cluster)))

  # ---- Merge ----

  integration$cell <- trimws(integration$cell)
  nt_cluster$cell  <- trimws(nt_cluster$cell)

  merged <- merge(integration, nt_cluster, by = "cell", all.x = TRUE)

  # ---- Summarize counts ----
  summary_counts <- aggregate(cell ~ integration_cluster + sample_type + cell_type,
                              data = merged, FUN = length)
  names(summary_counts)[4] <- "n_cells"

  # ---- Compute proportions per cluster ----
  cluster_totals <- aggregate(n_cells ~ integration_cluster, data = summary_counts, sum)
  names(cluster_totals)[2] <- "cluster_total"
  summary_counts <- merge(summary_counts, cluster_totals, by = "integration_cluster", all.x = TRUE)
  summary_counts$prop <- summary_counts$n_cells / summary_counts$cluster_total

  # ---- Sort for readability ----
  summary_counts <- summary_counts[order(summary_counts$integration_cluster, -summary_counts$n_cells), ]

  return(summary_counts)
}

