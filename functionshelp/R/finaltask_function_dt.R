#' final_task_dt — Associate cell types to integration clusters (data.table)
#'
#' This function merges Seurat integration results with cell type annotations
#' and summarizes the number and proportion of cells of each type
#' per integration cluster and tissue type (Normal/Tumor).
#'
#' @return data.table with columns:
#' integration_cluster, sample_type, cell_type, n_cells, prop
#' @export
final_task_dt <- function()
  {
     suppressPackageStartupMessages(library(data.table))

  # ---- Load data ----
  integration <- fread("~/workspace/annotated_GSM3516673_normal_annotated_GSM3516672_tumor_SeuratIntegration.csv")
  nt_cluster <- fread("~/workspace/nt_combined_clustering.output.csv")

  # ---- Standardize column names ----
  setnames(integration, tolower(trimws(names(integration))))
  setnames(nt_cluster, tolower(trimws(names(nt_cluster))))

  # ---- Clean Cell IDs ----
  integration[, cell := tolower(gsub("_X_", "", cell))]
  nt_cluster[, cell := tolower(cell)]  # nt_cluster already in correct form

  # ---- Merge ----
  merged <- merge(integration, nt_cluster, by = "cell", all.x = TRUE)

  # ---- Summarize counts ----
  summary_counts <- merged[, .N, by = .(integration_cluster, sample_type, cell_type)]
  setnames(summary_counts, "N", "n_cells")

  # ---- Add proportions per cluster ----
  summary_counts[, prop := n_cells / sum(n_cells), by = integration_cluster]

  # ---- Sort for readability ----
  setorder(summary_counts, integration_cluster, -n_cells)
  return(summary_counts)

  }
