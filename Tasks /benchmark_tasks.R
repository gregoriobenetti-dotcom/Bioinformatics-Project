############################################################
# benchmark_tasks.R
# Compare data.table vs data.frame functions (Tasks 1–12)
############################################################

############################################################
# Works automatically whether functions are in package or global env
############################################################

suppressPackageStartupMessages({
  library(microbenchmark)
  library(dplyr)
  library(tidyr)
})

# --- Auto-load all functions if not already available ---
if (!exists("task1_filter_counts")) {
  if (file.exists("~/functionshelp/R/Functionsdatatable.R")) source("~/functionshelp/R/Functionsdatatable.R")
}
if (!exists("task1_df_filter_counts")) {
  if (file.exists("~/functionshelp/R/Functionsdataframe.R")) source("~/functionshelp/R/Functionsdataframe.R")
}

# --- Helper: list available function names ---
dt_funcs <- ls(pattern = "^task[0-9]+_", envir = .GlobalEnv)
df_funcs <- grep("_df_", dt_funcs, value = TRUE)

# Match df and dt functions by task number
task_pairs <- tibble(
  task_num = gsub("^task([0-9]+).*", "\\1", df_funcs),
  df_fun = df_funcs
) %>%
  mutate(dt_fun = paste0("task", task_num, "_",
                         gsub("task[0-9]+_df_", "", df_fun)))

# --- Run benchmark for each pair ---
benchmark_results <- list()

for (i in seq_len(nrow(task_pairs))) {
  task_name <- paste0("task", task_pairs$task_num[i])
  dt_fun_name <- task_pairs$dt_fun[i]
  df_fun_name <- task_pairs$df_fun[i]

  cat("\n⏱️ Running", task_name, "...\n")

  if (!exists(dt_fun_name, mode = "function") ||
      !exists(df_fun_name, mode = "function")) {
    warning("⚠️ Skipping ", task_name, " — one of the functions not found.")
    next
  }

  dt_fun <- get(dt_fun_name, mode = "function")
  df_fun <- get(df_fun_name, mode = "function")

  mb <- microbenchmark(
    data.table = { suppressMessages(dt_fun()) },
    data.frame = { suppressMessages(df_fun()) },
    times = 3L
  )

  benchmark_results[[task_name]] <- data.frame(
    task   = task_name,
    method = c("data.table", "data.frame"),
    mean_ms = tapply(mb$time / 1e6, mb$expr, mean)
  )
}

# --- Combine all results ---
if (length(benchmark_results) == 0) {
  cat("\n⚠️ No benchmarks could be run — check that your functions are defined.\n")
} else

# ---- Combine results ----
benchmark_summary <- do.call(rbind, benchmark_results)
row.names(benchmark_summary) <- NULL

# ---- Compute speed ratio ----
library(dplyr)
library(tidyr)

# Pivot wide to have both methods in separate columns
benchmark_summary <- benchmark_summary |>
tidyr::pivot_wider(names_from = method, values_from = mean_ms) |>
dplyr::mutate(speed_ratio = data.frame / data.table)

# ---- Print summary ----
cat("\n✅ Benchmark Summary (mean time in ms):\n")
print(benchmark_summary[, c("task", "data.table", "data.frame", "speed_ratio")])

# ---- Optional: Plot automatically ----
library(ggplot2)
ggplot(benchmark_summary, aes(x = reorder(task, speed_ratio), y = speed_ratio)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    title = "Speed Comparison: data.frame vs data.table",
    x = "Task",
    y = "Speed ratio (data.frame / data.table)"
  ) +
  theme_minimal(base_size = 13)

