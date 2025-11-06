############################################################
# task7.R — Slice genomics windows efficiently
# Data: atac_peaks.bed.csv
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
peaks <- fread("~/workspace/atac_peaks.bed.csv")

# ---- Standardize column names ----
setnames(peaks, tolower(trimws(names(peaks))))

# Expect columns like: chr, start, end, score
# You can verify with: names(peaks)

# ---- Filter: peaks on chr2 with start between 2–4 Mb ----
filtered <- peaks[chr == "chr2" & start >= 2e6 & start <= 4e6]

# ---- Order by descending score and select top 50 ----
setorder(filtered, -score)
top50 <- filtered[1:50]

# ---- Output preview ----
print(head(top50))
print(paste("✅ Selected", nrow(top50), "peaks on chr2 between 2–4 Mb (top 50 by score)."))

