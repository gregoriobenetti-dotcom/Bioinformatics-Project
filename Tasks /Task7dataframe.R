############################################################
# task7_df.R — Slice genomics windows efficiently
# Data: atac_peaks.bed.csv
# Base R (data.frame) version using fread() for input
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
peaks <- as.data.frame(fread("~/workspace/atac_peaks.bed.csv"))

# ---- Standardize column names ----
names(peaks) <- tolower(trimws(names(peaks)))

# Expect: chr, start, end, score

# ---- 1️⃣ Filter peaks on chr2 with start between 2–4 Mb ----
filtered <- subset(peaks, chr == "chr2" & start >= 2e6 & start <= 4e6)

# ---- 2️⃣ Order by descending score and select top 50 ----
filtered <- filtered[order(-filtered$score), ]
top50 <- head(filtered, 50)

# ---- 3️⃣ Output preview ----
cat("\n✅ Top 50 peaks on chr2 between 2–4 Mb (by score):\n")
print(head(top50))
cat(paste("\nSelected", nrow(top50), "peaks.\n"))
