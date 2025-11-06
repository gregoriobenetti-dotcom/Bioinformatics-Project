############################################################
# task10_df.R — ATAC-to-gene mapping
# Data: atac_peaks.bed.csv, gene_annotation.bed.csv
# Base R (data.frame) version using fread() for input
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
peaks <- as.data.frame(fread("~/workspace/atac_peaks.bed.csv"))
genes <- as.data.frame(fread("~/workspace/gene_annotation.bed.csv"))

# ---- Standardize column names ----
names(peaks) <- tolower(trimws(names(peaks)))
names(genes) <- tolower(trimws(names(genes)))

# Expect:
# peaks: chr, start, end, score
# genes: chr, start, end, gene

# ---- 1️⃣ Find overlaps manually (base R loop version) ----
# Condition: peak.start <= gene.end  AND  peak.end >= gene.start

overlaps <- do.call(rbind, lapply(split(peaks, peaks$chr), function(pk_chr) {
  gn_chr <- subset(genes, chr == unique(pk_chr$chr))
  if (nrow(gn_chr) == 0) return(NULL)
  res <- list()
  for (i in seq_len(nrow(pk_chr))) {
    p <- pk_chr[i, ]
    overlap_idx <- which(gn_chr$start <= p$end & gn_chr$end >= p$start)
    if (length(overlap_idx) > 0) {
      tmp <- gn_chr[overlap_idx, ]
      tmp$peak_start <- p$start
      tmp$peak_end <- p$end
      tmp$peak_score <- p$score
      res[[length(res) + 1]] <- tmp
    }
  }
  if (length(res) > 0) do.call(rbind, res) else NULL
}))

# ---- 2️⃣ Compute overlap length ----
if (!is.null(overlaps)) {
  overlaps$overlap_bp <- pmax(0, pmin(overlaps$end, overlaps$peak_end) -
                                pmax(overlaps$start, overlaps$peak_start))
}

# ---- 3️⃣ Count peaks and sum overlap per gene ----
peaks_per_gene <- aggregate(peak_start ~ gene, data = overlaps, FUN = function(x) length(unique(x)))
names(peaks_per_gene)[2] <- "n_peaks"

bp_per_gene <- aggregate(overlap_bp ~ gene, data = overlaps, sum)
names(bp_per_gene)[2] <- "total_overlap_bp"

# ---- 4️⃣ Combine and rank ----
summary_df <- merge(peaks_per_gene, bp_per_gene, by = "gene", all = TRUE)
summary_df <- summary_df[order(-summary_df$total_overlap_bp), ]

# ---- 5️⃣ Select top 20 genes ----
top20 <- head(summary_df, 20)

# ---- Output preview ----
cat("\n✅ Top 20 genes by total overlapped base pairs:\n")
print(top20)
