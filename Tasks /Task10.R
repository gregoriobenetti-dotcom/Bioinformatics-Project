############################################################
# task10.R — ATAC-to-gene mapping
# Data: atac_peaks.bed.csv, gene_annotation.bed.csv
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
peaks <- fread("~/workspace/atac_peaks.bed.csv")
genes <- fread("~/workspace/gene_annotation.bed.csv")

# ---- Standardize column names ----
setnames(peaks, tolower(trimws(names(peaks))))
setnames(genes, tolower(trimws(names(genes))))

# Expect columns:
# peaks: chr, start, end, score (optionally others)
# genes: chr, start, end, gene

# ---- Key by chromosome and coordinates ----
setkey(peaks, chr, start, end)
setkey(genes, chr, start, end)

# ---- Non-equi join to find overlaps ----
# Condition: peak.start <= gene.end  AND  peak.end >= gene.start
overlaps <- foverlaps(
  peaks,
  genes,
  by.x = c("chr", "start", "end"),
  by.y = c("chr", "start", "end"),
  type = "any",
  nomatch = 0L
)

# ---- Compute overlap length in bp ----
# overlap = min(peak.end, gene.end) - max(peak.start, gene.start)
overlaps[, overlap_bp :=
           pmax(0, pmin(end, i.end) - pmax(start, i.start))
]

# ---- Count peaks per gene ----
peaks_per_gene <- overlaps[, .(n_peaks = uniqueN(.I)), by = gene]

# ---- Sum overlap length per gene ----
bp_per_gene <- overlaps[, .(total_overlap_bp = sum(overlap_bp)), by = gene]

# ---- Combine and rank ----
summary <- merge(peaks_per_gene, bp_per_gene, by = "gene")
setorder(summary, -total_overlap_bp)

# ---- Select top 20 genes ----
top20 <- summary[1:20]

# ---- Output preview ----
print(top20)
cat("\n✅ Top 20 genes by total overlapped base-pairs shown above.\n")
