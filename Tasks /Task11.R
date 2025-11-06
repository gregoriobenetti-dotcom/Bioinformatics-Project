############################################################
# task11.R — SNP-to-gene mapping
# Data: variants.csv, gene_annotation.bed.csv
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
variants <- fread("~/workspace/variants.csv")
genes    <- fread("~/workspace/gene_annotation.bed.csv")

# ---- Standardize column names ----
setnames(variants, tolower(trimws(names(variants))))
setnames(genes,    tolower(trimws(names(genes))))

# Expect:
# variants: chr, pos, sample_id, impact
# genes:    chr, start, end, gene

# ---- Convert variant positions to 1-bp intervals ----
variants[, start := pos]
variants[, end   := pos]

# ---- Set keys for overlap join ----
setkey(variants, chr, start, end)
setkey(genes,    chr, start, end)

# ---- Find overlaps (which gene each variant falls in) ----
overlaps <- foverlaps(
  variants,
  genes,
  by.x = c("chr", "start", "end"),
  by.y = c("chr", "start", "end"),
  type = "any",
  nomatch = 0L
)

# ---- Filter for HIGH-impact variants ----
high_imp <- overlaps[tolower(impact) == "high"]

# ---- Count HIGH variants by gene and sample ----
counts_gene_sample <- high_imp[, .N, by = .(gene, sample_id)]
setnames(counts_gene_sample, "N", "high_variant_count")

# ---- Summarize genes with HIGH-impact variants in all samples ----
samples_with_high <- unique(high_imp[, .(gene, sample_id)])
gene_sample_counts <- samples_with_high[, .N, by = gene]
total_samples <- uniqueN(variants$sample_id)

genes_all_samples <- gene_sample_counts[N == total_samples, gene]

# ---- Output ----
cat("✅ Count of HIGH-impact variants by gene and sample:\n")
print(head(counts_gene_sample))

cat("\n✅ Genes with HIGH-impact variants in all samples:\n")
print(genes_all_samples)
