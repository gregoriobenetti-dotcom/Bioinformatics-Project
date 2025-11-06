############################################################
# task11_df.R — SNP-to-gene mapping
# Data: variants.csv, gene_annotation.bed.csv
# Base R (data.frame) version using fread() for input
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
variants <- as.data.frame(fread("~/workspace/variants.csv"))
genes    <- as.data.frame(fread("~/workspace/gene_annotation.bed.csv"))

# ---- Standardize column names ----
names(variants) <- tolower(trimws(names(variants)))
names(genes)    <- tolower(trimws(names(genes)))

# Expect:
# variants: chr, pos, sample_id, impact
# genes: chr, start, end, gene

# ---- 1️⃣ Convert variant positions to 1-bp intervals ----
variants$start <- variants$pos
variants$end   <- variants$pos

# ---- 2️⃣ Find overlaps manually ----
# Condition: variant.start <= gene.end  AND  variant.end >= gene.start
overlaps <- do.call(rbind, lapply(split(variants, variants$chr), function(v_chr) {
  g_chr <- subset(genes, chr == unique(v_chr$chr))
  if (nrow(g_chr) == 0) return(NULL)
  res <- list()
  for (i in seq_len(nrow(v_chr))) {
    v <- v_chr[i, ]
    overlap_idx <- which(g_chr$start <= v$end & g_chr$end >= v$start)
    if (length(overlap_idx) > 0) {
      tmp <- g_chr[overlap_idx, ]
      tmp$sample_id <- v$sample_id
      tmp$impact <- v$impact
      res[[length(res) + 1]] <- tmp
    }
  }
  if (length(res) > 0) do.call(rbind, res) else NULL
}))

# ---- 3️⃣ Filter for HIGH-impact variants ----
high_imp <- subset(overlaps, tolower(impact) == "high")

# ---- 4️⃣ Count HIGH variants by gene and sample ----
counts_gene_sample <- aggregate(impact ~ gene + sample_id, data = high_imp, FUN = length)
names(counts_gene_sample)[3] <- "high_variant_count"

# ---- 5️⃣ Find genes with HIGH-impact variants in all samples ----
samples_with_high <- unique(high_imp[, c("gene", "sample_id")])
gene_sample_counts <- aggregate(sample_id ~ gene, data = samples_with_high, FUN = length)

total_samples <- length(unique(variants$sample_id))
genes_all_samples <- subset(gene_sample_counts, sample_id == total_samples)$gene

# ---- Output preview ----
cat("\n✅ Count of HIGH-impact variants by gene and sample:\n")
print(head(counts_gene_sample))

cat("\n✅ Genes with HIGH-impact variants in all samples:\n")
print(genes_all_samples)
