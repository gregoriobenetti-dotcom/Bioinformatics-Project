# task3_minimal.R — data.table only version
suppressPackageStartupMessages(library(data.table))

# Load data (assuming you're inside /home/rstudio/workspace)
meta   <- fread("~/workspace/sample_metadata.csv")
counts <- fread("~/workspace/bulk_counts_long.csv")

# 1. Key metadata by sample_id
setkey(meta, sample_id)

# 2. Join metadata into counts
joined <- meta[counts, on = "sample_id"]

# 3. Add secondary index on (gene, sample_id)
setindex(counts, gene, sample_id)

# 4. Benchmark subset query before and after index
query <- function() counts[gene == "GENE_005" & sample_id == "SAMPLE_01"]

system.time(query())   # before index (first run)
system.time(query())   # after index (second run, uses index)

# Preview joined data
head(joined)

