############################################################
# task5_df.R — Classify labs as normal or out_of_range
# Data: clinical_labs.csv, lab_reference_ranges.csv
# Base R (data.frame) version using fread() for input
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
labs   <- as.data.frame(fread("~/workspace/clinical_labs.csv"))
ranges <- as.data.frame(fread("~/workspace/lab_reference_ranges.csv"))

# ---- Standardize column names ----
names(labs)   <- tolower(trimws(names(labs)))
names(ranges) <- tolower(trimws(names(ranges)))

# Expect:
# labs: patient_id, lab, value
# ranges: lab, lower, upper

# ---- 1️⃣ Merge by lab ----
joined <- merge(labs, ranges, by = "lab", all.x = TRUE)

# ---- 2️⃣ Classify values ----
joined$status <- ifelse(
  joined$value >= joined$lower & joined$value <= joined$upper,
  "normal",
  "out_of_range"
)

# ---- 3️⃣ Compute abnormal rates ----
# Per patient
abnormal_rates_patient <- aggregate(status == "out_of_range" ~ patient_id, data = joined, mean)
names(abnormal_rates_patient)[2] <- "abnormal_rate"

# Per lab
abnormal_rates_lab <- aggregate(status == "out_of_range" ~ lab, data = joined, mean)
names(abnormal_rates_lab)[2] <- "abnormal_rate"

# ---- 4️⃣ Output preview ----
cat("\n✅ Abnormal rates per patient:\n")
print(head(abnormal_rates_patient))

cat("\n✅ Abnormal rates per lab:\n")
print(head(abnormal_rates_lab))
