############################################################
# task6.R — Nearest-time matching of vitals to lab draws
# Final minimal version (handles time_iso correctly)
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
labs   <- fread("~/workspace/clinical_labs.csv")
vitals <- fread("~/workspace/vitals_time_series.csv")

# ---- Standardize column names ----
setnames(labs,   tolower(trimws(names(labs))))
setnames(vitals, tolower(trimws(names(vitals))))

# Rename 'time_iso' to 'time'
setnames(labs,   "time_iso", "time")
setnames(vitals, "time_iso", "time")

# ---- Convert to POSIXct ----
labs[,   time := as.POSIXct(time,   tz = "UTC", tryFormats = c(
  "%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S", "%Y/%m/%d %H:%M:%S"
))]
vitals[, time := as.POSIXct(time,   tz = "UTC", tryFormats = c(
  "%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S", "%Y/%m/%d %H:%M:%S"
))]

# ---- Keep only HR and SBP ----
vitals_sub <- vitals[vital %in% c("HR", "SBP")]

# Split into HR and SBP tables
hr  <- vitals_sub[vital == "HR",  .(patient_id, time, HR  = value)]
sbp <- vitals_sub[vital == "SBP", .(patient_id, time, SBP = value)]

# ---- Rolling joins ----
setkey(hr,  patient_id, time)
setkey(sbp, patient_id, time)
setkey(labs, patient_id, time)

# Find nearest HR and SBP for each lab draw
labs_hr  <- hr[labs,  roll = "nearest"]
labs_sbp <- sbp[labs, roll = "nearest"]

# Merge HR and SBP results
merged <- merge(
  labs_hr[, .(patient_id, time, lab, value, HR,  time_HR  = time)],
  labs_sbp[, .(patient_id, time,        SBP, time_SBP = time)],
  by = c("patient_id", "time"),
  all.x = TRUE
)

# ---- Compute lag (minutes) ----
merged[, lag_HR  := abs(as.numeric(difftime(time, time_HR,  units = "mins")))]
merged[, lag_SBP := abs(as.numeric(difftime(time, time_SBP, units = "mins")))]

# ---- Correlation between CRP and vitals per patient ----
corr_summary <- merged[lab == "CRP",
                       .(
                         cor_HR  = cor(value, HR,  use = "pairwise.complete.obs"),
                         cor_SBP = cor(value, SBP, use = "pairwise.complete.obs")
                       ),
                       by = patient_id]

