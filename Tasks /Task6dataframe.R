############################################################
# task6_df.R — Nearest-time matching of vitals to lab draws
# Data: clinical_labs.csv, vitals_time_series.csv
# Base R (data.frame) version using fread() for input
############################################################

suppressPackageStartupMessages(library(data.table))

# ---- Load data ----
labs   <- as.data.frame(fread("~/workspace/clinical_labs.csv"))
vitals <- as.data.frame(fread("~/workspace/vitals_time_series.csv"))

# ---- Standardize column names ----
names(labs)   <- tolower(trimws(names(labs)))
names(vitals) <- tolower(trimws(names(vitals)))

# Rename time column if needed
if ("time_iso" %in% names(labs))   names(labs)[names(labs) == "time_iso"] <- "time"
if ("time_iso" %in% names(vitals)) names(vitals)[names(vitals) == "time_iso"] <- "time"

# ---- Convert time columns ----
labs$time   <- as.POSIXct(labs$time,   tz = "UTC", tryFormats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S"))
vitals$time <- as.POSIXct(vitals$time, tz = "UTC", tryFormats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S"))

# ---- Keep only HR and SBP ----
vitals_sub <- subset(vitals, vital %in% c("HR", "SBP"))

# Split into HR and SBP data.frames
hr  <- subset(vitals_sub, vital == "HR",  select = c("patient_id", "time", "value"))
names(hr)[3] <- "HR"

sbp <- subset(vitals_sub, vital == "SBP", select = c("patient_id", "time", "value"))
names(sbp)[3] <- "SBP"

# ---- Find nearest HR and SBP for each lab draw ----
# (Manual nearest match per patient)
nearest_match <- function(lab_df, vital_df, vital_name) {
  res <- lab_df
  res[[vital_name]] <- NA
  res[[paste0("lag_", vital_name)]] <- NA
  for (i in seq_len(nrow(res))) {
    pid <- res$patient_id[i]
    t   <- res$time[i]
    sub_vital <- subset(vital_df, patient_id == pid)
    if (nrow(sub_vital) > 0) {
      idx <- which.min(abs(difftime(sub_vital$time, t, units = "mins")))
      res[[vital_name]][i] <- sub_vital[[vital_name]][idx]
      res[[paste0("lag_", vital_name)]][i] <- abs(as.numeric(difftime(sub_vital$time[idx], t, units = "mins")))
    }
  }
  res
}

# Match HR and SBP
labs_hr  <- nearest_match(labs, hr,  "HR")
labs_sbp <- nearest_match(labs, sbp, "SBP")

# Combine HR and SBP
merged <- merge(labs_hr, labs_sbp[, c("patient_id", "time", "SBP", "lag_SBP")],
                by = c("patient_id", "time"), all.x = TRUE)

# ---- Correlation between CRP and vitals per patient ----
crp_data <- subset(merged, lab == "CRP")
corr_summary <- aggregate(
  cbind(HR, SBP) ~ patient_id, data = crp_data,
  FUN = function(x) cor(crp_data$value[!is.na(x)], x[!is.na(x)], use = "pairwise.complete.obs")
)

# ---- Output preview ----
cat("\n✅ Merged labs with nearest HR and SBP values:\n")
print(head(merged))

cat("\n✅ Correlation summary (CRP vs HR/SBP):\n")
print(head(corr_summary))
