# task5.R — Classify labs as normal/out_of_range
suppressPackageStartupMessages(library(data.table))

labs <- fread("~/workspace/clinical_labs.csv")
ranges <- fread("~/workspace/lab_reference_ranges.csv")

#--- data.table version ---
dt_labs <- copy(labs)
dt_ranges <- copy(ranges)

# Non-equi join: value between lower and upper
setkey(dt_ranges, lab)
classified <- dt_ranges[dt_labs, on = .(lab, lower <= value, upper >= value), nomatch = 0L]
classified[, status := ifelse(!is.na(lower), "normal", "out_of_range")]

# Fill missing labs as out_of_range
out_all <- merge(dt_labs, classified, all.x = TRUE)
out_all[is.na(status), status := "out_of_range"]

# Count abnormal rates
abnormal_rates_patient <- out_all[, .(abnormal_rate = mean(status == "out_of_range")), by = patient_id]
abnormal_rates_lab <- out_all[, .(abnormal_rate = mean(status == "out_of_range")), by = lab]

#--- data.frame version ---
df_labs <- as.data.frame(labs)
df_ranges <- as.data.frame(ranges)

df_joined <- merge(df_labs, df_ranges, by = "lab", all.x = TRUE)
df_joined$status <- ifelse(df_joined$value >= df_joined$lower & df_joined$value <= df_joined$upper, "normal", "out_of_range")

df_pat <- aggregate(status == "out_of_range" ~ patient_id, df_joined, mean)
df_lab <- aggregate(status == "out_of_range" ~ lab, df_joined, mean)

#--- Output check ---
print(head(abnormal_rates_patient))
print(head(abnormal_rates_lab))
