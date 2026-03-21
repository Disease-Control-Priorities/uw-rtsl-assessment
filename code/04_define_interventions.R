# Required inputs----

#...........................................................
# HTN Impact----
#...........................................................

# Import program data 2017-2025 from HTN global tracker

htn_tracker_file <- paste0(wd_raw, "HTN global tracker 2017-2025_15Mar2026.xlsx")

# -- 1. Import population enrolled in hypertension control programs by location and year
dt_enrolled_wide <- as.data.table(
  read_excel(htn_tracker_file,
             sheet = "Cumulative registrations",
             skip  = 1)
)

# Clean column names: year columns read as numeric -> force to integer character
yr_cols_enr <- setdiff(names(dt_enrolled_wide), c("Country","start_year"))
setnames(dt_enrolled_wide, yr_cols_enr, as.character(as.integer(as.numeric(yr_cols_enr))))
yr_cols_enr <- as.character(2017:2025)

# Drop non-country rows (blanks, notes, URLs at the bottom)
dt_enrolled_wide <- dt_enrolled_wide[!is.na(Country) & Country != "" &
                                       !grepl("^(PAHO|http|\\s*$)", Country, ignore.case = TRUE)]

# -- 2. Reshape enrolled wide to long format
dt_enrolled <- melt(
  dt_enrolled_wide,
  id.vars       = c("Country","start_year"),
  measure.vars  = yr_cols_enr,
  variable.name = "year",
  value.name    = "enrolled"
)
dt_enrolled[, year := as.integer(as.character(year))]
setnames(dt_enrolled, "Country", "location")

# -- 3. Import population with BP control by location and year
dt_controlled_wide <- as.data.table(
  read_excel(htn_tracker_file,
             sheet = "BP Controlled",
             skip  = 1)
)

yr_cols_ctrl <- setdiff(names(dt_controlled_wide), "Country")
setnames(dt_controlled_wide, yr_cols_ctrl, as.character(as.integer(as.numeric(yr_cols_ctrl))))
yr_cols_ctrl <- as.character(2017:2025)

# Drop non-country rows
dt_controlled_wide <- dt_controlled_wide[!is.na(Country) & Country != "" &
                                           !grepl("^(PAHO|http|\\s*$)", Country, ignore.case = TRUE)]

# -- 4. Reshape BP controlled wide to long format
dt_controlled <- melt(
  dt_controlled_wide,
  id.vars       = "Country",
  measure.vars  = yr_cols_ctrl,
  variable.name = "year",
  value.name    = "bp_controlled"
)
dt_controlled[, year := as.integer(as.character(year))]
setnames(dt_controlled, "Country", "location")

# -- 5. Merge datasets to get control rates by location and year
dt_htn_program <- merge(
  dt_enrolled,
  dt_controlled,
  by = c("location", "year"),
  all = TRUE
)


# Control rate = BP controlled / enrolled (NA where either is missing or enrolled == 0)

dt_htn_program[, `:=`(
  bp_controlled = as.numeric(bp_controlled),
  enrolled      = as.numeric(enrolled)
)]

dt_htn_program[, control_rate := fifelse(
  !is.na(bp_controlled) & !is.na(enrolled) & enrolled > 0,
  bp_controlled / enrolled,
  NA_real_
)]

setorder(dt_htn_program, location, year)

## Bound control rates to [0,1] and create a dummy to indicate where control rates are observed vs. imputed

dt_htn_program[, control_rate_observed := 0]
dt_htn_program[!is.na(control_rate), control_rate_observed := 1]
dt_htn_program[control_rate>1,control_rate_observed := 0]

dt_htn_program[control_rate>1,enrolled := bp_controlled]
dt_htn_program[, control_rate := pmin(1, pmax(0, control_rate))]

# -- 6. Computing growth rates in enrolled population and control rates by location
# Year-over-year growth rate for enrolled: (enrolled_t / enrolled_{t-1}) - 1
# Year-over-year change in control rate:   control_rate_t - control_rate_{t-1}
dt_htn_program[, `:=`(
  enrolled_prev     = shift(enrolled,     n = 1, type = "lag"),
  control_rate_prev = shift(control_rate, n = 1, type = "lag")
), by = location]

dt_htn_program[, `:=`(
  enrolled_growth_rate  = fifelse(
    !is.na(enrolled) & !is.na(enrolled_prev) & enrolled_prev > 0,
    enrolled / enrolled_prev - 1,
    NA_real_
  ),
  control_rate_change = fifelse(
    !is.na(control_rate) & !is.na(control_rate_prev),
    control_rate - control_rate_prev,
    NA_real_
  )
)]

# Median growth rates by location (across available years)
dt_htn_growth <- dt_htn_program[
  !is.na(enrolled_growth_rate),
  .(median_enrolled_growth = median(enrolled_growth_rate, na.rm = TRUE)),
  by = location
]

dt_ctrl_growth <- dt_htn_program[
  !is.na(control_rate_change),
  .(median_ctrl_change = median(control_rate_change, na.rm = TRUE)),
  by = location
]

dt_htn_program <- merge(dt_htn_program, dt_htn_growth, by = "location", all.x = TRUE)
dt_htn_program <- merge(dt_htn_program, dt_ctrl_growth, by = "location", all.x = TRUE)

# -- 7. Impute missing data based on growth rates
# Forward-fill enrolled using median growth rate where gaps exist
# Forward-fill control rate using median change in control rate
setorder(dt_htn_program, location, year)

dt_htn_program[, `:=`(
  enrolled_imputed     = as.numeric(enrolled),
  control_rate_imputed = as.numeric(control_rate)
)]

# Forward imputation: carry forward using location-specific growth rates
locations_htn <- unique(dt_htn_program$location)

for (loc in locations_htn) {
  idx <- which(dt_htn_program$location == loc)
  for (i in seq_along(idx)[-1]) {
    row_cur  <- idx[i]
    row_prev <- idx[i - 1]
    
    # Impute enrolled if missing but previous value and growth rate available
    if (is.na(dt_htn_program$enrolled_imputed[row_cur]) &
        !is.na(dt_htn_program$enrolled_imputed[row_prev]) &
        !is.na(dt_htn_program$median_enrolled_growth[row_cur])) {
      dt_htn_program$enrolled_imputed[row_cur] <-
        dt_htn_program$enrolled_imputed[row_prev] *
        (1 + dt_htn_program$median_enrolled_growth[row_cur])
    }
    
    # Impute control rate if missing but previous value and change rate available
    if (is.na(dt_htn_program$control_rate_imputed[row_cur]) &
        !is.na(dt_htn_program$control_rate_imputed[row_prev]) &
        !is.na(dt_htn_program$median_ctrl_change[row_cur])) {
      dt_htn_program$control_rate_imputed[row_cur] <- pmin(1, pmax(0,
                                                                   dt_htn_program$control_rate_imputed[row_prev] +
                                                                     dt_htn_program$median_ctrl_change[row_cur]
      ))
    }
  }
}

# Recompute imputed BP controlled from imputed enrolled and control rate
dt_htn_program[, bp_controlled_imputed := enrolled_imputed * control_rate_imputed]

# Drop working columns
dt_htn_program[, c("enrolled_prev", "control_rate_prev") := NULL]

# Diagnostics
cat("=== HTN Program Impact: Diagnostics ===\n")
cat("  Unique countries:              ", uniqueN(dt_htn_program$location), "\n")
cat("  Years covered:                 ", paste(range(dt_htn_program$year), collapse = "-"), "\n")
cat("  Rows with observed enrolled:   ", sum(!is.na(dt_htn_program$enrolled)), "\n")
cat("  Rows with imputed enrolled:    ", sum(!is.na(dt_htn_program$enrolled_imputed) &
                                               is.na(dt_htn_program$enrolled)), "\n")
cat("  Rows with observed control:    ", sum(!is.na(dt_htn_program$bp_controlled)), "\n")
cat("  Rows with imputed control:     ", sum(!is.na(dt_htn_program$bp_controlled_imputed) &
                                               is.na(dt_htn_program$bp_controlled)), "\n")

# Save output
fwrite(dt_htn_program, paste0(wd_data, "htn_program_impact_by_loc_year.csv"))
cat("Written: htn_program_impact_by_loc_year.csv\n")

# Clean up
rm(dt_enrolled_wide, dt_enrolled, dt_controlled_wide, dt_controlled,
   dt_htn_growth, dt_ctrl_growth, yr_cols_enr, yr_cols_ctrl,
   locations_htn, htn_tracker_file)


## Temp blood pressure coverage/targets

## Import targets excel
aim2_loc <- as.data.table(
  read_excel("C:/Users/wrgar/OneDrive - UW/02Work/WHO-CVD/Scenarios.xlsx",
             sheet = "Sheet1",range = "A6:O197")
)

# keep target control rates baseline htncov2, _aspirational,_ambitious, and _progress
aim2_loc <- aim2_loc[, .(location,htncov2, htncov2_aspirational, htncov2_ambitious, htncov2_progress,htn_ctrl_diabetes)]

# save .csv for input in running the model
fwrite(aim2_loc, paste0(wd_data,"htn_control_targets_by_loc.csv"))

# TFA Impact----
