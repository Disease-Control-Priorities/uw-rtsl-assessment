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

# Imputation to not misss countries

dt_enrolled[, enrolled_org := as.numeric(enrolled)]
dt_enrolled[is.na(start_year) & is.na(enrolled), enrolled := 1]
dt_enrolled[is.na(start_year) & is.na(enrolled_org), start_year := 2026]

# remove temp variable
dt_enrolled[, enrolled_org := NULL]

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
  enrolled_growth_rate = fifelse(
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

# Median growth in enrolled by location
dt_htn_growth <- dt_htn_program[
  !is.na(enrolled_growth_rate),
  .(median_enrolled_growth = median(enrolled_growth_rate, na.rm = TRUE)),
  by = location
]

# Control-rate history summary by location
dt_ctrl_info <- dt_htn_program[
  !is.na(control_rate),
  .(
    n_obs_control      = .N,
    first_obs_year     = min(year),
    first_obs_control  = control_rate[which.min(year)],
    median_ctrl_change = if (.N >= 2) median(diff(control_rate[order(year)]), na.rm = TRUE) else NA_real_
  ),
  by = location
]

dt_htn_program <- merge(dt_htn_program, dt_htn_growth,   by = "location", all.x = TRUE)
dt_htn_program <- merge(dt_htn_program, dt_ctrl_info,    by = "location", all.x = TRUE)

# -- 7. Impute missing data based on growth rates / single-observation fallback
setorder(dt_htn_program, location, year)

dt_htn_program[, `:=`(
  enrolled_imputed     = as.numeric(enrolled),
  control_rate_imputed = as.numeric(control_rate)
)]

locations_htn <- unique(dt_htn_program$location)

for (loc in locations_htn) {
  idx <- which(dt_htn_program$location == loc)
  
  for (i in seq_along(idx)) {
    row_cur <- idx[i]
    
    # ---------------------------
    # ENROLLED: forward imputation
    # ---------------------------
    if (i > 1) {
      row_prev <- idx[i - 1]
      
      if (is.na(dt_htn_program$enrolled_imputed[row_cur]) &&
          !is.na(dt_htn_program$enrolled_imputed[row_prev]) &&
          !is.na(dt_htn_program$median_enrolled_growth[row_cur])) {
        dt_htn_program$enrolled_imputed[row_cur] <-
          dt_htn_program$enrolled_imputed[row_prev] *
          (1 + dt_htn_program$median_enrolled_growth[row_cur])
      }
    }
    
    # -----------------------------------------
    # CONTROL RATE: two cases
    # 1) >=2 observations -> forward using median annual change
    # 2) exactly 1 observation -> carry that value forward to 2025
    # -----------------------------------------
    if (is.na(dt_htn_program$control_rate_imputed[row_cur])) {
      
      # Case 1: enough history to estimate own trend
      if (!is.na(dt_htn_program$n_obs_control[row_cur]) &&
          dt_htn_program$n_obs_control[row_cur] >= 2 &&
          i > 1) {
        
        row_prev <- idx[i - 1]
        
        if (!is.na(dt_htn_program$control_rate_imputed[row_prev]) &&
            !is.na(dt_htn_program$median_ctrl_change[row_cur])) {
          dt_htn_program$control_rate_imputed[row_cur] <- pmin(
            1, pmax(0,
                    dt_htn_program$control_rate_imputed[row_prev] +
                      dt_htn_program$median_ctrl_change[row_cur])
          )
        }
      }
      
      # Case 2: only one observed value -> carry forward constant after first observed year
      if (!is.na(dt_htn_program$n_obs_control[row_cur]) &&
          dt_htn_program$n_obs_control[row_cur] == 1 &&
          !is.na(dt_htn_program$first_obs_year[row_cur]) &&
          !is.na(dt_htn_program$first_obs_control[row_cur]) &&
          dt_htn_program$year[row_cur] > dt_htn_program$first_obs_year[row_cur] &&
          dt_htn_program$year[row_cur] <= 2025) {
        
        dt_htn_program$control_rate_imputed[row_cur] <- pmin(
          1, pmax(0, dt_htn_program$first_obs_control[row_cur])
        )
      }
    }
  }
}

# Recompute imputed BP controlled from imputed enrolled and control rate
dt_htn_program[, bp_controlled_imputed := enrolled_imputed * control_rate_imputed]

# Drop working columns
dt_htn_program[, c(
  "enrolled_prev", "control_rate_prev",
  "n_obs_control", "first_obs_year", "first_obs_control"
) := NULL]

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

# select only control rates for input into model and save as .csv
dt_htn_control_rates <- dt_htn_program[, .(location, year, start_year, control_rate,control_rate_imputed, control_rate_observed)]

## Temp blood pressure coverage/targets

## Import bussiness as usual hypertension control targets by location from Scenarios.xlsx (same as used in 03_clean_inputs.R)
dt_htn_bau <- as.data.table(
  read_excel("C:/Users/wrgar/OneDrive - UW/02Work/WHO-CVD/Scenarios.xlsx",
             sheet = "Sheet1",range = "A6:O197")
)

# keep target control rates baseline htncov2
dt_htn_bau <- dt_htn_bau[, .(location,htncov2)]

#identify locations with missing control rates in the HTN program dataset and add them to the BAU dataset with NA control rates
missing_locs <- setdiff(dt_htn_program$location, dt_htn_bau$location)

# Rename location in dt_htn_program to match dt_htn_bau (if needed)
dt_htn_program[location == "Antigua Barbuda", location := "Antigua and Barbuda"]
dt_htn_program[location == "Vietnam", location := "Viet Nam"]

# Locations definitely missing: c("Anguilla", "British Virgin Islands", "St Kitts and Nevis","Turks and Caicos") 

# Keep only countries in the HTN program dataset
dt_htn_bau <- dt_htn_bau[location %in% dt_htn_program$location,]

# Merge the two datasets to create a comprehensive list of locations with control rates (observed or imputed) and BAU targets

dt_htn_control_rates <- merge(
  dt_htn_program[, .(location, start_year,year, enrolled, control_rate, control_rate_imputed, control_rate_observed)],
  dt_htn_bau,
  by = "location",
  all.x = TRUE
)
# Filter NA htncov2 to keep only locations with no observed or imputed control rates in the HTN program datasetv2
dt_htn_control_rates <- dt_htn_control_rates[!(is.na(htncov2)),]

# Omit NA start year
dt_htn_control_rates <- dt_htn_control_rates[!is.na(start_year),]

# Count number of non NA observations in control rates by location and the filter locations with at least 1 non NA observation
dt_htn_control_rates[, non_na_control_rates := sum(!is.na(control_rate) | !is.na(control_rate_imputed)), by = location]
#dt_htn_control_rates <- dt_htn_control_rates[non_na_control_rates > 0 ,]

# Altered to not missed countries
dt_htn_control_rates <- dt_htn_control_rates[non_na_control_rates > 0 | start_year==2026,]

dt_htn_control_rates[, non_na_control_rates := NULL]
# Keep most recent observed data (2025)
dt_mean_25 <- dt_htn_control_rates[year %in% c(2025), ]
dt_mean_25 <- unique(dt_mean_25, by = "location")

# Extend constant values from 2026 to 2050
years <- 2026:2050
dt_htn_control_scenarios <- dt_mean_25[, {
  yearly_data <- data.table(year = years)
  yearly_data[, start_year := start_year]
  yearly_data[, enrolled := enrolled]
  yearly_data[, control_rate := control_rate]
  yearly_data[, control_rate_imputed := control_rate_imputed]
  yearly_data[, control_rate_observed := control_rate_observed]
  yearly_data[, htncov2 := htncov2]
  yearly_data
}, by = location]

# Combine observed and projected values
dt_htn_control_scenarios <- rbind(
  dt_htn_control_rates,
  dt_htn_control_scenarios,
  fill = TRUE
)

# Impute 1 to enrolled if NA. 0 could cause numeric issues, and since we are intered in difference
# in risk 1 vs. 0, imputing 1 to missing enrolled will allow us to keep the control rate as 0 where there is no data, which is a conservative assumption for the impact of the program (i.e., no impact where no data). This also allows us to avoid issues with division by zero when calculating control rates.

dt_htn_control_scenarios[is.na(enrolled), enrolled := 1]

dt_htn_control_scenarios[, control_rate_program := fifelse(
  !is.na(control_rate_imputed),
  control_rate_imputed,
  0
)]

# create two bau scenario vectors to compare with the program observed data:

# 1. Upper : Observed/imputed control rates from HTN program dataset (2025 values extended forward) - 0

dt_htn_control_scenarios[, control_rate_bau_upper := 0]

# 2. Conservative : Observed/imputed control rates from HTN program dataset (2025 values extended forward) - BAU control rates from Scenarios.xlsx (htncov2) - 0

dt_htn_control_scenarios[, control_rate_bau_conservative := fifelse(
  !is.na(htncov2),
  control_rate_program - htncov2,
  control_rate_program
)]

# Negative control rate reductions don't make sense, so set any negative values to 0
dt_htn_control_scenarios[control_rate_bau_conservative < 0, control_rate_bau_conservative := control_rate_program]

# Sort observations by location and year
setorder(dt_htn_control_scenarios, location, year)

fwrite(dt_htn_control_scenarios, paste0(wd_data, "htn_control_scenarios_by_loc_year.csv"))


# TFA Impact----


##? pending to re check IHME pop groups

dt_tfa_scenarios <- readRDS(file = paste0(wd_data,"tfa_policy_scenarios.rds"))

dt_tfa_bpp <- as.data.table(read_excel(paste0(wd_raw,"List of Countries with Policies Passed-updated March 2026.xlsx"), 
                                       sheet = "Sheet1", range = "A3:D56")) 

dt_tfa_bpp[, tfa_rtsl_include := "Yes"]

## Add countries for complete list including not involved in RTSL calculations

dt_tfa_bpp_nortsl <- as.data.table(read_excel(paste0(wd_raw,"List of Countries with Policies Passed-updated March 2026.xlsx"), 
                                       sheet = "Sheet1", range = "E3:G18")) 

# Create dummy before replace (not RTSL)
dt_tfa_bpp_nortsl[, tfa_rtsl_include := "No"]

# Rbind both files and keep unique locations

dt_tfa_bpp <- rbind(dt_tfa_bpp,dt_tfa_bpp_nortsl,fill = TRUE)

# Remove missing location
dt_tfa_bpp <- dt_tfa_bpp[!is.na(Country),]

# rename locations to match gbd locations (the baseline)
name_map <- c(
  "Czech Republic"                    = "Czechia",
  "Macedonia"                         = "North Macedonia",
  "Moldova"                           = "Republic of Moldova",
  "Republic of Cyprus"                = "Cyprus")

dt_tfa_bpp[, Country := fcoalesce(name_map[Country], Country)]

# Vector of N=67 locs of interest

locs_tfa <- unique(dt_tfa_bpp$Country)

# Check all dt_tfa_bpp location included in scenarios

check_locs <- setdiff(dt_tfa_bpp$Country, dt_tfa_scenarios$location)

#? Leichtein has no tfa data but policy
# So replicate Luxembourg policy data for Liechtenstein 
# since they are similar in many ways and Liechtenstein has no TFA data
# This wont be included in ccalculations, but in totalizing countries with policy

dt_tfa_scenarios_lei <- dt_tfa_scenarios[location == "Luxembourg",]

dt_tfa_scenarios_lei[, location := "Liechtenstein"]

dt_tfa_scenarios <- rbind(dt_tfa_scenarios, dt_tfa_scenarios_lei, fill = TRUE)

dt_tfa_scenarios <- merge(dt_tfa_scenarios, dt_tfa_bpp, by.x = "location", by.y= "Country", all.x = TRUE,all.y = F)

setnames(dt_tfa_scenarios, c("Date passed", "Date in Effect","RTSL Direct engagement"),
         c("tfa_bpp_date_passed", "tfa_bpp_date_in_effect","tfa_rtsl_engagement"))

# encode engagement as binary variable
dt_tfa_scenarios[, tfa_rtsl_engagement := ifelse(is.na(tfa_rtsl_engagement), "No", "Yes")]

# subset from year 2016 for counterfactual   
#dt_tfa_scenarios <- dt_tfa_scenarios[year>=2016,]

# subset countries of interest TFA bp passed 2016-2024
#dt_tfa_scenarios_counterfactual <- dt_tfa_scenarios[!is.na(tfa_bpp_date_passed) & year >= 2016 & year <= 2024 & tfa_current>0,]
dt_tfa_scenarios_counterfactual <- dt_tfa_scenarios[!is.na(tfa_bpp_date_passed) & year <= 2024 & tfa_current>0,]

# order most recent non zero observation by location and year
dt_tfa_scenarios_counterfactual <- dt_tfa_scenarios_counterfactual[order(location, year, decreasing = TRUE),]

# keep only one observation by location and only location year_recent and tfa_current columns
dt_tfa_scenarios_counterfactual <- dt_tfa_scenarios_counterfactual[, .SD[1], by = location, .SDcols = c("year","tfa_current")]

# rename tfa counterfacual and merge with dt_tfa scenarios
setnames(dt_tfa_scenarios_counterfactual, c("year","tfa_current"), c("year_counterfactual","tfa_counterfactual"))

dt_tfa_scenarios <- merge(dt_tfa_scenarios, dt_tfa_scenarios_counterfactual, by = "location", all.x = TRUE)

# create tfa_observed as tfa_current
dt_tfa_scenarios[, tfa_observed := tfa_current]

# impute tfa_current to counterfectual if year >= year_counterfactual
dt_tfa_scenarios[!is.na(year_counterfactual) & year >= year_counterfactual, tfa_current := tfa_counterfactual]

# Convert to percent scale
dt_tfa_scenarios[, tfa_current := tfa_current * 100]

# for this assessment target equals current
dt_tfa_scenarios[, tfa_target  := tfa_current]

# if year >= bpp in effect then tfa_target==0
dt_tfa_scenarios[!is.na(tfa_bpp_date_passed) & year>=tfa_bpp_date_in_effect, tfa_target := 0]

# Cases to highlight in presentation slideck

# 1) India with BPP passed and RTSL engagement and IHME data had errors
# 2) Colombia which passed and just started
# 3) Viet Nam no RTSL engagement and not passed
# 4) Thailand no RTSL engagement and passed

## Keep only countries of interest N= 67

dt_tfa_scenarios <- dt_tfa_scenarios[location %in% locs_tfa,]

# Save the data
saveRDS(dt_tfa_scenarios, file = paste0(wd_data,"tfa_policy_scenarios_assessment.rds"))

# cleaning: remove objects that are no longer needed
rm(dt_tfa_scenarios, dt_tfa_bpp, dt_tfa_scenarios_counterfactual,
   dt_tfa_scenarios_lei, check_locs, locs_tfa,dt_tfa_bpp_nortsl)


# #DT AI not included
# #dt_tfa <- fread(paste0(wd_raw,"TFAPolicy/","tfa_intake.csv"))
# # Extract first numeric value (including decimals)
# # dt_tfa[, tfa_current := as.numeric(stringr::str_extract(`Mean TFA Intake (%E)`, "\\d+\\.?\\d*"))]
# 
# # Input from IHME
# dt_tfa <- fread(paste0(wd_raw,"IHME_GBD_2023_DIET_RISK_1990_2024_TRANSFAT_Y2024M11D05",".csv"))
# 
# dt_tfa <- dt_tfa[,c("year_id","location_name","sex_name","age_group_name","val","upper","lower","age_group_id"),with = F]
# 
# setnames(dt_tfa,c("year_id","location_name","sex_name","age_group_name","val","upper","lower"),
#          c("year","location","sex","age","tfa_current","tfa_current_upper","tfa_current_lower"))
# 
# dt_tfa <- dt_tfa[location!="Global",]
# 
# # rename locations
# dt_tfa[location == "Türkiye", location := "Turkey"]
# dt_tfa[location == "Côte d'Ivoire", location := "Ivory Coast"]
# 
# # Only covers from 25 years and older, so imputing for 20-24 yearsas 25-29 years
# 
# dt_tfa_20 <- dt_tfa[age == "25 to 29",]
# dt_tfa_20[, age := "20 to 24"]
# dt_tfa <- rbind(dt_tfa, dt_tfa_20)
# 
# rm(dt_tfa_20)
# 
# # Countries with BPP
# 
# dt_tfa_bpp <- as.data.table(read_excel(paste0(wd_raw,"List of Countries with Policies Passed-updated March 2026.xlsx"), 
#                                        sheet = "Sheet1", range = "A3:D56")) 
# 
# # rename locations to match gbd locations (the baseline)
# name_map <- c(
#   "Czech Republic"                    = "Czechia",
#   "Macedonia"                         = "North Macedonia",
#   "Moldova"                           = "Republic of Moldova",
#   "Cyprus"                            = "Republic of Cyprus")
# 
# dt_tfa_bpp[, Country := fcoalesce(name_map[Country], Country)]
# 
# #? Leichtein has no tfa data but policy
# dt_tfa <- merge(dt_tfa, dt_tfa_bpp, by.x = "location", by.y= "Country", all.x = TRUE,all.y = F)
# 
# # check locations
# dt_check <- dt_tfa[is.na(sex),]
# 
# setnames(dt_tfa, c("Date passed", "Date in Effect","RTSL Direct engagement"),
#          c("tfa_bpp_date_passed", "tfa_bpp_date_in_effect","tfa_rtsl_engagement"))
# 
# #dt_tfa <- dt_tfa[, .(location, tfa_current, tfa_bpp_date_passed, tfa_bpp_date_in_effect)]
# 
# dt_tfa <- dt_tfa[sex!="Both",]
# 
# # encode engagement as binary variable
# dt_tfa[, tfa_rtsl_engagement := ifelse(is.na(tfa_rtsl_engagement), "No", "Yes")]
# 
# # Save the data
# saveRDS(dt_tfa, file = paste0(wd_data,"tfa_data.rds"))
# 
# # cleaning: remove objects that are no longer needed
# rm(dt_tfa, dt_tfa_bpp,dt_check)

#...........................................................
## TFA Policy----
#...........................................................

# dt_tfa_scenarios <- dt_tfa[, {
#   yearly_data <- data.table(year = years)
#   yearly_data[, tfa_current := tfa_current]
#   yearly_data[, tfa_target := ifelse(year<2025 | tfa_current< 0.5,tfa_current,0.5)] # Default value
#   yearly_data[year >= tfa_bpp_date_in_effect & year < 2025, tfa_target := tfa_current]
#   yearly_data
# }, by = location]
# 
# # TFA data
# dt_tfa <- readRDS(file = paste0(wd_data,"tfa_data.rds"))
# 
# name_map <- c(
#   "United States of America"          = "United States",
#   "Taiwan"                            = "Taiwan (Province of China)",
#   "Urkraine (without Crimea & Sevastopol)" = "Ukraine")
# 
# # P target
# p_tfa_target <- 0.00
# dt_tfa[, location := fcoalesce(name_map[location], location)]
# 
# # Computing mean intake per location
# 
# dt_pop <- readRDS(file = paste0("C:/Users/wrgar/OneDrive - UW/02Work/ResolveToSaveLives/100MLives/data/raw/GBD/","totalpop_ihme.rds"))
# 
# dt_pop <- dt_pop[, .(location, year_id, sex_name, age_group, val),with=T]
# 
# setnames(dt_pop, c("val", "year_id","sex_name"),
#          c("population","year","sex"))
# 
# name_map <- c(
#   "United States of America"          = "United States",
#   "Taiwan"                            = "Taiwan (Province of China)",
#   "Urkraine (without Crimea & Sevastopol)" = "Ukraine")
# 
# dt_pop[, location := fcoalesce(name_map[location], location)]
# 
# 
# age_match<-data.frame(age=20:95)%>%
#   mutate(age.group = ifelse(age<25, "20 to 24", NA),
#          age.group = ifelse(age>=25 & age<30, "25 to 29", age.group),
#          age.group = ifelse(age>=30 & age<35, "30 to 34", age.group),
#          age.group = ifelse(age>=35 & age<40, "35 to 39", age.group),
#          age.group = ifelse(age>=40 & age<45, "40 to 44", age.group),
#          age.group = ifelse(age>=45 & age<50, "45 to 49", age.group),
#          age.group = ifelse(age>=50 & age<55, "50 to 54", age.group),
#          age.group = ifelse(age>=55 & age<60, "55 to 59", age.group),
#          age.group = ifelse(age>=60 & age<65, "60 to 64", age.group),
#          age.group = ifelse(age>=65 & age<70, "65 to 69", age.group),
#          age.group = ifelse(age>=70 & age<75, "70 to 74", age.group),
#          age.group = ifelse(age>=75 & age<80, "75 to 79", age.group),
#          age.group = ifelse(age>=80 & age<85, "80 to 84", age.group),
#          age.group = ifelse(age>=85 & age<90, "85 to 89", age.group),
#          age.group = ifelse(age>=90 & age<95, "90 to 94", age.group),
#          age.group = ifelse(age==95, "95 plus", age.group))
# 
# age_match$age <- as.character(age_match$age)
# 
# dt_pop[age_group_name == "<1 year", age_group_name := "0"]
# dt_pop[age_group_name == "95 plus", age_group_name := "95"]
# 
# dt_pop <- merge(dt_pop, age_match, by.x = "age_group_name", by.y = "age", all.x = TRUE)
# 
# # Average population by age group (not year because POp comes up to 2019)
# dt_pop <- dt_pop[,list(population=mean(population)),by=list(location,sex,age.group)]
# 
# dt_tfa[,age.group:=age]
# 
# # Merge population data with TFA data
# dt_tfa <- merge(dt_tfa,dt_pop,
#                 by = c("location","sex","age.group"), all.x = TRUE,all.y = F)
# 
# # check not missing population . There are only 1 territories not included in our analysis
# #dt_check <- unique(dt_sodium[is.na(sodium_current),],by="location")
# dt_check <- unique(dt_tfa[is.na(population),],by="location")
# 
# 
# # compute mean TFA intake per location
# dt_tfa_mean <- dt_tfa[, .(tfa_current = weighted.mean(tfa_current,population, na.rm = TRUE),
#                           tfa_bpp_date_in_effect = min(tfa_bpp_date_in_effect, na.rm = TRUE),
#                           tfa_bpp_date_passed = min(tfa_bpp_date_passed, na.rm = TRUE)),
#                       by = list(year,location)]
# 
# #dt_tfa_mean_24 <- dt_tfa_mean[year==2024,]
# 
# #Most recent data from 2024:2019
# dt_tfa_mean_24 <- dt_tfa_mean[year %in% c(2024),]
# dt_tfa_mean_24 <- dt_tfa_mean_24[order(dt_tfa_mean_24$tfa_current,decreasing = T),]
# dt_tfa_mean_24 <- unique(dt_tfa_mean_24,by="location")
# 
# years <- 2025:2050
# dt_tfa_scenarios <- dt_tfa_mean_24[, {
#   yearly_data <- data.table(year = years)
#   yearly_data[, tfa_current := tfa_current]
#   yearly_data
# }, by = location]
# 
# dt_tfa_scenarios <- rbind(dt_tfa_mean,dt_tfa_scenarios,fill=T)
# 
# # Add the tfa_target to the main data table
# 
# dt_tfa_scenarios[,c("tfa_bpp_date_in_effect","tfa_bpp_date_passed"):=NULL]
# dt_tfa_scenarios <- merge(dt_tfa_scenarios, 
#                           dt_tfa_mean_24[, .(location, tfa_bpp_date_in_effect,tfa_bpp_date_passed)],
#                           by = "location", all.x = TRUE)
# 
# # Fix bpp date for specific locations
# dt_tfa_scenarios[ , tfa_target := tfa_current]
# dt_tfa_scenarios[ tfa_bpp_date_passed<=2023 & year>2023 & year<2025 , tfa_target := ifelse(tfa_current<p_tfa_target,tfa_current,p_tfa_target)]
# #dt_tfa_scenarios[ year>=2025 , tfa_target := 0.005]
# #dt_tfa_scenarios[ year>=2027 , tfa_target := ifelse(tfa_current<p_tfa_target,tfa_current,p_tfa_target)]
# 
# # ensure that once a country reaches 0, its tfa_target remain 0 thereafter.
# 
# # Set initial targets to current values
# dt_tfa_scenarios[, tfa_target := tfa_current]
# 
# # If TFA policy is passed before or in 2023, drop to target (0) by 2024
# dt_tfa_scenarios[
#   tfa_bpp_date_passed <= 2023 & year >= 2024, 
#   tfa_target := p_tfa_target
# ]
# 
# # For all countries, once TFA reaches 0, it cannot go up again
# dt_tfa_scenarios[, tfa_target := 
#                    ifelse(
#                      year >= min(year[tfa_target == p_tfa_target], na.rm = TRUE),
#                      p_tfa_target,
#                      tfa_target
#                    ),
#                  by = location
# ]
# 
# # Optional: also keep tfa_current at 0 after a country hits 0
# dt_tfa_scenarios[, tfa_current := 
#                    ifelse(
#                      year >= min(year[tfa_target == p_tfa_target], na.rm = TRUE),
#                      p_tfa_target,
#                      tfa_current
#                    ),
#                  by = location
# ]
# 
# dt_tfa_scenarios[ year>=2027 , tfa_target := ifelse(tfa_current<p_tfa_target,tfa_current,p_tfa_target)]
# 
# dt_tfa_scenarios <- dt_tfa_scenarios[, .(location, year, tfa_current,tfa_target)]
# 
# # Taiwan, Ukraine and other territories are NA, so assigning the minimum value
# dt_tfa_scenarios[is.na(tfa_target),  tfa_target := p_tfa_target]
# dt_tfa_scenarios[is.na(tfa_current), tfa_current := p_tfa_target]
# 
# # Check country names
# dt_tfa_scenarios[, location := gsub("United States of America", "United States", location)]
# dt_tfa_scenarios[location=="Bolivia (Plurinational State of)", location := "Bolivia"]
# 
# dt_tfa_scenarios[, location := gsub("Taiwan", "Taiwan (Province of China)", location)]
# 
# dt_tfa_scenarios[, location := gsub("United Republic of Tanzania", "Tanzania", location)]
# 
# # Fix countries names
# dt_tfa_scenarios[location == "Türkiye", location := "Turkey"]
# dt_tfa_scenarios[location == "Côte d'Ivoire", location := "Ivory Coast"]
# 
# #print length(locs_statins)
# locs_tfa <- unique(dt_tfa_scenarios$location)
# print(paste0("Number of locations with tfa data: ", length(locs_tfa)))
# 
# # Save the data
# saveRDS(dt_tfa_scenarios,file = paste0(wd_data,"TFAPolicy/", "tfa_policy_scenarios.rds"))
# 
# # clean up
# rm(dt_tfa, dt_pop, dt_tfa_mean, dt_tfa_mean_24,
#    years, dt_tfa_scenarios, name_map, dt_check, age_match)

