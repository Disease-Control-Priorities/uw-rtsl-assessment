# Required inputs----

#...........................................................
# HTN Impact----
#...........................................................

# Import program data 2017 -2025

# HTN global tracker 2017-2025_15Mar2026.xlsx

# Import population enrolled in hypertension control programs by location and year

# Reshape wide to long format


# Import population with BP control by location and year

# Reshape wide to long format

# Merge datasets to get control rates by location and year

# Computing growth rates in enrolled population and control rates by location

# Impute missing data based on growth rates






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

#...........................................................
# TFA Impact----
#...........................................................
