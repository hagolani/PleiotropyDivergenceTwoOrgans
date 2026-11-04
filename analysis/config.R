# Shared paths for all analysis scripts. Paths are relative to the repository
# root (found with the 'here' package), so scripts can be run from anywhere.
data_raw  <- here::here("data", "raw")         # simulation output
data_proc <- here::here("data", "processed")   # files built from the raw output
