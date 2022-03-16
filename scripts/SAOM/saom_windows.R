# SLIDING WINDOWS SAOM MODEL
#
# This file fits a sliding windows multi-layer SAOM model for the fixed time period
# of 1995:2018
#
# Network variables:
# - arms trade
# - conventional
#
# Covariates:
# - GDP (log)
# - polity differences
# - material capability
# - capital distances
# - alliances
# - lagged weapon trade (1,2,3)

library(RSiena)
source("utils/utils.R")
source("scripts/SAOM/fit_saom_sliding_windows.R")


# Fit model --------------------------------------------------------------------
dpath <- data_path.get()
fit_saom_sliding_windows(win_size = 4,
                         data_obj_path = path(dpath, "out/saom_data_objects", ext = "RData"),
                         save_dir = path(dpath, "models", "SAOM"),
                         model_id = paste0("test_tmp_", as.integer(Sys.time())),
                         test_mode = TRUE)

# TODO there are goodness of fits test for SAOMS
