library(RSiena)
library(fs)
source("utils/utils.R")
source("scripts/SAOM/fit_saom_pooled.R")

# Fit model --------------------------------------------------------------------
dpath <- data_path.get()
# fit_saom_pooled(model_id = paste0("test_tmp_", as.integer(Sys.time())),
#                 data_obj_path = path(dpath, "out/saom_data_objects", ext = "RData"),
#                 save_dir = path(dpath, "models", "SAOM"),
#                 test_mode = TRUE)

fit_saom_pooled(model_id = paste0("saom_pooled_220321"),
                data_obj_path = path(dpath, "out/saom_data_objects", ext = "RData"),
                save_dir = path(dpath, "models", "SAOM"),
                test_mode = FALSE)

# TODO there are goodness of fits test for SAOMS
# TODO error warning: only increases/decreases in certain years
# TODO decrease convergence ratio requirement