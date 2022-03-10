# POOLED SAOM MODEL
#
# This file implements a multi-layer SAOM model for the fixed time period
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
source("utils/custom_trade_to_binary.R")


# Load data --------------------------------------------------------------------
dpath <- data_path.get()
load(file = file.path(dpath, "out/saom_data_objects.RData"))





# Pooled RSiena/SAOM model -----------------------------------------------------
# define RSiena variables
dep_arm <- sienaNet(arm_main)
dep_trd <- sienaNet(trd)
cov_gdp_log <- varCovar(gdp_log)
cov_nmc <- varCovar(nmc_std)
cov_cdist <- coDyadCovar(cdist_std)
cov_pol <- varDyadCovar(pol_diff_std)
cov_allied <- varDyadCovar(allied)
cov_arm_lag1 <- varDyadCovar(arm_lag1)
cov_arm_lag2 <- varDyadCovar(arm_lag2)
cov_arm_lag3 <- varDyadCovar(arm_lag3)

data <- sienaDataCreate(dep_arm, dep_trd,
                        cov_gdp_log, cov_nmc, cov_cdist,
                        cov_pol, cov_allied,
                        cov_arm_lag1, cov_arm_lag2, cov_arm_lag3)



# define effects
eff <- getEffects(data)
## within networks effects
# standard effects (outdegree, reciprocity automatically added)
eff <- includeEffects(eff, gwespFF, name = "dep_arm")
eff <- includeEffects(eff, gwespFF, name = "dep_trd")
# covariate effects
eff <- includeEffects(eff, egoX, altX, name = "dep_arm", interaction1 = "cov_gdp_log")
eff <- includeEffects(eff, egoX, altX, name = "dep_arm", interaction1 = "cov_nmc")
eff <- includeEffects(eff, X, name = "dep_arm", interaction1 = "cov_pol")
eff <- includeEffects(eff, X, name = "dep_arm", interaction1 = "cov_cdist")
eff <- includeEffects(eff, X, name = "dep_arm", interaction1 = "cov_allied")
eff <- includeEffects(eff, X, name = "dep_arm", interaction1 = "cov_arm_lag1")
eff <- includeEffects(eff, X, name = "dep_arm", interaction1 = "cov_arm_lag2")
eff <- includeEffects(eff, X, name = "dep_arm", interaction1 = "cov_arm_lag3")


eff <- includeEffects(eff, egoX, altX, name = "dep_trd", interaction1 = "cov_gdp_log")
eff <- includeEffects(eff, egoX, altX, name = "dep_trd", interaction1 = "cov_nmc")
eff <- includeEffects(eff, X, name = "dep_trd", interaction1 = "cov_pol")
eff <- includeEffects(eff, X, name = "dep_trd", interaction1 = "cov_cdist")
eff <- includeEffects(eff, X, name = "dep_trd", interaction1 = "cov_allied")
eff <- includeEffects(eff, X, name = "dep_trd", interaction1 = "cov_arm_lag1")
eff <- includeEffects(eff, X, name = "dep_trd", interaction1 = "cov_arm_lag2")
eff <- includeEffects(eff, X, name = "dep_trd", interaction1 = "cov_arm_lag3")

## between networks effects
eff <- includeEffects(eff, crprod, name = "dep_arm", interaction1 = "dep_trd")
eff <- includeEffects(eff, crprod, name = "dep_trd", interaction1 = "dep_arm")
eff <- includeEffects(eff, crprodRecip, name = "dep_arm", interaction1 = "dep_trd")
eff <- includeEffects(eff, crprodRecip, name = "dep_trd", interaction1 = "dep_arm")

# run model
model_name <- "saom_pooled_220228"
alg <- sienaAlgorithmCreate(projname = model_name)
res <- siena07(alg, data = data, effects = eff)
# save SAOM
saveRDS(res, file.path(dpath, paste("models/SAOM/", model_name,".RDS", sep = "")))
