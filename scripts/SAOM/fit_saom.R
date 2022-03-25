require(RSiena)
require(fs)
source("utils/utils.R")
source("scripts/SAOM/siena07ToConverge.R")
source("scripts/SAOM/siena07ToConvergeMulticore.R")


# TODO check the density of last 3 years conventional trade

# Fit a multilayer SAOM model to the arms and trade networks
#
# Network variables:
# - arms trade
# - conventional trade
#
# Covariates:
# - GDP (log)
# - absolut polity differences
# - Military expenditure (log)
# - capital distances
# - alliances
# - weapon trade in last 3 years
# - conventional trade in last 3 years

# meta settings
dpath <- data_path.get()
model_id <- "Test"
data_obj_path <- path(dpath, "out/saom_data_objects", ext = "RData")
save_dir <- path(dpath, "models", "SAOM")

# test mode
test_mode <- TRUE
model_id <- paste0("test_tmp_", as.integer(Sys.time()))

# safety check to not override existing models
if (dir.exists(path(save_dir, model_id))) stop("Model already exsists")

# Setup
load(data_obj_path)
model_dir <- path(save_dir, model_id)
dir.create(model_dir)
actors <- dimnames(arm[[1]])[[1]] # actors
obs <- length(arm) # observations
n <- length(actors)
act <- 1:n

# test_mode, actor subselection
if (test_mode) {
  message("TEST MODE")
  act_sub <- c(sample(act, 10), c(4, 8, 21, 31, 36, 39, 47, 89, 91, 108, 109, 110))
  actors <- actors[act_sub]
  act <- act_sub
}

# save model meta information
model_meta <- list(
  actors = actors,
  years = names(arm)
)
saveRDS(model_meta, path(model_dir, "meta_info", ext = "rds"))

RSiena_vars <- list(
  # Dependent variables
  arm = sienaNet(lapply(arm, function(x) x[act, act])),
  trd = sienaNet(lapply(trd, function(x) x[act, act])),
  
  # Constant covariates
  # -
  
  # Varying covariates
  gdp_log = varCovar(gdp_log[act,]),
  milit_exp_log = varCovar(milit_exp_log[act,]),
  
  # Constant dyadic covariates
  cdist_std = coDyadCovar(cdist_std[act, act]),
  
  # Varying dyadic covariates
  pol_diff_std = varDyadCovar(lapply(pol_diff_std, function(x) x[act, act])),
  allied = varDyadCovar(lapply(allied, function(x) x[act, act])),
  arm_last3 = varDyadCovar(lapply(arm_last3, function(x) x[act, act])),
  trd_last3 = varDyadCovar(lapply(trd_last3, function(x) x[act, act]))
)
dat <- do.call(sienaDataCreate, RSiena_vars)

# define effects
eff <- getEffects(dat)
## within networks effects
# standard effects (outdegree, reciprocity automatically added)
eff <- includeEffects(eff, transTies, name = "trd", verbose = FALSE)
eff <- includeEffects(eff, gwespFF, name = "arm", verbose = FALSE)
eff <- includeEffects(eff, inPopSqrt, name = "arm", verbose = FALSE)
eff <- includeEffects(eff, inPopSqrt, name = "trd", verbose = FALSE)
eff <- includeEffects(eff, outPopSqrt, name = "arm", verbose = FALSE)
# covariate effects
eff <- includeEffects(eff, egoX, altX, name = "arm", interaction1 = "gdp_log", verbose = FALSE)
eff <- includeEffects(eff, egoX, altX, name = "arm", interaction1 = "milit_exp_log", verbose = FALSE)
eff <- includeEffects(eff, X, name = "arm", interaction1 = "pol_diff_std", verbose = FALSE)
eff <- includeEffects(eff, X, name = "arm", interaction1 = "cdist_std", verbose = FALSE)
eff <- includeEffects(eff, X, name = "arm", interaction1 = "allied", verbose = FALSE)
eff <- includeEffects(eff, X, name = "arm", interaction1 = "arm_last3", verbose = FALSE)
eff <- includeEffects(eff, X, name = "arm", interaction1 = "trd_last3", verbose = FALSE)
eff <- includeEffects(eff, egoX, altX, name = "trd", interaction1 = "gdp_log", verbose = FALSE)
eff <- includeEffects(eff, egoX, altX, name = "trd", interaction1 = "milit_exp_log", verbose = FALSE)
eff <- includeEffects(eff, X, name = "trd", interaction1 = "pol_diff_std", verbose = FALSE)
eff <- includeEffects(eff, X, name = "trd", interaction1 = "cdist_std", verbose = FALSE)
eff <- includeEffects(eff, X, name = "trd", interaction1 = "allied", verbose = FALSE)
eff <- includeEffects(eff, X, name = "trd", interaction1 = "arm_last3", verbose = FALSE)
eff <- includeEffects(eff, X, name = "trd", interaction1 = "trd_last3", verbose = FALSE)
## between networks effects
eff <- includeEffects(eff, crprod, name = "arm", interaction1 = "trd", verbose = FALSE)
eff <- includeEffects(eff, crprod, name = "trd", interaction1 = "arm", verbose = FALSE)
#eff <- includeEffects(eff, crprodRecip, name = "arm", interaction1 = "trd", verbose = FALSE)
#eff <- includeEffects(eff, crprodRecip, name = "trd", interaction1 = "arm", verbose = FALSE)


# # test mode effects
# if (test_mode) {
#   eff <- getEffects(dat)
# }

# run model
alg <- sienaAlgorithmCreate(projname = model_id)
ans <- siena07ToConvergenceMulticore(alg = alg, dat = dat, eff = eff,
                                     save_dir = model_dir,
                                     ans_id = model_id,
                                     threshold = 0.25)
saveRDS(ans, file=path(model_dir, paste0("final_fit_", model_id), ext = "rds"))
