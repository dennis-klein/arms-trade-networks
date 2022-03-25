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
# - absolute polity differences
# - Military expenditure (log)
# - capital distances
# - alliances
# - weapon trade in last 3 years
# - conventional trade in last 3 years

# meta settings
dpath <- data_path.get()
model_id <- paste0("test_tmp_", as.integer(Sys.time()))
data_obj_path <- path(dpath, "out/saom_data_objects", ext = "RData")
save_dir <- path(dpath, "models", "SAOM")

# safety check to not override existing models
if (dir.exists(path(save_dir, model_id))) stop("Model already exsists")

# setup
load(data_obj_path)
model_dir <- path(save_dir, model_id)
dir.create(model_dir)
actors <- dimnames(arm[[1]])[[1]] # actors
obs <- length(arm) # observations
act <- 1:length(actors)

# save model meta information
model_meta <- list(
  actors = actors,
  years = names(arm)
)
saveRDS(model_meta, path(model_dir, "meta_info", ext = "rds"))

# test mode
test_mode <- TRUE
if (test_mode) {
  act <- c(sample(act, 10), c(4, 8, 21, 31, 36, 39, 47, 89, 91, 108, 109, 110))
  
  arm <- lapply(arm, function(x) x[act, act])
  trd <- lapply(trd, function(x) x[act, act])
  gdp_log <- gdp_log[act,]
  milit_exp_log <- milit_exp_log[act,]
  cdist_std <- cdist_std[act, act]
  pol_diff_std <- lapply(pol_diff_std, function(x) x[act, act])
  allied <- lapply(allied, function(x) x[act, act])
  arm_last3 <- lapply(arm_last3, function(x) x[act, act])
  trd_last3 <- lapply(trd_last3, function(x) x[act, act])
}

RSiena_vars <- list(
  # Dependent variables
  arm = sienaNet(arm),
  trd = sienaNet(trd),
  
  # Constant covariates
  # -
  
  # Varying covariates
  gdp_log = varCovar(gdp_log),
  milit_exp_log = varCovar(milit_exp_log),
  
  # Constant dyadic covariates
  cdist_std = coDyadCovar(cdist_std),
  
  # Varying dyadic covariates
  pol_diff_std = varDyadCovar(pol_diff_std),
  allied = varDyadCovar(allied),
  arm_last3 = varDyadCovar(arm_last3),
  trd_last3 = varDyadCovar(trd_last3)
)
dat <- do.call(sienaDataCreate, RSiena_vars)

# define effects
eff <- getEffects(dat)
## within networks effects
# density, recip included (recip not possible for trd)
eff <- includeEffects(eff, transTies, name = "trd", verbose = FALSE)
eff <- includeEffects(eff, gwespFF, name = "arm", verbose = FALSE)
eff <- includeEffects(eff, inPopSqrt, name = "arm", verbose = FALSE)
eff <- includeEffects(eff, inPopSqrt, name = "trd", verbose = FALSE)
eff <- includeEffects(eff, outPopSqrt, name = "arm", verbose = FALSE)
# covariate effects
eff <- includeEffects(eff, egoX, altX, name = "arm", interaction1 = "gdp_log", verbose = FALSE)
# eff <- includeEffects(eff, egoX, altX, name = "arm", interaction1 = "milit_exp_log", verbose = FALSE)
eff <- includeEffects(eff, X, name = "arm", interaction1 = "pol_diff_std", verbose = FALSE)
eff <- includeEffects(eff, X, name = "arm", interaction1 = "cdist_std", verbose = FALSE)
eff <- includeEffects(eff, X, name = "arm", interaction1 = "allied", verbose = FALSE)
eff <- includeEffects(eff, X, name = "arm", interaction1 = "arm_last3", verbose = FALSE)
eff <- includeEffects(eff, X, name = "arm", interaction1 = "trd_last3", verbose = FALSE)
eff <- includeEffects(eff, egoX, altX, name = "trd", interaction1 = "gdp_log", verbose = FALSE)
# eff <- includeEffects(eff, egoX, altX, name = "trd", interaction1 = "milit_exp_log", verbose = FALSE)
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


# run model
alg <- sienaAlgorithmCreate(projname = model_id)
ans <- siena07ToConvergenceMulticore(alg = alg, dat = dat, eff = eff,
                                     save_dir = model_dir,
                                     ans_id = model_id,
                                     threshold = 0.25)
saveRDS(ans, file=path(model_dir, paste0("final_fit_", model_id), ext = "rds"))
