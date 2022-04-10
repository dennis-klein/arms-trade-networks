require(RSiena)
require(fs)
source("scripts/SAOM/siena07ToConverge.R")
source("scripts/SAOM/siena07ToConvergeMulticore.R")
source("scripts/SAOM/format_data_saom.R")
source("utils/utils.R")



# TODO check the density of last 3 years conventional trade
# TODO create non-symmetric trade network -> otherwise problems with the coefficients
# TODO "only increase" error warning: check section 4.2.4 (outdegree might be dropped)
# check if this is the reason for SE and conv ratio NAs

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
model_id <- "model_220326"
# model_id <- paste0("test_", as.integer(Sys.time()))
save_dir <- path(dpath, "models", "SAOM")

# safety check to not override existing models
if (dir.exists(path(save_dir, model_id))) stop("Model already exsists")

# setup
saom_data <- format_data_saom()
remove(list = names(saom_data))
attach(saom_data)
model_dir <- path(save_dir, model_id)
dir.create(model_dir)
actors <- dimnames(arm)[[1]] # actors
obs <- length(arm) # observations
act <- 1:length(actors)


# test mode
test_mode <- FALSE
if (test_mode) {
  # act <- c(4, 8, 21, 31, 36, 39, 47, 89, 91, 108, 109, 110)
  act <- c(sample(act, 10), c(4, 8, 21, 31, 36, 39, 47, 89, 91, 108, 109, 110))
  per <- paste(1998:2007)
  
  arm <- arm[act, act, per]
  trd <- trd[act, act, per]
  gdp_log <- gdp_log[act, head(per, -1)]
  mil_log <- mil_log[act, head(per, -1)]
  cdist <- cdist[act, act]
  pol_diff <- pol_diff[act, act, head(per, -1)]
  allied <- allied[act, act, head(per, -1)]
  arm_last3 <- arm_last3[act, act, head(per, -1)]
  trd_last3 <- trd_last3[act, act, head(per, -1)]
  
  # arm <- lapply(arm, function(x) x[act, act])
  # trd <- lapply(trd, function(x) x[act, act])
  # gdp_log <- gdp_log[act,]
  # mil_log <- mil_log[act,]
  # cdist <- cdist[act, act]
  # pol_diff <- lapply(pol_diff, function(x) x[act, act])
  # allied <- lapply(allied, function(x) x[act, act])
  # arm_last3 <- lapply(arm_last3, function(x) x[act, act])
  # trd_last3 <- lapply(trd_last3, function(x) x[act, act])
}

# create RSiena data
RSiena_vars <- list(
  # Dependent variables
  arm = sienaNet(arm),
  trd = sienaNet(trd),
  # arm = sienaNet(arm, allowOnly=FALSE),
  # trd = sienaNet(trd, allowOnly=FALSE),
  
  # Constant covariates
  # -
  
  # Varying covariates
  gdp_log = varCovar(gdp_log),
  mil_log = varCovar(mil_log),
  
  # Constant dyadic covariates
  cdist = coDyadCovar(cdist),
  
  # Varying dyadic covariates
  pol_diff = varDyadCovar(pol_diff),
  allied = varDyadCovar(allied)
  # arm_last3 = varDyadCovar(arm_last3),
  # trd_last3 = varDyadCovar(trd_last3)
)
dat <- do.call(sienaDataCreate, RSiena_vars)

# save model meta information
model_meta <- list(
  actors = actors,
  years = names(arm),
  RSiena_data = dat
)
saveRDS(model_meta, path(model_dir, "meta_info", ext = "rds"))

# define effects
eff <- getEffects(dat)
## within networks effects
# density, recip included (recip not possible for trd)
eff <- setEffect(eff, recip, name = "arm", include = FALSE, verbose = FALSE)
eff <- includeEffects(eff, transTies, name = "trd", verbose = FALSE)
# eff <- includeEffects(eff, gwespFF, name = "arm", verbose = FALSE)
# eff <- includeEffects(eff, inPopSqrt, name = "arm", verbose = FALSE)
# eff <- includeEffects(eff, inPopSqrt, name = "trd", verbose = FALSE)
# eff <- includeEffects(eff, outPopSqrt, name = "arm", verbose = FALSE)
# covariate effects
eff <- includeEffects(eff, egoX, altX, name = "arm", interaction1 = "gdp_log", verbose = FALSE)
eff <- includeEffects(eff, egoX, altX, name = "arm", interaction1 = "mil_log", verbose = FALSE)
eff <- includeEffects(eff, X, name = "arm", interaction1 = "pol_diff", verbose = FALSE)
eff <- includeEffects(eff, X, name = "arm", interaction1 = "cdist", verbose = FALSE)
eff <- includeEffects(eff, X, name = "arm", interaction1 = "allied", verbose = FALSE)
eff <- includeEffects(eff, X, name = "arm", interaction1 = "arm_last3", verbose = FALSE)
eff <- includeEffects(eff, X, name = "arm", interaction1 = "trd_last3", verbose = FALSE)
eff <- includeEffects(eff, egoX, name = "trd", interaction1 = "gdp_log", verbose = FALSE)
# eff <- includeEffects(eff, egoX, altX, name = "trd", interaction1 = "gdp_log", verbose = FALSE)
eff <- includeEffects(eff, altX, name = "trd", interaction1 = "mil_log", verbose = FALSE)
# eff <- includeEffects(eff, egoX, altX, name = "trd", interaction1 = "mil_log", verbose = FALSE)
eff <- includeEffects(eff, X, name = "trd", interaction1 = "pol_diff", verbose = FALSE)
eff <- includeEffects(eff, X, name = "trd", interaction1 = "cdist", verbose = FALSE)
eff <- includeEffects(eff, X, name = "trd", interaction1 = "allied", verbose = FALSE)
# eff <- includeEffects(eff, X, name = "trd", interaction1 = "arm_last3", verbose = FALSE)
# eff <- includeEffects(eff, X, name = "trd", interaction1 = "trd_last3", verbose = FALSE)
## between networks effects
eff <- includeEffects(eff, crprod, name = "arm", interaction1 = "trd", verbose = FALSE)
eff <- includeEffects(eff, crprod, name = "trd", interaction1 = "arm", verbose = FALSE)
#eff <- includeEffects(eff, crprodRecip, name = "arm", interaction1 = "trd", verbose = FALSE)
#eff <- includeEffects(eff, crprodRecip, name = "trd", interaction1 = "arm", verbose = FALSE)

# TODO: non-symmertric trade will also probably solve the problem with the collinearity of some effects on trade
# NOTE: excluded some effects, put them back in if possible
# TODO: check if it is actually ok to include the lag term (because of the decision interpretation)
# the lags are kind of considered because we model changes from the lags
# after many tries, there seems definitely some collinearity here...
# 6. eval arm: outdegree (density)            -74.7163 ( NA       )   -0.1021   
# 7. eval arm: reciprocity                     25.1812 ( NA       )    0.0000   
# 8. eval arm: GWESP I -> K -> J (69)          16.6241 ( NA       )    0.1498   
# 9. eval arm: indegree - popularity (sqrt)
# IDEA: indegree pop 
# found potential week link: reciprocity effect for arms
# gwespFF also does not seem to work for arms network
# TODO: mention in report why the endowment function does exactly what we want
# the design of the evaluation also means that the lags are not necessary
# NOTE: between refitting of the model, the effects are widely different
# NOTE: if trade is not same direction also some other effects like alter gdp on trade might work
# Try time dummy effects, test before

# run model
# alg <- sienaAlgorithmCreate(projname = model_id)
alg <- sienaAlgorithmCreate(projname = model_id, firstg = 0.01)
ans <- siena07ToConvergenceMulticore(alg = alg, dat = dat, eff = eff,
                                     save_dir = model_dir,
                                     ans_id = model_id,
                                     threshold = 0.25)
saveRDS(ans, file=path(model_dir, paste0("final_fit_", model_id), ext = "rds"))
