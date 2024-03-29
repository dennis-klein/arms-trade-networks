# Fit the main multilevel SAOM model with between networks effects


require(RSiena)
require(fs)
source("scripts/SAOM/siena07ToConverge.R")
source("scripts/SAOM/siena07ToConvergeMulticore.R")
source("utils/utils.R")


##### LOAD DATA ----------------------------------------------------------------
dpath <- data_path.get()
load(file = path(dpath, "out/saom_data.RData"))


##### FIT SAOM MODEL -----------------------------------------------------------

# meta settings
dpath <- data_path.get()
model_id <- "saom_multilevel_220511_import"
save_dir <- path(dpath, "models", "SAOM")

# setup
act <- dimnames(arm)[[1]] # actors
obs <- length(arm) # observations

# test mode
test_mode <- FALSE
if (test_mode) {
  model_id <- paste0("test_", as.integer(Sys.time()))
  act <- c(sample(act, 15), c("United States", "Russia", "France", "China", "Germany", "Italy",
              "United Kingdom", "South Korea",
              "Saudi Arabia", "India", "Egypt", "Australia",
              "China", "Algeria", "South Korea"))
  # act <- c(sample(act, 10), c(4, 8, 21, 31, 36, 39, 47, 89, 91, 108, 109, 110))
  per <- paste(1998:2007)
  
  arm <- arm[act, act, per]
  trd <- trd[act, act, per]
  gdp_log <- gdp_log[act, head(per, -1)]
  mil_log <- mil_log[act, head(per, -1)]
  cdist <- cdist[act, act]
  pol_diff <- pol_diff[act, act, head(per, -1)]
  allied <- allied[act, act, head(per, -1)]
}

# safety check to not override existing models
if (dir.exists(path(save_dir, model_id))) stop("Model already exsists")

# create directory
model_dir <- path(save_dir, model_id)
dir.create(model_dir)

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
)
dat <- do.call(sienaDataCreate, RSiena_vars)

# save model meta information
model_meta <- list(
  actors = act,
  years = names(arm),
  RSiena_data = dat
)
saveRDS(model_meta, path(model_dir, "meta_info", ext = "rds"))

# define effects
eff <- getEffects(dat)
## within networks effects
eff <- includeEffects(eff, density, name = "arm", verbose = FALSE)
eff <- includeEffects(eff, density, name = "trd", verbose = FALSE)
eff <- includeEffects(eff, recip, name = "arm", verbose = FALSE)
eff <- includeEffects(eff, recip, name = "trd", verbose = FALSE)
eff <- includeEffects(eff, gwespFF, name = "arm", verbose = FALSE)
eff <- includeEffects(eff, gwespFF, name = "trd", verbose = FALSE)
eff <- includeEffects(eff, transRecTrip, name = "arm", verbose = FALSE)
eff <- includeEffects(eff, transRecTrip, name = "trd", verbose = FALSE)
eff <- includeEffects(eff, cycle3, name = "arm", verbose = FALSE)
eff <- includeEffects(eff, cycle3, name = "trd", verbose = FALSE)
eff <- includeEffects(eff, inPopSqrt, name = "arm", verbose = FALSE)
eff <- includeEffects(eff, inPopSqrt, name = "trd", verbose = FALSE)
eff <- includeEffects(eff, outActSqrt, name = "arm", verbose = FALSE)
eff <- includeEffects(eff, outActSqrt, name = "trd", verbose = FALSE)
eff <- includeEffects(eff, inActSqrt, name = "arm", verbose = FALSE)
eff <- includeEffects(eff, inActSqrt, name = "trd", verbose = FALSE)

# covariate effects
eff <- includeEffects(eff, egoX, altX, name = "arm", interaction1 = "gdp_log", verbose = FALSE)
eff <- includeEffects(eff, egoX, altX, name = "arm", interaction1 = "mil_log", verbose = FALSE)
eff <- includeEffects(eff, X, name = "arm", interaction1 = "pol_diff", verbose = FALSE)
eff <- includeEffects(eff, X, name = "arm", interaction1 = "cdist", verbose = FALSE)
eff <- includeEffects(eff, X, name = "arm", interaction1 = "allied", verbose = FALSE)
eff <- includeEffects(eff, egoX, altX, name = "trd", interaction1 = "gdp_log", verbose = FALSE)
eff <- includeEffects(eff, egoX, altX, name = "trd", interaction1 = "mil_log", verbose = FALSE)
eff <- includeEffects(eff, X, name = "trd", interaction1 = "pol_diff", verbose = FALSE)
eff <- includeEffects(eff, X, name = "trd", interaction1 = "cdist", verbose = FALSE)
eff <- includeEffects(eff, X, name = "trd", interaction1 = "allied", verbose = FALSE)

## between networks effects
eff <- includeEffects(eff, crprod, name = "arm", interaction1 = "trd", verbose = FALSE)
eff <- includeEffects(eff, crprod, name = "trd", interaction1 = "arm", verbose = FALSE)
eff <- includeEffects(eff, crprodRecip, name = "arm", interaction1 = "trd", verbose = FALSE)
eff <- includeEffects(eff, crprodRecip, name = "trd", interaction1 = "arm", verbose = FALSE)

# run model
alg <- sienaAlgorithmCreate(projname = model_id,
                            n3 = 5000)
ans <- siena07ToConvergenceMulticore(alg = alg, dat = dat, eff = eff,
                                     save_dir = model_dir,
                                     ans_id = model_id,
                                     threshold = 0.30,
                                     batch = TRUE)
