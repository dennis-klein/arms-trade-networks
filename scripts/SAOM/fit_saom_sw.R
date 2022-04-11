# Fit the extended SAOM model with between networks effects and
# time heterogeneities


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
model_id <- "model_sw_220411"
# model_id <- paste0(as.integer(Sys.time()))
save_dir <- path(dpath, "models", "SAOM")

# safety check to not override existing models
if (dir.exists(path(save_dir, model_id))) stop("Model already exsists")

# setup
model_dir <- path(save_dir, model_id)
dir.create(model_dir)
act <- dimnames(arm)[[1]] # actors
obs <- length(arm) # observations

# test mode
test_mode <- FALSE
if (test_mode) {
  act <- c(sample(act, 15), c("United States", "Russia", "France", "China", "Germany", "Italy",
                              "United Kingdom", "South Korea",
                              "Saudi Arabia", "India", "Egypt", "Australia",
                              "China", "Algeria", "South Korea"))
  # act <- c(sample(act, 10), c(4, 8, 21, 31, 36, 39, 47, 89, 91, 108, 109, 110))
  per <- paste(1998:2005)
  
  arm <- arm[act, act, per]
  trd <- trd[act, act, per]
  gdp_log <- gdp_log[act, per]
  mil_log <- mil_log[act, per]
  cdist <- cdist[act, act]
  pol_diff <- pol_diff[act, act, per]
  allied <- allied[act, act, per]
}


# Loop: fit sliding windows SAOMs ----------------------------------------------

periods <- dimnames(arm)[[3]]
win_size <- 4
for (t in 1:(length(periods)-win_size+1)) {
  sub_model_id <- sprintf("win_%02d", t)
  win_per <- periods[t:(t+win_size-1)]
  vars <- list(
    # Dependent variables
    arm = sienaNet(arm[,,win_per]),
    trd = sienaNet(trd[,,win_per]),
    
    # Varying covariates
    gdp_log = varCovar(gdp_log[,win_per]),
    mil_log = varCovar(mil_log[,win_per]),
    
    # Constant dyadic covariates
    cdist = coDyadCovar(cdist),
    
    # Varying dyadic covariates
    pol_diff = varDyadCovar(pol_diff[,,head(win_per, -1)]),
    allied = varDyadCovar(allied[,,head(win_per, -1)])
  )
  dat <- do.call(sienaDataCreate, vars)
  
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
  
  ## covariate effects
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
  
  ## fitting
  alg <- sienaAlgorithmCreate()
  ans <- siena07ToConvergenceMulticore(alg = alg, dat = dat, eff = eff,
                                       save_dir = model_dir,
                                       ans_id = paste(model_id, sub_model_id, sep = "_"),
                                       threshold = 0.40)
}
