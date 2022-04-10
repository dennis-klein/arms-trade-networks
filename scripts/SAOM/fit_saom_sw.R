# Script to fit the SAOM model


require(RSiena)
require(fs)
require(abind)
source("scripts/preprocessing/trade_gdp_relevance_unidir.R")
source("scripts/SAOM/siena07ToConverge.R")
source("scripts/SAOM/siena07ToConvergeMulticore.R")
source("utils/utils.R")

##### FORMAT DATA --------------------------------------------------------------

# Preliminaries --------------------------------------------------------------
dpath <- data_path.get()
model_years <- 1998:2017


# Country list ---------------------------------------------------------------
load(file.path(dpath, "out/country_list.RData"))


# existence of countries over time -------------------------------------------
load(file.path(dpath, "out/EX.RData"))
ex <- rowSums(EX[, paste(model_years)]) == length(model_years)


# GDP data -------------------------------------------------------------------
gdp <- readRDS(file.path(dpath, "out/gdp.rds"))
gdp <- gdp[ex,]
gdp_list <- asplit(gdp, MARGIN = 2)
gdp_log <- log1p(gdp)


# Military expenditure -------------------------------------------------------
load(file.path(dpath, "out/milit_exp.RData"))
mil_log <- log1p(milit_exp[ex, paste(model_years)])


# BACI trade data ------------------------------------------------------------
baci <- readRDS(file.path(dpath, "out/baci_aggregated.rds"))
baci <- lapply(baci, function(x) x[ex, ex])
names(baci) <- paste(1995:2019)
trd <- mapply(function(x,y) trade_gdp_relevance_unidir(trd = x, gdp = y,
                                                       threshold = 0.01),
              baci[paste(1995:2018)], gdp_list[paste(1995:2018)],
              SIMPLIFY = FALSE)
tmp <- sum(trd[["2000"]])/length(trd[["2000"]])


# SIPRI arms trade data ------------------------------------------------------
sipri <- readRDS(file.path(dpath, "out/sipri_tiv.rds"))
sipri <- lapply(sipri, function(x) x[ex, ex])
arm <- lapply(sipri, function(x) 1*(x>1))
sipri_periods <- 1950:2018
names(arm) <- paste(sipri_periods)


# Alliances ------------------------------------------------------------------
# TODO: there are 68 periods, not 69 like for the SIPRI data
# assumed that it is 1951:2018 BUT CHECK
load(file.path(dpath, "out/atop_alliance.RData"))
allied <- lapply(atop_alliance, function(x) x[ex,ex])
allied_periods <- 1951:2018
names(allied) <- allied_periods


# Polity ---------------------------------------------------------------------
pol <- readRDS(file.path(dpath, "out/polity.rds"))
pol <- pol[ex, paste(model_years)]

pol_diff <- lapply(asplit(pol, MARGIN = 2),
                   FUN = function(x) abs(outer(x, x, "-")))
pol_diff <- lapply(pol_diff, function(x) x/20)


# National material capability -----------------------------------------------
nmc <- readRDS(file.path(dpath, "out/nmc_cinc.rds"))
nmc <- nmc[ex,]
nmc_std <- apply(nmc, MARGIN = 2,
                 function(x) x/sd(x, na.rm = T))


# Capital distance -----------------------------------------------------------
load(file.path(dpath, "out/cdist.RData"))
cdist <- cdist[ex, ex]
cdist_std <- cdist/sd(cdist)


# Calculate lags -------------------------------------------------------------
list_lag <- function(x) {
  r <- head(x, -1)
  names(r) <- tail(names(x), -1)
  return(r)
}

arm_lag1 <- list_lag(arm)
arm_lag2 <- list_lag(arm_lag1)
arm_lag3 <- list_lag(arm_lag2)

trd_lag1 <- list_lag(trd)
trd_lag2 <- list_lag(trd_lag1)
trd_lag3 <- list_lag(trd_lag2)


# filter for years -----------------------------------------------------------
year_filter <- paste(model_years)
year_filter_dyad <- paste(head(model_years, -1))
arm <- arm[year_filter]
arm_lag1 <- arm_lag1[year_filter_dyad]
arm_lag2 <- arm_lag2[year_filter_dyad]
arm_lag3 <- arm_lag3[year_filter_dyad]
trd <- trd[year_filter]
trd_lag1 <- trd_lag1[year_filter_dyad]
trd_lag2 <- trd_lag2[year_filter_dyad]
trd_lag3 <- trd_lag3[year_filter_dyad]
gdp_log <- gdp_log[, year_filter]
nmc_std <- nmc_std[, year_filter]
allied <- allied[year_filter_dyad]
pol_diff <- pol_diff[year_filter_dyad]


# any arms / conventional trade in the last 3 years --------------------------
arm_last3 <- mapply(function(l1, l2, l3) 1*((l1 == 1) | (l2 == 1) | (l3 == 1)),
                    arm_lag1, arm_lag2, arm_lag3, SIMPLIFY = FALSE)

trd_last3 <- mapply(function(l1, l2, l3) 1*((l1 == 1) | (l2 == 1) | (l3 == 1)),
                    trd_lag1, trd_lag2, trd_lag3, SIMPLIFY = FALSE)


# Convert to arrays ----------------------------------------------------------
arm <- abind(arm, along = 3)
trd <- abind(trd, along = 3)
arm_last3 <- abind(arm_last3, along = 3)
trd_last3<- abind(trd_last3, along = 3)
pol_diff <- abind(pol_diff, along = 3)
allied <- abind(allied, along = 3)


# robustness checks density --------------------------------------------------
all(arm[[3]] == arm_lag1[[4]])
any(arm[[3]] != arm_lag1[[5]])
sum(arm_lag2[[18]]) > 10



##### FIT SAOM MODEL -----------------------------------------------------------

# TODO check the density of last 3 years conventional trade
# TODO create non-symmetric trade network -> otherwise problems with the coefficients
# TODO "only increase" error warning: check section 4.2.4 (outdegree might be dropped)
# check if this is the reason for SE and conv ratio NAs

# meta settings
dpath <- data_path.get()
model_id <- "model_220407"
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

# TODO: maybe the SE get smaller when using window of 4
