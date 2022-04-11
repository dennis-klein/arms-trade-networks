# Create all data objects from raw / preprocessed data
# necessary for the different SAOM models

require(fs)
require(abind)
source("scripts/preprocessing/trade_gdp_relevance_unidir.R")
source("utils/utils.R")


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
nmc <- apply(nmc, MARGIN = 2,
                 function(x) x/sd(x, na.rm = T))


# Capital distance -----------------------------------------------------------
load(file.path(dpath, "out/cdist.RData"))
cdist <- cdist[ex, ex]
cdist <- cdist/sd(cdist)


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
nmc <- nmc[, year_filter]
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


# Save to a single object ------------------------------------------------------
save(arm, trd,
     gdp_log, mil_log,
     cdist, pol_diff, allied,
     file = path(dpath, "out/saom_data.RData"))
