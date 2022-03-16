# FORMAT DATA FOR RSIENA/SAOM MODELS
#
# This file produces the specific data objects used in the RSiena/SAOM models
#

library(Matrix)
source("utils/utils.R")
source("scripts/preprocessing/trade_flow_gdp_cutoff.R")

# Preliminaries ----------------------------------------------------------------
dpath <- data_path.get()
model_years <- 1998:2018


# Country list -----------------------------------------------------------------
load(file.path(dpath, "out/country_list.RData"))


# existence of countries over time ---------------------------------------------
load(file.path(dpath, "out/EX.RData"))
ex <- rowSums(EX[, paste(model_years)]) == length(model_years)


# GDP data ---------------------------------------------------------------------
gdp <- readRDS(file.path(dpath, "out/gdp.rds"))
gdp <- gdp[ex,]
gdp_list <- asplit(gdp, MARGIN = 2)
gdp_log <- log1p(gdp)


# BACI trade data --------------------------------------------------------------
baci <- readRDS(file.path(dpath, "out/baci_aggregated.rds"))
baci <- lapply(baci, function(x) x[ex, ex])
names(baci) <- paste(1995:2019)
trd <- mapply(function(x, y) trade_flow_gdp_cutoff(trd = x, gdp = y,
                                                   threshold = 0.01),
                baci[paste(1995:2018)],
                gdp_list[paste(1995:2018)],
                SIMPLIFY = FALSE)


# SIPRI arms trade data --------------------------------------------------------
sipri <- readRDS(file.path(dpath, "out/sipri_tiv.rds"))
sipri <- lapply(sipri, function(x) x[ex, ex])
arm <- lapply(sipri, function(x) 1*(x>1))
sipri_periods <- 1950:2018
names(arm) <- paste(sipri_periods)


# Alliances --------------------------------------------------------------------
# TODO: there are 68 periods, not 69 like for the SIPRI data
# assumed that it is 1951:2018 BUT CHECK
load(file.path(dpath, "out/atop_alliance.RData"))
allied <- lapply(atop_alliance, function(x) x[ex,ex])
allied_periods <- 1951:2018
names(allied) <- allied_periods


# Polity -----------------------------------------------------------------------
pol <- readRDS(file.path(dpath, "out/polity.rds"))
pol <- pol[ex, paste(model_years)]

pol_diff <- lapply(asplit(pol, MARGIN = 2),
             FUN = function(x) abs(outer(x, x, "-")))
pol_diff_std <- lapply(pol_diff, function(x) x/20)


# National material capability -------------------------------------------------
nmc <- readRDS(file.path(dpath, "out/nmc_cinc.rds"))
nmc_std <- apply(nmc, MARGIN = 2,
             function(x) x/sd(x, na.rm = T))


# Capital distance -------------------------------------------------------------
load(file.path(dpath, "out/cdist.RData"))
cdist <- cdist[ex, ex]
cdist_std <- cdist/sd(cdist)


# Calculate lags ---------------------------------------------------------------
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


# filter for years -------------------------------------------------------------
year_filter <- paste(model_years)
arm <- arm[year_filter]
arm_lag1 <- arm_lag1[year_filter]
arm_lag2 <- arm_lag2[year_filter]
arm_lag3 <- arm_lag3[year_filter]
trd <- trd[year_filter]
trd_lag1 <- trd_lag1[year_filter]
trd_lag2 <- trd_lag2[year_filter]
trd_lag3 <- trd_lag3[year_filter]


# Sparse matrices --------------------------------------------------------------
arm <- lapply(arm, function(x) as(x, "dgTMatrix"))
arm_lag1 <- lapply(arm_lag1, function(x) as(x, "dgTMatrix"))
arm_lag2 <- lapply(arm_lag2, function(x) as(x, "dgTMatrix"))
arm_lag3 <- lapply(arm_lag3, function(x) as(x, "dgTMatrix"))

trd <- lapply(trd, function(x) as(x, "dgTMatrix"))
trd_lag1 <- lapply(trd_lag1, function(x) as(x, "dgTMatrix"))
trd_lag2 <- lapply(trd_lag2, function(x) as(x, "dgTMatrix"))
trd_lag3 <- lapply(trd_lag3, function(x) as(x, "dgTMatrix"))

pol_diff_std <- lapply(pol_diff_std, function(x) as(x, "dgTMatrix"))
allied <- lapply(allied, function(x) as(x, "dgTMatrix"))



# robustness checks density ----------------------------------------------------
all(arm[[3]] == arm_lag1[[4]])
any(arm[[3]] != arm_lag1[[5]])
sum(arm_lag2[[18]]) > 10
lapply(trd, function(x) sum(x)/length(x))


# Save data --------------------------------------------------------------------
save(trd, trd_lag1, trd_lag2, trd_lag3,
     arm, arm_lag1, arm_lag2, arm_lag3,
     gdp_log, pol_diff_std, nmc_std, cdist_std, allied,
     file = file.path(dpath, "out/saom_data_objects.RData"))
