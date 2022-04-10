# FORMAT DATA FOR RSIENA/SAOM MODELS
#
# This file produces the specific data objects used in the RSiena/SAOM models
#


require(Matrix)
require(abind)
source("utils/utils.R")
source("scripts/preprocessing/trade_flow_gdp_cutoff.R")


format_data_saom <- function() {
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
  milit_exp_log <- log1p(milit_exp[ex, paste(model_years)])

  
  # BACI trade data ------------------------------------------------------------
  baci <- readRDS(file.path(dpath, "out/baci_aggregated.rds"))
  baci <- lapply(baci, function(x) x[ex, ex])
  names(baci) <- paste(1995:2019)
  trd <- mapply(function(x, y) trade_flow_gdp_cutoff(trd = x, gdp = y,
                                                     threshold = 0.01),
                baci[paste(1995:2018)],
                gdp_list[paste(1995:2018)],
                SIMPLIFY = FALSE)

  
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
  pol_diff_std <- lapply(pol_diff, function(x) x/20)

  
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
  pol_diff_std <- pol_diff_std[year_filter_dyad]
  
  
  # any arms / conventional trade in the last 3 years --------------------------
  arm_last3 <- mapply(function(l1, l2, l3) 1*((l1 == 1) | (l2 == 1) | (l3 == 1)),
                      arm_lag1, arm_lag2, arm_lag3, SIMPLIFY = FALSE)

  trd_last3 <- mapply(function(l1, l2, l3) 1*((l1 == 1) | (l2 == 1) | (l3 == 1)),
                      trd_lag1, trd_lag2, trd_lag3, SIMPLIFY = FALSE)

  
  # # Sparse matrices ------------------------------------------------------------
  # arm <- lapply(arm, function(x) as(x, "dgTMatrix"))
  # arm_lag1 <- lapply(arm_lag1, function(x) as(x, "dgTMatrix"))
  # arm_lag2 <- lapply(arm_lag2, function(x) as(x, "dgTMatrix"))
  # arm_lag3 <- lapply(arm_lag3, function(x) as(x, "dgTMatrix"))
  # arm_last3 <- lapply(arm_last3, function(x) as(x, "dgTMatrix"))
  # 
  # trd <- lapply(trd, function(x) as(x, "dgTMatrix"))
  # trd_lag1 <- lapply(trd_lag1, function(x) as(x, "dgTMatrix"))
  # trd_lag2 <- lapply(trd_lag2, function(x) as(x, "dgTMatrix"))
  # trd_lag3 <- lapply(trd_lag3, function(x) as(x, "dgTMatrix"))
  # trd_last3 <- lapply(trd_last3, function(x) as(x, "dgTMatrix"))
  # 
  # pol_diff_std <- lapply(pol_diff_std, function(x) as(x, "dgTMatrix"))
  # allied <- lapply(allied, function(x) as(x, "dgTMatrix"))
  
  # Convert to arrays ----------------------------------------------------------
  arm <- abind(arm, along = 3)
  trd <- abind(trd, along = 3)
  arm_last3 <- abind(arm_last3, along = 3)
  trd_last3<- abind(trd_last3, along = 3)
  pol_diff_std <- abind(pol_diff_std, along = 3)
  allied <- abind(allied, along = 3)

  
  # robustness checks density --------------------------------------------------
  all(arm[[3]] == arm_lag1[[4]])
  any(arm[[3]] != arm_lag1[[5]])
  sum(arm_lag2[[18]]) > 10
  lapply(trd, function(x) sum(x)/length(x))
  
  
  
  
  
  # collect data in object -----------------------------------------------------
  saom_data <- list()
  saom_data[["arm"]] <- arm
  saom_data[["trd"]] <- trd
  saom_data[["gdp_log"]] <- gdp_log
  saom_data[["mil_log"]] <- milit_exp_log
  saom_data[["allied"]] <- allied
  saom_data[["pol_diff"]] <- pol_diff_std
  saom_data[["nmc"]] <- nmc_std
  saom_data[["cdist"]] <- cdist_std
  saom_data[["arm_last3"]] <- arm_last3
  saom_data[["trd_last3"]] <- trd_last3
  
  return(saom_data)
}
