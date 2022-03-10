# FORMAT DATA FOR RSIENA/SAOM MODELS
#
# This file produces the specific data objects used in the RSiena/SAOM models
#


library(assertthat)
source("utils/utils.R")
source("utils/custom_trade_to_binary.R")

# TODO shift function for variables

tmp <- array()

shift_matrix <- function(mat, dim) {
  
}

# TODO implement different cutoff rule (from the paper)
# TODO put cut off rule in citations

# Load data --------------------------------------------------------------------
dpath <- data_path.get()
# BACI conventional trade data
baci <- readRDS(file.path(dpath, "out/baci_aggregated.rds"))
# SIPRI arms trade data
sipri <- readRDS(file.path(dpath, "out/sipri_tiv.rds"))
# existence of countries over time
load(file.path(dpath, "out/EX.RData"))
# country list
load(file.path(dpath, "out/country_list.RData"))
# polity score data
pol <- readRDS(file.path(dpath, "out/polity.rds"))
# GDP data
gdp <- readRDS(file.path(dpath, "out/gdp.rds"))
# alliances
load(file.path(dpath, "out/atop_alliance.RData"))
# national material capability index
nmc <- readRDS(file.path(dpath, "out/nmc_cinc.rds"))
load(file.path(dpath, "out/cdist.RData"))


# Meta specifications for SAOM data --------------------------------------------
periods <- 1995:2018
n_all <- nrow(country_list)


# Country existence ------------------------------------------------------------
exist <- rowSums(EX[, paste(periods)]) == length(periods)



# BACI trade data --------------------------------------------------------------
# to 3 dim array
baci_periods <- 1995:2019
trd <- array(
  0,
  dim = c(n_all, n_all, length(baci_periods)),
  dimnames = list(dimnames(baci[[1]])[[1]], dimnames(baci[[1]])[[2]],
                  paste(baci_periods))
)
for (t in 1:length(baci)) {
  trd[,,t] <- baci[[t]]
}
# filter existing countries over time frame
trd <- trd[exist, exist, paste(periods)]
for (year in dimnames(trd)[[3]]) {
  trd[,,year] <- custom_trade_to_binary(trd[,,year], type = "B", threshold = 0.01)
}


# SIPRI arms trade data --------------------------------------------------------
# to 3 dim array
sipri_periods <- 1950:2018
arm <- array(
  0,
  dim = c(n_all, n_all, length(sipri_periods)),
  dimnames = list(dimnames(sipri[[1]])[[1]], dimnames(sipri[[1]])[[2]],
                  paste(sipri_periods))
)
for (t in 1:length(sipri)) {
  arm[,,t] <- sipri[[t]]
}
# filter existing countries over time frame
arm <- 1*(arm>1)
arm_ex <- arm[exist, exist,]
arm_main <- arm_ex[,,paste(periods)]

arm_lag1 <- arm_ex[,,paste(periods-1)]
dimnames(arm_lag1)[[3]] <- paste(periods)

arm_lag2 <- arm_ex[,,paste(periods-2)]
dimnames(arm_lag2)[[3]] <- paste(periods)

arm_lag3 <- arm_ex[,,paste(periods-3)]
dimnames(arm_lag3)[[3]] <- paste(periods)

# drop last observation for lagged arms deals (RSiena specification)

arm_lag1 <- arm_lag1[,,1:(length(periods)-1)]
arm_lag2 <- arm_lag2[,,1:(length(periods)-1)]
arm_lag3 <- arm_lag3[,,1:(length(periods)-1)]

assert_that(all(arm_main[,,1:(length(periods)-1)] == arm_lag1[,,2:length(periods)]))
assert_that(all(arm_main[,,1:(length(periods)-2)] == arm_lag2[,,3:length(periods)]))
assert_that(all(arm_main[,,1:(length(periods)-3)] == arm_lag3[,,4:length(periods)]))


# Alliances --------------------------------------------------------------------
# TODO: there are 68 periods, not 69 like for the SIPRI data
# assumed that it is 1951:2018 BUT CHECK
allied_periods <- 1951:2018
allied <- array(
  0,
  dim = c(n_all, n_all, length(allied_periods)),
  dimnames = list(dimnames(atop_alliance[[1]])[[1]], dimnames(atop_alliance[[1]])[[2]],
                  paste(allied_periods))
)
for (t in 1:length(atop_alliance)) {
  allied[,,t] <- atop_alliance[[t]]
}
allied <- allied[exist, exist, paste(periods)]
# drop last observation (RSiena specification)
allied <- allied[,,1:(length(periods)-1)]

# GDP --------------------------------------------------------------------------
gdp_log <- log1p(gdp[exist, paste(periods)])


# Polity -----------------------------------------------------------------------
# drop last observation (RSiena specification)
pol <- pol[exist, paste(periods)]
abs_pol_diff_year <- function(year) abs(outer(pol[, paste(year)], pol[,paste(year)], "-"))
pol_diff <- sapply(periods, abs_pol_diff_year, simplify = "array")
dimnames(pol_diff)[[3]] <- paste(periods)
pol_diff <- pol_diff[,,1:(length(periods)-1)]
pol_diff_std <- pol_diff / 20




# National material capability -------------------------------------------------
nmc <- nmc[exist, paste(periods)]
nmc_std <- nmc/sd(nmc, na.rm = T)


# Capital distance -------------------------------------------------------------
cdist <- cdist[exist, exist]
cdist_std <- cdist/sd(cdist)


# Save data --------------------------------------------------------------------
save(trd, arm_main, arm_lag1, arm_lag2, arm_lag3,
     gdp_log, pol_diff_std, nmc_std, cdist_std, allied,
     file = file.path(dpath, "out/saom_data_objects.RData"))
