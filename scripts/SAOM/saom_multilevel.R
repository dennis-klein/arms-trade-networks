library(RSiena)
library(dplyr)
library(assertthat)

# specify path to data folder
source("utils/utils.R")
dpath <- data_path.get()


# read in data

# NOTE: code was written when baci_aggregated.rds was set to 0 if below treshold
sip <- readRDS(file.path(dpath, "out/sipri_tiv.rds")) # arms deals data 1950-2018
baci <- readRDS(file.path(dpath, "out/baci_aggregated.rds")) # trade data 1995-2019
load(file.path(dpath, "out/country_list.RData")) # country_list
load(file.path(dpath, "out/EX.RData")) # EX, existence of countries over time

# transform networks into arrays
N <- dim(sip[[1]])[1]
sip.net <- array(
  0, dim = c(N, N, length(sip)),
  dimnames = list(country_list$iso3, country_list$iso3, paste(1950:2018))
)
for (t in 1:dim(sip.net)[3]) {
  sip.net[,,t] <- sip[[t]]
}
baci.net <- array(
  0, dim = c(N, N, length(baci)),
  dimnames = list(country_list$iso3, country_list$iso3, paste(1995:2019))
)
for (t in 1:dim(baci.net)[3]) {
  baci.net[,,t] <- baci[[t]]
}

# truncate sip and baci to common years and existing countries
periods <- 1995:2005
exist <- apply(EX[,paste(periods)], MARGIN = 1, all)
sip.net <- sip.net[exist, exist, paste(periods)]
baci.net <- baci.net[exist, exist, paste(periods)]
sip.net.bin <- ifelse(sip.net > 0, 1, 0)
baci.net.bin <- ifelse(baci.net > 0, 1, 0)
mode(sip.net.bin) <- "integer"
mode(baci.net.bin) <- "integer"

# Siena model
# dependent variables / arms and trade networks
model.dep.arms <- sienaNet(sip.net.bin)
model.dep.trd <- sienaNet(baci.net.bin)
# model.dep.arms <- sienaDependent(netarray = sip.net.bin)
# model.dep.trd <- sienaDependent(netarray = baci.net.bin)

# Siena data object
model.data <- sienaDataCreate(
  model.dep.arms,
  model.dep.trd
)
#print01Report(model.data)

# Siena effects
model.eff <- getEffects(model.data)
#effectsDocumentation(model.eff)
#model.eff

# on-networks effects
model.eff <- includeEffects(model.eff, gwespFF, name = "model.dep.arms")
model.eff <- includeEffects(model.eff, gwespFF, name = "model.dep.trd")

# between-networks effects
model.eff <- includeEffects(model.eff, crprod, name = "model.dep.arms", interaction1 = "model.dep.trd")
model.eff <- includeEffects(model.eff, crprod, name = "model.dep.trd", interaction1 = "model.dep.arms")

# Siena algorithm
model.alg <- sienaAlgorithmCreate(projname = "saom_multi_base")

# run model
model.results <- siena07(model.alg, data = model.data, effects = model.eff)
saveRDS(model.results, file.path(dpath, "models/SAOM/saom_multi_base.RDS"))

# second run
#model.results <- siena07(model.alg, data = model.data, effects = model.eff, prevAns = model.results)
