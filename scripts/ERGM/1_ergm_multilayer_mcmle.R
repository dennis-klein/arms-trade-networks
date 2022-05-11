#------------------------------------------------------------------------------#
# Replication code for the ERGM Analysis
# MCMLE Estimation
#------------------------------------------------------------------------------#

library(network)
library(ergm)
library(multilayer.ergm)
library(texreg)
library(kableExtra)
library(stargazer)
library(coda)
library(plyr)
library(dplyr)
library(tidyr)
library(scales)


rm(list = ls(all.names = TRUE))


# setup
source("utils/utils.R")
source("utils/custom_trade_to_binary.R")
path <- data_path.get()
set.seed(1234)



#------------------------------------------------------------------------------#
# Data
#------------------------------------------------------------------------------#

# load data 
load(file.path(path, "out/EX.RData"))
load(file.path(path, "out/colony.RData"))
load(file.path(path, "out/country_list.RData"))
load(file.path(path, "out/cdist.RData"))
load(file.path(path, "out/atop_alliance.RData"))
load(file.path(path, "out/milit_exp.RData"))
nmc_cinc <- readRDS(file.path(path, "out/nmc_cinc.rds"))
polity <- readRDS(file.path(path, "out/polity.rds"))
gdp <- readRDS(file.path(path, "out/gdp.rds"))
gdppc <- readRDS(file.path(path, "out/gdppc.rds"))
sipri_tiv <- readRDS(file.path(path, "out/sipri_tiv.rds"))
trade <- readRDS(file.path(path, "out/baci_aggregated.rds"))


# select years of analysis
year <- 2003
start <- 1995
end <- 2018
maxlag <- 1


# selection of countries included: present over the complete time horizon (cf. soam)
included <- rowSums(EX[, (start:end) - 1949]) == length(start:end)
n <- sum(included)


# some data transformations
# sipri military expenditure is in constant USD$m
# cinc is index which sums up to 0 -> better interpretability using pct points
milit_exp <- milit_exp * 1000000
nmc_cinc100 <- nmc_cinc * 100


# setup constraints
free <- matrix(c(rep(c(rep(1, n), rep(0, n)), n), rep(c(rep(0, n), rep(1, n)), n)), nrow = 2 * n, byrow = T)
diag(free) <- 0
free <- network(free, directed = FALSE)


# prepare data 
# set correct indices
i1 = year - 1949 # Cornelius Replication Data 1950:2018
i2 = year - 1994 # CEPII BACI 1995:2018


# dependent multi-layer network
netC <- to.multiplex(
  1 * (sipri_tiv[[i1]][included, included] > 0),
  custom_trade_to_binary(trade[[i2]][included, included], type = "C", threshold = 0.01),
  output = "network"
)

netD <- to.multiplex(
  1 * (sipri_tiv[[i1]][included, included] > 0),
  custom_trade_to_binary(trade[[i2]][included, included], type = "D", threshold = 0.01),
  output = "network"
)


## include exogenous covariates lagged by 1 time period for each net

# gdppc out
log_gdppc_out <- matrix(log(gdppc[included, i1 - 1]), n, n, byrow = FALSE)

# gdppc in
log_gdppc_in <- matrix(log(gdppc[included, i1 - 1]), n, n, byrow = TRUE)

# gdp out
log_gdp_out <- matrix(log(gdp[included, i1 - 1]), n, n, byrow = FALSE)

# gdp in
log_gdp_in <- matrix(log(gdp[included, i1 - 1]), n, n, byrow = TRUE)

# atop alliance
alliance <- atop_alliance[[i1 - 1]][included, included]

# absolute polity diff
polity_diff <- abs(outer(polity[included, i1 - 1], polity[included, i1 - 1], "-"))

# military expenditure out
log_milit_exp_out <- matrix(log(milit_exp[included, i1 - 1]), n, n, byrow = FALSE)

# military expenditure in
log_milit_exp_in <- matrix(log(milit_exp[included, i1 - 1]), n, n, byrow = TRUE)

# define path dependency for each layer and once for import once export dependency
pathdep_layer1 <- 1 * (
  sipri_tiv[[i1 - 1]][included, included] +
    sipri_tiv[[i1 - 2]][included, included] +
    sipri_tiv[[i1 - 3]][included, included] > 0)
pathdep_layer2C <- 1 * (
  custom_trade_to_binary(trade[[i2 - 1]][included, included], type = "C", threshold = 0.01) +
    custom_trade_to_binary(trade[[i2 - 2]][included, included], type = "C", threshold = 0.01) +
    custom_trade_to_binary(trade[[i2 - 3]][included, included], type = "C", threshold = 0.01))
pathdep_layer2D <- 1 * (
  custom_trade_to_binary(trade[[i2 - 1]][included, included], type = "D", threshold = 0.01) +
    custom_trade_to_binary(trade[[i2 - 2]][included, included], type = "D", threshold = 0.01) +
    custom_trade_to_binary(trade[[i2 - 3]][included, included], type = "D", threshold = 0.01))

# distance
log_cdist <- log(cdist[included, included] + diag(n))



#------------------------------------------------------------------------------#
# Model Specification
#------------------------------------------------------------------------------#

# _nn indicates no interlayer network effects

modelC <- netC ~ edges_layer(layer = 1) +
  edges_layer(layer = 2) +
  mutual(same = "layer.mem", diff = TRUE) +
  gwidegree(decay = 1.5, fixed = TRUE, attr = "layer.mem") +
  gwodegree(decay = 1.5, fixed = TRUE, attr = "layer.mem") +
  gwesp_layer(decay = 0.69, fixed = TRUE, layer = 1) +
  gwesp_layer(decay = 0.69, fixed = TRUE, layer = 2) +
  duplexdyad(c("e", "f", "h"), layers = list(1, 2)) +
  edgecov_layer(log_cdist, layer = 1) +
  edgecov_layer(log_cdist, layer = 2) +
  edgecov_layer(log_gdp_out, layer = 1) +
  edgecov_layer(log_gdp_out, layer = 2) +
  edgecov_layer(log_gdp_in, layer = 1) +
  edgecov_layer(log_gdp_in, layer = 2) +
  edgecov_layer(alliance, layer = 1) +
  edgecov_layer(alliance, layer = 2) +
  edgecov_layer(polity_diff, layer = 1) +
  edgecov_layer(polity_diff, layer = 2) +
  edgecov_layer(log_milit_exp_out, layer = 1) +
  edgecov_layer(log_milit_exp_out, layer = 2) +
  edgecov_layer(log_milit_exp_in, layer = 1) +
  edgecov_layer(log_milit_exp_in, layer = 2) +
  edgecov_layer(pathdep_layer1, layer = 1) +
  edgecov_layer(pathdep_layer2C, layer = 2) 

modelC_nn <- netC ~
  edges_layer(layer = 1) +
  edges_layer(layer = 2) +
  mutual(same = "layer.mem", diff = TRUE) +
  gwidegree(decay = 1.5, fixed = TRUE, attr = "layer.mem") +
  gwodegree(decay = 1.5, fixed = TRUE, attr = "layer.mem") +
  gwesp_layer(decay = 0.69, fixed = TRUE, layer = 1) +
  gwesp_layer(decay = 0.69, fixed = TRUE, layer = 2) +
  edgecov_layer(log_cdist, layer = 1) +
  edgecov_layer(log_cdist, layer = 2) +
  edgecov_layer(log_gdp_out, layer = 1) +
  edgecov_layer(log_gdp_out, layer = 2) +
  edgecov_layer(log_gdp_in, layer = 1) +
  edgecov_layer(log_gdp_in, layer = 2) +
  edgecov_layer(alliance, layer = 1) +
  edgecov_layer(alliance, layer = 2) +
  edgecov_layer(polity_diff, layer = 1) +
  edgecov_layer(polity_diff, layer = 2) +
  edgecov_layer(log_milit_exp_out, layer = 1) +
  edgecov_layer(log_milit_exp_out, layer = 2) +
  edgecov_layer(log_milit_exp_in, layer = 1) +
  edgecov_layer(log_milit_exp_in, layer = 2) +
  edgecov_layer(pathdep_layer1, layer = 1) +
  edgecov_layer(pathdep_layer2C, layer = 2) 

modelD <- netD ~
  edges_layer(layer = 1) +
  edges_layer(layer = 2) +
  mutual(same = "layer.mem", diff = TRUE) +
  gwidegree(decay = 1.5, fixed = TRUE, attr = "layer.mem") +
  gwodegree(decay = 1.5, fixed = TRUE, attr = "layer.mem") +
  gwesp_layer(decay = 0.69, fixed = TRUE, layer = 1) +
  gwesp_layer(decay = 0.69, fixed = TRUE, layer = 2) +
  duplexdyad(c("e", "f", "h"), layers = list(1, 2)) +
  edgecov_layer(log_cdist, layer = 1) +
  edgecov_layer(log_cdist, layer = 2) +
  edgecov_layer(log_gdp_out, layer = 1) +
  edgecov_layer(log_gdp_out, layer = 2) +
  edgecov_layer(log_gdp_in, layer = 1) +
  edgecov_layer(log_gdp_in, layer = 2) +
  edgecov_layer(alliance, layer = 1) +
  edgecov_layer(alliance, layer = 2) +
  edgecov_layer(polity_diff, layer = 1) +
  edgecov_layer(polity_diff, layer = 2) +
  edgecov_layer(log_milit_exp_out, layer = 1) +
  edgecov_layer(log_milit_exp_out, layer = 2) +
  edgecov_layer(log_milit_exp_in, layer = 1) +
  edgecov_layer(log_milit_exp_in, layer = 2) +
  edgecov_layer(pathdep_layer1, layer = 1) +
  edgecov_layer(pathdep_layer2D, layer = 2) 

modelD_nn <- netD ~
  edges_layer(layer = 1) +
  edges_layer(layer = 2) +
  mutual(same = "layer.mem", diff = TRUE) +
  gwidegree(decay = 1.5, fixed = TRUE, attr = "layer.mem") +
  gwodegree(decay = 1.5, fixed = TRUE, attr = "layer.mem") +
  gwesp_layer(decay = 0.69, fixed = TRUE, layer = 1) +
  gwesp_layer(decay = 0.69, fixed = TRUE, layer = 2) +
  edgecov_layer(log_cdist, layer = 1) +
  edgecov_layer(log_cdist, layer = 2) +
  edgecov_layer(log_gdp_out, layer = 1) +
  edgecov_layer(log_gdp_out, layer = 2) +
  edgecov_layer(log_gdp_in, layer = 1) +
  edgecov_layer(log_gdp_in, layer = 2) +
  edgecov_layer(alliance, layer = 1) +
  edgecov_layer(alliance, layer = 2) +
  edgecov_layer(polity_diff, layer = 1) +
  edgecov_layer(polity_diff, layer = 2) +
  edgecov_layer(log_milit_exp_out, layer = 1) +
  edgecov_layer(log_milit_exp_out, layer = 2) +
  edgecov_layer(log_milit_exp_in, layer = 1) +
  edgecov_layer(log_milit_exp_in, layer = 2) +
  edgecov_layer(pathdep_layer1, layer = 1) +
  edgecov_layer(pathdep_layer2D, layer = 2) 



#------------------------------------------------------------------------------#
# Model Estimation
#------------------------------------------------------------------------------#

# using the values estimated in the stochastic approximation estimation before
init_val <- readRDS(file = paste0(path, "/models/ERGM/ergm_results_2003_init_val.rds"))

fit_C <- ergm(
  modelC, 
  constraints = ~fixallbut(free), 
  estimate = "MLE", 
  eval.loglik = TRUE, 
  check.degeneracy = TRUE, 
  verbose = TRUE,
  control = control.ergm(seed = 1234, 
                         init = init_val$fit_C,
                         parallel = 6, 
                         main.method = "MCMLE", 
                         MCMLE.maxit = 30)
)

fit_C_nn <- ergm(
  modelC_nn, 
  constraints = ~fixallbut(free), 
  estimate = "MLE", 
  eval.loglik = TRUE, 
  check.degeneracy = TRUE, 
  verbose = TRUE,
  control = control.ergm(seed = 1234, 
                         init = init_val$fit_C_nn,
                         parallel = 6, 
                         main.method = "MCMLE", 
                         MCMLE.maxit = 30)
)

fit_D <- ergm(
  modelD, 
  constraints = ~fixallbut(free), 
  estimate = "MLE", 
  eval.loglik = TRUE, 
  check.degeneracy = TRUE, 
  verbose = TRUE,
  control = control.ergm(seed = 1234, 
                         init = init_val$fit_D,
                         parallel = 6, 
                         main.method = "MCMLE", 
                         MCMLE.maxit = 30)
)

fit_D_nn <- ergm(
  modelD_nn, 
  constraints = ~fixallbut(free), 
  estimate = "MLE", 
  eval.loglik = TRUE, 
  check.degeneracy = TRUE, 
  verbose = TRUE,
  control = control.ergm(seed = 1234, 
                         init = init_val$fit_D_nn,
                         parallel = 6, 
                         main.method = "MCMLE", 
                         MCMLE.maxit = 30)
)

cat("Estimation of models finished. \n")



#------------------------------------------------------------------------------#
# Simulations for Goodness of Fit Assessment
#------------------------------------------------------------------------------#

# helper function to simulate, split and compute network stats - layer specific
sim_split_stat <- function(fit, n_sim = 1000){
  
  list_of_nets <- simulate(fit,
                           nsim = n_sim,
                           seed = 42L,
                           output = "network",
                           control = control.simulate.ergm(parallel = 6))
  
  list_of_mat <- lapply(list_of_nets, as.matrix.network)
  
  list_dt <- list()
  
  for (i in seq_along(list_of_mat)){
    
    mat <- list_of_mat[[i]]
    
    arms <- mat[1:114, 1:114]
    arms <- as.network.matrix(arms)
    
    trade <- mat[115:228, 115:228]
    trade <- as.network.matrix(trade)
    
    geodist_arms <- ergm.geodistdist(arms)
    names(geodist_arms) <- paste0("geodist.", names(geodist_arms))
   
    geodist_trade <- ergm.geodistdist(trade)
    names(geodist_trade) <- paste0("geodist.", names(geodist_trade))
    
    summary_arms <- summary(arms ~ edges + mutual + idegree(0:50) + odegree(0:50) + esp(0:50) + + dsp(0:50))
    summary_trade <- summary(trade ~ edges + mutual + idegree(0:50) + odegree(0:50) + esp(0:50) + + dsp(0:50))
    
    vec_arms <- c("network" = 1, geodist_arms, summary_arms)
    vec_trade <- c("network" = 2, geodist_trade, summary_trade)
    
    list_dt[[i]] <- rbind(c("network" = 1, geodist_arms, summary_arms), 
                            c("network" = 2, geodist_trade, summary_trade))
  }
  
  dt_cat <- do.call(rbind, list_dt)
  
  return(dt_cat)
}


sim_C <- sim_split_stat(fit_C, n_sim = 1000)
sim_C_nn <- sim_split_stat(fit_C_nn, n_sim = 1000)
sim_D <- sim_split_stat(fit_D, n_sim = 1000)
sim_D_nn <- sim_split_stat(fit_D_nn, n_sim = 1000)


# simulate cross layer stats
sim_cross_C <- simulate(fit_C, 
                        nsim = 1000,
                        seed = 42L,
                        output = "stats",
                        monitor = NULL,
                        control = control.simulate.ergm(parallel = 6))

sim_cross_C_nn <- simulate(fit_C_nn, 
                            nsim = 1000,
                            seed = 42L,
                            output = "stats",
                            monitor = ~ duplexdyad(c("e", "f", "h"), layers = list(1, 2)),
                            control = control.simulate.ergm(parallel = 6))

sim_cross_D <- simulate(fit_D, 
                         nsim = 1000,
                         seed = 42L,
                         output = "stats",
                         monitor = NULL,
                         control = control.simulate.ergm(parallel = 6))

sim_cross_D_nn <- simulate(fit_D_nn, 
                            nsim = 1000,
                            seed = 42L,
                            output = "stats",
                            monitor = ~ duplexdyad(c("e", "f", "h"), layers = list(1, 2)),
                            control = control.simulate.ergm(parallel = 6))



#------------------------------------------------------------------------------#
# Save results
#------------------------------------------------------------------------------#

# save
save.image(file = paste0(path, "/models/ERGM/ergm_results_2003_mcmle.RData"))


