#------------------------------------------------------------------------------#
# Replication code for the ERGM Analysis
# Starting Values for MCMLE: Stochastic Approximation
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


#------------------------------------------------------------------------------#
# Construct Covariates
#------------------------------------------------------------------------------#

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

# chose parallel = 1 only because function seems broken and does not aggregate 
# with multiple chains

control_config_SA <- control.ergm(seed = 42L, 
                                  main.method = "Stochastic-Approximation", 
                                  SA.phase1_n = NULL,
                                  SA.initial_gain = NULL,
                                  SA.nsubphases = 4,
                                  SA.niterations = NULL,
                                  SA.phase3_n = 2500, # default is 1000 
                                  SA.interval = 1024 * 96,
                                  SA.burnin = 1024 * 12,
                                  SA.samplesize = 1024 * 4)

control_config_RM <- control.ergm(seed = 42L, 
                                  parallel = FALSE, 
                                  main.method = "Robbins-Monro")

control_config_MCMLE <- control.ergm(seed = 1234, 
                                     parallel = 8, 
                                     main.method = "MCMLE", 
                                     MCMLE.maxit = 30)


fit_C <- ergm(
  modelC, 
  constraints = ~fixallbut(free), 
  estimate = "MLE", 
  eval.loglik = TRUE, 
  check.degeneracy = TRUE, 
  verbose = TRUE,
  control = control_config_SA
)

fit_C_nn <- ergm(
  modelC_nn, 
  constraints = ~fixallbut(free), 
  estimate = "MLE", 
  eval.loglik = TRUE, 
  check.degeneracy = TRUE, 
  verbose = TRUE,
  control = control_config_SA
)

fit_D <- ergm(
  modelD, 
  constraints = ~fixallbut(free), 
  estimate = "MLE", 
  eval.loglik = TRUE, 
  check.degeneracy = TRUE, 
  verbose = TRUE,
  control = control_config_SA
)

fit_D_nn <- ergm(
  modelD_nn, 
  constraints = ~fixallbut(free), 
  estimate = "MLE", 
  eval.loglik = TRUE, 
  check.degeneracy = TRUE, 
  verbose = TRUE,
  control = control_config_SA
)

cat("Estimation of models finished. \n")



#------------------------------------------------------------------------------#
# Save results
#------------------------------------------------------------------------------#

# save
save.image(file = paste0(path, "/models/ERGM/ergm_results_2003_SA.RData"))

# load
#load(file = paste0(path, "/models/ERGM/ergm_results_2003_SA.RData"))



#------------------------------------------------------------------------------#
# Save Coefficients for MCMLE Estimation
#------------------------------------------------------------------------------#

init_val <- list()

init_val$fit_C <- coef(fit_C)
init_val$fit_C_nn <- coef(fit_C_nn)
init_val$fit_D <- coef(fit_D)
init_val$fit_D_nn <- coef(fit_D_nn)

saveRDS(init_val, paste0(path, "/models/ERGM/ergm_results_2003_init_val.rds"))



#------------------------------------------------------------------------------#
# MCMC Assessment
#------------------------------------------------------------------------------#

# output mcmc statistics
mcmc.diagnostics(fit_C, vars.per.page = 8, which = c("plots"))
mcmc.diagnostics(fit_C_nn, vars.per.page = 8, which = c("plots"))
mcmc.diagnostics(fit_D, vars.per.page = 8, which = c("plots"))
mcmc.diagnostics(fit_D_nn, vars.per.page = 8, which = c("plots"))


# Gweke Diagnostic?
mcmc.diagnostics(fit_C, which = c("burnin"))
mcmc.diagnostics(fit_C_nn, which = c("burnin"))
mcmc.diagnostics(fit_D, which = c("burnin"))
mcmc.diagnostics(fit_D_nn, which = c("burnin"))


# Effective Sample Size
neff_C <- data.frame("neff" = coda::effectiveSize(as.mcmc.list(fit_C$sample)), "Name" = names(coda::effectiveSize(as.mcmc.list(fit_C$sample))))
neff_C_nn <- data.frame("neff" = coda::effectiveSize(as.mcmc.list(fit_C_nn$sample)), "Name" = names(coda::effectiveSize(as.mcmc.list(fit_C_nn$sample))))
neff_D <- data.frame("neff" = coda::effectiveSize(as.mcmc.list(fit_D$sample)), "Name" = names(coda::effectiveSize(as.mcmc.list(fit_D$sample))))
neff_D_nn <- data.frame("neff" = coda::effectiveSize(as.mcmc.list(fit_D_nn$sample)), "Name" = names(coda::effectiveSize(as.mcmc.list(fit_D_nn$sample))))


neff_C %>% 
  rename("Import Dep. with" = neff) %>%
  full_join(neff_C_nn, by = "Name") %>%
  rename("Import Dep. without" = neff) %>%
  full_join(neff_D, by = "Name") %>%
  rename("Export Dep. with" = neff) %>%
  full_join(neff_D_nn, by = "Name") %>%
  rename("Export Dep. without" = neff)%>%
  relocate(Name) %>%
  pivot_longer(2:5, names_to = "Model", values_to = "eff") %>%
  group_by(Model) %>%
  summarise(Mean = mean(eff, na.rm = TRUE), SD = sd(eff, na.rm = TRUE)) -> neff



#------------------------------------------------------------------------------#
# Plot Results
#------------------------------------------------------------------------------#

custom.coef.names <- c(
  "Edges",
  "Edges",
  "Reciprocity",
  "Reciprocity",
  "Gw Indegree (d = 1.5)",
  "Gw Indegree (d = 1.5)",
  "Gw Outdegree (d = 1.5)",
  "Gw Outdegree (d = 1.5)",
  "GWESP Outgoing Two-path (d = 0.69)",
  "GWESP Outgoing Two-path (d = 0.69)",
  "E",
  "F",
  "H",
  "Distance (log)",
  "Distance (log)",
  "GDP in (log)",
  "GDP in (log)",
  "GDP out (log)",
  "GDP out (log)",
  "Alliance",
  "Alliance",
  "Polity Diff. (abs)",
  "Polity Diff. (abs)",
  "Military Expenditure in (log)",
  "Military Expenditure in (log)",
  "Military Expenditure out (log)",
  "Military Expenditure out (log)",
  "Path Dependency",
  "Path Dependency", 
  "Path Dependency"
)

custom.coef.names2 <- c(
  "Edges",
  "Edges",
  "Reciprocity",
  "Reciprocity",
  "Gw Indegree (d = 1.5)",
  "Gw Indegree (d = 1.5)",
  "Gw Outdegree (d = 1.5)",
  "Gw Outdegree (d = 1.5)",
  "GWESP Outgoing Two-path (d = 0.69)",
  "GWESP Outgoing Two-path (d = 0.69)",
  "Distance (log)",
  "Distance (log)",
  "GDP in (log)",
  "GDP in (log)",
  "GDP out (log)",
  "GDP out (log)",
  "Alliance",
  "Alliance",
  "Polity Diff. (abs)",
  "Polity Diff. (abs)",
  "Military Expenditure in (log)",
  "Military Expenditure in (log)",
  "Military Expenditure out (log)",
  "Military Expenditure out (log)",
  "Path Dependency",
  "Path Dependency", 
  "Path Dependency"
)


# output 
screenreg(list(fit_C, fit_D),
       single.row = T,
       float.pos = "H",
       use.packages = FALSE,
       custom.coef.names = custom.coef.names, 
       custom.model.names = c("Import Dep.", "Export Dep."),
       booktabs = T, dcolumn = T, include.nobs = F,
       reorder.coef = c(1, 3, 5, 7, 9, 14, 16, 18, 20, 22, 24, 26, 28,
                        2, 4, 6, 8, 10, 15, 17, 19, 21, 23, 25, 27, 29, 
                        11, 12, 13), 
       groups = list("Layer 1: Arms Trade" = 1:13,
                     "Layer 2: Conventional Trade" = 14:26, 
                     "Cross Layer Network Effects" = 27:29),
       caption = "MERGM results for two-layer network of weapons and import (left) or export (right) trade dependency in the year 2003.", 
       label = "tab:ergm_estimates_model1C", 
       custom.note = "Estimates based on Stochastic Approximation. Standard Errors in parenthesis.\\newline%stars. ")

# output 
texreg(list(fit_C_nn, fit_D_nn),
       single.row = T, 
       use.packages = FALSE,
       float.pos = "H",
       custom.coef.names = custom.coef.names2, 
       custom.model.names = c("Import Dep.", "Export Dep."),
       booktabs = T, dcolumn = T, include.nobs = F,
       reorder.coef = c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25,
                        2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26),
       groups = list("Layer 1: Arms Trade" = 1:13,
                     "Layer 2: Conventional Trade" = 14:26),
       caption = "MERGM results for two-layer network of weapons and import (left) or export (right) trade dependency in the year 2003 - estimated without cross-layer effects.", 
       label = "tab:ergm_estimates_model1CnnDnn", 
       custom.note = "Estimates based on Stochastic Approximation. Standard Errors in parenthesis.\\newline%stars. ")
