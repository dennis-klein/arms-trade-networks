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
source("utils/construct_header.R")
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

# using the values estimated in the Stochastic Approx before
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
# MCMC Diagnostics and Goodness of Fit
#------------------------------------------------------------------------------#

# helper function to simulate, split and compute network stats - network specific
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

# load
#load(file = paste0(path, "/models/ERGM/ergm_results_2003_mcmle.RData"))



#------------------------------------------------------------------------------#
# MCMC Assessment
#------------------------------------------------------------------------------#

# output mcmc statistics
pdf(file = "figures/ergm_mcmc_1C_2003.pdf", width = 10, height = 0.9*sqrt(2)*10)
par(mfrow = c(6,2))
plot(fit_C$sample, smooth = TRUE, auto.layout = FALSE)
dev.off()

pdf(file = "figures/ergm_mcmc_1Cnn_2003.pdf", width = 10, height = 0.9*sqrt(2)*10)
par(mfrow = c(6,2))
plot(fit_C_nn$sample, smooth = TRUE, auto.layout = FALSE)
dev.off()

pdf(file = "figures/ergm_mcmc_1D_2003.pdf", width = 10, height = 0.9*sqrt(2)*10)
par(mfrow = c(6,2))
plot(fit_D$sample, smooth = TRUE, auto.layout = FALSE)
dev.off()

pdf(file = "figures/ergm_mcmc_1Dnn_2003.pdf", width = 10, height = 0.9*sqrt(2)*10)
par(mfrow = c(6,2))
plot(fit_D_nn$sample, smooth = TRUE, auto.layout = FALSE)
dev.off()


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
# Goodness of Fit Assessment
#------------------------------------------------------------------------------#

# recover observed cross layer statistics
obs_stats_C <- as.data.frame(as.list(summary(netC ~ duplexdyad(c("e", "f",  "h"), layers = list(1, 2)))))
obs_stats_D <- as.data.frame(as.list(summary(netD ~ duplexdyad(c("e", "f",  "h"), layers = list(1, 2)))))


# plot raster of distributions with observed and simulated cross-layer statistics
pdf(file = "figures/ergm_gof_duplex_2003.pdf", width = 10, height = 0.3*sqrt(2)*10)
par(mfrow = c(2,3))

plot(density(sim_cross_C[, "duplexdyad.e"]) , main = "Import Dependency Model: E", lty = 3)
lines(density(sim_cross_C_nn[, "duplexdyad.e"]), lty = 2)
abline(v = obs_stats_C$duplexdyad.e, lty = 1)
legend("topright", legend = c("with", "without", "observed"), lty = c(3,2,1))

plot(density(sim_cross_C[, "duplexdyad.f"]) , main = "Import Dependency Model: F", lty = 3)
lines(density(sim_cross_C_nn[, "duplexdyad.f"]), lty = 2)
abline(v = obs_stats_C$duplexdyad.f, lty = 1)
legend("topright", legend = c("with", "without", "observed"), lty = c(3,2,1))

plot(density(sim_cross_C[, "duplexdyad.h"]) , main = "Import Dependency Model: H", lty = 3)
lines(density(sim_cross_C_nn[, "duplexdyad.h"]), lty = 2)
abline(v = obs_stats_C$duplexdyad.h, lty = 1)
legend("topright", legend = c("with", "without", "observed"), lty = c(3,2,1))


plot(density(sim_cross_D[, "duplexdyad.e"]) , main = "Export Dependency Model: E", lty = 3)
lines(density(sim_cross_D_nn[, "duplexdyad.e"]), lty = 2)
abline(v = obs_stats_D$duplexdyad.e, lty = 1)
legend("topright", legend = c("with", "without", "observed"), lty = c(3,2,1))

plot(density(sim_cross_D[, "duplexdyad.f"]) , main = "Export Dependency Model: F", lty = 3)
lines(density(sim_cross_D_nn[, "duplexdyad.f"]), lty = 2)
abline(v = obs_stats_D$duplexdyad.f, lty = 1)
legend("topright", legend = c("with", "without", "observed"), lty = c(3,2,1))

plot(density(sim_cross_D[, "duplexdyad.h"]) , main = "Export Dependency Model: H", lty = 3)
lines(density(sim_cross_D_nn[, "duplexdyad.h"]), lty = 2)
abline(v = obs_stats_D$duplexdyad.h, lty = 1)
legend("topright", legend = c("with", "without", "observed"), lty = c(3,2,1))

dev.off()


#------------------------------------------------------------------------------#
# Plot comparison of degree esp desp  
#------------------------------------------------------------------------------#

# C vs C_nn comparison
net_arms <- as.network(as.matrix.network(netC)[1:114, 1:114])
net_trade_C <- as.network(as.matrix.network(netC)[115:228, 115:228])
net_trade_D <- as.network(as.matrix.network(netD)[115:228, 115:228])

obs_stats_arms <- as.data.frame(as.list(summary(net_arms ~ edges + mutual + idegree(0:50) + odegree(0:50) + esp(0:50) + + dsp(0:50))))
obs_stats_C <- as.data.frame(as.list(summary(net_trade_C ~ edges + mutual + idegree(0:50) + odegree(0:50) + esp(0:50) + + dsp(0:50))))
obs_stats_D <- as.data.frame(as.list(summary(net_trade_D ~ edges + mutual + idegree(0:50) + odegree(0:50) + esp(0:50) + + dsp(0:50))))

sim_cross_C <- as.data.frame(sim_cross_C)
sim_cross_C_nn <- as.data.frame(sim_cross_C_nn)
sim_cross_D <- as.data.frame(sim_cross_D)
sim_cross_D_nn <- as.data.frame(sim_cross_D_nn)

sim_C <- as.data.frame(sim_C)
sim_C_nn <- as.data.frame(sim_C_nn)
sim_D <- as.data.frame(sim_D)
sim_D_nn <- as.data.frame(sim_D_nn)

geo_dist_arms <- ergm.geodistdist(net_arms)
names(geo_dist_arms) <- paste0("geodist.", names(geo_dist_arms))
geo_dist_C <- ergm.geodistdist(net_trade_C)
names(geo_dist_C) <- paste0("geodist.", names(geo_dist_C))
geo_dist_D <- ergm.geodistdist(net_trade_D)
names(geo_dist_D) <- paste0("geodist.", names(geo_dist_D))

range(sim_cross_C$edges_layer.1)
range(sim_cross_C$edges_layer.1)
range(sim_cross_C$edges_layer.1)
range(sim_cross_C$edges_layer.1)


## prepare data

# edges + reciprocity
bind_rows(list(sim_cross_C, sim_cross_C_nn)) %>%
  mutate(model = as.factor(c(rep("Import - with", 1000), rep("Import - without", 1000)))) -> pl_edges

# geodist distribution
bind_rows(list(sim_C, sim_C_nn)) %>%
  mutate(model = as.factor(c(rep("Import - with", 2000), rep("Import - without", 2000)))) -> pl_sim




pdf(file = "figures/ergm_gof_import_models_pt1_2003.pdf",  width = 10, height = 0.9*sqrt(2)*10)
layout(matrix(c(1,2,3,4,
                5,5,5,6,
                7,7,7,8, 
                9,9,9,10,
                11,11,11,12), nrow = 5, ncol = 4, byrow = TRUE))

# Plot Edges
par(mar = c(5.1, 7.1, 4.1, 2.1))
boxplot(edges_layer.1 ~ model,  data = pl_edges, xlab = "Edges - MCW", ylab = "", horizontal=TRUE, las = 1, main = 1)
points(rep(obs_stats_arms$edges, 2), 1:2, pch = 19, col = "red", lwd = 4)

par(mar = c(5.1, 7.1, 4.1, 2.1))
boxplot(edges_layer.2 ~ model,  
        data = pl_edges, 
        xlab = "Edges - Conventional Trade", ylab = "", 
        horizontal=TRUE, 
        las = 1, main = 2)
points(c(obs_stats_C$edges, obs_stats_C$edges), 1:2, pch = 19, col = "red", lwd = 4)


# Plot Reciprocity
par(mar = c(5.1, 7.1, 4.1, 2.1))
boxplot(mutual.same.layer.mem.1 ~ model,  data = pl_edges, xlab = "Reciprocity - MCW", ylab = "", horizontal=TRUE, las = 1, main = 3)
points(rep(obs_stats_arms$mutual, 4), 4:1, pch = 19, col = "red", lwd = 4)

par(mar = c(5.1, 7.1, 4.1, 2.1))
boxplot(mutual.same.layer.mem.2 ~ model,  
        data = pl_edges, 
        xlab = "Reciprocity - Conventional Trade", ylab = "", 
        horizontal=TRUE, 
        las = 1, main = 4)
points(c(obs_stats_C$mutual, obs_stats_C$mutual), 1:2, pch = 19, col = "red", lwd = 4)


# Plot Geodist for ARMS
tmp <- filter(pl_sim, model == "Import - with" & network == 1)
take <- c("geodist.1", "geodist.2", "geodist.3", "geodist.4", "geodist.5", "geodist.6")
boxplot(tmp[, take], 
        xlab = "Geodesic Distance - MCW - with cross-layer effects", ylab = "",
        las = 1, main = 5)
points(1:length(take), geo_dist_arms[take], pch = 19, col = "red", lwd = 4)

boxplot(tmp[, c("geodist.Inf")], 
        xlab = "", ylab = "",
        las = 1, main = 6)
axis(side = 1, at = 1, label = "geodist.Inf")
points(1, geo_dist_arms[c("geodist.Inf")], pch = 19, col = "red", lwd = 4)


tmp <- filter(pl_sim, model == "Import - without" & network == 1)
boxplot(tmp[, take], 
        xlab = "Geodesic Distance - MCW - without cross-layer effects", ylab = "",
        las = 1, main = 7)
points(1:length(take), geo_dist_arms[take], pch = 19, col = "red", lwd = 4)

boxplot(tmp[, c("geodist.Inf")], 
        xlab = "", ylab = "",
        las = 1, main = 8)
axis(side = 1, at = 1, label = "geodist.Inf")
points(1, geo_dist_arms[c("geodist.Inf")], pch = 19, col = "red", lwd = 4)


# Plot Geodist for TRADE
tmp <- filter(pl_sim, model == "Import - with" & network == 2)
boxplot(tmp[, take], 
        xlab = "Geodesic Distance - Conventional Trade - with cross-layer effects", ylab = "",
        las = 1, main = 9)
points(1:length(take), geo_dist_C[take], pch = 19, col = "red", lwd = 4)

boxplot(tmp[, c("geodist.Inf")], 
        xlab = "", ylab = "",
        las = 1, main = 10)
axis(side = 1, at = 1, label = "geodist.Inf")
points(1, geo_dist_C[c("geodist.Inf")], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 2)
boxplot(tmp[, take], 
        xlab = "Geodesic Distance - Conventional Trade - without cross-layer effects", ylab = "",
        las = 1, main = 11)
points(1:length(take), geo_dist_C[take], pch = 19, col = "red", lwd = 4)

boxplot(tmp[, c("geodist.Inf")], 
        xlab = "", ylab = "",
        las = 1, main = 12)
axis(side = 1, at = 1, label = "geodist.Inf")
points(1, geo_dist_C[c("geodist.Inf")], pch = 19, col = "red", lwd = 4)


dev.off()





## Plot Edgewise Shared Partners

pdf(file = "figures/ergm_gof_import_models_pt2_2003.pdf",  width = 10, height = 0.9*sqrt(2)*10)
layout(matrix(c(1,1,2,2,
                3,3,4,4,
                5,5,6,6,
                7,7,8,8), nrow = 4, ncol = 4, byrow = TRUE))
# Indegree ARMS
tmp <- filter(pl_sim, model == "Import - with" & network == 1)
take <- c("idegree0", "idegree1", "idegree2", "idegree3", "idegree4", "idegree5", "idegree6", "idegree7", "idegree8", "idegree9", "idegree10")
boxplot(tmp[, take], 
        xlab = "Indegree - MCW - with cross-layer effects", 
        ylab = "", las = 1, main = 1)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 1)
boxplot(tmp[, take], 
        xlab = "Indegree - MCW - without cross-layer effects", 
        ylab = "", las = 1, main = 2)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)


# Indegree TRADE
tmp <- filter(pl_sim, model == "Import - with" & network == 2)
boxplot(tmp[, take], 
        xlab = "Indegree - Conventional Trade - with cross-layer effects", 
        ylab = "", las = 1, main = 3)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 2)
boxplot(tmp[, take], 
        xlab = "Indegree - Conventional Trade - without cross-layer effects", 
        ylab = "", las = 1, main = 4)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)



# Outdegree ARMS
tmp <- filter(pl_sim, model == "Import - with" & network == 1)
take <- c("odegree0", "odegree1", "odegree2", "odegree3", "odegree4", "odegree5", "odegree6", "odegree7", "odegree8", "odegree9", "odegree10")
boxplot(tmp[, take], 
        xlab = "Outdegree - MCW - with cross-layer effects", 
        ylab = "", las = 1, main = 5)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 1)
boxplot(tmp[, take], 
        xlab = "Outdegree - MCW - without cross-layer effects", 
        ylab = "", las = 1, main = 6)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)


# Outdegree TRADE
tmp <- filter(pl_sim, model == "Import - with" & network == 2)
boxplot(tmp[, take], 
        xlab = "Outdegree - Conventional Trade - with cross-layer effects", 
        ylab = "", las = 1, main = 7)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 2)
boxplot(tmp[, take], 
        xlab = "Outdegree - Conventional Trade - without cross-layer effects", 
        ylab = "", las = 1, main = 8)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)

dev.off()



pdf(file = "figures/ergm_gof_import_models_sl1_2003.pdf", width = 10, height = 6)
layout(matrix(c(1,1,2,2,
                3,3,4,4
                ), nrow = 2, ncol = 4, byrow = TRUE))

# Outdegree ARMS
tmp <- filter(pl_sim, model == "Import - with" & network == 1)
take <- c("odegree0", "odegree1", "odegree2", "odegree3", "odegree4", "odegree5", "odegree6", "odegree7", "odegree8", "odegree9", "odegree10")
boxplot(tmp[, take], 
        xlab = "Outdegree - MCW - with cross-layer effects", 
        ylab = "", las = 1, main = 1)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 1)
boxplot(tmp[, take], 
        xlab = "Outdegree - MCW - without cross-layer effects", 
        ylab = "", las = 1, main = 2)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)


# Outdegree TRADE
tmp <- filter(pl_sim, model == "Import - with" & network == 2)
boxplot(tmp[, take], 
        xlab = "Outdegree - Conventional Trade - with cross-layer effects", 
        ylab = "", las = 1, main = 3)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 2)
boxplot(tmp[, take], 
        xlab = "Outdegree - Conventional Trade - without cross-layer effects", 
        ylab = "", las = 1, main = 4)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)

dev.off()



# Plot ESP

pdf(file = "figures/ergm_gof_import_models_pt3_2003.pdf",  width = 10, height = 0.9*sqrt(2)*10)
layout(matrix(c(1,1,2,2,
                3,3,4,4,
                5,5,6,6,
                7,7,8,8), nrow = 4, ncol = 4, byrow = TRUE))

# ESP ARMS
tmp <- filter(pl_sim, model == "Import - with" & network == 1)
take <- c("esp0", "esp1", "esp2", "esp3" ,"esp4" ,"esp5", "esp6" ,"esp7", "esp8", "esp9" ,"esp10" ,"esp11" ,"esp12" ,"esp13", "esp14", "esp15", "esp16", "esp17" ,"esp18" ,"esp19", "esp20")
boxplot(tmp[, take], 
        xlab = "ESP - MCW - with cross-layer effects", 
        ylab = "", las = 1, main = 1)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 1)
boxplot(tmp[, take], 
        xlab = "ESP - MCW - without cross-layer effects", 
        ylab = "", las = 1, main = 2)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)


# ESP TRADE
tmp <- filter(pl_sim, model == "Import - with" & network == 2)
boxplot(tmp[, take], 
        xlab = "ESP - Conventional Trade - with cross-layer effects", 
        ylab = "", las = 1, main = 3)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 2)
boxplot(tmp[, take], 
        xlab = "ESP - Conventional Trade - without cross-layer effects", 
        ylab = "", las = 1, main = 4)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)


# DSP ARMS
tmp <- filter(pl_sim, model == "Import - with" & network == 1)
take <- c("dsp0", "dsp1", "dsp2", "dsp3" ,"dsp4" ,"dsp5", "dsp6" ,"dsp7", "dsp8", "dsp9" ,"dsp10" ,"dsp11" ,"dsp12" ,"dsp13", "dsp14", "dsp15", "dsp16", "dsp17" ,"dsp18" ,"dsp19", "dsp20")
boxplot(tmp[, take], 
        xlab = "DSP - MCW - with cross-layer effects", 
        ylab = "", las = 1, main = 5)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 1)
boxplot(tmp[, take], 
        xlab = "DSP - MCW - without cross-layer effects", 
        ylab = "", las = 1, main = 6)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)


# DSP TRADE
tmp <- filter(pl_sim, model == "Import - with" & network == 2)
boxplot(tmp[, take], 
        xlab = "DSP - Conventional Trade - with cross-layer effects", 
        ylab = "", las = 1, main = 7)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 2)
boxplot(tmp[, take], 
        xlab = "DSP - Conventional Trade - without cross-layer effects", 
        ylab = "", las = 1, main = 8)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)

dev.off()



pdf(file = "figures/ergm_gof_import_models_sl2_2003.pdf",  width = 10, height = 6)
layout(matrix(c(1,1,2,2,
                3,3,4,4), nrow = 2, ncol = 4, byrow = TRUE))

# ESP ARMS
tmp <- filter(pl_sim, model == "Import - with" & network == 1)
take <- c("esp0", "esp1", "esp2", "esp3" ,"esp4" ,"esp5", "esp6" ,"esp7", "esp8", "esp9" ,"esp10" ,"esp11" ,"esp12" ,"esp13", "esp14", "esp15", "esp16", "esp17" ,"esp18" ,"esp19", "esp20")
boxplot(tmp[, take], 
        xlab = "ESP - MCW - with cross-layer effects", 
        ylab = "", las = 1, main = 1)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 1)
boxplot(tmp[, take], 
        xlab = "ESP - MCW - without cross-layer effects", 
        ylab = "", las = 1, main = 2)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)


# ESP TRADE
tmp <- filter(pl_sim, model == "Import - with" & network == 2)
boxplot(tmp[, take], 
        xlab = "ESP - Conventional Trade - with cross-layer effects", 
        ylab = "", las = 1, main = 3)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 2)
boxplot(tmp[, take], 
        xlab = "ESP - Conventional Trade - without cross-layer effects", 
        ylab = "", las = 1, main = 4)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)

dev.off()


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


sink("figures/ergm_effectivesize_2003.txt")
kbl(neff, booktabs = T,  format = "latex", 
    digits = 0,  escape = F, linesep = "",
    caption = "Effective Sample Size of the MCMC Sample Statistics", label = "ergm_effectivesize_2003") %>%
  kable_styling(font_size = 11, full_width = FALSE) %>%
  kable_styling(latex_options = "HOLD_position")
sink()


# output 
sink(file = "figures/ergm_estimates_1CD_2003.txt")
texreg(list(fit_C, fit_D),
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
       custom.note = "Estimates based on Monte Carlo MLE. Standard Errors in parenthesis.\\newline%stars. ")
sink()


# output 
sink(file = "figures/ergm_estimates_1CnnDnn_2003.txt")
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
       custom.note = "Estimates based on Monte Carlo MLE. Standard Errors in parenthesis.\\newline%stars. ")
sink()






