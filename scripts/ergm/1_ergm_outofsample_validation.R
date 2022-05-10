#------------------------------------------------------------------------------#
# out of Sample Prediction
#------------------------------------------------------------------------------#

library(network)
library(ergm)
library(multilayer.ergm)
library(texreg)
library(stargazer)
library(coda)
library(PRROC)
library(dplyr)
library(xtable)

rm(list = ls(all.names = TRUE))


# setup
source("utils/utils.R")
source("utils/construct_header.R")
source("utils/custom_trade_to_binary.R")
path <- data_path.get()
set.seed(1234)


# load output from results file including model, fits, original networks
load(file = paste0(path, "/models/ERGM/ergm_results_2003_mcmle.RData"))


# define year for out-of-sample prediction
year_next <- 2004


## update coefficients 
# such that these are from time point t (because still lagged)

# gdppc out
log_gdppc_out <- matrix(log(gdppc[included, i1]), n, n, byrow = FALSE)

# gdppc in
log_gdppc_in <- matrix(log(gdppc[included, i1]), n, n, byrow = TRUE)

# gdp out
log_gdp_out <- matrix(log(gdp[included, i1]), n, n, byrow = FALSE)

# gdp in
log_gdp_in <- matrix(log(gdp[included, i1]), n, n, byrow = TRUE)

# atop alliance
alliance <- atop_alliance[[i1]][included, included]

# absolute polity diff
polity_diff <- abs(outer(polity[included, i1], polity[included, i1], "-"))

# military expenditure out
log_milit_exp_out <- matrix(log(milit_exp[included, i1]), n, n, byrow = FALSE)

# military expenditure in
log_milit_exp_in <- matrix(log(milit_exp[included, i1]), n, n, byrow = TRUE)

# define path dependency for each layer and once for import once export dependency
pathdep_layer1 <- 1 * (
  sipri_tiv[[i1]][included, included] +
    sipri_tiv[[i1 - 1]][included, included] +
    sipri_tiv[[i1 - 2]][included, included] > 0)

pathdep_layer2C <- 1 * (
  custom_trade_to_binary(trade[[i2]][included, included], type = "C", threshold = 0.01) +
    custom_trade_to_binary(trade[[i2 - 1]][included, included], type = "C", threshold = 0.01) +
    custom_trade_to_binary(trade[[i2 - 2]][included, included], type = "C", threshold = 0.01))

pathdep_layer2D <- 1 * (
  custom_trade_to_binary(trade[[i2]][included, included], type = "D", threshold = 0.01) +
    custom_trade_to_binary(trade[[i2 - 1]][included, included], type = "D", threshold = 0.01) +
    custom_trade_to_binary(trade[[i2 - 2]][included, included], type = "D", threshold = 0.01))

# dependent multi-layer network
netC_future <- to.multiplex(
  1 * (sipri_tiv[[i1 + 1]][included, included] > 0),
  custom_trade_to_binary(trade[[i2 + 1]][included, included], type = "C", threshold = 0.01),
  output = "network"
)

netD_future <- to.multiplex(
  1 * (sipri_tiv[[i1 + 1]][included, included] > 0),
  custom_trade_to_binary(trade[[i2 + 1]][included, included], type = "D", threshold = 0.01),
  output = "network"
)

# matrix to select paramters which are under consideration
free_mat <- as.matrix.network(free)

# true observed values as vector
obs_C <- as.matrix.network(netC_future)[which(free_mat == 1, arr.ind = T)]
obs_D <- as.matrix.network(netD_future)[which(free_mat == 1, arr.ind = T)]




# config
control_simulate_config <- control.simulate(MCMC.interval = 1024 * 4, parallel = 6)


# For each of the four we simulate 1000 models with the old formula (but
# updated covariates, fitted coefs)


## C
pred_Sim_C <- simulate(modelC, 
                       coef = coef(fit_C), 
                       nsim = 1000,
                       seed = 42L,
                       output = "network",
                       constraints = ~fixallbut(free), 
                       control = control_simulate_config)

pred_Sim_C <- lapply(pred_Sim_C, as.matrix.network)
pred_Sim_C <- Reduce("+", pred_Sim_C) / length(pred_Sim_C)
pred_Sim_C <- pred_Sim_C[which(free_mat == 1, arr.ind = T)]

fg <- pred_Sim_C[obs_C == 1]
bg <- pred_Sim_C[obs_C == 0]

pr_auc_C <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)$auc.integral
bs_C <- mean((pred_Sim_C - obs_C)**2)


## C_nn
pred_Sim_C_nn <- simulate(modelC_nn, 
                       coef = coef(fit_C_nn), 
                       nsim = 1000,
                       seed = 42L,
                       output = "network",
                       constraints = ~fixallbut(free), 
                       control = control_simulate_config)

pred_Sim_C_nn <- lapply(pred_Sim_C_nn, as.matrix.network)
pred_Sim_C_nn <- Reduce("+", pred_Sim_C_nn) / length(pred_Sim_C_nn)
pred_Sim_C_nn <- pred_Sim_C_nn[which(free_mat == 1, arr.ind = T)]

fg <- pred_Sim_C_nn[obs_C == 1]
bg <- pred_Sim_C_nn[obs_C == 0]

pr_auc_C_nn <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)$auc.integral
bs_C_nn <- mean((pred_Sim_C_nn - obs_C)**2)


## D

pred_Sim_D <- simulate(modelD, 
                       coef = coef(fit_D), 
                       nsim = 1000,
                       seed = 42L,
                       output = "network",
                       constraints = ~fixallbut(free), 
                       control = control_simulate_config)

pred_Sim_D <- lapply(pred_Sim_D, as.matrix.network)
pred_Sim_D <- Reduce("+", pred_Sim_D) / length(pred_Sim_D)
pred_Sim_D <- pred_Sim_D[which(free_mat == 1, arr.ind = T)]

fg <- pred_Sim_D[obs_D == 1]
bg <- pred_Sim_D[obs_D == 0]

pr_auc_D <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)$auc.integral
bs_D <- mean((pred_Sim_D - obs_D)**2)


## D_nn
pred_Sim_D_nn <- simulate(modelD_nn, 
                           coef = coef(fit_D_nn), 
                           nsim = 1000,
                           seed = 42L,
                           output = "network",
                           constraints = ~fixallbut(free), 
                           control = control_simulate_config)

pred_Sim_D_nn <- lapply(pred_Sim_D_nn, as.matrix.network)
pred_Sim_D_nn <- Reduce("+", pred_Sim_D_nn) / length(pred_Sim_D_nn)
pred_Sim_D_nn <- pred_Sim_D_nn[which(free_mat == 1, arr.ind = T)]

fg <- pred_Sim_D_nn[obs_D == 1]
bg <- pred_Sim_D_nn[obs_D == 0]

pr_auc_D_nn <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)$auc.integral
bs_D_nn <- mean((pred_Sim_D_nn - obs_D)**2)



## Output

out <- data.frame(
  "Model" = c("Import Dependency", "Import Dependency", "Export Dependency", "Export Dependency"),
  "Metric" = c("PR AUC", "Brier Score", "PR AUC", "Brier Score"),
  "without" = c(pr_auc_C_nn, bs_C_nn, pr_auc_D_nn, bs_D_nn),
  "with" = c(pr_auc_C, bs_C, pr_auc_D, bs_D)
  )


sink("figures/ergm_outofsample_scores_2003.txt")

xtable(out, booktabs = T,  format = "latex", 
    digits = 4,  escape = F, linesep = "",
    caption = "Out-of-sample validation for the ERGM \\newline \\footnotesize Comparison of out-of-sample validation with and without cross-layer effects.", 
    label = "tab:ergm_outofsample_scores") 

sink()
