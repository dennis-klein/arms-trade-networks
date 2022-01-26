library(network)
library(sna)
library(data.table)
library(ergm)
library(multilayer.ergm)
library(xtable)


# help 
?'multilayer_terms'
?'ergm-terms'


# setup
rm(list = ls(all.names = TRUE))
set.seed(1234)
source("utils/utils.R")
source("utils/construct_header.R")
source("utils/custom_trade_to_binary.R")
path = data_path.get()


# load data 
load(file.path(path, "out/EX.RData"))
load(file.path(path, "out/colony.RData"))
load(file.path(path, "out/country_list.RData"))
load(file.path(path, "out/cdist.RData"))

nmc_cinc = readRDS(file.path(path, "out/nmc_cinc.rds"))
polity = readRDS(file.path(path, "out/polity.rds"))
gdp = readRDS(file.path(path, "out/gdp.rds"))
gdppc = readRDS(file.path(path, "out/gdppc.rds"))
sipri_tiv = readRDS(file.path(path, "out/sipri_tiv.rds"))
trade = readRDS(file.path(path, "out/baci_aggregated.rds"))
load(file.path(path, "out/atop_alliance.RData"))


# set parameters of analyses
# trade data starting 1995
start = 1995
end = 2018
year = 2003
n_sim = 100


# set correct indices corresponding data source
i1 = year - 1949 # Cornelius Replication Data 1950:2018
i2 = year - 1994 # CEPII BACI


# selection of countries included
included = rowSums(EX[, (start:end)-1949]) == length(start:end)
n = sum(included) 


# construct dependent multi-layer net
# layer 1: arms, layer 2: trade + apply thresholds
tmp1 = (sipri_tiv[[i1]][included, included] > 0) * 1
tmp2 = custom_trade_to_binary(trade[[i2]][included, included], type = "D", threshold = 0.01)
net = to.multiplex(tmp1, tmp2, output = "network")
check.multilayer(net)


# define sample space constraint
free = matrix(c(rep(c(rep(1, n), rep(0, n)), n), rep(c(rep(0, n), rep(1, n)), n)), nrow = 2*n, byrow = T)
diag(free) = 0
free = network(free, directed = T)


# nodal and dyadic attributes, lagged by t = 1
log_cdist = log(cdist[included, included]+1)
log_gdp_ocov = matrix(log(gdp[included, i1-1]), n, n, byrow = FALSE)
log_gdppc_icov = matrix(log(gdppc[included, i1-1]), n, n, byrow = TRUE)
absdiff_polity = abs(outer(polity[included, i1-1], polity[included, i1-1],'-'))
cinc100_ocov = matrix(100 * nmc_cinc[included, i1-1], n, n, byrow = FALSE)
alliance = atop_alliance[[i1-1]][included, included]

pathdep_arms = (sipri_tiv[[i1-1]][included, included] > 0) * 1
pathdep_trade = custom_trade_to_binary(trade[[i2-1]][included, included], type = "D", threshold = 0.01)


# measure time
start_time = Sys.time()


# layer independence model
if(FALSE){
fit1 <- ergm(net ~ mutual(same = "layer.mem", diff = TRUE) +
               edges_layer(layer = 1) +
               edges_layer(layer = 2) +
               gwesp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
               gwesp_layer(decay = 1.5, fixed = TRUE, layer = 2) +
               gwidegree(decay = 1, fixed = TRUE, attr = "layer.mem") +
               gwodegree(decay = 1, fixed = TRUE, attr = "layer.mem") +
               edgecov_layer(log_cdist, layer = 1) +
               edgecov_layer(log_gdppc_icov, layer = 1) +
               edgecov_layer(log_gdp_ocov, layer = 1) +
               edgecov_layer(log_cdist, layer = 2) +
               edgecov_layer(log_gdppc_icov, layer = 2) +
               edgecov_layer(log_gdp_ocov, layer = 2),
             eval.loglik = TRUE, 
             check.degeneracy = TRUE,
             verbose = TRUE, 
             estimate = c("MLE"),
             control = control.ergm(parallel = 4, seed = 0000),
             constraints = ~ fixallbut(free)
)}


# layer dependence model
fit2 <- ergm(net ~ mutual(same = "layer.mem", diff = TRUE) +
               edges_layer(layer = 1) +
               edges_layer(layer = 2) +
               gwesp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
               gwesp_layer(decay = 2, fixed = TRUE, layer = 2) +
               gwidegree(decay = 1, fixed = TRUE, attr = "layer.mem") +
               gwodegree(decay = 1, fixed = TRUE, attr = "layer.mem") +
               edgecov_layer(log_cdist, layer = 1) +
               edgecov_layer(log_gdppc_icov, layer = 1) +
               edgecov_layer(log_gdp_ocov, layer = 1) +
               edgecov_layer(log_cdist, layer = 2) +
               edgecov_layer(log_gdppc_icov, layer = 2) +
               edgecov_layer(log_gdp_ocov, layer = 2) +
               duplexdyad(c("e", "f", "h"), layers = list(1, 2)),
             eval.loglik = TRUE, 
             check.degeneracy = TRUE,
             verbose = TRUE, 
             estimate = c("MLE"),
             control = control.ergm(parallel = 4, seed = 0000),
             constraints = ~ fixallbut(free)
)


# Goodness of fit simulations
gof1 = gof(fit1,  control = control.gof.ergm(nsim = n_sim, seed = 1234), verbose = T)
#gof2 = gof(fit2,  control = control.gof.ergm(nsim = n_sim, seed = 1234), verbose = T)


# check duration
end_time = Sys.time()
duration = end_time-start_time


# Saving outputs
save(fit1, gof1, duration, file = paste0(path, "/models/ERGM/estimation_mle_", year,".RData"))
#save(fit1, fit2, gof1, gof2, duration, file = paste0(path, "/models/ERGM/estimation_mle_", year,".RData"))
#load(file = paste0(path, "/models/ERGM/estimation_mle_", year,".RData"))


# Save Summary 
sink(file = paste0("scripts/ERGM/1 summary mle ", year, ".txt"))
summary(fit1)
cat("\n \n \n")
summary(fit2)
sink()


# Summary Statistics
sink(file = paste0("scripts/ERGM/1 gof statistics mle ", year, ".txt"))
options(width = 160)
print(gof1)
cat("\n \n Model Layer Dependence: \n \n")
print(gof2)
sink()


# Save GOF Plots
pdf(paste0("scripts/ERGM/1 gof statistics mle indep ", year, ".pdf"), width = 5, height = 7)
par(mfrow = c(3,2))
plot(gof1, main = "")
dev.off()
pdf(paste0("scripts/ERGM/1 gof statistics mle dep ", year, ".pdf"), width = 5, height = 7)
par(mfrow = c(3,2))
plot(gof2, main = "")
dev.off()


# Save MCMC Diagnostics
pdf(paste0("scripts/ERGM/1 mcmc statistics mle indep ", year, ".pdf"), paper = "a4", width = 5, height = 7)
mcmc.diagnostics(fit1, which = "plots", vars.per.page = 9)
dev.off()
pdf(paste0("scripts/ERGM/1 mcmc statistics mle dep ", year, ".pdf"), paper = "a4", width = 5, height = 7)
mcmc.diagnostics(fit2, which = "plots", vars.per.page = 9)
dev.off()


# Assessing fit of cross-layer terms
dyadfit_observed = summary(net ~ duplexdyad(type = c("e", "f","h"), layers = list(1, 2)))
dyadfit_reduced = simulate(fit1, nsim = 100, monitor = ~duplexdyad(type = c("e","f", "h"), layers = list(1, 2)), output = "stats", seed = 1234)
dyadfit_full = simulate(fit2, nsim = 100, output = "stats", seed = 1234)

out = cbind(
  names(dyadfit_observed),
  dyadfit_observed,
  colMeans(dyadfit_reduced[, (ncol(dyadfit_reduced)-2) : ncol(dyadfit_reduced)]), 
  round(apply(dyadfit_reduced[, (ncol(dyadfit_reduced)-2) : ncol(dyadfit_reduced)], 2, sd),2),
  colMeans(dyadfit_full[, (ncol(dyadfit_full)-2) : ncol(dyadfit_full)]), 
  round(apply(dyadfit_full[, (ncol(dyadfit_full)-2) : ncol(dyadfit_full)], 2, sd),2)
)

out.table = xtable(out, auto = TRUE, label = "tab:dyadfit", caption = "Estimates based on 100 simulated networks.")
colnames(out.table) = c("cross-layer effect", "observed", "mean", "sd", "mean", "sd")
a_header <- construct_header(
  out,
  grp_names = c("", "layer independence", "layer dependence"), 
  span = c(2, 2, 2), 
  align = "c"
)

print(out.table,
      booktabs = T,
      include.rownames = F,
      add.to.row = a_header,
      hline.after = F,
      table.placement = "htp",
      file = paste0("scripts/ERGM/1 latex crosslayer fit mle ", year, ".txt")
)


# Save Latex Coefficients Table
out = matrix(NA, length(coefficients(fit2)), 7)
out = data.frame(out)

out = cbind(
  names(coefficients(fit2)),
  c(round(coefficients(fit1), 2), c(NA, NA, NA)),
  c(round(confint(fit1), 2), c(NA, NA, NA)),
  round(c(coefficients(fit2)), 2),
  round(confint(fit2), 2)
)

#out = rbind(out[1:10, ], rep(NA, 7), out[11:24, ], rep(NA, 7), out[25:26,])
out.table = xtable(out, auto = TRUE, label = "tab:mergm_mle_2003", caption = "Results for the year 2003. Multilayer exponential random graph model estimated with mcmc mle.")
names(out.table) = c("variable", "estimate", "lower ci", "upper ci", "estimate", "lower ci", "upper ci")
align(out.table) = "llrrrrrr"

a_header <- construct_header(
  # the data.frame or matrix that should be plotted  
  out,
  # the labels of the groups that we want to insert
  grp_names = c("", "layer independence", "layer dependence"), 
  # the number of columns each group spans
  span = c(1, 3, 3), 
  # the alignment of each group, can be a single character (lcr) or a vector
  align = "c"
)

print(out.table,
      booktabs = TRUE,
      include.rownames = FALSE,
      add.to.row = a_header,
      hline.after = F,
      table.placement = "htp",
      file = paste0("scripts/ERGM/1 latex coeffs mle ", year, ".txt")
)


Sys.time()

