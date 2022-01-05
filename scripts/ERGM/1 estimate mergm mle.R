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
path = data_path.get()


# load data 
load(file.path(path, "out/EX.RData"))
load(file.path(path, "out/colony.RData"))
load(file.path(path, "out/country_list.RData"))
load(file.path(path, "out/cdist.RData"))

conflict = readRDS(file.path(path, "out/conflict_intrastate.rds"))
nmc_cinc = readRDS(file.path(path, "out/nmc_cinc.rds")) * 100
polity = readRDS(file.path(path, "out/polity.rds"))
gdp = readRDS(file.path(path, "out/gdp.rds"))
sipri_tiv = readRDS(file.path(path, "out/sipri_tiv.rds"))
trade = readRDS(file.path(path, "out/baci_aggregated.rds"))
load(file.path(path, "out/atop_alliance.RData"))


# set parameters of analyses
# trade data starting 1995, sipri ends 2017
start = 1995
end = 2017
year = 2003
n_sim = 200

thrshld1 = 0
thrshld2 = 0 

i1 = year - 1949 # Cornelius Replication Data
i2 = year - 1994 # CEPII BACI
i3 = year - 1999

included = rowSums(EX[, (start:end)-1949]) == length(start:end)
n = sum(included) 


# construct dependent nets 
# layer 1: arms, layer 2: trade + apply thresholds
tmp1 = (sipri_tiv[[i1]][included, included] > thrshld1) * 1
tmp2 = (trade[[i2]][included, included] > thrshld2) * 1
net = to.multiplex(tmp1, tmp2, output = "network")
check.multilayer(net)


# define constraint
free = matrix(c(rep(c(rep(1, n), rep(0, n)), n), rep(c(rep(0, n), rep(1, n)), n)), nrow = 2*n, byrow = T)
diag(free) = 0
free = network(free, directed = T)


# nodal attributes, lagged by t = 1
log_cdist = log(cdist[included, included] + 1)
alliance = atop_alliance[[i1 - 1]][included, included]
nmc_icov = matrix(nmc_cinc[included, i1-1], length(nmc_cinc[included, i1-1]), n, byrow = TRUE)
nmc_ocov = matrix(nmc_cinc[included, i1-1], length(nmc_cinc[included, i1-1]), n, byrow = FALSE)
lgdp_icov = matrix(log(gdp[included, i1-1]), length(gdp[included, i1-1]), n, byrow = TRUE)
lgdp_ocov = matrix(log(gdp[included, i1-1]), length(gdp[included, i1-1]), n, byrow = FALSE)
conflict_icov = matrix(conflict[included, i1-1], length(conflict[included, i1-1]), n, byrow = TRUE)
conflict_ocov = matrix(conflict[included, i1-1], length(conflict[included, i1-1]), n, byrow = FALSE)
pol_absdiff = abs(outer(polity[included, i1-1], polity[included, i1-1],'-'))
pathdep_arms = (sipri_tiv[[i1-1]][included, included] > thrshld1) * 1
pathdep_trade = (trade[[i2-1]][included, included] > thrshld2) * 1


# measure time
start_time = Sys.time()


# contrastive divergence estimation
# layer independence model
fit1 <- ergm(net ~ mutual(same = "layer.mem", diff = TRUE) +
               edges_layer(layer = 1) +
               edges_layer(layer = 2) +
               gwesp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
               gwesp_layer(decay = 1.5, fixed = TRUE, layer = 2) +
               gwidegree(decay = 1, fixed = TRUE, attr = "layer.mem") +
               gwodegree(decay = 1, fixed = TRUE, attr = "layer.mem") +
               edgecov_layer(pol_absdiff, layer = 1) +
               edgecov_layer(alliance, layer = 1) +
               edgecov_layer(nmc_ocov, layer = 1) +
               edgecov_layer(log_cdist, layer = 1) +
               edgecov_layer(lgdp_icov, layer = 1) +
               edgecov_layer(lgdp_ocov, layer = 1) +
               edgecov_layer(pathdep_arms, layer = 1) +
               edgecov_layer(pathdep_trade, layer = 1) +
               edgecov_layer(pol_absdiff, layer = 2) +
               edgecov_layer(alliance, layer = 2) +
               edgecov_layer(nmc_ocov, layer = 2) +
               edgecov_layer(log_cdist, layer = 2) +
               edgecov_layer(lgdp_icov, layer = 2) +
               edgecov_layer(lgdp_ocov, layer = 2) +
               edgecov_layer(pathdep_arms, layer = 2) +
               edgecov_layer(pathdep_trade, layer = 2),
             eval.loglik = TRUE, 
             check.degeneracy = TRUE,
             verbose = TRUE, 
             estimate = c("MLE"),
             control = control.ergm(parallel = 4, seed = 1234),
             constraints = ~ fixallbut(free)
)


# layer dependence model
fit2 <- ergm(net ~ mutual(same = "layer.mem", diff = TRUE) +
               edges_layer(layer = 1) +
               edges_layer(layer = 2) +
               gwesp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
               gwesp_layer(decay = 1.5, fixed = TRUE, layer = 2) +
               gwidegree(decay = 1, fixed = TRUE, attr = "layer.mem") +
               gwodegree(decay = 1, fixed = TRUE, attr = "layer.mem") +
               edgecov_layer(pol_absdiff, layer = 1) +
               edgecov_layer(alliance, layer = 1) +
               edgecov_layer(nmc_ocov, layer = 1) +
               edgecov_layer(log_cdist, layer = 1) +
               edgecov_layer(lgdp_icov, layer = 1) +
               edgecov_layer(lgdp_ocov, layer = 1) +
               edgecov_layer(pathdep_arms, layer = 1) +
               edgecov_layer(pathdep_trade, layer = 1) +
               edgecov_layer(pol_absdiff, layer = 2) +
               edgecov_layer(alliance, layer = 2) +
               edgecov_layer(nmc_ocov, layer = 2) +
               edgecov_layer(log_cdist, layer = 2) +
               edgecov_layer(lgdp_icov, layer = 2) +
               edgecov_layer(lgdp_ocov, layer = 2) +
               edgecov_layer(pathdep_arms, layer = 2) +
               edgecov_layer(pathdep_trade, layer = 2) +
               duplexdyad(c("e", "h"), layers = list(1, 2)),
             eval.loglik = TRUE, 
             check.degeneracy = TRUE,
             verbose = TRUE, 
             estimate = c("MLE"),
             control = control.ergm(parallel = 4, seed = 1234),
             constraints = ~ fixallbut(free)
)


end_time = Sys.time()
duration = end_time-start_time


# Goodness of fit simulations
gof1 = gof(fit1,  control = control.gof.ergm(nsim = n_sim, seed = 1234), verbose = T)
gof2 = gof(fit2,  control = control.gof.ergm(nsim = n_sim, seed = 1234), verbose = T)


# Saving outputs
save(fit1, fit2, gof1, gof2, duration, file = paste0(path, "/models/ERGM/estimation_mle_", year,".RData"))
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
pdf(paste0("scripts/ERGM/1 gof statistics mle indep ", year, ".pdf"), width=5, height=7)
par(mfrow = c(3,2))
plot(gof1, main = "")
dev.off()
pdf(paste0("scripts/ERGM/1 gof statistics mle dep ", year, ".pdf"), width=5, height=7)
par(mfrow = c(3,2))
plot(gof2, main = "")
dev.off()


# Save MCMC Diagnostics
pdf(paste0("scripts/ERGM/1 mcmc statistics mle ", year, ".pdf"), paper = "a4", width=8, height=11)
#mcmc.diagnostics(fit1, which = "plots")
mcmc.diagnostics(fit2, which = "plots", vars.per.page = 9)
dev.off()


# Save Latex Coefficients Table
out = matrix(NA, 28, 7)
out = data.frame(out)

out[, 1] = names(coefficients(fit2))
out[1:26, 2] = round(coefficients(fit1), 2)
out[1:26, 3:4] = round(confint(fit1), 2)
out[, 5]= round(c(coefficients(fit2)), 2)
out[, 6:7] = round(confint(fit2), 2)

out = rbind(out[1:10, ], rep(NA, 7), out[11:26, ], rep(NA, 7), out[27:28,])
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




