library(network)
library(sna)
library(data.table)
library(ergm)
library(multilayer.ergm)
library(stargazer)


# help 
?'multilayer_terms'
?'ergm-terms'


# setup
rm(list = ls(all.names = TRUE))
set.seed(1234)
source("utils/utils.R")
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
nsim = 200

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
               edgecov_layer(pol_absdiff, layer = 2) +
               edgecov_layer(alliance, layer = 1) +
               edgecov_layer(alliance, layer = 2) +
               edgecov_layer(nmc_ocov, layer = 1) +
               edgecov_layer(nmc_ocov, layer = 2) +
               edgecov_layer(log_cdist, layer = 1) +
               edgecov_layer(log_cdist, layer = 2) +
               edgecov_layer(lgdp_icov, layer = 1) +
               edgecov_layer(lgdp_icov, layer = 2) +
               edgecov_layer(lgdp_ocov, layer = 1) +
               edgecov_layer(lgdp_ocov, layer = 2) +
               edgecov_layer(pathdep_arms, layer = 1) +
               edgecov_layer(pathdep_arms, layer = 2) +
               edgecov_layer(pathdep_trade, layer = 1) +
               edgecov_layer(pathdep_trade, layer = 2),
             eval.loglik = TRUE, 
             check.degeneracy = TRUE,
             verbose = TRUE, 
             estimate = c("MLE"),
             control = control.ergm(parallel = 4),
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
               edgecov_layer(pol_absdiff, layer = 2) +
               edgecov_layer(alliance, layer = 1) +
               edgecov_layer(alliance, layer = 2) +
               edgecov_layer(nmc_ocov, layer = 1) +
               edgecov_layer(nmc_ocov, layer = 2) +
               edgecov_layer(log_cdist, layer = 1) +
               edgecov_layer(log_cdist, layer = 2) +
               edgecov_layer(lgdp_icov, layer = 1) +
               edgecov_layer(lgdp_icov, layer = 2) +
               edgecov_layer(lgdp_ocov, layer = 1) +
               edgecov_layer(lgdp_ocov, layer = 2) +
               edgecov_layer(pathdep_arms, layer = 1) +
               edgecov_layer(pathdep_arms, layer = 2) +
               edgecov_layer(pathdep_trade, layer = 1) +
               edgecov_layer(pathdep_trade, layer = 2) +
               duplexdyad(c("e", "h"), layers = list(1, 2)),
             eval.loglik = TRUE, 
             check.degeneracy = TRUE,
             verbose = TRUE, 
             estimate = c("MLE"),
             control = control.ergm(parallel = 4),
             constraints = ~ fixallbut(free)
)


end_time <- Sys.time()
duration = end_time-start_time

# add aic and bic if method cd
#fit1 <-logLik(fit1, add=TRUE)
#fit2 <-logLik(fit2, add=TRUE)


# Goodness of fit simulations
gof1 = gof(fit1,  control = control.gof.ergm(nsim = nsim, seed = 1234), verbose = T)
gof2 = gof(fit2,  control = control.gof.ergm(nsim = nsim, seed = 1234), verbose = T)


# Saving outputs
save(fit1, fit2, gof1, gof2, duration, file = paste0(path, "/models/ERGM/estimation_mle_", year,".RData"))
#load(file = paste0(path, "/models/ERGM/estimation_", year,".RData"))

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
cat("\n \n Model Layer Independence: \n \n")
print(gof2)
sink()


# Save GOF Plots
pdf(paste0("scripts/ERGM/1 gof statistics mle ", year, ".pdf"), paper = "a4r", width=11, height=8)
par(mfrow = c(2,3))
plot(gof1, main = paste0("Goodness-of-fit diagnostics: Independence Model ", year))
plot.new()
plot(gof2, main = paste0("Goodness-of-fit diagnostics: Dependence Model ", year))
plot.new()
dev.off()


# Save MCMC Diagnostics
pdf(paste0("scripts/ERGM/1 mcmc statistics mle ", year, ".pdf"), paper = "a4r", width=11, height=8)
par(mfrow = c(3,4))
mcmc.diagnostics(fit1, which = "plots")
mcmc.diagnostics(fit2, which = "plots")
dev.off()


# Save Coefficients Table for Latex
sink(file = paste0("scripts/ERGM/1 latex coeffs mle ", year, ".txt"))
stargazer(fit1, fit2, 
          title="Multilayer ERGM - Results for the Year 2003", 
          align = TRUE, 
          column.labels = c("Layer \n Independence", "Layer \n Dependence"),
          star.cutoffs = c(NA, NA, NA), 
          model.numbers = FALSE,
          single.row = TRUE,
          no.space = TRUE,
          digits = 2,
          label = "table:results2003",
          notes = "Estimated with MCMC MLE.",  
          notes.append = TRUE
          )
sink()
