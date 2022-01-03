library(network)
library(sna)
library(data.table)
library(ergm)
library(multilayer.ergm)
library(stargazer)
library(foreach)
library(parallel)
library(xtable)


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


# set parameters of analysis
# trade data starting 1995, sipri ends 2017
start = 1995
end = 2017
year = 2003
n_sim = 200
n_bootstrap = 100

thrshld1 = 0
thrshld2 = 0 

i1 = year - 1949 # Cornelius Replication Data
i2 = year - 1994 # CEPII BACI
i3 = year - 1999


# countries included in analysis
included = rowSums(EX[, (start:end)-1949]) == length(start:end)
n = sum(included) 


# construct dependent multilayer network
# layer 1: arms, layer 2: trade + apply thresholds
tmp1 = (sipri_tiv[[i1]][included, included] > thrshld1) * 1
tmp2 = (trade[[i2]][included, included] > thrshld2) * 1
net = to.multiplex(tmp1, tmp2, output = "network")
check.multilayer(net)


# define sample space constraint
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
             check.degeneracy = TRUE,
             verbose = FALSE, 
             estimate = c("CD"),
             control = control.ergm(
               CD.nsteps = 16,
               CD.multiplicity = 2
             ),
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
             check.degeneracy = TRUE,
             verbose = FALSE, 
             estimate = c("CD"),
             control = control.ergm(
               CD.nsteps = 16,
               CD.multiplicity = 2
             ),
             constraints = ~ fixallbut(free)
)


# add aic and bic 
fit1 <-logLik(fit1, add=TRUE)
fit2 <-logLik(fit2, add=TRUE)


# Goodness of fit simulations
gof1 = gof(fit1,  control = control.gof.ergm(nsim = n_sim, seed = 1234), verbose = T)
gof2 = gof(fit2,  control = control.gof.ergm(nsim = n_sim, seed = 1234), verbose = T)


# Save Summary 
sink(file = paste0("scripts/ERGM/1 summary cd ", year, ".txt"))
summary(fit1)
cat("\n \n \n")
summary(fit2)
sink()


# Summary Statistics
sink(file = paste0("scripts/ERGM/1 gof statistics cd ", year, ".txt"))
options(width = 160)
print(gof1)
cat("\n \n Model Layer Independence: \n \n")
print(gof2)
sink()


# Save GOF Plots
pdf(paste0("scripts/ERGM/1 gof statistics cd ", year, ".pdf"), paper = "a4r", width=11, height=8)
par(mfrow = c(2,3))
plot(gof1, main = paste0("Goodness-of-fit diagnostics: Independence Model ", year))
plot.new()
plot(gof2, main = paste0("Goodness-of-fit diagnostics: Dependence Model ", year))
plot.new()
dev.off()


# Saving outputs
save(fit1, fit2, gof1, gof2, file = paste0(path, "/models/ERGM/estimation_cd_", year,".RData"))


# Assessing fit of cross-layer terms
#dyadfit_observed = summary(net ~ duplexdyad(type = c("e", "f", "g", "h"), layers = list(1, 2)))
#dyadfit_full = simulate(fit2, nsim = 500, output = "stats", seed = 1234)[, 11:13]
#dyadfit_reduced = simulate(fit1, nsim = 500, monitor = ~duplexdyad(type = c("e", "f", "g", "h"), layers = list(1, 2)),output = "stats", seed = 1234)[, 11:13]


# Output Table 2
#ans = cbind(dyadfit_observed,
#      colMeans(dyadfit_reduced), apply(dyadfit_reduced, 2, sd),
#      colMeans(dyadfit_full), apply(dyadfit_full, 2, sd)
#)

#colnames(ans) = c("observed", 
#                 "reduced mean", "reduced sd", 
#                "full mean", "full sd")

#sink(file = paste0("scripts/ERGM/1 latex coeffs cd ", year, ".txt"))
#ans
#sink()


# bootstrap ci's if estimated with contrastive divergence
#define quantiles of interest
q = c(.05, .95)


# layer independence model
cat("Simulating ", n_bootstrap , " networks\n")
fit1_sim = simulate(fit1, constraints = ~ fixallbut(free), nsim = n_bootstrap, seed = 1234, verbose = T)
n_variables = length(coef(fit1))


# parallel not working - memory issue ?
# use %dopar%
#n.cores = parallel::detectCores() - 1
#n.cores = 3
#my.cluster = parallel::makeCluster(n.cores, type = "PSOCK")
#doParallel::registerDoParallel(cl = my.cluster)
#foreach::getDoParRegistered()


cat("Estimating simulated networks\n")
bootstrap <- foreach(
  i = 1:n_bootstrap, 
  .export = 'ergm',
  .packages = c('ergm', 'multilayer.ergm'), 
  .combine = rbind
) %do% {
  
  # replace response network in formula
  tt = terms(formula(fit1))
  sim_net = fit1_sim[[i]]
  new_formula = reformulate(attr(tt, "term.labels"), "sim_net")
  
  # estimate simulated
  tmp <- suppressMessages(
    ergm(new_formula, 
         estimate = c("CD"),
         control = control.ergm(
           CD.nsteps = 8,
           CD.multiplicity = 1
         ),
         constraints = ~ fixallbut(free)
    )
  )
  
  ans = as.data.table(coef(tmp),  keep.rownames = T)
  
} # end parallel

fit1_bconf = bootstrap[,  .(q5 = quantile(V2, probs = q[1]), q95 = quantile(V2, probs = q[2])), by = .(V1)]


# layer dependence model
cat("Simulating ", n_bootstrap , " networks (dependence)\n")
fit2_sim = simulate(fit2, constraints = ~ fixallbut(free), nsim = n_bootstrap, seed = 1234, verbose = T)
n_variables = length(coef(fit2))


cat("Estimating simulated networks (dependence)\n")
bootstrap <- foreach(
  i = 1:n_bootstrap, 
  .export = 'ergm', 
  .packages = c('ergm', 'multilayer.ergm'),  
  .combine = rbind
) %do% {
  # replace response network in formula
  tt = terms(formula(fit2))
  sim_net = fit2_sim[[i]]
  new_formula = reformulate(attr(tt, "term.labels"), "sim_net")
  
  tmp <- suppressMessages(
    ergm(new_formula, 
         estimate = c("CD"),
         control = control.ergm(
           CD.nsteps = 8,
           CD.multiplicity = 1
         ),
         constraints = ~ fixallbut(free)
    )
  )
  
  ans = as.data.table(coef(tmp),  keep.rownames = T)
  
} # end parallel


#parallel::stopCluster(cl = my.cluster)

fit2_bconf = bootstrap[,  .(q5 = quantile(V2, probs = q[1]), q95 = quantile(V2, probs = q[2])), by = .(V1)]


# Saving outputs
save(fit1, fit2, fit1_bconf, fit2_bconf, file = paste0(path, "/models/ERGM/estimation_cd_boot_", year,".RData"))
#load(file = paste0(path, "/models/ERGM/estimation_cd_boot_", year,".RData"))


# Save Latex Coefficients Table
out = matrix(NA, 28, 7)
out = data.frame(out)

out[, 1] = names(coefficients(fit2))
out[1:26, 2] = coefficients(fit1)
out[1:26, 3] = fit1_bconf$q5
out[1:26, 4] = fit1_bconf$q95
out[, 5]= c(coefficients(fit2))
out[, 6] = fit2_bconf$q5
out[, 7] = fit2_bconf$q95

out = rbind(out[1:10, ], rep(NA, 7), out[11:26, ], rep(NA, 7), out[27:28,])
out.table = xtable(out, auto = TRUE, digits=2, label = "tab:mergm_cd_2003", caption = "Estimated with Contrastive Divergence. Confidence intervals calculated based on 100 bootstrap iterations. 95pct Confidence Intervals provided.")
names(out.table) = c("Effect", "Layer Indep.", "LCI", "UCI", "Layer Dep.", "LCI", "UCI")
align(out.table) = "llrrrrrr"

print(out.table,
  booktabs = TRUE,
  include.rownames = FALSE,
  file = "scripts/ERGM/1 latex coeffs cd 2003.txt"
)

