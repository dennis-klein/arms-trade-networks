library(network)
library(sna)
library(data.table)
library(ergm)
library(multilayer.ergm)
library(foreach)
library(parallel)
library(xtable)
library(texreg)


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
n_bootstrap = 100


# set correct indices depending on data source
i1 = year - 1949 # Cornelius Replication Data
i2 = year - 1994 # CEPII BACI


# selection of countries included
included = rowSums(EX[, (start:end)-1949]) == length(start:end)
n = sum(included) 


# construct dependent multi-layer net
# layer 1: arms, layer 2: trade + apply thresholds
tmp1 = (sipri_tiv[[i1]][included, included] > 0) * 1
tmp2 = custom_trade_to_binary(trade[[i2]][included, included], type = "C", threshold = 0.01)
net = to.multiplex(tmp1, tmp2, output = "network")
check.multilayer(net)


# define sample space constraint
free = matrix(c(rep(c(rep(1, n), rep(0, n)), n), rep(c(rep(0, n), rep(1, n)), n)), nrow = 2*n, byrow = T)
diag(free) = 0
free = network(free, directed = T)


# nodal and dyadic attributes, lagged by t = 1
log_cdist = log(cdist[included, included] + 1)
log_gdp_out = matrix(log(gdp[included, i1-1]), n, n, byrow = FALSE)
log_gdppc_in = matrix(log(gdppc[included, i1-1]), n, n, byrow = TRUE)
absdiff_polity = abs(outer(polity[included, i1-1], polity[included, i1-1],'-'))
cinc100_out = matrix(100 * nmc_cinc[included, i1-1], n, n, byrow = F)
alliance = atop_alliance[[i1 - 1]][included, included]

pathdep_arms = (sipri_tiv[[i1-1]][included, included] > 0) * 1
pathdep_trade = custom_trade_to_binary(trade[[i2-1]][included, included], type = "C", threshold = 0.01)


# contrastive divergence estimation
# layer independence model
fit1 <- ergm(net ~ mutual(same = "layer.mem", diff = TRUE) +
               gwidegree(decay = 1, fixed = TRUE, attr = "layer.mem") +
               gwodegree(decay = 1, fixed = TRUE, attr = "layer.mem") +
               edges_layer(layer = 1) +
               gwesp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
               edgecov_layer(log_cdist, layer = 1) +
               edgecov_layer(log_gdppc_in, layer = 1) +
               edgecov_layer(log_gdp_out, layer = 1) +
               edgecov_layer(cinc100_out, layer = 1) +
               edgecov_layer(alliance, layer = 1) + 
               edgecov_layer(absdiff_polity, layer = 1) +
               
               edges_layer(layer = 2) +
               gwesp_layer(decay = 1.5, fixed = TRUE, layer = 2) +
               edgecov_layer(log_cdist, layer = 2) +
               edgecov_layer(log_gdppc_in, layer = 2) +
               edgecov_layer(log_gdp_out, layer = 2) +
               edgecov_layer(cinc100_out, layer = 2) +
               edgecov_layer(alliance, layer = 2) +
               edgecov_layer(absdiff_polity, layer = 2),
             check.degeneracy = TRUE,
             verbose = F, 
             estimate = c("CD"),
             control = control.ergm(
               parallel = 4,
               CD.nsteps = 4096,
               CD.multiplicity = 1
             ),
             constraints = ~ fixallbut(free)
)


# layer dependence model
fit2 <- ergm(net ~ mutual(same = "layer.mem", diff = TRUE) +
               gwidegree(decay = 1, fixed = TRUE, attr = "layer.mem") +
               gwodegree(decay = 1, fixed = TRUE, attr = "layer.mem") +
               edges_layer(layer = 1) +
               gwesp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
               edgecov_layer(log_cdist, layer = 1) +
               edgecov_layer(log_gdppc_in, layer = 1) +
               edgecov_layer(log_gdp_out, layer = 1) +
               edgecov_layer(cinc100_out, layer = 1) +
               edgecov_layer(alliance, layer = 1) + 
               edgecov_layer(absdiff_polity, layer = 1) +
               
               edges_layer(layer = 2) +
               gwesp_layer(decay = 1.5, fixed = TRUE, layer = 2) +
               edgecov_layer(log_cdist, layer = 2) +
               edgecov_layer(log_gdppc_in, layer = 2) +
               edgecov_layer(log_gdp_out, layer = 2) +
               edgecov_layer(cinc100_out, layer = 2) +
               edgecov_layer(alliance, layer = 2) +
               edgecov_layer(absdiff_polity, layer = 2) +
               duplexdyad(c("e", "f", "h"), layers = list(1, 2)),
             check.degeneracy = TRUE,
             verbose = FALSE, 
             estimate = c("CD"),
             control = control.ergm(
               parallel = 4,
               CD.nsteps = 4096,
               CD.multiplicity = 2
             ),
             constraints = ~ fixallbut(free)
)


fit2_mple <- ergm(net ~ mutual(same = "layer.mem", diff = TRUE) +
               gwidegree(decay = 1, fixed = TRUE, attr = "layer.mem") +
               gwodegree(decay = 1, fixed = TRUE, attr = "layer.mem") +
               edges_layer(layer = 1) +
               gwesp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
               edgecov_layer(log_cdist, layer = 1) +
               edgecov_layer(log_gdppc_in, layer = 1) +
               edgecov_layer(log_gdp_out, layer = 1) +
               edgecov_layer(cinc100_out, layer = 1) +
               edgecov_layer(alliance, layer = 1) + 
               edgecov_layer(absdiff_polity, layer = 1) +
               
               edges_layer(layer = 2) +
               gwesp_layer(decay = 1.5, fixed = TRUE, layer = 2) +
               edgecov_layer(log_cdist, layer = 2) +
               edgecov_layer(log_gdppc_in, layer = 2) +
               edgecov_layer(log_gdp_out, layer = 2) +
               edgecov_layer(cinc100_out, layer = 2) +
               edgecov_layer(alliance, layer = 2) +
               edgecov_layer(absdiff_polity, layer = 2) +
               duplexdyad(c("e", "f", "h"), layers = list(1, 2)),
             check.degeneracy = TRUE,
             verbose = FALSE, 
             estimate = c("MPLE"),
             constraints = ~ fixallbut(free)
)


# add aic and bic 
#fit1 <-logLik(fit1, add=TRUE)
#fit2 <-logLik(fit2, add=TRUE)


# Goodness of fit Simulations
gof1 = gof(fit1, verbose = T)
gof2 = gof(fit2, verbose = T)
gof2_mple = gof(fit2_mple, verbose = T)


# Saving outputs
save(fit1, fit2, gof1, gof2, fit2_mple, gof2_mple, file = paste0(path, "/models/ERGM/estimation_cd_", year,".RData"))
load(file = paste0(path, "/models/ERGM/estimation_cd_", year,".RData"))


# Save Summary 
sink(file = paste0("scripts/ERGM/1 screenreg cd ", year, ".txt"))
screenreg(list(fit1, fit2, fit2_mple))
sink()

Sys.time()

# Summary Statistics
#sink(file = paste0("scripts/ERGM/1 gof statistics cd ", year, ".txt"))
#options(width = 160)
#print(gof1)
#cat("\n \n Model Layer Independence: \n \n")
#print(gof2)
#sink()


# Save GOF Plots
pdf(paste0("scripts/ERGM/1 gof statistics cd ", year, ".pdf"), paper = "a4r", width=11, height=8)
par(mfrow = c(2,3))
plot(gof1, main = paste0("Goodness-of-fit diagnostics: Independence Model ", year))
plot.new()
plot(gof2, main = paste0("Goodness-of-fit diagnostics: Dependence Model cd ", year))
plot.new()
plot(gof2, main = paste0("Goodness-of-fit diagnostics: Dependence Model mple ", year))
plot.new()
dev.off()



if(FALSE){
# Bootstrap ci's if estimated with contrastive divergence
# Define quantiles of interest
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
out = data.frame(matrix(NA, length(coef(fit2)), 7))
out[, 1] = names(coefficients(fit2))
out[1:length(coef(fit1)), 2] = round(coefficients(fit1), 2)
out[1:length(coef(fit1)), 3] = round(fit1_bconf$q5, 2)
out[1:length(coef(fit1)), 4] = round(fit1_bconf$q95, 2)
out[, 5]= round(c(coefficients(fit2)), 2)
out[, 6] = round(fit2_bconf$q5, 2)
out[, 7] = round(fit2_bconf$q95, 2)

#out = rbind(out[1:10, ], rep(NA, 7), out[11:24, ], rep(NA, 7), out[25:26,])
out.table = xtable(out, auto = TRUE, label = "tab:mergm_cd_2003", caption = "Results for the year 2003. Multilayer exponential random graph model estimated with contrastive divergence. Confidence intervals based on 100 bootstrap iterations, 95pct confidence intervals provided.")
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
  table.placement = "!h",
  file = paste0("scripts/ERGM/1 latex coeffs cd ", year, ".txt")
)



}