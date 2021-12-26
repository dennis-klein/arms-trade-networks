library(network)
library(sna)
library(data.table)
library(ergm)
library(multilayer.ergm)
library(stargazer)


# help 
?'multilayer_terms'
?'ergm-terms'


set.seed(1234)
rm(list = ls(all.names = TRUE))
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

start = 1995
end = 2017
year = 2003

thrshld1 = 0
thrshld2 = 0 

i1 = year - 1949
i2 = year - 1994
i3 = year - 1999

included = rowSums(EX[, (start:end)-1949]) == length(start:end)


# construct dependent nets; 
# layer 1: arms, layer 2: trade + apply thresholds
n = sum(included)

mat1 = (sipri_tiv[[i1]][included, included] > thrshld1) * 1
mat2 = (trade[[i2]][included, included] > thrshld2) * 1
net = to.multiplex(mat1, mat2, output = "network"); check.multilayer(net)

free = matrix(c(rep(c(rep(1, n), rep(0, n)), n), rep(c(rep(0, n), rep(1, n)), n)), nrow = 2*n, byrow = T)
diag(free) = 0
free = network(free, directed = T)


# other edge co-variates, lagged by 1
log_cdist = log(cdist[included, included] + 1)
alliance = atop_alliance[[i1 - 1]][included, included]
nmc_icov = matrix(nmc_cinc[included, i1-1], length(nmc_cinc[included, i1-1]), n, byrow = TRUE)
nmc_ocov = matrix(nmc_cinc[included, i1-1], length(nmc_cinc[included, i1-1]), n, byrow = FALSE)
gdp_icov = matrix(log(gdp[included, i1-1]), length(gdp[included, i1-1]), n, byrow = TRUE)
gdp_ocov = matrix(log(gdp[included, i1-1]), length(gdp[included, i1-1]), n, byrow = FALSE)
conflict_icov = matrix(conflict[included, i1-1], length(conflict[included, i1-1]), n, byrow = TRUE)
conflict_ocov = matrix(conflict[included, i1-1], length(conflict[included, i1-1]), n, byrow = FALSE)
pol_absdiff = abs(outer(polity[included, i1-1], polity[included, i1-1],'-'))
pathdep_arms = (sipri_tiv[[i1 - 1]][included, included] > thrshld1) * 1
pathdep_trade = (trade[[i2-1]][included, included] > thrshld2) * 1


# contrastive divergence estimation
# layer independence model
fit1 <- ergm(net ~ mutual(same = "layer.mem", diff = TRUE) +
  edges_layer(layer = 1) +
  edges_layer(layer = 2) +
  gwesp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
  gwesp_layer(decay = 1.5, fixed = TRUE, layer = 2) +
  gwdsp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
  gwdsp_layer(decay = 1.5, fixed = TRUE, layer = 2) +
  edgecov_layer(pol_absdiff, layer = 1) +
  edgecov_layer(pol_absdiff, layer = 2) +
  edgecov_layer(alliance, layer = 1) +
  edgecov_layer(alliance, layer = 2) +
  edgecov_layer(nmc_ocov, layer = 1) +
  edgecov_layer(nmc_ocov, layer = 2) +
  edgecov_layer(log_cdist, layer = 1) +
  edgecov_layer(log_cdist, layer = 2) +
  edgecov_layer(gdp_icov, layer = 1) +
  edgecov_layer(gdp_icov, layer = 2) +
  edgecov_layer(gdp_ocov, layer = 1) +
  edgecov_layer(gdp_ocov, layer = 2) +
  edgecov_layer(pathdep_arms, layer = 1) +
  edgecov_layer(pathdep_arms, layer = 2) +
  edgecov_layer(pathdep_trade, layer = 1) +
  edgecov_layer(pathdep_trade, layer = 2),
eval.loglik = TRUE, check.degeneracy = TRUE,
verbose = TRUE, estimate = c("CD"),
constraints = ~ fixallbut(free)
)


# layer dependence model
fit2 <- ergm(net ~ mutual(same = "layer.mem", diff = TRUE) +
  edges_layer(layer = 1) +
  edges_layer(layer = 2) +
  gwesp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
  gwesp_layer(decay = 1.5, fixed = TRUE, layer = 2) +
  gwdsp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
  gwdsp_layer(decay = 1.5, fixed = TRUE, layer = 2) +
  edgecov_layer(pol_absdiff, layer = 1) +
  edgecov_layer(pol_absdiff, layer = 2) +
  edgecov_layer(alliance, layer = 1) +
  edgecov_layer(alliance, layer = 2) +
  edgecov_layer(nmc_ocov, layer = 1) +
  edgecov_layer(nmc_ocov, layer = 2) +
  edgecov_layer(log_cdist, layer = 1) +
  edgecov_layer(log_cdist, layer = 2) +
  edgecov_layer(gdp_icov, layer = 1) +
  edgecov_layer(gdp_icov, layer = 2) +
  edgecov_layer(gdp_ocov, layer = 1) +
  edgecov_layer(gdp_ocov, layer = 2) +
  edgecov_layer(pathdep_arms, layer = 1) +
  edgecov_layer(pathdep_arms, layer = 2) +
  edgecov_layer(pathdep_trade, layer = 1) +
  edgecov_layer(pathdep_trade, layer = 2) +
  duplexdyad(c("e", "f", "g", "h"), layers = list(1, 2)),
eval.loglik = TRUE, check.degeneracy = TRUE,
verbose = TRUE, estimate = c("CD"),
constraints = ~ fixallbut(free)
)

fit1 <-logLik(fit1, add=TRUE)
fit2 <-logLik(fit2, add=TRUE)

gof1 = gof(fit1,  control = control.gof.ergm(nsim = 200, seed = 1234), verbose = TRUE)
gof2 = gof(fit2,  control = control.gof.ergm(nsim = 200, seed = 1234), verbose = TRUE)


# Saving outputs
save(fit1, fit2, gof1, gof2, file = paste0(path, "/models/ERGM/estimation_", year,".RData"))



# Save Summary 
sink(file = paste0("scripts/ERGM/1 summary ", year, ".txt"))
summary(fit1)
cat("\n \n \n")
summary(fit2)
sink()


# Save GOF Plots
pdf(paste0("figures/1 gof statistics ", year, ".pdf"), paper = "a4r", width=11, height=8)
par(mfrow = c(2,3))
plot(gof1, main = paste0("Goodness-of-fit diagnostics: Independence Model ", year))
plot.new()
plot(gof2, main = paste0("Goodness-of-fit diagnostics: Dependence Model ", year))
plot.new()
dev.off()


# Save Coefficients Table for Latex
sink(file = paste0("scripts/ERGM/1 latex coeffs ", year, ".txt"))
stargazer(fit1, fit2, 
          title="Multilayer ERGM - Results for the Year 2003", 
          align = TRUE, 
          column.labels = c("Independence", "Dependence"),
          star.cutoffs = c(NA, NA, NA), 
          model.numbers = FALSE,
          single.row = TRUE,
          no.space = TRUE,
          label = "table:results2003",
          keep.stat=c("n"),
          notes = "Estimated with Contrastive Divergence.",  
          notes.append = FALSE
          )
sink()


# Assessing fit of cross-layer terms
dyadfit_observed = summary(net ~ duplexdyad(type = c("e", "f", "g", "h"), layers = list(1, 2)))

dyadfit_reduced = simulate(fit1, nsim = 500, 
                           monitor = ~duplexdyad(type = c("e", "f", "h"), layers = list(1, 2)),
                           output = "stats", seed = 1234)[, 11:13]

dyadfit_full = simulate(fit2, nsim = 500, output = "stats", seed = 1234)[, 11:13]


# Output Table 2
ans = cbind(dyadfit_observed,
      colMeans(dyadfit_reduced), apply(dyadfit_reduced, 2, sd),
      colMeans(dyadfit_full), apply(dyadfit_full, 2, sd)
)
colnames(ans) = c("dyadfit_observed", 
                  "reduced mean", "reduced sd", 
                  "full mean", "full sd")

sink(file = "scripts/ERGM/1 table crosslayer stats.txt")
ans
sink()

