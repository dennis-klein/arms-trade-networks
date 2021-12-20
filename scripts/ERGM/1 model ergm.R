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
trade = readRDS(file.path(path, "out/itpd_aggregated.rds"))

start = 2000
end = 2016
year = 2001

thrshld1 = 0
thrshld2 = 0 

i1 = year - 1949
i2 = year - 1994
i3 = year - 1999

included = rowSums(EX[, (2000:2016)-1949]) == 17


# construct dependent nets; 
# layer 1: arms, layer 2: trade + apply thresholds
n = sum(included)

mat1 = (sipri_tiv[[i1]][included, included] > thrshld1) * 1
mat2 = (trade[[i3]][included, included] > thrshld2) * 1
net = to.multiplex(mat1, mat2, output = "network"); check.multilayer(net)
free = matrix(c(rep(c(rep(1, n), rep(0, n)), n), rep(c(rep(0, n), rep(1, n)), n)), nrow = 2*n, byrow = T)
diag(free) = 0
free = network(free, directed = T)


# other edge co-variates, lagged by 1
tmp_cdist = log(cdist[included, included] + 1)
nmc_icov = matrix(nmc_cinc[included, i1-1], length(nmc_cinc[included, i1-1]), n, byrow = TRUE)
nmc_ocov = matrix(nmc_cinc[included, i1-1], length(nmc_cinc[included, i1-1]), n, byrow = FALSE)
gdp_icov = matrix(log(gdp[included, i1-1]), length(gdp[included, i1-1]), n, byrow = TRUE)
gdp_ocov = matrix(log(gdp[included, i1-1]), length(gdp[included, i1-1]), n, byrow = FALSE)
conflict_icov = matrix(conflict[included, i1-1], length(conflict[included, i1-1]), n, byrow = TRUE)
conflict_ocov = matrix(conflict[included, i1-1], length(conflict[included, i1-1]), n, byrow = FALSE)
pol_absdiff = abs(outer(polity[included, i1-1], polity[included, i1-1],'-'))


# load hyper parameters
hyper = readRDS("scripts/ERGM/output/hyper_sipri_trade.rds")
d1 = hyper[which.min(hyper$aic), 1]
d2 = hyper[which.min(hyper$aic), 2]


# contrastive divergence estimation
# layer independence model
fit1 <- ergm(net ~ mutual(same = "layer.mem", diff = TRUE) +
  edges_layer(layer = 1) +
  edges_layer(layer = 2) +
  gwesp_layer(decay = d1, fixed = TRUE, layer = 1) +
  gwdsp_layer(decay = d2, fixed = TRUE, layer = 2) +
  edgecov_layer(pol_absdiff, layer = 1) +
  edgecov_layer(nmc_ocov, layer = 1) +
  edgecov_layer(tmp_cdist, layer = 2) +
  edgecov_layer(gdp_icov, layer = 2),
eval.loglik = TRUE, check.degeneracy = TRUE,
verbose = TRUE, estimate = c("MPLE"),
constraints = ~ fixallbut(free)
)


# layer dependence model
fit2 <- ergm(net ~ mutual(same = "layer.mem", diff = TRUE) +
  edges_layer(layer = 1) +
  edges_layer(layer = 2) +
  gwesp_layer(decay = d1, fixed = TRUE, layer = 1) +
  gwdsp_layer(decay = d2, fixed = TRUE, layer = 2) +
  edgecov_layer(pol_absdiff, layer = 1) +
  edgecov_layer(nmc_ocov, layer = 1) +
  edgecov_layer(tmp_cdist, layer = 2) +
  edgecov_layer(gdp_icov, layer = 2) +
  duplexdyad(c("e", "f", "h"), layers = list(1, 2)),
eval.loglik = TRUE, check.degeneracy = TRUE,
verbose = TRUE, estimate = c("MPLE"),
constraints = ~ fixallbut(free)
)


gof1 = gof(fit1,  control = control.gof.ergm(nsim = 500, seed = 1234), verbose = TRUE)
gof2 = gof(fit2,  control = control.gof.ergm(nsim = 500, seed = 1234), verbose = TRUE)


# Saving ERGM outputs
save(fit1, fit2, gof1, gof2, file = paste0("scripts/ERGM/models/estimation_", year,".RData"))



# Output 
sink(file = paste0("scripts/ERGM/models/summary_", year, ".txt"))
summary(fit1)
cat("\n \n \n")
summary(fit2)
sink()


pdf(paste0("scripts/ERGM/models/gof_statistics_", year, ".pdf"), paper = "a4r", width=11, height=8)
par(mfrow = c(2,3))
plot(gof1, main = paste0("Goodness-of-fit diagnostics: Independence Model ", year))
plot.new()
plot(gof2, main = paste0("Goodness-of-fit diagnostics: Dependence Model ", year))
plot.new()
dev.off()


stargazer(fit1, fit2, title="ERGM Estimation Results Year 2001", align = TRUE, 
          column.labels = c("Independence", "Dependence"),
          star.cutoffs = c(NA, NA, NA), 
          notes = "",  notes.append = FALSE)


# Assessing fit of cross-layer terms
dyadfit_observed = summary(net ~ duplexdyad(type = c("e", "f", "h"), layers = list(1, 2)))

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
sink(file = "scripts/ERGM/models/gof_crsslyr_table.txt")
ans
sink()

