library(network)
library(sna)
library(data.table)
library(ergm)
library(multilayer.ergm)


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

year = 2001
start = 2000
end = 2016

thrshld1 = 0
thrshld2 = 0 

# hyper = expand.grid(seq(2, 5, 1), seq(2, 5, 1), KEEP.OUT.ATTRS = TRUE)
hyper = expand.grid(seq(0, 2, 0.1), seq(0, 2, 0.1), KEEP.OUT.ATTRS = T)
hyper$aic = NA

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



# grid search for hyper param. layer 2 mining
for (i in 1:nrow(hyper)) {

  # layer dependence model
  fit <- ergm(net ~ mutual(same = "layer.mem", diff = TRUE) +
    edges_layer(layer = 1) +
    edges_layer(layer = 2) +
    gwesp_layer(decay = hyper[i, 1], fixed = TRUE, layer = 1) +
    gwdsp_layer(decay = hyper[i, 1], fixed = TRUE, layer = 2) +
    edgecov_layer(pol_absdiff, layer = 1) +
    edgecov_layer(nmc_ocov, layer = 1) +
    edgecov_layer(tmp_cdist, layer = 2) +
    edgecov_layer(gdp_icov, layer = 2) +
    duplexdyad(c("e", "f", "h"), layers = list(1, 2)),
  eval.loglik = TRUE, check.degeneracy = TRUE,
  verbose = TRUE, estimate = c("MPLE"),
  constraints = ~ fixallbut(free)
  )

  hyper[i, "aic"] <- summary(fit)$aic[1]
}

hyper[which.min(hyper$aic), 1]
hyper[which.min(hyper$aic), 2]

saveRDS(hyper, "scripts/ERGM/models/hyper_sipri_trade.rds")



# grid search for hyper param. layer 2 trade
hyper$aic = NA






