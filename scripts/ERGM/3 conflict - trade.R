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

midb_conflict = readRDS(file.path(path, "out/midb_conflict.rds"))
nmc_cinc = readRDS(file.path(path, "out/nmc_cinc.rds")) * 100
polity = readRDS(file.path(path, "out/polity.rds"))
gdp = readRDS(file.path(path, "out/gdp.rds"))
sipri_tiv = readRDS(file.path(path, "out/sipri_tiv.rds"))
trade = readRDS(file.path(path, "out/baci_aggregated.rds"))

year = 2003
included = EX[, 2003-1949] == 1
n = sum(included)


mat1 = (sipri_tiv[[2003-1949]][included, included] > 0) * 1
mat1 = ((mat1 + t(mat1)) > 0) * 1
mat2 = (midb_conflict[[2003-1994]][included, included] > 0) * 1
mat3 = (trade[[2003-1994]][included, included] > 0) * 1
mat4 = t((trade[[2003-1994]][included, included] > 0) * 1)


net = rbind(cbind(mat1, mat3), cbind(mat4, mat2))
net = network(net, directed = FALSE)
set.vertex.attribute(net, "layer.mem", c(rep(1, n), rep(2, n)))
check.multilayer(net)


free = matrix(c(rep(c(rep(1, n), rep(0, n)), n), rep(c(rep(0, n), rep(1, n)), n)), nrow = 2*n, byrow = T)
diag(free) = 0
free = network(free, directed = FALSE)


log_cdist = log(cdist[included, included] + 1)
pol_absdiff = abs(outer(polity[included, 2003-1949], polity[included, 2003-1949],'-'))


# layer dependence model
fit <- ergm(net ~ degree_layer(0, layer = 1) +
  degree_layer(0, layer = 2) +
  altkstar.fixed_layer(layer = 1) +
  altkstar.fixed_layer(layer = 2) +
  edges_layer(layer = 1) +
  edges_layer(layer = 2) +
  gwesp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
  gwesp_layer(decay = 1.5, fixed = TRUE, layer = 2) +
  edgecov_layer(log_cdist, layer = 1) +
  edgecov_layer(log_cdist, layer = 2) +
  edgecov_layer(pol_absdiff, layer = 1) +
  edgecov_layer(pol_absdiff, layer = 2) +
  altkstar.fixed_crosslayer(lambda = 1, layers = list(1, c(1, 2))) +
  altkstar.fixed_crosslayer(lambda = 1, layers = list(c(1, 2), 2)) +
  threetrail_crosslayer(layers = list(1, c(1, 2), 2), incident = 1),
eval.loglik = TRUE, check.degeneracy = TRUE,
verbose = TRUE, estimate = c("CD"),
constraints = ~ fixallbut(free)
)

gof = gof(fit, control = control.gof.ergm(nsim = 500, seed = 1234), verbose = TRUE)


# Saving ERGM outputs
save(fit, gof, file = paste0("scripts/ERGM/models/conflict-trade_", year,".RData"))





