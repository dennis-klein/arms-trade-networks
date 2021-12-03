library(network)
library(sna)
library(data.table)
library(ergm)
library(multilayer.ergm)
library(mgcv)


# help 
?'multilayer_terms'
?'ergm-terms'


set.seed(1234)
rm(list = ls(all.names = TRUE))


# load data 
load("data/out/EX.RData")
load("data/out/colony.RData")
load("data/out/country_list.RData")
load("data/out/sipri_tiv.RData")
load("data/out/polity.RData")
load("data/out/cdist.RData")
load("data/out/real_gdp_p_c.RData")

nmc_cinc = readRDS("data/out/nmc_cinc.rds")
gdp = readRDS("data/out/gdp.rds")
baci_aggregated = readRDS("data/out/baci_aggregated.rds")


period = 1996:2018 # cepii baci availability + EX limit
thrshld1 = 0
thrshld2 = 100 

info <- matrix(NA, 2, length(1996:2018), dimnames = list(c("AIC Indep.", "AIC Dep."), c(period)))
change = list()


for (year in period){
  year =
  i1 = year - 1949
  i2 = year - 1994
  i3 = year - 1995
  present = (EX[, i1] == 1)
  
  
  # construct dependent net, layer 1: arms, layer 2: trade + apply thresholds
  n = sum(present)
  mat1 = (sipri_tiv[[i1]][present, present] > thrshld1) * 1
  mat2 = (baci_aggregated[[i2]][present, present] > thrshld2) * 1
  
  net = to.multiplex(mat1, mat2, output = "network"); check.multilayer(net)
  free = to.multiplex(matrix(1, n, n), matrix(1, n, n), output = "network", offzeros = TRUE)
  
  
  # other edge co-variates, lagged by 1
  tmp_cdist = log(cdist[present, present] + 1)
  nmc_icov = matrix(nmc_cinc[present, i1-1], length(nmc_cinc[present, i1-1]), n, byrow = TRUE)
  nmc_ocov = matrix(nmc_cinc[present, i1-1], length(nmc_cinc[present, i1-1]), n, byrow = FALSE)
  gdp_icov = matrix(gdp[present, i1-1], length(gdp[present, i1-1]), n, byrow = TRUE)
  gdp_ocov = matrix(gdp[present, i1-1], length(gdp[present, i1-1]), n, byrow = FALSE)
  pol_absdiff = abs(outer(polity[present, i1-1], polity[present, i1-1],'-'))
  
  
  # layer independence model
  fit = ergm(net ~ mutual(same="layer.mem", diff = TRUE) +
               edges_layer(layer = 1) + 
               edges_layer(layer = 2) +
               gwesp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
               gwesp_layer(decay = 1.5, fixed = TRUE, layer = 2) +
               edgecov_layer(tmp_cdist, layer = 1) +
               edgecov_layer(tmp_cdist, layer = 2) +
               edgecov_layer(nmc_ocov, layer = 1) +
               edgecov_layer(gdp_icov, layer = 2) +
               edgecov_layer(pol_absdiff, layer = 1),
             eval.loglik = TRUE, check.degeneracy = TRUE, verbose = TRUE,
             control = control.ergm(seed = 1234, parallel = 4), constraints = ~fixallbut(free)
  )
  
  info[1, i3] = summary(fit)$aic[1]
  saveRDS(fit, file = paste("scripts/ERGM/output/model_indep_", year,".rds", sep = ""))
  
  
  fit = ergmMPLE(net ~ mutual(same="layer.mem", diff = TRUE) +
               edges_layer(layer = 1) + 
               edges_layer(layer = 2) +
               gwesp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
               gwesp_layer(decay = 1.5, fixed = TRUE, layer = 2) +
               edgecov_layer(tmp_cdist, layer = 1) +
               edgecov_layer(tmp_cdist, layer = 2) +
               edgecov_layer(nmc_ocov, layer = 1) +
               edgecov_layer(gdp_icov, layer = 2) +
               edgecov_layer(pol_absdiff, layer = 1),
             constraints = ~fixallbut(free)
  )
  
  change[[i3]] = data.table(response = fit$response, predictor = fit$predictor, 
                                   weights = fit$weights, year = year)
  
  
  # layer dependence model
  fit = ergm(net ~ mutual(same="layer.mem", diff = TRUE) +
               edges_layer(layer = 1) + 
               edges_layer(layer = 2) +
               gwesp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
               gwesp_layer(decay = 1.5, fixed = TRUE, layer = 2) +
               edgecov_layer(tmp_cdist, layer = 1) +
               edgecov_layer(tmp_cdist, layer = 2) +
               edgecov_layer(nmc_ocov, layer = 1) +
               edgecov_layer(gdp_icov, layer = 2) +
               edgecov_layer(pol_absdiff, layer = 1) +
               duplexdyad(c("e", "f", "h"), layers = list(1, 2)),
             eval.loglik = TRUE, check.degeneracy = TRUE, verbose = TRUE, estimate = c("MPLE"),
             control = control.ergm(seed = 1234, parallel = 4), constraints = ~fixallbut(free)
  )
  
  info[2, i3] = summary(fit)$aic[1]
  saveRDS(fit, file = paste("scripts/ERGM/output/model_dep_", year,".rds", sep = ""))
  
  
  if(year == 2018){
    saveRDS(info, "scripts/ERGM/output/aic_comparison.rds")
    saveRDS(change, "scripts/ERGM/output/change_list.rds")
  }
}





# gof analysis



# binomial glm + splines
# try_mod =  glm.fit(y = trying$response,x = trying$predictor,weights = trying$weights, family = binomial())

