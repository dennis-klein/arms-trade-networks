library(network)
library(sna)
library(data.table)
library(ergm)
library(multilayer.ergm)
library(mgcv)
library(stargazer)


set.seed(1234)
rm(list = ls(all.names = TRUE))

source("utils/utils.R")
path = data_path.get()


# load data 
load(file.path(path, "out/EX.RData"))
load(file.path(path, "out/colony.RData"))
load(file.path(path, "out/country_list.RData"))


# results table in latex year 2002
fit1 = readRDS("scripts/ERGM/output/model_indep_2002.rds")
fit2 = readRDS("scripts/ERGM/output/model_dep_2002.rds")
stargazer(fit1, fit2, title="ERGM Estimation Results Year 2002", align = TRUE, 
          column.labels = c("Independence", "Dependence"),
          star.cutoffs = c(NA, NA, NA),
          notes = "MPLE Estimates. Std. Error in Brackets.",  notes.append = FALSE)



# read in change statistics
data = list(), start = 2000; end = 2016

for (year in start:end){
  model = readRDS(file = paste("scripts/ERGM/output/model_change_", year,".rds", sep = ""))
  ans = data.table(y = model$response, x = model$predictor, w = model$weights)
  ans[, t := year]
  data[[year-1999]] = ans
}

data = rbindlist(data)



# binomial glm + splines
fit1 = glm.fit(y = data$y, x = data[, 2:11], weights = data$w, family = binomial())
summary.glm(fit1)

frml = y ~ x.mutual.same.layer.mem.1 + 
  x.mutual.same.layer.mem.2 + 
  x.edges_layer.1 +
  x.edges_layer.2 +
  x.gwesp.layer.1.fixed.3 +
  x.gwesp.layer.2.fixed.3 +
  x.edgecov.layer.2.tmp_cdist +
  x.edgecov.layer.1.nmc_ocov +
  x.edgecov.layer.2.gdp_icov +
  x.edgecov.layer.1.pol_absdiff

fit2 = gam(frml, weights = w, data = data, family = binomial())


