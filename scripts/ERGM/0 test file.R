# read in change statistics
data = list(); start = 2000; end = 2016

for (year in start:end){
  model = readRDS(file = paste("scripts/ERGM/models/model_change_", year,".rds", sep = ""))
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


