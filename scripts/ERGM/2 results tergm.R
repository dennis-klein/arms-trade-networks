library(network)
library(data.table)
library(xtable)
library(btergm)
library(multilayer.ergm)


rm(list = ls(all.names = TRUE))
source("utils/utils.R")
path = data_path.get()


# load model
vs = "E"
load(file = file.path(path, paste0("models/ERGM/model_estimated_", vs, ".RData")))


# gather estimates and confint
Year = (start + maxlag):(end-2)
confidence =  cbind(confint(results[[1]]), rep(start+maxlag, nrow(confint(results[[1]]))))
confidence = as.data.frame(confidence)
confidence$name = rownames(confidence); rownames(confidence) = NULL

for (i in 2:length(results)){
  yr = start + maxlag + i - 1
  tmp = cbind(confint(results[[i]]), rep(yr, nrow(confint(results[[i]]))))
  tmp = as.data.frame(tmp)
  tmp$name = rownames(tmp); rownames(tmp) = NULL
  confidence = rbind(confidence, tmp)
}

colnames(confidence) = c("Estimate", "LCI", "UCI", "Year", "name")
names = unique(confidence$name)


# plot estimates
pdf(paste0("figures/btergm_estimates_model_", vs, ".pdf"), paper = "a4", height = 10.5)
par(mfrow = c(4, 2))

for (object in names){
  tmp = confidence[confidence$name == object, ]
  plot(tmp$Year, tmp$Estimate, type = "b", pch = 20, ylab = "Estimate", xlab = "Year", 
       main = object, ylim = c(min(tmp$LCI)-0.1, max(tmp$UCI)+0.1))
  arrows(x0 = tmp$Year, y0 = tmp$LCI, x1 = tmp$Year, y1 = tmp$UCI, code=3, angle=90, length=0.025)
  if ((min(tmp$LCI)-0.2 <= 0) & (0 <= max(tmp$UCI)+0.2)) {abline(h = 0, col="grey")}
}

dev.off()



# test output latex table
# add bootstrapping sample size as note
print(xtable(confint(fit), auto = TRUE, type = "latex"), file = "output_latex.txt")




