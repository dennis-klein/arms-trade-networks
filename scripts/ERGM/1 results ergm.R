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
fit1 = readRDS("scripts/ERGM/models/model_indep_mple_2003.rds")
fit2 = readRDS("scripts/ERGM/models/model_dep_mple_2003.rds")



pdf("scripts/ERGM/models/gof_statistics.pdf", paper = "a4r", width = 11, height = 8)
par(mfrow = c(2,3))
sim = gof(fit1, control = control.gof.ergm(parallel = 4), verbose = TRUE)
plot(sim, main = paste0("Goodness-of-fit diagnostics: Independence Model ", year))
plot.new()
sim = gof(fit2, control = control.gof.ergm(parallel = 4), verbose = TRUE)
plot(sim, main = paste0("Goodness-of-fit diagnostics: Dependence Model ", year))
plot.new()
dev.off()


stargazer(fit1, fit2, title="ERGM Estimation Results Year 2002", align = TRUE, 
          column.labels = c("Independence", "Dependence"),
          star.cutoffs = c(NA, NA, NA), 
          notes = "",  notes.append = FALSE)



