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
load(file.path(path, "out/milit_exp.RData"))

conflict = readRDS(file.path(path, "out/conflict_intrastate.rds"))
nmc_cinc = readRDS(file.path(path, "out/nmc_cinc.rds")) * 100
polity = readRDS(file.path(path, "out/polity.rds"))
gdp = readRDS(file.path(path, "out/gdp.rds"))
sipri_tiv = readRDS(file.path(path, "out/sipri_tiv.rds"))
mining = readRDS(file.path(path, "out/itpd_mining-energy.rds"))


start = 2000
end = 2017




# test for NA's in the data
for (year in start:end){
  
  i1 = year - 1949
  i2 = year - 1994
  i3 = year - 1999
  present = (EX[, i1] == 1)
  n = sum(present)
  
  tmp_cdist = log(cdist[present, present] + 1)
  nmc_icov = matrix(nmc_cinc[present, i1-1], length(nmc_cinc[present, i1-1]), n, byrow = TRUE)
  nmc_ocov = matrix(nmc_cinc[present, i1-1], length(nmc_cinc[present, i1-1]), n, byrow = FALSE)
  gdp_icov = matrix(log(gdp[present, i1-1]), length(gdp[present, i1-1]), n, byrow = TRUE)
  gdp_ocov = matrix(log(gdp[present, i1-1]), length(gdp[present, i1-1]), n, byrow = FALSE)
  conflict_icov = matrix(conflict[present, i1-1], length(conflict[present, i1-1]), n, byrow = TRUE)
  conflict_ocov = matrix(conflict[present, i1-1], length(conflict[present, i1-1]), n, byrow = FALSE)
  pol_absdiff = abs(outer(polity[present, i1-1], polity[present, i1-1],'-'))
  
  print(year)
  print(any(is.na(tmp_cdist)))
  print(any(is.na(nmc_icov)))
  print(any(is.na(gdp_icov)))
  print(any(is.na(conflict_icov)))
  print(any(is.na(pol_absdiff)))
  
  
}


for (year in start:end){
  i1 = year - 1949
  present = (EX[, i1] == 1)
  n = sum(present)
  
  tmp1 = matrix(milit_exp[present, i1-1], n, n, byrow = TRUE)
  tmp2 = matrix(milit_exp[present, i1-2], n, n, byrow = TRUE)
  tmp3 = matrix(milit_exp[present, i1-3], n, n, byrow = TRUE)
  tmp4 = matrix(milit_exp[present, i1-4], n, n, byrow = TRUE)
  
  print(any(is.na(tmp1)))
  print(any(is.na(tmp2)))
  print(any(is.na(tmp3)))
  print(any(is.na(tmp4)))
}
