library(data.table)
library(countrycode)


rm(list = ls(all.names = TRUE))


# https://correlatesofwar.org/data-sets/national-material-capabilities  
load("data/out/country_list.RData")
nmc = fread("data/raw/NMC-60-abridged.csv")

nmc = nmc[year %in% 1950:2018]
nmc[stateabb == "GFR", ccode := 255]
nmc[, id := match(ccode, country_list$COW)]

tmp = matrix(0, 257, 69)
rownames(tmp) = country_list$V1
colnames(tmp) = 1950:2018

tmp[cbind(nmc$id, nmc$year-1949)] = nmc$cinc
saveRDS(tmp, file = "data/out/nmc_cinc.rds")



# https://correlatesofwar.org/data-sets/COW-war/intra-state-wars-v5-1.zip/view
intrastate = fread("data/raw/INTRA-STATE_State_participants v5.1 CSV.csv")
intrastate[EndYr1 == -7, EndYr1 := 2014]

preserve <- intrastate
help_func <- function(i) return(c(intrastate$StartYr1[i]:intrastate$EndYr1[i]))
list_years <- lapply(1:nrow(intrastate), FUN = help_func)

intrastate[, year := list_years]
intrastate = intrastate[,.(year = unlist(year)), by = setdiff(names(intrastate), 'year')]
intrastate = intrastate[year %in% c(1950:2018)]
intrastate = intrastate[, .(WarName, CcodeA, SideA, SideB, year)]

intrastate[CcodeA == -8] # what to do with them? 

intrastate = intrastate[CcodeA != -8][,conflict := 1, by = .(CcodeA, SideA, SideB, year)]
intrastate[, id := match(CcodeA, country_list$COW)]

tmp = matrix(0, 257, 69)
rownames(tmp) = country_list$V1
colnames(tmp) = 1950:2018

tmp[cbind(intrastate$id, intrastate$year-1949)] = intrastate$conflict 
saveRDS(tmp, file = "data/out/conflict_intrastate.rds")
