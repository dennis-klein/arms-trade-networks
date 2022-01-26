library(data.table)
library(countrycode)
library(readxl)
library(imputeTS)


rm(list = ls(all.names = TRUE))

source("utils/utils.R")
path = data_path.get()


# load data
load(file.path(path, "out/EX.RData"))
load(file.path(path, "out/colony.RData"))
load(file.path(path, "out/country_list.RData"))


# https://correlatesofwar.org/data-sets/national-material-capabilities  
nmc = fread(file.path(path, "raw/NMC-60-abridged.csv"))

nmc = nmc[year %in% 1950:2018]
nmc[stateabb == "GFR", ccode := 255]
nmc[, id := match(ccode, country_list$COW)]

tmp = matrix(0, 257, 69)
rownames(tmp) = country_list$V1
colnames(tmp) = 1950:2018

tmp[cbind(nmc$id, nmc$year-1949)] = nmc$cinc
tmp[, 2017-1949] <- tmp[, 2016-1949]
tmp[, 2018-1949] <- tmp[, 2016-1949]

saveRDS(tmp, file = file.path(path, "out/nmc_cinc.rds"))



# http://www.systemicpeace.org/inscrdata.html
polity = as.data.table(read_xls(file.path(path, "raw/p5v2018.xls")))
polity = polity[year>1949 & year<2019]

polity[, id := match(ccode, country_list$COW)]; ind_na = is.na(polity$id)
polity[ind_na, id := match(country, country_list$SIPRI)] 
unique(polity[is.na(id),"country"]) # Germany, USSR, Vietnam

polity$id[polity$ccode == 364] = 176 # Sovietunion is 365 not 364 (id 176)
polity$id[polity$ccode == 818] = 214 # Viet Nam is 816 not 818  (id 214)
polity$id[which(polity$ccode == 260)] = 72 # Germany has two cow ids 

tmp = matrix(NA, 257, 69)
rownames(tmp) = country_list$V1
colnames(tmp) = 1950:2018
tmp[cbind(polity$id, polity$year-1949)] = polity$polity2


# imputation strategy: impute if at least 3 observations given in period 1995:2018
ind_na = is.na(rowSums(tmp))
#tmp[ind_na,]

for (i in 1:nrow(tmp)){
  if(sum(!is.na(tmp[i, (1995:2018)-1949])) >= 3){
    tmp[i,]=na_interpolation(tmp[i,])
  }
}


# check if countries still have NA observations
ind_incl = rowSums(EX[, (1995:2018)-1949]) == length(1995:2018)
any(is.na(tmp[ind_incl, (1995:2018)-1949]))

saveRDS(tmp, file = file.path(path, "out/polity.rds"))



# https://correlatesofwar.org/data-sets/COW-war/intra-state-wars-v5-1.zip/view
intrastate = fread(file.path(path, "raw/INTRA-STATE_State_participants v5.1 CSV.csv"))
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
saveRDS(tmp, file = file.path(path, "out/conflict_intrastate.rds"))



# rebase time series with cpi us 2010
cpi_us = fread(file.path(path, "raw/united-states.index.cpi-u (statbureau.org).csv"))
setnames(cpi_us, "Year", "year")
cpi_us = cpi_us[year %in% c(1919:2020)]

months = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
cpi_us = melt(cpi_us, c("year"), measure.vars = months, variable.name = "month", value.name = "cpi")
cpi_us = cpi_us[, .(cpi = mean(cpi)), by = .(year)]

cpi_factor = as.numeric(cpi_us[year == 2010, 2]) / as.numeric(cpi_us[year == 2015, 2])



# https://data.worldbank.org/indicator/NY.GDP.MKTP.KD (GDP constant 2015USD)
gdp = fread(file.path(path, "raw/API_NY.GDP.MKTP.KD_DS2_en_csv_v2_3469458/API_NY.GDP.MKTP.KD_DS2_en_csv_v2_3469458.csv"), header = TRUE)
gdp$'Country Name'[gdp$'Country Name'== "Congo, Dem. Rep."] = "DR Congo"

gdp$id = match(gdp$'Country Code', country_list$iso3)
tmp_ind_missing = which(is.na(gdp$id))

gdp$id[tmp_ind_missing] = match(gdp$'Country Name'[tmp_ind_missing], country_list$V1)
tmp_ind_missing = which(is.na(gdp$id)) 
gdp$'Country Name'[tmp_ind_missing] # All unmatched countries  are not independent states 

gdp = gdp[-tmp_ind_missing]

tmp = matrix(NA, 257, 69)
rownames(tmp) = country_list$V1
colnames(tmp) = 1950:2018

gdp_mat =as.matrix(gdp[, 5:63])
colnames(gdp_mat) = 1960:2018
rownames(gdp_mat) = gdp$id 


# fill table and rebase to 2010 values
tmp[as.numeric(rownames(gdp_mat)), 11:69] = gdp_mat * cpi_factor

for (i in 1:nrow(tmp)){
  if(sum(!is.na(tmp[i, (1995:2018)-1949])) >= 12){
    tmp[i,]=na_interpolation(tmp[i,])
  }
}

any(is.na(tmp[ind_incl, (1995:2018)-1949]))

saveRDS(tmp, file = file.path(path, "out/gdp.rds"))



# https://data.worldbank.org/indicator/NY.GDP.PCAP.KD (GDP per capita constant 2015USD)
gdppc = fread(file.path(path, "raw/API_NY.GDP.PCAP.KD_DS2_en_csv_v2_3469368/API_NY.GDP.PCAP.KD_DS2_en_csv_v2_3469368.csv"), header = TRUE)
gdppc$'Country Name'[gdppc$'Country Name'== "Congo, Dem. Rep."] = "DR Congo"

gdppc$id = match(gdppc$'Country Code', country_list$iso3)
tmp_ind_missing = which(is.na(gdppc$id))

gdppc$id[tmp_ind_missing] = match(gdppc$'Country Name'[tmp_ind_missing], country_list$V1)
tmp_ind_missing = which(is.na(gdppc$id)) 
gdppc$'Country Name'[tmp_ind_missing] # All unmatched countries  are not independent states 

gdppc = gdppc[-tmp_ind_missing]

tmp = matrix(NA, 257, 69)
rownames(tmp) = country_list$V1
colnames(tmp) = 1950:2018tmp[ind_na,]

gdppc_mat =as.matrix(gdppc[, 5:63])
colnames(gdppc_mat) = 1960:2018
rownames(gdppc_mat) = gdppc$id 


# fill table and rebase to 2010 values
tmp[as.numeric(rownames(gdppc_mat)), 11:69] = gdppc_mat * cpi_factor

for (i in 1:nrow(tmp)){
  if(sum(!is.na(tmp[i, (1995:2018)-1949])) >= 12){
    tmp[i,]=na_interpolation(tmp[i,])
  }
}

any(is.na(tmp[ind_incl, (1995:2018)-1949]))

saveRDS(tmp, file = file.path(path, "out/gdppc.rds"))
