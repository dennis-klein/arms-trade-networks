library(data.table)
library(countrycode)
library(readxl)
library(imputeTS)



rm(list = ls(all.names = TRUE))
load("data/out/country_list.RData")
load("data/out/EX.RData")


# https://correlatesofwar.org/data-sets/national-material-capabilities  
nmc = fread("data/raw/NMC-60-abridged.csv")

nmc = nmc[year %in% 1950:2018]
nmc[stateabb == "GFR", ccode := 255]
nmc[, id := match(ccode, country_list$COW)]

tmp = matrix(0, 257, 69)
rownames(tmp) = country_list$V1
colnames(tmp) = 1950:2018

tmp[cbind(nmc$id, nmc$year-1949)] = nmc$cinc
saveRDS(tmp, file = "data/out/nmc_cinc.rds")



# http://www.systemicpeace.org/inscrdata.html
polity = as.data.table(read_xls("data/raw/p5v2018.xls"))
polity = polity[year>1949 & year<2019]

polity[, id := match(ccode, country_list$COW)]; ind_na = is.na(polity$id)
polity[ind_na, id := match(country, country_list$SIPRI)] 
unique(polity[is.na(id),"country"]) # Germany, USSR, Vietnam

polity$id[polity$ccode == 364] = 176 # Sovietunion is 365 not 364 (id 176)
polity$id[polity$ccode == 818] = 214 # Viet Nam is 816 not 818  (id 214)
polity$id[which(polity$ccode == 260)] = 72 # Germany has two cow ids 

tmp = matrix(0, 257, 69)
rownames(tmp) = country_list$V1
colnames(tmp) = 1950:2018
tmp[cbind(polity$id, polity$year-1949)] = polity$polity2

ind_na = is.na(rowSums(tmp))
tmp[ind_na,]

tmp[country_list$V1 == "Hungary", 7] = -7
tmp[country_list$V1 == "India", 1:2] = -7
tmp[country_list$V1 == "Bosnia and Herzegovina", ] = 0
tmp[country_list$V1 == "Cambodia", 30:39] = -7
tmp[country_list$V1 == "Lebanon", 41:55] = 3
tmp[country_list$V1 == "German Democratic Republic", 40:42] = 0
tmp[country_list$V1 == "Japan", 1:2] = 10
tmp[country_list$V1 == "Afghanistan", 30:39] = -7
tmp[country_list$V1 == "Afghanistan", 52:64] = -4
tmp[country_list$V1 == "Czechoslovakia", 19] = -7
tmp[country_list$V1 == "Czechoslovakia", 30:39] = -7
tmp[country_list$V1 == "Iraq", 30:39] = -7
tmp[country_list$V1 == "Iraq", 54:60] = 0
tmp[country_list$V1 == "South Viet Nam", 16:23] = -2
tmp[country_list$V1 == "Kuwait", 41] = -9
tmp[country_list$V1 == "Solomon Islands", 54] = 4
tmp[country_list$V1 == "Somalia", 62] = 3
tmp[country_list$V1 == "Syria", 9:11] = 5
tmp[country_list$V1 == "Uganda", 30] = 0

saveRDS(tmp, file = "data/out/polity.rds")



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



# https://data.worldbank.org/indicator/NY.GDP.MKTP.KD (GDP constant 2015USD)
gdp = fread("data/raw/API_NY.GDP.MKTP.KD_DS2_en_csv_v2_3358328/API_NY.GDP.MKTP.KD_DS2_en_csv_v2_3358328.csv", header = TRUE)
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

# gdp_mat[is.na(gdp_mat)] = 0

tmp[as.numeric(rownames(gdp_mat)), 11:69] = gdp_mat



# https://www.rug.nl/ggdc/historicaldevelopment/maddison/releases/maddison-project-database-2020
maddison = as.data.table(read_xlsx("data/raw/mpd2020.xlsx", sheet = "Full data"))
maddison = maddison[year > 1949 & year < 2019]  
maddison[, gdp := gdppc * pop] 

maddison[, id := match(countrycode, country_list$iso3)]
unique(maddison$country[is.na(maddison$id)])

maddison[country == "Former USSR", id := 154]
maddison[country == "Czechoslovakia", id := 52]
maddison = maddison[!is.na(id)]

tmp2 = matrix(NA, 257, 69)
rownames(tmp2) = country_list$V1
colnames(tmp2) = 1950:2018

tmp2[cbind(maddison$id, maddison$year-1949)] = maddison$gdp
tmp[is.na(tmp)] = tmp2[is.na(tmp)] # before imputing correct rebase of values

which(rowSums(is.na(tmp))==69) 


# linear interpolation of time series 
gdp_na = which(is.na(tmp) , arr.ind = T)
existing = which(matrix((rowSums(EX) != 0)*1, 257, 69) == 1 , arr.ind = T)
to_impute = rbind(gdp_na, existing)[duplicated(rbind(gdp_na, existing)), ]


# only impute if at least 40% of the time series is observed
for (i in unique(to_impute[, 1])){
  if(sum(is.na(tmp[i, ])) <  sum(EX[i,])*0.4){
    # plotNA.distribution(tmp[i,])
    tmp[i,] = na_interpolation(tmp[i,])
  }
}

saveRDS(tmp, file = "data/out/gdp.rds")