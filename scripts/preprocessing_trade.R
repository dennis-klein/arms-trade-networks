library(data.table)
library(countrycode)
library(readxl)
library(imputeTS)

rm(list = ls(all.names = TRUE))

source("utils/utils.R")
path = data_path.get()


# load data
load(file.path(path, "out/country_list.RData"))
load(file.path(path, "out/EX.RData"))


# prepare cpi usd series
cpi_us = fread(file.path(path, "raw/united-states.index.cpi-u (statbureau.org).csv"))
setnames(cpi_us, "Year", "year")
cpi_us = cpi_us[year %in% c(1919:2020)]

months = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
cpi_us = melt(cpi_us, c("year"), measure.vars = months, variable.name = "month", value.name = "cpi")
cpi_us = cpi_us[, .(cpi = mean(cpi)), by = .(year)]

fctr = as.numeric(cpi_us[year == 2010, 2])
cpi_us[, cpi_rebased := cpi / fctr]



# http://www.cepii.fr/cepii/en/bdd_modele/download.asp?id=37
# - HS92 (1995-2019, 2.18 Go)

#  t	Year  k	Product category (HS 6-digit code)
#  i	Exporter (ISO 3-digit country code) j	Importer (ISO 3-digit country code)
#  v	Value of the trade flow (in thousands current USD) q	Quantity (in metric tons)

codes = fread(file.path(path,"raw/CEPII BACI/country_codes_V202102.csv"))
file_lst = paste0(file.path(path,"raw/CEPII BACI/BACI_HS92_V202102/"), "/", list.files(file.path(path, "raw/CEPII BACI/BACI_HS92_V202102/"))) 
file_lst = file_lst[-26]
baci = data.table()

for (string in file_lst){
  tmp = fread(string, select = c("t", "i", "j", "k", "v"))
  tmp = tmp[, .(v = sum(v)), by = .(t, i, j)]
  baci = rbind(baci, tmp)
}


# rebase to 2010 values
baci = merge(baci, cpi_us, by.x = "t", by.y = "year", all.x = T)
baci = baci[, v := v * 1000 / cpi_rebased][, .(t, i, j, v)]


# apply relevant exporter / importer thresholds
exp = baci[, .(total = sum(v), count = .N), by = .(t, i)]
imp = baci[, .(total = sum(v), count = .N), by = .(t, j)]
thrshld = 0.01

baci = merge(baci, exp, by.x = c("t", "i"), by.y = c("t", "i"), all.x = T)
baci[, dependency := v / total]
baci[, v := ifelse(dependency > thrshld, v, 0)]


baci = merge(baci, codes[, .(country_code, iso_3digit_alpha)], by.x = "i", by.y = "country_code", all.x = T)
setnames(baci, "iso_3digit_alpha", "i_iso3")
baci = merge(baci, codes[, c("country_code", "iso_3digit_alpha")], by.x = "j", by.y = "country_code", all.x = T)
setnames(baci, "iso_3digit_alpha", "j_iso3")

baci[, i_iso3 := as.factor(i_iso3)]
baci[, j_iso3 := as.factor(j_iso3)]

baci[, from_id:=match(i_iso3, country_list$iso3)]
baci[, to_id:=match(j_iso3, country_list$iso3)]

#baci$from_id[is.na(baci$from_id)] = match(baci$from[is.na(baci$from_id)],country_list$iso3)
#baci$to_id[is.na(baci$to_id)] = match(baci$to[is.na(baci$to_id)],country_list$iso3)

unique(baci$i_iso3[is.na(baci$from_id)]) 
unique(baci$j_iso3[is.na(baci$to_id)]) 

baci = baci[!(is.na(baci$from_id) | is.na(baci$to_id))] # Delete all observations where the id is not known 
tmp = list(); 


for (i in 1:25){
  tmp[[i]] = matrix(0, 257, 257)
  colnames(tmp[[i]]) = country_list$V1
  rownames(tmp[[i]]) = country_list$V1
  
  yr = i + 1994
  tmp_dfl = baci[t == yr]; 
  if (yr != unique(tmp_dfl$t)) {print("Error")}
  tmp[[i]][cbind(tmp_dfl$from_id, tmp_dfl$to_id)] <- tmp_dfl$v
}

saveRDS(tmp, file = file.path(path, "out/baci_aggregated.rds"))
rm(tmp, tmp_dfl, baci, codes)




# https://www.usitc.gov/data/gravity/itpde.htm
# 17 years, 2000-2016, 243 countries, 170 industries
# includes 26 industries in agriculture, 7 in mining and energy,
# 120 in manufacturing, and 17 in services.
# Trade values are expressed in millions of current United States dollars, rebased to 2010$values
# https://www.usitc.gov/data/gravity/itpde_guide/industries_list/

itpd = fread(file.path(path, "raw/itpd_e_r01/ITPD_E_R01.csv"))
itpd[, broad_sector := as.factor(broad_sector)]; levels(itpd$broad_sector)

itpd = itpd[, .(trade = sum(trade)), by = .(exporter_iso3, importer_iso3, year)]
itpd[, trade := trade * 1000000]

itpd = merge(itpd, cpi_us, by.x = "year", by.y = "year", all.x = TRUE)
itpd = itpd[, trade := trade / cpi_rebased][, .(exporter_iso3, importer_iso3, year, trade)]


# apply relevant exporter / importer thresholds
exp = itpd[, .(total = sum(trade), count = .N), by = .(year, exporter_iso3)]
imp = itpd[, .(total = sum(trade), count = .N), by = .(year, importer_iso3)]
thrshld = 0.01

itpd = merge(itpd, exp, by.x = c("year", "exporter_iso3"), by.y = c("year", "exporter_iso3"), all.x = T)
itpd[, dependency := trade / total]
itpd[, trade := ifelse(dependency > thrshld, trade, 0)]


# impute values if at least five are given
#itpd = dcast(itpd, exporter_iso3 + importer_iso3 ~ year, value.var = c("trade"))
#tmp = as.matrix(itpd[, -c(1,2)])

#for(i in 1:nrow(tmp)){
#  if(sum(is.na(tmp[i, ])) <  17-5){
#    tmp[i, ] = na_interpolation(tmp[i, ])
#  }
#}

#tmp = as.data.frame(tmp)
#tmp$exporter_iso3 = itpd$exporter_iso3; tmp$importer_iso3 = itpd$importer_iso3
#setDT(tmp)

#itpd = melt(tmp, id.vars = c("exporter_iso3", "importer_iso3"), measure.vars = 1:17, variable.name = "year", value.name = "trade")
#sum(is.na(itpd$trade))

itpd[, from_id:=match(exporter_iso3, country_list$iso3)]
itpd[, to_id:=match(importer_iso3, country_list$iso3)]

unique(itpd$importer_iso3[is.na(itpd$to_id)]) 
unique(itpd$exporter_iso3[is.na(itpd$from_id)]) 

itpd = itpd[!(is.na(itpd$from_id) | is.na(itpd$to_id))] # Delete all observations where the id is not known 
tmp = list()

for (i in 2000:2016){
  tmp[[i-1999]] = matrix(0, 257, 257)
  colnames(tmp[[i-1999]]) = country_list$V1
  rownames(tmp[[i-1999]]) = country_list$V1
  
  tmp_dfl = itpd[year == i]
  tmp[[i- 1999]][cbind(tmp_dfl$from_id, tmp_dfl$to_id)] = tmp_dfl$trade
}

saveRDS(tmp, file = file.path(path, "out/itpd_aggregated.rds"))
rm(tmp, tmp_dfl, itpd)




# CEPII TRADHIST <http://www.cepii.fr/cepii/en/bdd_modele/presentation.asp?id=32>
# Bilateral Trade and Bilateral Tariffs
# Exchange Rates <https://www.measuringworth.com/datasets/exchangeglobal/>
# CPI US <https://www.statbureau.org/en/united-states/cpi-u>

TRADHIST_BITRADE_BITARIFF_1 = read_excel(file.path(path, "raw/CEPII TRADEHIST/TRADHIST_BITRADE_BITARIFF_1.xlsx"))
TRADHIST_BITRADE_BITARIFF_2 = read_excel(file.path(path, "raw/CEPII TRADEHIST/TRADHIST_BITRADE_BITARIFF_2.xlsx"))
TRADHIST_BITRADE_BITARIFF_3 = read_excel(file.path(path, "raw/CEPII TRADEHIST/TRADHIST_BITRADE_BITARIFF_3.xlsx"))

tradhist = rbindlist(list(
  TRADHIST_BITRADE_BITARIFF_1,
  TRADHIST_BITRADE_BITARIFF_2,
  TRADHIST_BITRADE_BITARIFF_3
), fill = TRUE)

rm(TRADHIST_BITRADE_BITARIFF_1,TRADHIST_BITRADE_BITARIFF_2,TRADHIST_BITRADE_BITARIFF_3)


# currency is British Pound Sterling. 
# nominal trade flows, rebase to constant 2010$ USD.

exchange_usd_gbp = fread(file.path(path, "raw/EXCHANGEGLOBAL_1945-2020.csv"))
setnames(exchange_usd_gbp, c("V1", "V2"), c("year", "GBP_USD"))
exchange_usd_gbp[, GBP_USD := as.numeric(gsub(" British Pound", "", GBP_USD))]

cpi_us[, cpi_rebased := cpi / fctr]


# merge data, filter for out time-frame. 
# do conversions and replace flow_0 values. 

tradhist = tradhist[year %in% 1950:2018]
tradhist = tradhist[FLOW_0 == 0, FLOW := 0]

tradhist = merge(tradhist, exchange_usd_gbp, by = c("year"), all.x = TRUE)
tradhist = merge(tradhist, cpi_us, by = c("year"), all.x = TRUE)

tradhist = tradhist[, FLOW := (FLOW / GBP_USD) / cpi_rebased][, .(iso_o, iso_d, year, FLOW)]


# impute values if at least five are given
tradhist = dcast(tradhist, iso_o + iso_d ~ year, value.var = c("FLOW"))
tmp = as.matrix(tradhist[, -c(1,2)])

for(i in 1:nrow(tmp)){
  if(sum(is.na(tmp[i, ])) <  65-22){
    tmp[i, ] = na_interpolation(tmp[i, ])
  }
}

tmp = as.data.frame(tmp)
tmp$iso_o = tradhist$iso_o; tmp$iso_d = tradhist$iso_d
setDT(tmp)

tradhist = melt(tmp, id.vars = c("iso_o", "iso_d"), measure.vars = 1:65, variable.name = "year", value.name = "FLOW")
sum(is.na(tradhist$FLOW))


# Match the country names with the ids in country_list
tradhist[, from_id := match(iso_o, country_list$cepii)]
tradhist[, to_id := match(iso_d, country_list$cepii)]


# check not matched countries
unique(tradhist$iso_o[is.na(tradhist$from_id)])
unique(tradhist$iso_d[is.na(tradhist$to_id)])

  # some are micro states and will be dropped in the analysis
  # how to handle East / West Germany, Russia - has to be adapted
  # check these again!

tradhist$from_id[tradhist$iso_o == "USSR"] <- 154
tradhist$to_id[tradhist$iso_d == "USSR"] <- 154

tradhist$from_id[tradhist$iso_o == "EDEU"] <- 71
tradhist$to_id[tradhist$iso_d == "EDEU"] <- 71

tradhist$from_id[tradhist$iso_o == "WDEU"] <- 72
tradhist$to_id[tradhist$iso_d == "WDEU"] <- 72

tradhist$from_id[tradhist$iso_o == "CZSK"] <- 52
tradhist$to_id[tradhist$iso_d == "CZSK"] <- 52

tradhist$from_id[tradhist$iso_o == "ROM"] <- 153
tradhist$to_id[tradhist$iso_d == "ROM"] <- 153


# exclude countries which are not part of the analysis
tradhist = tradhist[!(is.na(from_id) | is.na(to_id))]
tmp = list()


# create a list of 65 adjacency matrices, one for each year from 1950:2014
for (i in 1950:2014) {
  tmp[[i-1949]] = matrix(0, 257, 257)
  colnames(tmp[[i-1949]]) = country_list$V1
  rownames(tmp[[i-1949]]) = country_list$V1
  
  tmp_dfl = tradhist[year == i]
  tmp[[i - 1949]][cbind(tmp_dfl$from_id, tmp_dfl$to_id)] = tmp_dfl$FLOW
}


# save network list 
saveRDS(tmp, file = file.path(path, "out/tradhist_aggregated.rds"))
rm(tmp, tmp_dfl, tradhist)


