library(data.table)
library(countrycode)
library(readxl)


rm(list = ls(all.names = TRUE))


load("data/out/country_list.RData")


# http://www.cepii.fr/cepii/en/bdd_modele/download.asp?id=37
# - HS92 (1995-2019, 2.18 Go)

#  t	Year  k	Product category (HS 6-digit code)
#  i	Exporter (ISO 3-digit country code) j	Importer (ISO 3-digit country code)
#  v	Value of the trade flow (in thousands current USD) q	Quantity (in metric tons)

codes = fread("data/raw/CEPII BACI/country_codes_V202102.csv")
file_lst = paste0("data/raw/CEPII BACI/BACI_HS92_V202102/", list.files(path = "data/raw/CEPII BACI/BACI_HS92_V202102/")) 
baci = data.table()

for (string in file_lst){
  tmp = fread(string, select = c("t", "i", "j", "k", "v"))
  tmp = tmp[, .(v = sum(v)), by = .(t, i, j)]
  baci = rbind(baci, tmp)
}

baci = merge(baci, codes[, .(country_code, iso_3digit_alpha)], by.x = "i", by.y = "country_code", all.x = TRUE)
setnames(baci, "iso_3digit_alpha", "i_iso3")
baci = merge(baci, codes[, c("country_code", "iso_3digit_alpha")], by.x = "j", by.y = "country_code", all.x = TRUE)
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


# rebase to 2005 values
cpi_us <- fread("data/raw/united-states.index.cpi-u (statbureau.org).csv")
setnames(cpi_us, "Year", "year")
cpi_us = cpi_us[year %in% c(1919:2020)]

months <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
cpi_us = melt(cpi_us, c("year"), measure.vars = months, variable.name = "month", value.name = "cpi")
cpi_us = cpi_us[, .(cpi = mean(cpi)), by = .(year)]
cpi_us[, cpi_rebased := cpi / as.numeric(filter(cpi_us, year == 2005)[,2])]

baci = merge(baci, cpi_us, by.x = "t", by.y = "year", all.x = TRUE)
baci = baci[, v := v / cpi_rebased][, .(t, i, j, v)]
tmp = list(); 

yr = 1995
for (i in 1:25){
  tmp[[i]] = matrix(0, 257, 257)
  colnames(tmp[[i]]) = country_list$V1
  rownames(tmp[[i]]) = country_list$V1
  
  tmp_dfl = baci[t == yr]; yr = yr + 1 # don't ask me why
  tmp[[i]][cbind(tmp_dfl$from_id, tmp_dfl$to_id)] = tmp_dfl$v
}

saveRDS(tmp, file = "data/out/baci_aggregated.rds")
rm(tmp, tmp_dfl, baci, codes)




# https://www.usitc.gov/data/gravity/itpde.htm
# 17 years, 2000-2016, 243 countries, 170 industries
# includes 26 industries in agriculture, 7 in mining and energy,
# 120 in manufacturing, and 17 in services.
# https://www.usitc.gov/data/gravity/itpde_guide/industries_list/

itpd = fread("data/raw/itpd_e_r01/ITPD_E_R01.csv")
itpd = itpd[, .(trade = sum(trade)), by = .(exporter_iso3, importer_iso3, year)]

itpd[, from_id:=match(exporter_iso3, country_list$iso3)]
itpd[, to_id:=match(importer_iso3, country_list$iso3)]

unique(itpd$importer_iso3[is.na(itpd$to_id)]) 
unique(itpd$exporter_iso3[is.na(itpd$from_id)]) 

itpd = itpd[!(is.na(itpd$from_id) | is.na(itpd$to_id))] # Delete all observations where the id is not known 
tmp = list()

for (i in 2000:2016){
  tmp[[i-1999]] <- matrix(0, 257, 257)
  colnames(tmp[[i-1999]]) <- country_list$V1
  rownames(tmp[[i-1999]]) <- country_list$V1
  
  tmp_dfl = itpd[year == i]
  tmp[[i- 1999]][cbind(tmp_dfl$from_id, tmp_dfl$to_id)] = tmp_dfl$trade
}

saveRDS(tmp, file = "data/out/itpd_aggregated.rds")
rm(tmp, tmp_dfl, itpd)




# CEPII TRADEHIST <http://www.cepii.fr/cepii/en/bdd_modele/presentation.asp?id=32>
#   Bilateral Trade and Bilateral Tariffs
# Exchange Rates <https://www.measuringworth.com/datasets/exchangeglobal/>
#   CPI US <https://www.statbureau.org/en/united-states/cpi-u>

TRADHIST_BITRADE_BITARIFF_1 <- read_excel("data/raw/CEPII TRADEHIST/TRADHIST_BITRADE_BITARIFF_1.xlsx")
TRADHIST_BITRADE_BITARIFF_2 <- read_excel("data/raw/CEPII TRADEHIST/TRADHIST_BITRADE_BITARIFF_2.xlsx")
TRADHIST_BITRADE_BITARIFF_3 <- read_excel("data/raw/CEPII TRADEHIST/TRADHIST_BITRADE_BITARIFF_3.xlsx")

tradehist <- rbindlist(list(
  TRADHIST_BITRADE_BITARIFF_1,
  TRADHIST_BITRADE_BITARIFF_2,
  TRADHIST_BITRADE_BITARIFF_3
), fill = TRUE)

rm(TRADHIST_BITRADE_BITARIFF_1,TRADHIST_BITRADE_BITARIFF_2,TRADHIST_BITRADE_BITARIFF_3)


# currency is British Pound Sterling. 
# nominal trade flows, rebase to constant 2005$ USD.

exchange_usd_gbp <- fread("data/raw/EXCHANGEGLOBAL_1945-2020.csv") 
setnames(exchange_usd_gbp, c("V1", "V2"), c("year", "GBP_USD"))
exchange_usd_gbp[, GBP_USD := as.numeric(gsub(" British Pound", "", GBP_USD))]

cpi_us[, cpi_rebased := cpi / as.numeric(filter(cpi_us, year == 2005)[,2])]


# merge data, filter for out time-frame. 
# do conversions and replace flow_0 values. 

tradehist = tradehist[year %in% c(1950:2018)]
tradehist = tradehist[FLOW_0 == 0, FLOW := 0]

tradehist = merge(tradehist, exchange_usd_gbp, by = c("year"), all.x = TRUE)
tradehist = merge(tradehist, cpi_us, by = c("year"), all.x = TRUE)

tradehist = tradehist[, FLOW := (FLOW / GBP_USD) / cpi_rebased][, .(iso_o, iso_d, year, FLOW)]


# Match the country names with the ids in country_list
tradehist[, from_id := match(tradehist$iso_o, country_list$cepii)]
tradehist[, to_id := match(tradehist$iso_d, country_list$cepii)]


# check not matched countries
unique(tradehist$iso_o[is.na(tradehist$from_id)])
unique(tradehist$iso_d[is.na(tradehist$to_id)])

  # some are micro states and will be dropped in the analysis
  # how to handle East / West Germany, Russia - has to be adapted
  # check these again!

tradehist$from_id[tradehist$iso_o == "USSR"] <- 154
tradehist$to_id[tradehist$iso_d == "USSR"] <- 154

tradehist$from_id[tradehist$iso_o == "EDEU"] <- 71
tradehist$to_id[tradehist$iso_d == "EDEU"] <- 71

tradehist$from_id[tradehist$iso_o == "WDEU"] <- 72
tradehist$to_id[tradehist$iso_d == "WDEU"] <- 72

tradehist$from_id[tradehist$iso_o == "CZSK"] <- 52
tradehist$to_id[tradehist$iso_d == "CZSK"] <- 52

tradehist$from_id[tradehist$iso_o == "ROM"] <- 153
tradehist$to_id[tradehist$iso_d == "ROM"] <- 153

# exclude countries which are not part of the analysis
tradehist = tradehist[!(is.na(from_id) | is.na(to_id))]
tmp = list()

# create a list of 69 adjacency matrices, one for each year from 1950:2018
for (i in 1950:2018) {
  tmp[[i-1949]] = matrix(0, 257, 257)
  colnames(tmp[[i-1949]]) = country_list$V1
  rownames(tmp[[i-1949]]) = country_list$V1
  
  tmp_dfl = tradehist[year == i]
  tmp[[i - 1949]][cbind(tmp_dfl$from_id, tmp_dfl$to_id)] = tmp_dfl$FLOW
}

# save network list 
saveRDS(tmp, file = "data/out/cepii_trade.rds")
rm(tmp, tmp_dfl, tradehist)


