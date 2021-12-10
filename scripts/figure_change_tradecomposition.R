library(data.table)
library(countrycode)
library(areaplot)

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


# select countries included in analysis
include = rowSums(EX[, (2000:2016)-1949]) == 17
country_list = country_list[include, ]
sum(include)



# https://www.usitc.gov/data/gravity/itpde.htm
# 17 years, 2000-2016, 243 countries, 170 industries
# includes 26 industries in agriculture, 7 in mining and energy,
# 120 in manufacturing, and 17 in services.
# Trade values are expressed in millions of current United States dollars, rebased to 2010$values
# https://www.usitc.gov/data/gravity/itpde_guide/industries_list/

itpd = fread(file.path(path, "raw/itpd_e_r01/ITPD_E_R01.csv"))
itpd[, broad_sector := as.factor(broad_sector)]; levels(itpd$broad_sector)

itpd = itpd[, .(trade = sum(trade)), by = .(exporter_iso3, importer_iso3, year, broad_sector)]
itpd[, trade := trade / 1000] # in billions

itpd = merge(itpd, cpi_us, by.x = "year", by.y = "year", all.x = TRUE)
itpd = itpd[, trade := trade / cpi_rebased][, .(exporter_iso3, importer_iso3, year, trade, broad_sector)]

itpd[, from_id:=match(exporter_iso3, country_list$iso3)]
itpd[, to_id:=match(importer_iso3, country_list$iso3)]

itpd = itpd[!(is.na(itpd$from_id) | is.na(itpd$to_id))] # Delete all observations where the id is not known 
itpd = itpd[, .(trade = sum(trade)), by = .(year, broad_sector)]

itpd = dcast(itpd, year~broad_sector)

data = as.matrix(itpd)


pdf("figures/change_in_tradecomposition.pdf", family = "ArialMT",
    height = 17.35 / 2.54, width = 23.35 / 2.54, pointsize = 18)

areaplot(data[, 1], data[, 2:5], main = "ITPD: Changes in Trade Composition of Broad Sectors", 
         legend = TRUE, args.legend = list(x = "topleft", cex = 0.75),
         col = hcl.colors(4, palette = "viridis", alpha = 0.8),
         xlab = "Year", ylab = "")

dev.off()
