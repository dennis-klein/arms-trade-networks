#------------------------------------------------------------------------------#
# Summary table to see whats in SIPRI data
#------------------------------------------------------------------------------#

# This script reruns the SIPRI Data Aggregation and filtes for the included countries.

library(readr)
library(data.table)
library(dplyr)
library(tidyr)
library(kableExtra)


rm(list = ls(all.names = TRUE))


# setup
source("utils/utils.R")
source("utils/construct_header.R")
source("utils/custom_trade_to_binary.R")
path <- data_path.get()
set.seed(1234)


# load data 
load(file.path(path, "out/EX.RData"))
load(file.path(path, "out/colony.RData"))
load(file.path(path, "out/country_list.RData"))

# setup 
year <- 2003
start <- 1995
end <- 2018

# selection of countries included: present over the complete time horizon (cf. soam)
included <- rowSums(EX[, (start:end) - 1949]) == length(start:end)
n <- sum(included)


# load sipri data
paths <- list(
  file.path(path, "raw/SIPRI/1950_1970.txt"),
  file.path(path, "raw/SIPRI/1971_1980.txt"),
  file.path(path, "raw/SIPRI/1981_1990.txt"),
  file.path(path, "raw/SIPRI/1991_2000.txt"),
  file.path(path, "raw/SIPRI/2001_2019.txt")
)

dfl = lapply(paths, read.csv, sep = ";", skip = 4) # Read the files into a list
dfl = rbindlist(dfl) # Make one data.table out of the read files
changeCols<- names(dfl)


# Make factors to characters
dfl[, (changeCols) := lapply(.SD, as.character), .SDcols = changeCols]


# Unkown countries are deleted 
dfl = dfl[dfl$Seller != ""]


# There seems to be Problem with deals whose description includes in the raw data a ; in the Description 
missing_tmp = which(suppressWarnings(is.na(as.numeric(as.character(dfl$Delivery.year)))))
dfl[missing_tmp,6:16] = dfl[missing_tmp,7:17]

dfl$Delivery.year = suppressWarnings(as.numeric(as.character(dfl$Delivery.year)))
dfl$Order.date = suppressWarnings(as.numeric(as.character(dfl$Order.date)))

dfl$Seller[dfl$Seller == "Cote d'Ivoire"] = "Cote dIvoire"
dfl$Seller[dfl$Seller == "Czechia"] = "Czech Republic"
dfl$Seller[dfl$Seller == "Macedonia"] = "Macedonia (FYROM)"
dfl$Seller[dfl$Seller == "East Germany (GDR)"] = "German Democratic Republic"
dfl$Seller[dfl$Seller == "Bosnia-Herzegovina"] = "Bosnia and Herzegovina"

dfl$Buyer[dfl$Buyer == "Cote d'Ivoire"] = "Cote dIvoire"
dfl$Buyer[dfl$Buyer == "Czechia"] = "Czech Republic"
dfl$Buyer[dfl$Buyer == "Macedonia"] = "Macedonia (FYROM)"
dfl$Buyer[dfl$Buyer == "East Germany (GDR)"] = "German Democratic Republic"
dfl$Buyer[dfl$Buyer == "Bosnia-Herzegovina"] = "Bosnia and Herzegovina"


# Change all transactions with the Soviet Union to Russia 
dfl$Seller[dfl$Seller == "Soviet Union"] = "Russia" 
dfl$Buyer[dfl$Buyer == "Soviet Union"] = "Russia" 

dfl$SIPRI.estimate = as.numeric(dfl$SIPRI.estimate)
dfl$TIV.deal.unit = as.numeric(dfl$TIV.deal.unit)
dfl$TIV.delivery.values = as.numeric(dfl$TIV.delivery.values)


# select transfers in our time frame 
dfl = dfl[Order.date >= start & Order.date <= end]


# select countries included in our analysis 
dfl = dfl[Seller %in% country_list$V1 & Buyer %in% country_list$V1]

dfl_summary = dfl[, .(Count = .N, "TIV (Deal Unit)" = sum(TIV.deal.unit, na.rm = TRUE)), by=c("Armament.category")]
dfl_summary = dfl_summary[order(-Count)]
setnames(dfl_summary, "Armament.category", "Armament Category")


dfl_volume = dfl[, .("tiv_deal_unit" = sum(TIV.deal.unit, na.rm = TRUE)), by=c("Order.date")]
dfl_volume = dfl_volume[order(Order.date)]


## Output

sink("figures/table_overview_sipri.txt")

kbl(dfl_summary, booktabs = T,  format = "latex", 
    digits = 0,  escape = F, linesep = "",
    caption = "Major Conventional Weapons Transfers (1995-2018)", label = "overview_sipri") %>%
  kable_styling(full_width = FALSE) %>%
  kable_styling(latex_options = "HOLD_position")

sink()
