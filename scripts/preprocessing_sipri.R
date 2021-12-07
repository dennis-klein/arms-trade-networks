# Preprocessing of the SIPRI data

# Large Parts from the following code are taken from the 
# replication files of Fritz, C et al. 2021

# run to create usable dataset on ARMS DEALS and TOTAL ARMS TRADE PER YEAR
# data/processed/sipri_arms_deals_1950_2019.rds
# data/processed/sipri_arms_trade_1990_2019.rds


library(readr)
library(data.table)
library(dplyr)
library(tidyr)


rm(list = ls(all.names = TRUE))

source("utils/utils.R")
path = data_path.get()

# load countrylist
load(file.path(path, "out/country_list.RData"))

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


# should we select for different arms categories?
unique(dfl$Armament.category)
tiv_dfl = dfl

tiv_dfl = aggregate(tiv_dfl$TIV.delivery.values, 
                    by = list(paste(tiv_dfl$Seller, tiv_dfl$Buyer, tiv_dfl$Order.date, sep = "_")), 
                    sum) # Aggregate all TIV values by country combination and order year 

tmp = unlist(strsplit(tiv_dfl$Group.1,split = "_")) # 1,4,7 ... are sender 2,5,8 ... are recevier 3,6,9 ... is the year  

tiv_dfl$from = tmp[seq(1, to = length(tmp), by = 3)]
tiv_dfl$to = tmp[seq(2, to = length(tmp), by = 3)]
tiv_dfl$year = as.numeric(tmp[seq(3, to = length(tmp), by = 3)])
tiv_dfl$from_id = match(tiv_dfl$from, country_list$SIPRI) # Match the country names with the ids in country_list
tiv_dfl$to_id = match(tiv_dfl$to ,country_list$SIPRI)


# Try matching the unmatched names with the alternative name also given in country_list
tiv_dfl$from_id[is.na(tiv_dfl$from_id)] = match(tiv_dfl$from[is.na(tiv_dfl$from_id)],country_list$V1)
tiv_dfl$to_id[is.na(tiv_dfl$to_id)] = match(tiv_dfl$to[is.na(tiv_dfl$to_id)],country_list$V1)

unique(tiv_dfl$to[is.na(tiv_dfl$to_id)]) # Rest of the unmatched countries are not independent states and thus not part of the analysis 
unique(tiv_dfl$from[is.na(tiv_dfl$from_id)]) # Rest of the unmatched countries are not independent states and thus not part of the analysis 


# Delete all observations where the id is not known 
tiv_dfl = tiv_dfl[!(is.na(tiv_dfl$from_id) | is.na(tiv_dfl$to_id)),]
sipri_tiv<- list()


# create a list of 69 adjacency matrices, one for each year from 1950:2017. 
for (i in 1:69){
  sipri_tiv[[i]]<- matrix(0,257,257)
  colnames(sipri_tiv[[i]])<-country_list$V1
  rownames(sipri_tiv[[i]])<-country_list$V1
}


# fill the matrix per year 
for(i in 1950:2018){
  tmp <- tiv_dfl[tiv_dfl$year == i,]
  sipri_tiv[[i- 1949]][cbind(tmp$from_id, tmp$to_id)] = tmp$x
}


# save
saveRDS(sipri_tiv, file = file.path(path, "out/sipri_tiv.rds"))

