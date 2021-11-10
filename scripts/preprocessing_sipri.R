# clean the SIPRI data from Cornelius Fritz
# code is partly taken from Cornelius Fritz

# run to create usable dataset on ARMS DEALS and TOTAL ARMS TRADE PER YEAR
# data/processed/sipri_arms_deals_1950_2019.rds
# data/processed/sipri_arms_trade_1990_2019.rds


library(readr)
library(dplyr)
library(tidyr)
library(purr)


paths <- c(
  "data/raw/SIPRI/1950_1970.txt",
  "data/raw/SIPRI/1971_1980.txt",
  "data/raw/SIPRI/1981_1990.txt",
  "data/raw/SIPRI/1991_2000.txt",
  "data/raw/SIPRI/2001_2019.txt"
)

# read in raw SIPRI data
deals_raw <- read_delim(paths, delim = ";", skip = 5)


# converting variables to correct type
deals <- deals_raw %>%
  mutate(
    `Deal ID` = as.integer(`Deal ID`),
    Seller = factor(Seller),
    Buyer = factor(Buyer),
    `Armament category` = factor(`Armament category`),
    Designation = factor(Designation),
    `Order date` = as.integer(`Order date`),
    `Order date is estimate` = if_else(`Order date is estimate`=="Yes", 1, 0) %>% 
      as.integer(),
    `Numbers delivered` = as.integer(`Numbers delivered`),
    `Numbers delivered is estimate` = if_else(`Numbers delivered is estimate`=="Yes", 1, 0) %>% 
      as.integer(),
    Status = factor(Status),
    `Local production` = if_else(`Local production`=="Yes", 1, 0) %>% 
      as.integer()
  )

# NOTE: warning message, but data seems fine

# save cleaned deal data
saveRDS(deals, "data/processed/sipri_arms_deals_1950_2019.rds")

### build panel (seller, buyer, year)

# restrict to year 1990 to 2019
deals <- deals %>% 
  filter(between(`Order date`, 1990, 2019))

# group to panel (seller, buyer, year)
# total trade, number of orders
arms <- deals %>% 
  group_by(Seller, Buyer, `Order date`) %>% 
  summarise(
    TIV_total_delivered = sum(`TIV delivery values`),
    num_orders = n()
  ) %>% 
  arrange(desc(TIV_total_delivered))

# save data
saveRDS(arms, "data/processed/sipri_arms_trade_1990_2019.rds")


# BACKUP
# expand panel
# index <- arms %>% expand(Seller, Buyer, `Order date`)
