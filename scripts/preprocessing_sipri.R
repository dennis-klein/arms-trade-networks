# clean the SIPRI data from Cornelius Fritz
# code is partly taken from Cornelius Fritz

# run to create usable dataset on ARMS DEALS and TOTAL ARMS TRADE PER YEAR
# data/processed/sipri_arms_deals_1950_2019.rds
# data/processed/sipri_arms_trade_1990_2019.rds


library(readr)
library(dplyr)
library(tidyr)
library(purr)
library(stringr)


paths <- c(
  "data/raw/SIPRI/1950_1970.txt",
  "data/raw/SIPRI/1971_1980.txt",
  "data/raw/SIPRI/1981_1990.txt",
  "data/raw/SIPRI/1991_2000.txt",
  "data/raw/SIPRI/2001_2019.txt"
)

# read in raw SIPRI data
deals_raw <- read_delim(paths, delim = ";", skip = 5)
deals <- deals_raw
names(deals) <- str_replace_all(names(deals), " ", ".") 

# converting variables to correct type
deals <- deals %>%
  mutate(
    `Deal.ID` = as.integer(`Deal.ID`),
    Seller = factor(Seller),
    Buyer = factor(Buyer),
    `Armament.category` = factor(`Armament.category`),
    Designation = factor(Designation),
    `Order.date` = as.integer(`Order.date`),
    `Order.date.is.estimate` = if_else(`Order.date.is.estimate`=="Yes", 1, 0) %>% 
      as.integer(),
    `Numbers.delivered` = as.integer(`Numbers.delivered`),
    `Numbers.delivered is estimate` = if_else(`Numbers.delivered.is.estimate`=="Yes", 1, 0) %>% 
      as.integer(),
    Status = factor(Status),
    `Local.production` = if_else(`Local.production`=="Yes", 1, 0) %>% 
      as.integer()
  ) %>% 
  rename(Order.year = Order.date)

# NOTE: warning message, but data seems fine

# save cleaned deal data
saveRDS(deals, "data/out/sipri_arms_deals_1950_2019.rds")

### build panel (seller, buyer, year)

# restrict to year 1990 to 2019
deals <- deals %>% 
  filter(between(`Order.year`, 1990, 2019))

# group to panel (seller, buyer, year)
# total trade, number of orders
arms <- deals %>% 
  group_by(Seller, Buyer, `Order.year`) %>% 
  summarise(
    TIV.total.delivered = sum(`TIV.delivery.values`),
    num.orders = n()
  ) %>% 
  arrange(desc(TIV.total.delivered)) %>% 
  ungroup()

# save data
saveRDS(arms, "data/out/sipri_arms_trade_1990_2019.rds")


### create networks matrices for MCW trades
country.list <- union(arms$Seller, arms$Buyer) %>% sort()

net.frame <- arms %>%
  select(Seller, Buyer, Order.year, TIV.total.delivered) %>% 
  mutate(
    period = Order.year - min(Order.year) + 1,
    seller.id = match(arms$Seller, country.list),
    buyer.id = match(arms$Buyer, country.list)
  )

years <- unique(arms$Order.year) %>% sort()
T <- length(years)
N <- length(country.list)

# create continuous networks
net.MCW.cont <- array(rep(0.0, T*N*N), dim = c(T, N, N))
for (t in 1:T) {
  # frame for year t
  net.frame.t <- net.frame %>% 
    filter(period == t)
  # fill in year t values
  net.MCW.cont[t, net.frame.t$seller.id, net.frame.t$buyer.id] <- net.frame.t$TIV.total.delivered
}

# create discrete networks
net.MCW.discr <- ifelse(net.MCW.cont > 0, 1, 0) %>% as.integer()


# set attributes
attr(net.MCW.cont, "years") <- years
attr(net.MCW.discr, "years") <- years
attr(net.MCW.cont, "country.list") <- country.list
attr(net.MCW.discr, "country.list") <- country.list

# save networks
saveRDS(net.MCW.cont, "data/out/MCW_trade_net_cont.rds")
saveRDS(net.MCW.discr, "data/out/MCW_trade_net_discr.rds")



# BACKUP
# expand panel
# index <- arms %>% expand(Seller, Buyer, `Order date`)
