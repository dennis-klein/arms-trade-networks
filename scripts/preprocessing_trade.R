# Data Preprocessing 

# Sources

# CEPII TRADEHIST <http://www.cepii.fr/cepii/en/bdd_modele/presentation.asp?id=32>
#   Bilateral Trade and Bilateral Tariffs
# Exchange Rates <https://www.measuringworth.com/datasets/exchangeglobal/>
#   CPI US <https://www.statbureau.org/en/united-states/cpi-u>
# Replication files <https://www.cambridge.org/core/journals/network-science/article/separable-and-semiparametric-networkbased-counting-processes-applied-to-the-international-combat-aircraft-trades/0D57EC7B7E1775B0BEF72BDE101E507F#article>
# Peacesciencer <http://svmiller.com/peacesciencer/index.html>

library(countrycode)
library(readxl)
library(tidyverse)
library(peacesciencer)
library(statnet)


# define time frame of analysis
time_frame <- 1991:2014


# create state-year nodelist
create_stateyears(system = "cow", mry = FALSE) %>%
  add_democracy() %>%
  add_nmc() %>%
  add_sdp_gdp() %>%
  filter(year %in% time_frame) -> nodes


# create directed state-state-year nodelist
create_dyadyears(system = "cow", mry = FALSE) %>%
  add_capital_distance %>%
  add_contiguity() %>%
  add_minimum_distance() %>%
  add_atop_alliance() %>%
  filter(year %in% time_frame) -> dyads


# load CEPII TRADEHIST data
TRADHIST_BITRADE_BITARIFF_1 <- read_excel("data/raw/CEPII TRADEHIST/TRADHIST_BITRADE_BITARIFF_1.xlsx")
TRADHIST_BITRADE_BITARIFF_2 <- read_excel("data/raw/CEPII TRADEHIST/TRADHIST_BITRADE_BITARIFF_2.xlsx")
TRADHIST_BITRADE_BITARIFF_3 <- read_excel("data/raw/CEPII TRADEHIST/TRADHIST_BITRADE_BITARIFF_3.xlsx")


# currency is British Pound Sterling. 
# nominal trade flows, rebase to constant 2011$ USD.
exchange_usd_gbp <- read.csv(file = "data/raw/EXCHANGEGLOBAL_1945-2020.csv", 
                             header = FALSE) %>%
  rename(year = V1, GBP_USD = V2) %>%
  mutate(GBP_USD = as.numeric(gsub(" British Pound", "", GBP_USD)))

cpi_us <- read.csv(file = "data/raw/united-states.index.cpi-u (statbureau.org).csv") %>%
  rename(year = ï..Year) %>%
  filter(year %in% c(1950:2020)) %>%
  pivot_longer(cols = c("January", "February", "March", "April", "May", "June",
                        "July", "August", "September", "October", "November",
                        "December"),
               names_to = "month", values_to = "cpi") %>% 
  group_by(year) %>%
  summarize(cpi = mean(cpi)) %>%
  mutate(cpi_factor = cpi / as.numeric(filter(., year == 2011)[,2]))


# merge data, filter for out time-frame. 
# do conversions and replace flow_0 values. 
TRADHIST_BITRADE_BITARIFF_1 %>%
  bind_rows(TRADHIST_BITRADE_BITARIFF_2) %>%
  bind_rows(TRADHIST_BITRADE_BITARIFF_3) %>%
  filter(year %in% time_frame) %>%
  mutate(FLOW = replace(FLOW, FLOW_0 == 0, 0)) %>%
  left_join(exchange_usd_gbp, by = c("year" = "year"), na_matches = "never") %>%
  left_join(cpi_us, by = c("year" = "year"), na_matches = "never") %>%
  mutate(FLOW = (FLOW / GBP_USD) / cpi_factor) %>% 
  select(iso_o, iso_d, year, FLOW) -> bilateral_trade


# include COW codes using countrycode 
dict <- c(
  "ABW" = NA, 
  "ANT" = NA, 
  "BMU" = NA, 
  "CUW" = NA,
  "CZSK" = 315, 
  "EDEU" = 265, 
  "FLK" = NA, 
  "FRO" = NA, 
  "GIB" = NA, 
  "GLP" = NA, 
  "GRL" = NA, 
  "GUF" = NA, 
  "HKG" = NA, 
  "MAC" = NA,
  "MTQ" = NA, 
  "NCL" = NA, 
  "PYF" = NA,
  "REU" = NA,
  "ROM" = 360, 
  "SHN" = NA,
  "SPM" = NA,
  "USSR" = 365, 
  "WDEU" = 260, 
  "YAR" = NA,
  "YMD" = NA, 
  "YUG"= 345
)

bilateral_trade %>%
  mutate(ccode1 = countrycode(iso_o, "iso3c", "cown", custom_match = dict)) %>%
  mutate(ccode2 = countrycode(iso_d, "iso3c", "cown", custom_match = dict)) -> bilateral_trade


# merge to other covariates
dyads %>% left_join(
    bilateral_trade,
    by = c("ccode1" = "ccode1", "ccode2" = "ccode2", "year" = "year"),
    na_matches = "never"
  ) -> dyads


# what is the distribution of bilateral trade flows normed by country GDP? Thurner 2019?
# select appropiate cut-off to binarize the network

## to do ##


# create network object
network_trade <- vector("list", length(time_frame))
network_trade_binary <- vector("list", length(time_frame))

for (num in 1:length(time_frame)) {
  
  filter(dyads, year == time_frame[num])[, c("ccode1", "ccode2", "FLOW")] %>%
    pivot_wider(names_from = ccode2, names_sort = TRUE, values_from = FLOW) %>%
    column_to_rownames(var = "ccode1") %>%
    as.matrix() %>%
    as.network() -> network_trade[[num]]
  
}


# create dynamic network object
# issue of states entering and leaving the network

## to do ##

# trade_dynamic_network <- networkDynamic(network.list=network_trade)


# save network list 
saveRDS(network_trade, "data/out/network_trade.rds")
