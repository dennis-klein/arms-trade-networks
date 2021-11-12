# Data Preprocessing: CEPII Tradehist

# Sources

# CEPII TRADEHIST <http://www.cepii.fr/cepii/en/bdd_modele/presentation.asp?id=32>
#   Bilateral Trade and Bilateral Tariffs
# Exchange Rates <https://www.measuringworth.com/datasets/exchangeglobal/>
#   CPI US <https://www.statbureau.org/en/united-states/cpi-u>
# Replication files <https://www.cambridge.org/core/journals/network-science/article/separable-and-semiparametric-networkbased-counting-processes-applied-to-the-international-combat-aircraft-trades/0D57EC7B7E1775B0BEF72BDE101E507F#article>


library(data.table)
library(tidyverse)
library(readxl)


# load countrylist
load("data/out/country_list.RData")


# create tradematrix
cepii_trade <- list()


# load CEPII TRADEHIST data
TRADHIST_BITRADE_BITARIFF_1 <- read_excel("data/raw/CEPII TRADEHIST/TRADHIST_BITRADE_BITARIFF_1.xlsx")
TRADHIST_BITRADE_BITARIFF_2 <- read_excel("data/raw/CEPII TRADEHIST/TRADHIST_BITRADE_BITARIFF_2.xlsx")
TRADHIST_BITRADE_BITARIFF_3 <- read_excel("data/raw/CEPII TRADEHIST/TRADHIST_BITRADE_BITARIFF_3.xlsx")

tradehist <- rbindlist(list(
  TRADHIST_BITRADE_BITARIFF_1,
  TRADHIST_BITRADE_BITARIFF_2,
  TRADHIST_BITRADE_BITARIFF_3
), fill = TRUE)

rm(
  TRADHIST_BITRADE_BITARIFF_1,
  TRADHIST_BITRADE_BITARIFF_2,
  TRADHIST_BITRADE_BITARIFF_3
)


# currency is British Pound Sterling. 
# nominal trade flows, rebase to constant 2005$ USD.
exchange_usd_gbp <- read.csv(file = "data/raw/EXCHANGEGLOBAL_1945-2020.csv", header = FALSE) %>%
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
  mutate(cpi_factor = cpi / as.numeric(filter(., year == 2005)[,2]))


# merge data, filter for out time-frame. 
# do conversions and replace flow_0 values. 
tradehist %>%
  filter(year %in% 1950:2018) %>%
  mutate(FLOW = replace(FLOW, FLOW_0 == 0, 0)) %>%
  left_join(exchange_usd_gbp, by = c("year" = "year"), na_matches = "never") %>%
  left_join(cpi_us, by = c("year" = "year"), na_matches = "never") %>%
  mutate(FLOW = (FLOW / GBP_USD) / cpi_factor) %>% 
  select(iso_o, iso_d, year, FLOW) -> tradehist


# Match the country names with the ids in country_list
tradehist$from_id = match(tradehist$iso_o, country_list$cepii)
tradehist$to_id = match(tradehist$iso_d, country_list$cepii)


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
tradehist <- tradehist[!(is.na(tradehist$from_id) | is.na(tradehist$to_id)),]


# create a list of 69 adjacency matrices, one for each year from 1950:2018
for (i in 1:69){
  cepii_trade[[i]]<- matrix(NA,257,257)
  colnames(cepii_trade[[i]])<-country_list$V1
  rownames(cepii_trade[[i]])<-country_list$V1
}

for (i in 1950:2018){
  tmp_tradehist = tradehist[tradehist$year == i, ]
  cepii_trade[[i - 1949]][cbind(tmp_tradehist$from_id, tmp_tradehist$to_id)] <- tmp_tradehist$FLOW
}


# what is the distribution of bilateral trade flows normed by country GDP? 
# select appropriate cut-off to binarize the network - Thurner 2019?

## to do ##


# save network list 
save(cepii_trade, file = "data/out/cepii_trade.RData")
