# Data Sources 
## Folder: data/raw/
BACI: International Trade Database at the Product-Level. The 1995-2019 Version
http://www.cepii.fr/cepii/en/bdd_modele/presentation.asp?id=37


International Trade and Production Database for Estimation (ITPD-E)
17 years, 2000-2016, 243 countries, 170 industries
https://www.usitc.gov/data/gravity/itpde.htm


TRADHIST Bilateral Trade Historical Series: New Dataset 1827-2014
http://www.cepii.fr/CEPII/en/bdd_modele/presentation.asp?id=32


API_NY.GDP.MKTP.KD_DS2_en_csv_v2_3358328 - World Bank GDP (current US$)
https://data.worldbank.org/indicator/NY.GDP.MKTP.CD


INTRA-STATE_State_participants v5.1 CSV - Intra-State War Data (v5.1)
https://correlatesofwar.org/data-sets/COW-war/intra-state-wars-v5-1.zip/view


NMC-60-abridged - NMC: National military capabilities of states from 1816-2016
https://correlatesofwar.org/news/nmc-6-0-data-available


mpd2020 - Maddison Historical Statistics
https://www.rug.nl/ggdc/historicaldevelopment/maddison/?lang=en


p5v2018 - Polity5d Polity-Case Format, 1800-2018
https://www.systemicpeace.org/inscrdata.html


SIPRI Arms Transfers Data
Consumer Price Index CPI US
GBP / USD Exchange Rate


## Folder: data/out/
Note: different time-frames apply!

atop_alliance.RData: List of yearly (69) matrices (257 x 257) indicating if in 
year t, country i had an alliance with country j. (Check: defense alliance?)  

baci_aggregated.rds: List of yearly (25, 1995 to 2018) matrices (257 x 257); 
numerical, aggregated export flows. Source: cepii baci. Rebased to constant 2010 USD.

cdist.RData: Matrix with dimensions 257 x 257, numerical for capital distance 
between country i and j.

colony.RData: Matrix of dimensions 257 x 257, dummy 1 for colonial dependency

conflict_intrastate.rds: Matrix of dimensions 257 x 69, dummy 1 if intrastate 
conflict present in country i and year j.

country_list.RData: data.frame with 257 observations: countries with respective 
identification.

EX.RData: Matrix of dimensions 257 x 69 with dummys 1 if country i existed in year j

gdp.rds: Matrix of dimensions 257 x 69 with numerical constant USD GDP from the 
World Bank, Augumented with data from Maddison (check correct denomination!).

itpd_mining-energy.rds: List of yearly (17, 2000 to 2016) matrices (257 x 257); 
numerical, aggregated export flows in category mining/energy. 
Not rebased because we don't need it any more.

itpd_aggregated.rds: List of yearly (17, 2000 to 2016) matrices (257 x 257); 
numerical, aggregated export flows. Rebased to constant 2010 USD.

milit_exp.RData: Matrix of dimensions 257 x 69, numerical for military expenditure. Source: NMC

nmc_cinc.rds: Matrix of dimensions 257 x 69, numerical for nmc index (percentage), imputed with 0.

polity.RData: Matrix with dimensions 257 x 69, integer from -10 to 10 indicating 
country system, time series harmonized. Plausible imputed values.

pop.RData: Matrix with dimensions 257 x 69, numerical for total population. Source: ?

sipri_tiv.rds: List of yearly (69, 1950 to 2018) matrices (257 x 257); numerical, 
aggregated TIV Values for Order Dates (!) indicating the ordered transfer in year t, from country i to j.

tradhist_aggregated.rds.rds: List of yearly (69, 1950 to 2018) matrices (257 x 257); 
numerical, aggregated export flows. Source: CEPII TRADHIST. Constant 2010 USD

RData files taken from the replication files of https://www.cambridge.org/core/journals/network-science/article/separable-and-semiparametric-networkbased-counting-processes-applied-to-the-international-combat-aircraft-trades/0D57EC7B7E1775B0BEF72BDE101E507F


