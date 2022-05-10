#------------------------------------------------------------------------------#
# Figure: Summary Statistics for Appendix
#------------------------------------------------------------------------------#

library(network)
library(dplyr)
library(reshape2)
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
load(file.path(path, "out/cdist.RData"))
load(file.path(path, "out/atop_alliance.RData"))
load(file.path(path, "out/milit_exp.RData"))
nmc_cinc <- readRDS(file.path(path, "out/nmc_cinc.rds"))
polity <- readRDS(file.path(path, "out/polity.rds"))
gdp <- readRDS(file.path(path, "out/gdp.rds"))
gdppc <- readRDS(file.path(path, "out/gdppc.rds"))
sipri_tiv <- readRDS(file.path(path, "out/sipri_tiv.rds"))
trade <- readRDS(file.path(path, "out/baci_aggregated.rds"))

# select years of analysis
year <- 2003
start <- 1995
end <- 2017
period <- start:end

# selection of countries included: present over the complete time horizon (cf. soam)
ind <- rowSums(EX[, (start:end) - 1949]) == length(start:end)
n <- sum(ind)

# some data transformations
# sipri military expenditure is in constant USD$m
# cinc is index which sums up to 0 -> better interpretability using pct points
milit_exp <- milit_exp * 1000000
nmc_cinc100 <- nmc_cinc * 100

# prepare data 
# set correct indices
i1 = period - 1949 # Cornelius Replication Data 1950:2018
i2 = period - 1994 # CEPII BACI 1995:2018

# summary statistics of exogenous covariates (time-pooled) over complete period
df <- data.frame(
  "Name" = character(),
  "N" = integer(),
  "Mean" = numeric(), "SD" = numeric(),
  "Min" = numeric(), "Max" = numeric(),
  "Median" = numeric()
)



## Time - Varying

# MCW
sipri_tiv_pooled <- list()

for (t in i2) {
  tmp <- sipri_tiv[[t]][ind, ind] / 1000
  diag(tmp) <- NA
  tmp <- melt(tmp)
  sipri_tiv_pooled[[t]] <- tmp
}

sipri_tiv_pooled <- do.call(rbind, sipri_tiv_pooled)

df[1, "Name"] <- "MCW Trade Tie (TIV, '000s)"
df[1, "N"] <- sum(!is.na(sipri_tiv_pooled$value))
df[1, "Mean"] <- mean(sipri_tiv_pooled$value, na.rm = TRUE)
df[1, "SD"] <- sd(sipri_tiv_pooled$value, na.rm = TRUE)
df[1, "Median"] <- median(sipri_tiv_pooled$value, na.rm = TRUE)
df[1, "Min"] <- min(sipri_tiv_pooled$value, na.rm = TRUE)
df[1, "Max"] <- max(sipri_tiv_pooled$value, na.rm = TRUE)


# Transfers
trade_pooled <- list()

for (t in i2) {
  tmp <- trade[[t]][ind, ind] / 1000000
  diag(tmp) <- NA
  tmp <- melt(tmp)
  trade_pooled[[t]] <- tmp
}

trade_pooled <- do.call(rbind, trade_pooled)

df[2, "Name"] <- "Conventional Trade Tie (2010 USD '000 000s)"
df[2, "N"] <- sum(!is.na(trade_pooled$value))
df[2, "Mean"] <- mean(trade_pooled$value, na.rm = TRUE)
df[2, "SD"] <- sd(trade_pooled$value, na.rm = TRUE)
df[2, "Median"] <- median(trade_pooled$value, na.rm = TRUE)
df[2, "Min"] <- min(trade_pooled$value, na.rm = TRUE)
df[2, "Max"] <- max(trade_pooled$value, na.rm = TRUE)


# GDP per capita
tmp <- log(gdp[ind, i1])
df[3, "Name"] <- "GDP (log)"
df[3, "N"] <- sum(ind)*length(period)
df[3, "Mean"] <- mean(tmp)
df[3, "SD"] <- sd(tmp)
df[3, "Median"] <- median(tmp)
df[3, "Min"] <- min(tmp)
df[3, "Max"] <- max(tmp)


# Military Expenditure
tmp <- log(milit_exp[ind, i1])
df[4, "Name"] <- "Military Expenditure (log)"
df[4, "N"] <- sum(ind)*length(period)
df[4, "Mean"] <- mean(tmp)
df[4, "SD"] <- sd(tmp)
df[4, "Median"] <- median(tmp)
df[4, "Min"] <- min(tmp)
df[4, "Max"] <- max(tmp)



# Polity Index
polity_abs_diff_pooled <- list()

for (t in i1) {
  tmp <- polity[ind, t]
  tmp <- abs(outer(tmp, tmp, FUN = "-"))
  
  diag(tmp) <- NA
  
  polity_abs_diff_pooled[[t]] <- melt(tmp)
}

polity_abs_diff_pooled <- do.call(rbind, polity_abs_diff_pooled)

df[5, "Name"] <- "Difference in Polity (abs)"
df[5, "N"] <- sum(!is.na(polity_abs_diff_pooled$value))
df[5, "Mean"] <- mean(polity_abs_diff_pooled$value, na.rm = TRUE)
df[5, "SD"] <- sd(polity_abs_diff_pooled$value, na.rm = TRUE)
df[5, "Median"] <- median(polity_abs_diff_pooled$value, na.rm = TRUE)
df[5, "Min"] <- min(polity_abs_diff_pooled$value, na.rm = TRUE)
df[5, "Max"] <- max(polity_abs_diff_pooled$value, na.rm = TRUE)



# Defense Alliance
atop_defense_pooled = list()

for (t in i1) {
  tmp <- atop_alliance[[t]][ind, ind]
  diag(tmp) <- NA
  tmp <- melt(tmp)
  atop_defense_pooled[[t]] <- tmp
}

atop_defense_pooled <- do.call(rbind, atop_defense_pooled)

df[6, "Name"] <- "Alliance"
df[6, "N"] <- sum(!is.na(atop_defense_pooled$value))
df[6, "Mean"] <- mean(atop_defense_pooled$value, na.rm = TRUE)
df[6, "SD"] <- sd(atop_defense_pooled$value, na.rm = TRUE)
df[6, "Min"] <- min(atop_defense_pooled$value, na.rm = TRUE)
df[6, "Max"] <- max(atop_defense_pooled$value, na.rm = TRUE)


## Time Constant

# Distance
diag(cdist) <- NA
tmp <- log(cdist[ind, ind])
df[7, "Name"] <- "Distance in km (log)"
df[7, "N"] <- sum(!is.na(tmp))
df[7, "Mean"] <- mean(tmp, na.rm = TRUE)
df[7, "SD"] <- sd(tmp, na.rm = TRUE)
df[7, "Median"] <- median(tmp, na.rm = TRUE)
df[7, "Min"] <- min(tmp, na.rm = TRUE)
df[7, "Max"] <- max(tmp, na.rm = TRUE)

df



## Output

sink("figures/table_summary_statistics.txt")

kbl(df[, 1:6], booktabs = T,  format = "latex", 
    digits = 2,  escape = F, linesep = "",
    caption = "Summary statistics of exogenous covariates (1995--2017)", label = "summary_statistics") %>%
  kable_styling(font_size = 11, full_width = FALSE) %>%
  kable_styling(latex_options = "HOLD_position")

sink()
