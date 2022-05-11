library(testthat)
library(tidyverse)

source("scripts/SAOM/trade_import_gdp_relevance.R")
source("utils/utils.R")


dpath <- data_path.get()

gdp <- readRDS(file.path(dpath, "out/gdp.rds"))
baci <- readRDS(file.path(dpath, "out/baci_aggregated.rds"))
load(file.path(dpath, "out/EX.RData"))

names(baci) <- paste(1995:2019)

test1 <- trade_import_gdp_relevance(baci[["2000"]], gdp[,"2000"], threshold = NULL)
test2 <- trade_import_gdp_relevance(baci[["2000"]], gdp[,"2000"])
test3 <- trade_import_gdp_relevance(baci[["2000"]], gdp[,"2000"], threshold = 0.001)

t1 <- baci[["2000"]]["France", "Germany"]
t2 <- gdp["Germany","2000"]
t3 <- t1/t2
t4 <- 1*(t3 >= 0.01)

expect_equal(test1["Germany", "France"], t3)
expect_equal(test2["Germany", "France"], t4)
expect_equal(test2["United States", "Liechtenstein"], 0)
expect_equal(test2["Chad", "China"], 0)
expect_equal(test2["United States", "China"], 1)
expect_equal(test2["China", "United States"], 1)



sum(test2, na.rm = T)/length(test2)
mean(test2, na.rm = T)
mean(test3, na.rm = T)
test1["China", "United States"]
mean(test1["China",])
mean(test1["United States",])
mean(test1["Germany",])
mean(test1[,"Germany"], na.rm = T)
mean(test1[,"United States"], na.rm = T)
mean(test1[,"China"], na.rm = T)

x <- names(baci)
y <- dimnames(gdp)[[2]]
years_both <- intersect(x, y)
ex <- rowSums(EX[, paste(years_both)]) == length(years_both)
baci1 <- baci[years_both]
baci1 <- lapply(baci1, function(x) x[ex, ex])
gdp1 <- gdp[ex, years_both]

testing_thresholds <- seq(0.0005, 0.02, length.out = 100)

# identify adequate threshold
dat <- expand.grid(testing_thresholds, as.integer(years_both))
names(dat) <- c("tresh", "year")
dat <- tibble(dat)
head(dat)

net_density <- function(year, threshold) {
  trd_bin <- trade_import_gdp_relevance(trd = baci1[[paste(year)]], gdp = gdp1[,paste(year)],
                                        threshold = threshold)
  return(mean(trd_bin, na.rm = T))
}

dat$density <- mapply(function(year, thresh) net_density(year, thresh), dat$year, dat$tresh)

tmp <- trade_import_gdp_relevance(trd = baci1[["2000"]], gdp = gdp1[,"2000"],
                           threshold = 0.0035)
mean(tmp, na.rm = T)
sum(tmp, na.rm = T)/length(tmp)


# Result:
# We aim for about 10% density in the network and
# concluded that this happens at a threshold of around 0.004