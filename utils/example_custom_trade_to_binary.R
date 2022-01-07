library(network)
library(sna)

rm(list = ls(all.names = TRUE))
set.seed(1234)

source("utils/utils.R")
source("utils/custom_trade_to_binary.R")
path = data_path.get()


# load data
load(file.path(path, "out/EX.RData"))
trade = readRDS(file.path(path, "out/baci_aggregated.rds"))


# select existing countries (one way), include all which exist continously 
ind <- rowSums(EX[, (1995:2018)-1949]) == length(1995:2018)
for (i in 1:length(trade)) {
  trade[[i]] <- trade[[i]][ind, ind]
}


# apply thresholding to list
trade_binary = lapply(trade, FUN=custom_trade_to_binary, type = "B", threshold = 0.01)
trade_net = lapply(trade_binary, FUN=network, directed = T)

network.density(trade_net[[1]])


