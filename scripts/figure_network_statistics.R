library(scales)
library(network)
library(networkDynamic)
library(sna)
library(tsna)


rm(list = ls(all.names = TRUE))


# load data
load("data/out/EX.RData")
load("data/out/colony.RData")
load("data/out/country_list.RData")
load("data/out/sipri_tiv.RData")

baci_aggregated = readRDS("data/out/baci_aggregated.rds")

period = 1995:2018
networks_yearly_arms = list()
networks_yearly_trade = list()

# some data frames are in the time period 1950:2018 
# trade data is only available for 1995:2018


for(i in 1:length(period)){
  j = i + 1994 - 1949
  ind = (EX[, j] == 1) 
  
  net = network(sipri_tiv[[j]][ind, ind])
  set.edge.value(net, "TIV (in thousands)", sipri_tiv[[j]][ind, ind])
  set.vertex.attribute(net, "index", country_list$ID[ind])
  set.vertex.attribute(net, "iso3", country_list$iso3[ind])
  networks_yearly_arms[[i]] = net
  
  net = network(baci_aggregated[[i]][ind, ind])
  set.edge.value(net, "Trade Flow", baci_aggregated[[i]][ind, ind])
  set.vertex.attribute(net, "index", country_list$ID[ind])
  set.vertex.attribute(net, "iso3", country_list$iso3[ind])
  networks_yearly_trade[[i]] = net
  
}


# plot networks
# needs improvements, take in / outdegree for nodal size
# tune colors and labels, e.g. with ggraph

pdf("figures/network_plots.pdf", paper = "a4")
par(mfrow = c(2,2))

for (i in c(1, 6, 11, 16)) {
  plot(networks_yearly_arms[[i]],
       displayisolates = FALSE,  
       label.pos=5,  label.cex = 0.6, 
       vertex.cex = 2.2, pad = 0.5,
       label = networks_yearly_arms[[i]] %v% "iso3", 
       main = paste("International Arms Transfers in", i + 1994))
}

for (i in c(1, 6, 11, 16)) {
  plot(networks_yearly_trade[[i]],
       displayisolates = FALSE,  
       label.pos=5,  label.cex = 0.6, 
       vertex.cex = 2.2, pad = 0.5,
       label = networks_yearly_trade[[i]] %v% "iso3", 
       main = paste("International Trade in", i + 1994))
}
dev.off()


# create networkDynamic object
# By default, each observation is assumed to span a unit interval, (so the 1st goes from 0 to 1, 2nd from 1-2, etc).
networks_dynamic_arms = networkDynamic(network.list = networks_yearly_arms, vertex.pid = 'index', start = 1995, end = 2018)
# networks_dynamic_trade = networkDynamic(network.list = networks_yearly_trade, vertex.pid = 'index', start = 1995, end = 2018)



# retrieve network stats and plot them using tsna 
pdf("figures/network_statistics.pdf", paper = "a4")
par(mfrow = c(3, 2))

#plot(tEdgeDissolution(networks_dynamic_arms, start = 1995, end = 2018), ylab = "Series 1", main="Edge dissolution counts")
#lines(tEdgeDissolution(networks_dynamic_trade, start = 1995, end = 2018))

#plot(tEdgeFormation(networks_dynamic_arms, start = 1995, end = 2018), ylab = "Series 1", main="Edge formation counts")
#lines(tEdgeFormation(networks_dynamic_trade, start = 1995, end = 2018))

plot(tErgmStats(networks_dynamic_arms, '~edges', start = 1995, end = 2018), main = "Arms Transfers: Edge Count")
plot(1995:2018, unlist(lapply(networks_yearly_trade, FUN = function(x){network.edgecount(x)})), type = "l", xlab = "Time", ylab = "Series 1", main = "Conventional Trade: Edge Count")

plot(tSnaStats(networks_dynamic_arms, 'gden', start = 1995, end = 2018), main = "Arms Transfers: Network Density")
plot(1995:2018, unlist(lapply(networks_yearly_trade, FUN = function(x){network.density(x)})), type = "l", xlab = "Time", ylab = "Series 1", main = "Conventional Trade: Network Density")

plot(tSnaStats(networks_dynamic_arms, 'grecip', measure = "dyadic.nonnull", start = 1995, end = 2018), main = "Arms Transfers: Network Reciprocity")
plot(1995:2018, unlist(lapply(networks_yearly_trade, FUN = function(x){grecip(x, measure = "dyadic.nonnull")})), type = "l", xlab = "Time", ylab = "Series 1", main = "Conventional Trade: Network Reciprocity")

plot(tSnaStats(networks_dynamic_arms, 'gtrans', start = 1995, end = 2018), main = "Arms Transfers: Network Transitivity")
plot(1995:2018, unlist(lapply(networks_yearly_trade, FUN = function(x){gtrans(x)})), type = "l", xlab = "Time", ylab = "Series 1", main = "Conventional Trade: Network Transitivity")

plot(tSnaStats(networks_dynamic_arms, "centralization", FUN = "betweenness", start = 1995, end = 2018), main = "Arms Transfers: Betweenness-centrality")
plot(1995:2018, unlist(lapply(networks_yearly_trade, FUN = function(x){centralization(x, FUN = "betweenness")})), type = "l", xlab = "Time", ylab = "Series 1", main = "Conventional Trade: Betweenness-centrality")

plot(tSnaStats(networks_dynamic_arms, "centralization", FUN = "degree", cmode = "indegree", start = 1995, end = 2018), main = "Arms Transfers: In-centrality")
plot(1995:2018, unlist(lapply(networks_yearly_trade, FUN = function(x){centralization(x, FUN = "degree", cmode = "indegree")})), type = "l", xlab = "Time", ylab = "Series 1", main = "Conventional Trade: In-centrality")

plot(tSnaStats(networks_dynamic_arms, "centralization", FUN = "degree", cmode = "outdegree", start = 1995, end = 2018), main = "Arms Transfers: Out-centrality")
plot(1995:2018, unlist(lapply(networks_yearly_trade, FUN = function(x){centralization(x, FUN = "degree", cmode = "outdegree")})), type = "l", xlab = "Time", ylab = "Series 1", main = "Conventional Trade: Out-centrality")

plot(tErgmStats(networks_dynamic_arms, '~idegree(d=1)', start = 1995, end = 2018), main = "Single Supplier")
plot(1995:2018, unlist(lapply(networks_yearly_arms, FUN = function(x){network.size(x)})),  type = "l", xlab = "Time", ylab = "Series 1", main="Network Size")
plot(1995:2018, unlist(lapply(networks_yearly_arms, FUN = function(x){sum(get.edge.value(x,"TIV (in thousands)"))})),  type = "l", xlab = "Time", ylab = "Series 1", main="Sum TIV (in thousands)")
plot(1995:2018, unlist(lapply(networks_yearly_trade, FUN = function(x){sum(get.edge.value(x,"Trade Flow"))})),  type = "l", xlab = "Time", ylab = "Series 1", main="Sum of Conventional Trade Flow")



# In- and Outdegree Distributions
# to fix: degree = 1 should be 0, log transform issues. Simply new labels?

in_degrees1 = unlist(lapply(networks_yearly_arms[(1995:2005)-1994], FUN = function(x){degree(x,cmode="indegree")}))
in_degrees2 = unlist(lapply(networks_yearly_trade[(1995:2005)-1994], FUN = function(x){degree(x,cmode="indegree")}))
tmp1 = table(in_degrees1) / length(in_degrees1)
tmp2 = table(in_degrees2) / length(in_degrees2)
dat = list(as.matrix(tmp1), as.matrix(tmp2))

plot(unlist(dat),type="n", log = "x", xlim=c(1,max(sapply(dat,length))), ylab = "Proportion", xlab = "Degree", main = "In-Degree Distribution 1995-2005 (log-scale)")
mapply(lines, dat, col=seq_along(dat), type = "b", pch = 20)
legend("topright", c("Arms", "Trade"), lty=1, col=seq_along(dat))


out_degrees1 = unlist(lapply(networks_yearly_arms[(1995:2005)-1994], FUN = function(x){degree(x,cmode="outdegree")}))
out_degrees2 = unlist(lapply(networks_yearly_trade[(1995:2005)-1994], FUN = function(x){degree(x,cmode="outdegree")}))
tmp1 = table(out_degrees1) / length(out_degrees1)
tmp2 = table(out_degrees2) / length(out_degrees2)
dat = list(as.matrix(tmp1), as.matrix(tmp2))

plot(unlist(dat), type="n", log = "x", xlim=c(1,max(sapply(dat,length))), ylab = "Proportion", xlab = "Degree", main = "Out-Degree Distribution 1995-2005 (log-scale)")
mapply(lines,dat, col=seq_along(dat), type = "b", pch = 20)
legend("topright", c("Arms", "Trade"), lty=1, col=seq_along(dat))


in_degrees1 = unlist(lapply(networks_yearly_arms[(2005:2018)-1994], FUN = function(x){degree(x,cmode="indegree")}))
in_degrees2 = unlist(lapply(networks_yearly_trade[(2005:2018)-1994], FUN = function(x){degree(x,cmode="indegree")}))
tmp1 = table(in_degrees1) / length(in_degrees1)
tmp2 = table(in_degrees2) / length(in_degrees2)
dat = list(as.matrix(tmp1), as.matrix(tmp2))

plot(unlist(dat), type="n", log = "x", xlim=c(1,max(sapply(dat,length))), ylab = "Proportion", xlab = "Degree", main = "In-Degree Distribution 2005-2018 (log-scale)")
mapply(lines,dat, col=seq_along(dat), type = "b", pch = 20)
legend("topright", c("Arms", "Trade"), lty=1, col=seq_along(dat))


out_degrees1 = unlist(lapply(networks_yearly_arms[(2005:2018)-1994], FUN = function(x){degree(x,cmode="outdegree")}))
out_degrees2 = unlist(lapply(networks_yearly_trade[(2005:2018)-1994], FUN = function(x){degree(x,cmode="outdegree")}))
tmp1 = table(out_degrees1) / length(out_degrees1)
tmp2 = table(out_degrees2) / length(out_degrees2)
dat = list(as.matrix(tmp1), as.matrix(tmp2))

plot(unlist(dat), type="n", log = "x", xlim=c(1,max(sapply(dat,length))), ylab = "Proportion", xlab = "Degree", main = "Out-Degree Distribution 2005-2018 (log-scale)")
mapply(lines,dat, col=seq_along(dat), type = "b", pch = 20)
legend("topright", c("Arms", "Trade"), lty=1, col=seq_along(dat))

dev.off()
