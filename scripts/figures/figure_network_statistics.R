library(scales)
library(network)
library(networkDynamic)
library(sna)
library(tsna)


rm(list = ls(all.names = TRUE))

source("utils/utils.R")
path = data_path.get()


# load data
load(file.path(path, "out/EX.RData"))
load(file.path(path, "out/colony.RData"))
load(file.path(path, "out/country_list.RData"))

sipri_tiv = readRDS(file.path(path, "out/sipri_tiv.rds"))
baci_aggregated = readRDS(file.path(path, "out/baci_aggregated.rds"))
itpd_mining = readRDS(file.path(path, "out/itpd_mining-energy.rds"))

period = 1995:2018
networks_yearly_sipri = list()
networks_yearly_cepii = list()
networks_yearly_itpd = list()

# some data frames are in the time period 1950:2018 
# trade data is only available for 1995:2018
# itpd mining is only available for 2000:2016

for(i in 1:length(period)){
  j = i + 1994 - 1949
  ind = (EX[, j] == 1) 
  
  net = network(sipri_tiv[[j]][ind, ind])
  set.edge.value(net, "TIV (in thousands)", sipri_tiv[[j]][ind, ind])
  set.vertex.attribute(net, "index", country_list$ID[ind])
  set.vertex.attribute(net, "iso3", country_list$iso3[ind])
  networks_yearly_sipri[[i]] = net
  
  net = network(baci_aggregated[[i]][ind, ind])
  set.edge.value(net, "Trade Flow", baci_aggregated[[i]][ind, ind])
  set.vertex.attribute(net, "index", country_list$ID[ind])
  set.vertex.attribute(net, "iso3", country_list$iso3[ind])
  networks_yearly_cepii[[i]] = net
  
}

for(i in 1:length(2000:2016)){
  j = i + 1999 - 1949
  ind = (EX[, j] == 1) 
  
  net = network(itpd_mining[[i]][ind, ind])
  set.edge.value(net, "ITPD Mining/Energy", itpd_mining[[i]][ind, ind])
  set.vertex.attribute(net, "index", country_list$ID[ind])
  set.vertex.attribute(net, "iso3", country_list$iso3[ind])
  networks_yearly_itpd[[i]] = net
}

# plot networks
# needs improvements, take in / outdegree for nodal size
# tune colors and labels, e.g. with ggraph


par(mfrow = c(2,2))

for (i in c(1, 6, 11, 16)) {
  plot(networks_yearly_sipri[[i]],
       displayisolates = FALSE,  
       label.pos=5,  label.cex = 0.6, 
       vertex.cex = 2.2, pad = 0.5,
       label = networks_yearly_sipri[[i]] %v% "iso3", 
       main = paste("International Arms Transfers in", i + 1994))
}

for (i in c(1, 6, 11, 16)) {
  plot(networks_yearly_cepii[[i]],
       displayisolates = FALSE,  
       label.pos=5,  label.cex = 0.6, 
       vertex.cex = 2.2, pad = 0.5,
       label = networks_yearly_cepii[[i]] %v% "iso3", 
       main = paste("International Trade in", i + 1994))
}

for (i in c(1, 6, 11)) {
  plot(networks_yearly_itpd[[i]],
       displayisolates = FALSE,  
       label.pos=5,  label.cex = 0.6, 
       vertex.cex = 2.2, pad = 0.5,
       label = networks_yearly_itpd[[i]] %v% "iso3", 
       main = paste("International Trade Energy / Mining in", i + 1994))
}


# create networkDynamic object
# By default, each observation is assumed to span a unit interval, (so the 1st goes from 0 to 1, 2nd from 1-2, etc).
networks_dynamic_arms = networkDynamic(network.list = networks_yearly_sipri, vertex.pid = 'index', start = 1995, end = 2018)
# networks_dynamic_trade = networkDynamic(network.list = networks_yearly_cepii, vertex.pid = 'index', start = 1995, end = 2018)
# networks_dynamic_itpd = networkDynamic(network.list = networks_yearly_itpd, vertex.pid = 'index', start = 2000, end = 2016)


# retrieve network stats and plot them using tsna 
pdf("figures/network_statistics.pdf", paper = "a4")
par(mfrow = c(3, 2))

#plot(tEdgeDissolution(networks_dynamic_arms, start = 1995, end = 2018), ylab = "Series 1", main="Edge dissolution counts")
#lines(tEdgeDissolution(networks_dynamic_trade, start = 1995, end = 2018))

#plot(tEdgeFormation(networks_dynamic_arms, start = 1995, end = 2018), ylab = "Series 1", main="Edge formation counts")
#lines(tEdgeFormation(networks_dynamic_trade, start = 1995, end = 2018))

plot(1995:2018, unlist(lapply(networks_yearly_sipri, FUN = function(x){network.size(x)})),  type = "l", xlab = "Year", ylab = "", main="Network Size")
plot(tErgmStats(networks_dynamic_arms, '~edges', start = 1995, end = 2018), xlab = "Year", ylab = "", main = "SIPRI: Edge Count")
plot(1995:2018, unlist(lapply(networks_yearly_cepii, FUN = function(x){network.edgecount(x)})), type = "l", xlab = "Year", ylab = "", main = "CEPII BACI: Edge Count")
plot(2000:2016, unlist(lapply(networks_yearly_itpd, FUN = function(x){network.edgecount(x)})), type = "l", xlab = "Year", ylab = "", main = "ITPD Mining/Energy: Edge Count")

plot(tErgmStats(networks_dynamic_arms, '~idegree(d=1)', start = 1995, end = 2018), xlab = "Year", ylab = "", main = "SIPRI: Single Supplier")
plot.new()


plot(1995:2018, unlist(lapply(networks_yearly_sipri, FUN = function(x){sum(get.edge.value(x,"TIV (in thousands)"))})),  type = "l", xlab = "Year", ylab = "", main ="SIPRI: TIV (in thousands)")
plot(1995:2018, unlist(lapply(networks_yearly_cepii, FUN = function(x){sum(get.edge.value(x,"Trade Flow"))})),  type = "l", xlab = "Year", ylab = "", main ="CEPII BACI: Flow Value")
plot(2000:2016, unlist(lapply(networks_yearly_itpd, FUN = function(x){sum(get.edge.value(x,"ITPD Mining/Energy"))})),  type = "l", xlab = "Year", ylab = "", main ="ITPD Mining/Energy: Flow Value")

plot(tSnaStats(networks_dynamic_arms, 'gden', start = 1995, end = 2018), xlab = "Year", ylab = "", main = "SIPRI: Network Density")
plot(1995:2018, unlist(lapply(networks_yearly_cepii, FUN = function(x){network.density(x)})), type = "l", xlab = "Year", ylab = "", main = "CEPII BACI: Network Density")
plot(2000:2016, unlist(lapply(networks_yearly_itpd, FUN = function(x){network.density(x)})), type = "l", xlab = "Year", ylab = "", main = "ITPD Mining/Energy: Network Density")



plot(tSnaStats(networks_dynamic_arms, 'grecip', measure = "dyadic.nonnull", start = 1995, end = 2018), xlab = "Year", ylab = "", main = "SIPRI: Network Reciprocity")
plot(1995:2018, unlist(lapply(networks_yearly_cepii, FUN = function(x){grecip(x, measure = "dyadic.nonnull")})), type = "l", xlab = "Year", ylab = "", main = "CEPII BACI: Network Reciprocity")
plot(2000:2016, unlist(lapply(networks_yearly_itpd, FUN = function(x){grecip(x, measure = "dyadic.nonnull")})), type = "l", xlab = "Year", ylab = "", main = "ITPD Mining/Energy: Network Reciprocity")

plot(tSnaStats(networks_dynamic_arms, 'gtrans', start = 1995, end = 2018), xlab = "Year", ylab = "", main = "SIPRI: Network Transitivity")
plot(1995:2018, unlist(lapply(networks_yearly_cepii, FUN = function(x){gtrans(x)})), type = "l", xlab = "Year", ylab = "", main = "CEPII BACI: Network Transitivity")
plot(2000:2016, unlist(lapply(networks_yearly_itpd, FUN = function(x){gtrans(x)})),  type = "l", xlab = "Year", ylab = "", main = "ITPD Mining/Energy: Network Transitivity")



plot(tSnaStats(networks_dynamic_arms, "centralization", FUN = "betweenness", start = 1995, end = 2018), xlab = "Year", ylab = "", main = "SIPRI: Betweenness-centrality")
plot(1995:2018, unlist(lapply(networks_yearly_cepii, FUN = function(x){centralization(x, FUN = "betweenness")})), type = "l", xlab = "Year", ylab = "", main = "CEPII BACI: Betweenness-centrality")
plot(2000:2016, unlist(lapply(networks_yearly_itpd, FUN = function(x){centralization(x, FUN = "betweenness")})), type = "l", xlab = "Year", ylab = "", main = "ITPD Mining/Energy: Betweenness-centrality")

plot(tSnaStats(networks_dynamic_arms, "centralization", FUN = "degree", cmode = "indegree", start = 1995, end = 2018), xlab = "Year", ylab = "", main = "SIPRI: In-centrality")
plot(1995:2018, unlist(lapply(networks_yearly_cepii, FUN = function(x){centralization(x, FUN = "degree", cmode = "indegree")})), type = "l", xlab = "Year", ylab = "", main = "CEPII BACI: In-centrality")
plot(2000:2016, unlist(lapply(networks_yearly_itpd, FUN = function(x){centralization(x, FUN = "degree", cmode = "indegree")})), type = "l", xlab = "Year", ylab = "", main = "ITPD Mining/Energy: In-centrality")



plot(tSnaStats(networks_dynamic_arms, "centralization", FUN = "degree", cmode = "outdegree", start = 1995, end = 2018), xlab = "Year", ylab = "", main = "SIPRI: Out-centrality")
plot(1995:2018, unlist(lapply(networks_yearly_cepii, FUN = function(x){centralization(x, FUN = "degree", cmode = "outdegree")})), type = "l", xlab = "Year", ylab = "", main = "CEPII BACI: Out-centrality")
plot(1995:2018, unlist(lapply(networks_yearly_itpd, FUN = function(x){centralization(x, FUN = "degree", cmode = "outdegree")})), type = "l", xlab = "Year", ylab = "", main = "CITPD Mining/Energy: Out-centrality")

plot.new()
plot.new()
plot.new()


par(mfrow = c(2, 2))
# In- and Outdegree, SIPRI vs BACI
# to fix: degree = 1 should be 0, log transform issues. Simply new labels?
in_degrees1 = unlist(lapply(networks_yearly_sipri[(1995:2005)-1994], FUN = function(x){degree(x,cmode="indegree")}))
in_degrees2 = unlist(lapply(networks_yearly_cepii[(1995:2005)-1994], FUN = function(x){degree(x,cmode="indegree")}))
tmp1 = table(in_degrees1) / length(in_degrees1)
tmp2 = table(in_degrees2) / length(in_degrees2)
dat = list(as.matrix(tmp1), as.matrix(tmp2))

plot(unlist(dat),type="n", log = "x", xlim=c(1,max(sapply(dat,length))), ylab = "Proportion", xlab = "Degree", main = "In-Degree Distribution 1995-2005 \n (log-scale)")
mapply(lines, dat, col=seq_along(dat), type = "b", pch = 20)
legend("topright", c("SIPRI", "CEPII"), lty=1, col=seq_along(dat))


out_degrees1 = unlist(lapply(networks_yearly_sipri[(1995:2005)-1994], FUN = function(x){degree(x,cmode="outdegree")}))
out_degrees2 = unlist(lapply(networks_yearly_cepii[(1995:2005)-1994], FUN = function(x){degree(x,cmode="outdegree")}))
tmp1 = table(out_degrees1) / length(out_degrees1)
tmp2 = table(out_degrees2) / length(out_degrees2)
dat = list(as.matrix(tmp1), as.matrix(tmp2))

plot(unlist(dat), type="n", log = "x", xlim=c(1,max(sapply(dat,length))), ylab = "Proportion", xlab = "Degree", main = "Out-Degree Distribution 1995-2005 \n (log-scale)")
mapply(lines,dat, col=seq_along(dat), type = "b", pch = 20)
legend("topright", c("SIPRI", "CEPII"), lty=1, col=seq_along(dat))


in_degrees1 = unlist(lapply(networks_yearly_sipri[(2005:2018)-1994], FUN = function(x){degree(x,cmode="indegree")}))
in_degrees2 = unlist(lapply(networks_yearly_cepii[(2005:2018)-1994], FUN = function(x){degree(x,cmode="indegree")}))
tmp1 = table(in_degrees1) / length(in_degrees1)
tmp2 = table(in_degrees2) / length(in_degrees2)
dat = list(as.matrix(tmp1), as.matrix(tmp2))

plot(unlist(dat), type="n", log = "x", xlim=c(1,max(sapply(dat,length))), ylab = "Proportion", xlab = "Degree", main = "In-Degree Distribution 2005-2018 \n (log-scale)")
mapply(lines,dat, col=seq_along(dat), type = "b", pch = 20)
legend("topright", c("SIPRI", "CEPII"), lty=1, col=seq_along(dat))


out_degrees1 = unlist(lapply(networks_yearly_sipri[(2005:2018)-1994], FUN = function(x){degree(x,cmode="outdegree")}))
out_degrees2 = unlist(lapply(networks_yearly_cepii[(2005:2018)-1994], FUN = function(x){degree(x,cmode="outdegree")}))
tmp1 = table(out_degrees1) / length(out_degrees1)
tmp2 = table(out_degrees2) / length(out_degrees2)
dat = list(as.matrix(tmp1), as.matrix(tmp2))

plot(unlist(dat), type="n", log = "x", xlim=c(1,max(sapply(dat,length))), ylab = "Proportion", xlab = "Degree", main = "Out-Degree Distribution 2005-2018 \n (log-scale)")
mapply(lines,dat, col=seq_along(dat), type = "b", pch = 20)
legend("topright", c("SIPRI", "CEPII"), lty=1, col=seq_along(dat))




# In- and Outdegree, SIPRI vs Mining
# to fix: degree = 1 should be 0, log transform issues. Simply new labels?
in_degrees1 = unlist(lapply(networks_yearly_sipri[(2000:2008)-1994], FUN = function(x){degree(x,cmode="indegree")}))
in_degrees2 = unlist(lapply(networks_yearly_itpd[(2000:2008)-1999], FUN = function(x){degree(x,cmode="indegree")}))
tmp1 = table(in_degrees1) / length(in_degrees1)
tmp2 = table(in_degrees2) / length(in_degrees2)
dat = list(as.matrix(tmp1), as.matrix(tmp2))

plot(unlist(dat),type="n", log = "x", xlim=c(1,max(sapply(dat,length))), ylab = "Proportion", xlab = "Degree", main = "In-Degree Distribution 1995-2005 \n (log-scale)")
mapply(lines, dat, col=seq_along(dat), type = "b", pch = 20)
legend("topright", c("SIPRI", "ITPD Mining/Energy"), lty=1, col=seq_along(dat))


out_degrees1 = unlist(lapply(networks_yearly_sipri[(2000:2008)-1994], FUN = function(x){degree(x,cmode="outdegree")}))
out_degrees2 = unlist(lapply(networks_yearly_itpd[(2000:2008)-1999], FUN = function(x){degree(x,cmode="outdegree")}))
tmp1 = table(out_degrees1) / length(out_degrees1)
tmp2 = table(out_degrees2) / length(out_degrees2)
dat = list(as.matrix(tmp1), as.matrix(tmp2))

plot(unlist(dat), type="n", log = "x", xlim=c(1,max(sapply(dat,length))), ylab = "Proportion", xlab = "Degree", main = "Out-Degree Distribution 1995-2005 \n (log-scale)")
mapply(lines,dat, col=seq_along(dat), type = "b", pch = 20)
legend("topright", c("SIPRI", "ITPD Mining/Energy"), lty=1, col=seq_along(dat))


in_degrees1 = unlist(lapply(networks_yearly_sipri[(2009:2016)-1994], FUN = function(x){degree(x,cmode="indegree")}))
in_degrees2 = unlist(lapply(networks_yearly_itpd[(2009:2016)-1999], FUN = function(x){degree(x,cmode="indegree")}))
tmp1 = table(in_degrees1) / length(in_degrees1)
tmp2 = table(in_degrees2) / length(in_degrees2)
dat = list(as.matrix(tmp1), as.matrix(tmp2))

plot(unlist(dat), type="n", log = "x", xlim=c(1,max(sapply(dat,length))), ylab = "Proportion", xlab = "Degree", main = "In-Degree Distribution 2005-2018 \n (log-scale)")
mapply(lines,dat, col=seq_along(dat), type = "b", pch = 20)
legend("topright", c("SIPRI", "ITPD Mining/Energy"), lty=1, col=seq_along(dat))


out_degrees1 = unlist(lapply(networks_yearly_sipri[(2009:2016)-1994], FUN = function(x){degree(x,cmode="outdegree")}))
out_degrees2 = unlist(lapply(networks_yearly_itpd[(2009:2016)-1999], FUN = function(x){degree(x,cmode="outdegree")}))
tmp1 = table(out_degrees1) / length(out_degrees1)
tmp2 = table(out_degrees2) / length(out_degrees2)
dat = list(as.matrix(tmp1), as.matrix(tmp2))

plot(unlist(dat), type="n", log = "x", xlim=c(1,max(sapply(dat,length))), ylab = "Proportion", xlab = "Degree", main = "Out-Degree Distribution 2005-2018 \n (log-scale)")
mapply(lines,dat, col=seq_along(dat), type = "b", pch = 20)
legend("topright", c("SIPRI", "ITPD Mining/Energy"), lty=1, col=seq_along(dat))

dev.off()
