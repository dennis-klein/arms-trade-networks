# replication of the btergm model in Thurner et al. 2019, J. of Conflict Resolution

library(network)
library(btergm)


set.seed(1234)
rm(list = ls(all.names = TRUE))


# load data
load("data/out/EX.RData")
load("data/out/country_list.RData")
arms = readRDS("data/out/sipri_tiv.rds")
trade = readRDS("data/out/baci_aggregated.rds")
cinc = readRDS("data/out/nmc_cinc.rds")
conflict = readRDS("data/out/conflict_intrastate.rds")
load("data/out/polity.RData")
load("data/out/real_gdp.RData")
load("data/out/cdist.RData")
load("data/out/atop_alliance.RData")


# model year 2003, lag covariates of interest approprietly 
c(1950:2018)[54] # 2003
present <- EX[, 54] == 1 

for (i in 1:69) {
  arms[[i]] <- arms[[i]][present, present]
  if (i == 69){
    atop_alliance[[i]] <- atop_alliance[[68]]
  } else {
    atop_alliance[[i]] <- atop_alliance[[i]][present, present]
  }
}

for (i in 1:25) {
  trade[[i]] <- trade[[i]][present, present]
}

cinc = cinc[present, ]
conflict = conflict[present, ]
real_gdp = real_gdp[present, ]
polity = polity[present, ]
cdist = cdist[present, present]

networks = list()
alliance = list()
distance = list()
flow = list()
pathdep = list()


# create list of network objects
for (i in 1:5){
  net <- network(arms[[49 + i]], directed = TRUE)
  set.vertex.attribute(net, "gdp", real_gdp[, 49+i-2])
  set.vertex.attribute(net, "cinc", cinc[, 49+i-2])
  set.vertex.attribute(net, "polity2", polity[, 49+i-2])
  set.vertex.attribute(net, "conflict", conflict[, 49+i-2])
  
  networks[[i]] <- net
  
  # path dependency: summed last five arms transfers networks
  AR5 = arms[[49+i-1]] + arms[[49+i-2]] + arms[[49+i-3]] + arms[[49+i-4]] + arms[[49+i-5]]
  pathdep[[i]] <- AR5
  
  alliance[[i]] <- atop_alliance[[49+i-2]]
  distance[[i]] <- cdist
  flow[[i]] <- trade[[4+i-2]]
}


# estimate tergm for the year 2003 without exogenous effects
model_simple <- btergm(
  networks ~ 
    # endogneous effects
    edges +
    gwodegree(1, fixed = TRUE) +
    gwidegree(1, fixed = TRUE) + 
    mutual +
    gwesp(1.5, fixed = TRUE) +
    nodeicov("cinc") + 
    nodeocov("cinc") + 
    absdiff('polity2'),
  parallel = c("multicore"), 
  ncpus = 4
)


# model summary output
sink(file = "model_simple.txt")
summary(model_simple)
sink()


model_full <- btergm(
  networks ~ 
    # endogneous effects
    edges +
    gwodegree(1, fixed = TRUE) +
    gwidegree(1, fixed = TRUE) + 
    mutual +
    gwesp(1.5, fixed = TRUE) +
    
    # exogenous effects
    edgecov(alliance) + 
    nodeocov("gdp") + 
    nodeicov("gdp") +
    edgecov(distance) + 
    edgecov(flow) +
    absdiff('polity2') +
    nodeicov("conflict") +
    nodeicov("cinc") + 
    nodeocov("cinc") + 
    edgecov(pathdep), 
  parallel = c("multicore"), 
  ncpus = 4
  )


# model summary output
sink(file = "model_full.txt")
summary(model_full)
sink()



