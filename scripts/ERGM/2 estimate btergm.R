library(network)
library(ergm)
library(btergm)
library(multilayer.ergm)


# init setup
rm(list = ls(all.names = TRUE))
set.seed(1234)
source("utils/utils.R")
source("utils/construct_header.R")
source("utils/custom_trade_to_binary.R")
path = data_path.get()


# load data 
load(file.path(path, "out/EX.RData"))
load(file.path(path, "out/colony.RData"))
load(file.path(path, "out/country_list.RData"))
load(file.path(path, "out/cdist.RData"))
load(file.path(path, "out/atop_alliance.RData"))
load(file.path(path, "out/milit_exp.RData"))
nmc_cinc = readRDS(file.path(path, "out/nmc_cinc.rds"))
polity = readRDS(file.path(path, "out/polity.rds"))
gdp = readRDS(file.path(path, "out/gdp.rds"))
gdppc = readRDS(file.path(path, "out/gdppc.rds"))
sipri_tiv = readRDS(file.path(path, "out/sipri_tiv.rds"))
trade = readRDS(file.path(path, "out/baci_aggregated.rds"))


# set parameters of analyses
# own model version number for saving 
version = "A"


# trade data starting 1995
# "D" is export relevance, "C" is import relevance
start = 1995
end = 2018
maxlag = 4
n_boot = 1000
n_sim = 10000
type = "D"


# selection of countries included
included = rowSums(EX[, (start:end)-1949]) == length(start:end)
n = sum(included) 


# some data transformations 
# sipri military expenditure is in constant USD$m
# cinc is index which sums up to 0 -> better interpretability using pct points
milit_exp = milit_exp * 1000000
nmc_cinc100 = nmc_cinc * 100 


# setup constraints
free = matrix(c(rep(c(rep(1, n), rep(0, n)), n), rep(c(rep(0, n), rep(1, n)), n)), nrow = 2*n, byrow = T)
diag(free) = 0 
free = network(free, directed = T)


# matrices for kronecker products 
A = matrix(c(1,0,0,0), nrow = 2)
B = matrix(c(0,0,0,1), nrow = 2)


# measure computing time and results lists
start_time = Sys.time()
results_restricted = list()
results_full = list()
cov_restricted = list()
cov_full = list()
predictions_cond_restricted = list()
predictions_cond_full = list()
predictions_sim_restricted = list()
predictions_sim_full = list()
target = list()


# estimate tergms for each year
for (year in (start+maxlag):end){
  
  # set correct indices corresponding data source
  i1 = year - 1949 # Cornelius Replication Data 1950:2018
  i2 = year - 1994 # CEPII BACI 1995:2018


  # init covariates 
  net = list() 
  log_gdp_out_layer1 = list()
  log_gdp_out_layer2 = list()
  log_gdp_in_layer1 = list()
  log_gdp_in_layer2 = list()
  alliance_layer1 = list() 
  alliance_layer2 = list()
  polity_diff_layer1 = list() 
  polity_diff_layer2 = list()
  log_milit_exp_out_layer1 = list()
  log_milit_exp_out_layer2 = list()
  log_milit_exp_in_layer1 = list()
  log_milit_exp_in_layer2 = list()
  pathdep_layer1 = list() 
  pathdep_layer2 = list() 

  
  # dependent multi-layer network
  net[[4]] = to.multiplex(1* (sipri_tiv[[i1]][included, included] > 0), 
                          custom_trade_to_binary(trade[[i2]][included, included], type = type, threshold = 0.01), 
                          output = "network")
  net[[3]] = to.multiplex(1* (sipri_tiv[[i1-1]][included, included] > 0), 
                          custom_trade_to_binary(trade[[i2-1]][included, included], type = type, threshold = 0.01), 
                          output = "network")
  net[[2]] = to.multiplex(1* (sipri_tiv[[i1-2]][included, included] > 0), 
                          custom_trade_to_binary(trade[[i2-2]][included, included], type = type, threshold = 0.01), 
                          output = "network")
  net[[1]] = to.multiplex(1* (sipri_tiv[[i1-3]][included, included] > 0), 
                          custom_trade_to_binary(trade[[i2-3]][included, included], type = type, threshold = 0.01), 
                          output = "network")
  
  # distance
  log_cdist_layer1 = kronecker(A, log(cdist[included, included] + 1)) 
  log_cdist_layer2 = kronecker(B, log(cdist[included, included] + 1))
  colnames(log_cdist_layer1) <- rownames(log_cdist_layer1) <- 1:(2*n)
  colnames(log_cdist_layer2) <- rownames(log_cdist_layer2) <- 1:(2*n)
  
  
  # include exogenous covariates lagged by 1 time period for each net
  # gdp out
  tmp1 = matrix(log(gdp[included, i1-1]), n, n, byrow = FALSE)
  tmp2 = matrix(log(gdp[included, i1-1]), n, n, byrow = FALSE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_gdp_out_layer1[[4]] = tmp1 ; log_gdp_out_layer2[[4]] = tmp2
  
  tmp1 = matrix(log(gdp[included, i1-2]), n, n, byrow = FALSE)
  tmp2 = matrix(log(gdp[included, i1-2]), n, n, byrow = FALSE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_gdp_out_layer1[[3]] = tmp1 ; log_gdp_out_layer2[[3]] = tmp2
  
  tmp1 = matrix(log(gdp[included, i1-3]), n, n, byrow = FALSE)
  tmp2 = matrix(log(gdp[included, i1-3]), n, n, byrow = FALSE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_gdp_out_layer1[[2]] = tmp1 ; log_gdp_out_layer2[[2]] = tmp2
  
  tmp1 = matrix(log(gdp[included, i1-4]), n, n, byrow = FALSE)
  tmp2 = matrix(log(gdp[included, i1-4]), n, n, byrow = FALSE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_gdp_out_layer1[[1]] = tmp1 ; log_gdp_out_layer2[[1]] = tmp2
  
  # gdp in
  tmp1 = matrix(log(gdp[included, i1-1]), n, n, byrow = TRUE)
  tmp2 = matrix(log(gdp[included, i1-1]), n, n, byrow = TRUE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_gdp_in_layer1[[4]] = tmp1 ; log_gdp_in_layer2[[4]] = tmp2
  
  tmp1 = matrix(log(gdp[included, i1-2]), n, n, byrow = TRUE)
  tmp2 = matrix(log(gdp[included, i1-2]), n, n, byrow = TRUE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_gdp_in_layer1[[3]] = tmp1 ; log_gdp_in_layer2[[3]] = tmp2
  
  tmp1 = matrix(log(gdp[included, i1-3]), n, n, byrow = TRUE)
  tmp2 = matrix(log(gdp[included, i1-3]), n, n, byrow = TRUE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_gdp_in_layer1[[2]] = tmp1 ; log_gdp_in_layer2[[2]] = tmp2
  
  tmp1 = matrix(log(gdp[included, i1-4]), n, n, byrow = TRUE)
  tmp2 = matrix(log(gdp[included, i1-4]), n, n, byrow = TRUE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_gdp_in_layer1[[1]] = tmp1 ; log_gdp_in_layer2[[1]] = tmp2
  
  # atop alliance
  tmp1 = kronecker(A, atop_alliance[[i1 - 1]][included, included]) 
  tmp2 = kronecker(B, atop_alliance[[i1 - 1]][included, included])
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  alliance_layer1[[4]] = tmp1; alliance_layer2[[4]] = tmp2
  
  tmp1 = kronecker(A, atop_alliance[[i1 - 2]][included, included]) 
  tmp2 = kronecker(B, atop_alliance[[i1 - 2]][included, included])
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  alliance_layer1[[3]] = tmp1; alliance_layer2[[3]] = tmp2
  
  tmp1 = kronecker(A, atop_alliance[[i1 - 3]][included, included]) 
  tmp2 = kronecker(B, atop_alliance[[i1 - 3]][included, included])
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  alliance_layer1[[2]] = tmp1; alliance_layer2[[2]] = tmp2
  
  tmp1 = kronecker(A, atop_alliance[[i1 - 4]][included, included]) 
  tmp2 = kronecker(B, atop_alliance[[i1 - 4]][included, included])
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  alliance_layer1[[1]] = tmp1; alliance_layer2[[1]] = tmp2
  
  # absolute polity diff
  tmp1 = kronecker(A, abs(outer(polity[included, i1-1], polity[included, i1-1],'-')))
  tmp2 = kronecker(B, abs(outer(polity[included, i1-1], polity[included, i1-1],'-')))
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  polity_diff_layer1[[4]] = tmp1; polity_diff_layer2[[4]] = tmp2
  
  tmp1 = kronecker(A, abs(outer(polity[included, i1-2], polity[included, i1-2],'-')))
  tmp2 = kronecker(B, abs(outer(polity[included, i1-2], polity[included, i1-2],'-')))
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  polity_diff_layer1[[3]] = tmp1; polity_diff_layer2[[3]] = tmp2
  
  tmp1 = kronecker(A, abs(outer(polity[included, i1-3], polity[included, i1-3],'-')))
  tmp2 = kronecker(B, abs(outer(polity[included, i1-3], polity[included, i1-3],'-')))
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  polity_diff_layer1[[2]] = tmp1; polity_diff_layer2[[2]] = tmp2
  
  tmp1 = kronecker(A, abs(outer(polity[included, i1-4], polity[included, i1-4],'-')))
  tmp2 = kronecker(B, abs(outer(polity[included, i1-4], polity[included, i1-4],'-')))
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  polity_diff_layer1[[1]] = tmp1; polity_diff_layer2[[1]] = tmp2
  
  # military expenditure out
  tmp1 = matrix(log(milit_exp[included, i1-1]), n, n, byrow = FALSE)
  tmp2 = matrix(log(milit_exp[included, i1-1]), n, n, byrow = FALSE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_milit_exp_out_layer1[[4]] = tmp1 ; log_milit_exp_out_layer2[[4]] = tmp2
  
  tmp1 = matrix(log(milit_exp[included, i1-2]), n, n, byrow = FALSE)
  tmp2 = matrix(log(milit_exp[included, i1-2]), n, n, byrow = FALSE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_milit_exp_out_layer1[[3]] = tmp1 ; log_milit_exp_out_layer2[[3]] = tmp2
  
  tmp1 = matrix(log(milit_exp[included, i1-3]), n, n, byrow = FALSE)
  tmp2 = matrix(log(milit_exp[included, i1-3]), n, n, byrow = FALSE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_milit_exp_out_layer1[[2]] = tmp1 ; log_milit_exp_out_layer2[[2]] = tmp2
  
  tmp1 = matrix(log(milit_exp[included, i1-4]), n, n, byrow = FALSE)
  tmp2 = matrix(log(milit_exp[included, i1-4]), n, n, byrow = FALSE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_milit_exp_out_layer1[[1]] = tmp1 ; log_milit_exp_out_layer2[[1]] = tmp2
  
  # military expenditure in
  tmp1 = matrix(log(milit_exp[included, i1-1]), n, n, byrow = TRUE)
  tmp2 = matrix(log(milit_exp[included, i1-1]), n, n, byrow = TRUE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_milit_exp_in_layer1[[4]] = tmp1 ; log_milit_exp_in_layer2[[4]] = tmp2
  
  tmp1 = matrix(log(milit_exp[included, i1-2]), n, n, byrow = TRUE)
  tmp2 = matrix(log(milit_exp[included, i1-2]), n, n, byrow = TRUE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_milit_exp_in_layer1[[3]] = tmp1 ; log_milit_exp_in_layer2[[3]] = tmp2
  
  tmp1 = matrix(log(milit_exp[included, i1-3]), n, n, byrow = TRUE)
  tmp2 = matrix(log(milit_exp[included, i1-3]), n, n, byrow = TRUE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_milit_exp_in_layer1[[2]] = tmp1 ; log_milit_exp_in_layer2[[2]] = tmp2
  
  tmp1 = matrix(log(milit_exp[included, i1-4]), n, n, byrow = TRUE)
  tmp2 = matrix(log(milit_exp[included, i1-4]), n, n, byrow = TRUE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_milit_exp_in_layer1[[1]] = tmp1 ; log_milit_exp_in_layer2[[1]] = tmp2
  
  # define path dependency
  tmp1 = 1 * (sipri_tiv[[i1-1]][included, included] > 0)
  tmp2 = custom_trade_to_binary(trade[[i2-1]][included, included], type = type, threshold = 0.01)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  pathdep_layer1[[4]] = tmp1 ; pathdep_layer2[[4]] = tmp2
  
  tmp1 = 1 * (sipri_tiv[[i1-2]][included, included] > 0)
  tmp2 = custom_trade_to_binary(trade[[i2-2]][included, included], type = type, threshold = 0.01)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  pathdep_layer1[[3]] = tmp1 ; pathdep_layer2[[3]] = tmp2
  
  tmp1 = 1 * (sipri_tiv[[i1-3]][included, included] > 0)
  tmp2 = custom_trade_to_binary(trade[[i2-3]][included, included], type = type, threshold = 0.01)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  pathdep_layer1[[2]] = tmp1 ; pathdep_layer2[[2]] = tmp2
  
  tmp1 = 1* (sipri_tiv[[i1-4]][included, included] > 0)
  tmp2 = custom_trade_to_binary(trade[[i2-4]][included, included], type = type, threshold = 0.01)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  pathdep_layer1[[1]] = tmp1 ; pathdep_layer2[[1]] = tmp2
  
  
  # specify base model
  model = net ~
    edges_layer(layer = 1) +
    edges_layer(layer = 2) + 
    mutual(same = "layer.mem", diff = TRUE) +
    gwidegree(decay = 1, fixed = TRUE, attr = "layer.mem") +
    gwodegree(decay = 1, fixed = TRUE, attr = "layer.mem") +
    gwesp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
    gwesp_layer(decay = 1.5, fixed = TRUE, layer = 2) +
    edgecov(log_cdist_layer1) + 
    edgecov(log_cdist_layer2) +
    edgecov(log_gdp_in_layer1) +
    edgecov(log_gdp_in_layer2) +
    edgecov(log_gdp_out_layer1) +
    edgecov(log_gdp_out_layer2) +
    edgecov(alliance_layer1) +
    edgecov(alliance_layer2) +
    edgecov(polity_diff_layer1) +
    edgecov(polity_diff_layer2) +
    edgecov(log_milit_exp_in_layer1) +
    edgecov(log_milit_exp_in_layer2) +
    edgecov(log_milit_exp_out_layer1) +
    edgecov(log_milit_exp_out_layer2) +
    edgecov(pathdep_layer1) +
    edgecov(pathdep_layer2)
  
  
  # estimate restricted model without cross-layer interdependence, save and simulate
  fit = btergm(model, parallel = "snow", ncpus = 4, verbose = FALSE, R = n_boot)
  results_restricted[[year + 1 - start - maxlag]] = confint(fit)
  
  sim = simulate(fit, index = 4, nsim = n_sim, output = "stats", constraints = ~fixallbut(free), verbose = F)
  cov_restricted[[year + 1 - start - maxlag]] = cov(sim)
  

  # estimate full model with cross-layer interdependence, save and simulate 
  model = update(model, ~ . + duplexdyad(c("e", "h"), layers = list(1, 2)))
  fit = btergm(model, parallel = "snow", ncpus = 4, verbose = FALSE, R = n_boot)
  results_full[[year + 1 - start - maxlag]] = confint(fit)
  
  sim = simulate(fit, index = 4, nsim = n_sim, output = "stats", constraints = ~fixallbut(free), verbose = F)
  cov_full[[year + 1 - start - maxlag]] = cov(sim)
  
  
  # mcmle estimation -> constraints & cd method throw errors bc of object type
  if(FALSE){ 
    fit = mtergm(model, estimate = c("CD"), constraints = ~fixallbut(free), 
                 control = control.ergm(parallel = 4, CD.nsteps = 8, CD.multiplicity = 1))
  }
  
  
  # estimate benchmark logit model 
  # To Do
  
  
  # compute predictions for t+1
  if(year != end){
    
    # load corresponding t+1 covariates
    future = to.multiplex(1* (sipri_tiv[[i1+1]][included, included] > 0), 
                          custom_trade_to_binary(trade[[i2+1]][included, included], type = type, threshold = 0.01), 
                          output = "network")
    # gdp out
    tmp1 = matrix(log(gdp[included, i1]), n, n, byrow = FALSE)
    tmp2 = matrix(log(gdp[included, i1]), n, n, byrow = FALSE)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    log_gdp_out_layer1 = tmp1 ; log_gdp_out_layer2 = tmp2
    
    # gdp in
    tmp1 = matrix(log(gdp[included, i1]), n, n, byrow = TRUE)
    tmp2 = matrix(log(gdp[included, i1]), n, n, byrow = TRUE)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    log_gdp_in_layer1 = tmp1 ; log_gdp_in_layer2 = tmp2
    
    # atop alliance
    tmp1 = kronecker(A, atop_alliance[[i1]][included, included]) 
    tmp2 = kronecker(B, atop_alliance[[i1]][included, included])
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    alliance_layer1 = tmp1; alliance_layer2 = tmp2
    
    # absolute polity diff
    tmp1 = kronecker(A, abs(outer(polity[included, i1], polity[included, i1],'-')))
    tmp2 = kronecker(B, abs(outer(polity[included, i1], polity[included, i1],'-')))
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    polity_diff_layer1 = tmp1; polity_diff_layer2 = tmp2
    
    # military expenditure out
    tmp1 = matrix(log(milit_exp[included, i1]), n, n, byrow = FALSE)
    tmp2 = matrix(log(milit_exp[included, i1]), n, n, byrow = FALSE)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    log_milit_exp_out_layer1 = tmp1 ; log_milit_exp_out_layer2 = tmp2
    
    # military expenditure in
    tmp1 = matrix(log(milit_exp[included, i1]), n, n, byrow = TRUE)
    tmp2 = matrix(log(milit_exp[included, i1]), n, n, byrow = TRUE)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    log_milit_exp_in_layer1 = tmp1 ; log_milit_exp_in_layer2 = tmp2
    
    # define path dependency
    tmp1 = 1 * (sipri_tiv[[i1]][included, included] > 0)
    tmp2 = custom_trade_to_binary(trade[[i2]][included, included], type = type, threshold = 0.01)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    pathdep_layer1 = tmp1 ; pathdep_layer2 = tmp2
    
    
    # predict t+1 without cross-layer dependence
    model_future = update(model, future ~ . - duplexdyad(c("e", "h"), layers = list(1, 2)))
    coef = results_restricted[[year + 1 - start - maxlag]][, 1]
    names(coef) = gsub("[[i]]", "", names(coef), fixed = T)
    
    pred = predict(model_future, 
                   theta = coef,
                   conditional = T,
                   type = "response", 
                   output = "matrix", 
                   constraints = ~fixallbut(free)
    )
    
    predictions_cond_restricted[[year + 1 - start - maxlag]] = pred
    
    
    # same but with sampling using simulate function, non-conditional
    sim = simulate(model_future, nsim = n_sim, 
                   output = "network", coef = coef, 
                   constraints = ~fixallbut(free),
                   verbose = F
    )
    sim = lapply(sim, FUN = as.matrix.network)
    sim = Reduce('+', sim)/n_sim
    
    predictions_sim_restricted[[year + 1 - start - maxlag]] = sim
    
    
    # predict t+1 with cross-layer dependence
    model_future = update(model_future, ~ . + duplexdyad(c("e", "h"), layers = list(1, 2)))
    coef = results_full[[year + 1 - start - maxlag]][, 1]
    names(coef) = gsub("[[i]]", "", names(coef), fixed = T)
    
    pred = predict(model_future, 
                   theta = coef,
                   conditional = T,
                   type = "response", 
                   nsim = n_sim, 
                   output = "matrix"
    )
    
    predictions_cond_full[[year + 1 - start - maxlag]] = pred
    
    
    # same but with sampling using simulate function, non-conditional
    sim = simulate(model_future, nsim = n_sim, 
                   output = "network", coef = coef, 
                   constraints = ~fixallbut(free), 
                   verbose = F
    )
    sim = lapply(sim, FUN = as.matrix.network)
    sim = Reduce('+', sim)/n_sim
    
    predictions_sim_full[[year + 1 - start - maxlag]] = sim
    
    
    # save target function
    target[[year + 1 - start - maxlag]] = as.matrix.network(future)
    
  }
  
  print(paste0("Year ", year , " finished."))
}


end_time = Sys.time()
duration = end_time - start_time


# save results for evaluation function
save(start, end, maxlag, type, duration, 
     results_restricted, results_full,
     cov_restricted, cov_full, target, 
     predictions_cond_restricted, predictions_cond_ful,
     predictions_sim_restricted, predictions_sim_full,
     n_sim, n_boot, n, 
     file = paste0(path, "/models/ERGM/model_estimated_", version, ".RData"))

#load(file = paste0(path, "/models/ERGM/model_estimated_", version, ".RData")) 


